#include "stdafx.h"

#include "compression.h"

#include "twigg/linalg.h"
#include "twigg/random.h"
#include "twigg/ioutil.h"
#include "twigg/EulerAngles.h"

#ifdef NAG
#include <nag.h>
#include <nag_stdlib.h>
#include <nage04.h>
#endif

#include <numeric>
#include <algorithm>

//#define INSTRUMENT_COMPRESSION

namespace planning {

template <typename Pr>
class Compare1st
{
public:
	bool operator()( const Pr& left, const Pr& right ) const
	{
		return std::less<typename Pr::first_type>()( left.first, right.first );
	}
};

#ifdef NAG
class NAGOptions
{
public:
	NAGOptions()
	{
		nag_opt_init(&options_);

		options_.inf_bound = 1e20;
		options_.prob = Nag_LP;
		options_.print_level = Nag_NoPrint;
		options_.list = FALSE;
	}

	Nag_E04_Opt* get()
	{
		return &options_;
	}

	~NAGOptions()
	{
		NagError fail;
		INIT_FAIL(fail);

		nag_opt_free( &options_, "", &fail );

		if( fail.code != NE_NOERROR )
		{
			throw SolverException("Error freeing NAG_Options struct: "
				+ std::string(fail.message) );
		}
	}

private:
	Nag_E04_Opt options_;
};
#endif

#ifdef NAG
std::pair<double, std::vector<double> > leastInfinityNorm( 
	const fortran_matrix& A, const std::vector<double>& b )
{
	const double epsilon = 1e-7;

	std::vector<double> x;
	{
		fortran_matrix temp_A( A );
		x = leastSquares_QR( temp_A, &b[0] );
		assert( x.size() == A.ncols() );
		std::vector<double> err( b );
		prodAdd( A, scaled(x, -1.0), err );
		double norminf_err = norminf( err );
		if( norminf_err < epsilon || b.size() <= A.ncols() )
			return std::make_pair( norminf_err, x );

		x.push_back( norminf_err );
	}

	assert( A.nrows() == b.size() );
	Integer n = A.ncols() + 1;
	Integer nclin = 2*b.size();
	Integer tda = n;

	NAGOptions options;
	const double inf_bound = options.get()->inf_bound;

	c_matrix a( nclin, tda );
	std::vector<double> bl( nclin + n, -inf_bound );
	bl[ n-1 ] = 0.0;
	std::vector<double> bu( nclin + n, inf_bound );

	for( size_t i = 0; i < b.size(); ++i )
	{
		for( size_t j = 0; j < A.ncols(); ++j )
		{
			a(2*i, j) = A(i, j);
			a(2*i+1, j) = A(i, j);
		}

		a(2*i, n-1) = 1.0;
		bl[n + 2*i] = b[i];
		assert( bl[n + 2*i] < bu[n + 2*i] );

		a(2*i + 1, n-1) = -1.0;
		bu[n + 2*i + 1] = b[i];
		assert( bl[n + 2*i+1] < bu[n + 2*i+1] );
	}

	// this scalar should be strictly less than our error
	//   tolerance to ensure that it doesn't affect the outcome
	std::vector<double> cvec( n, 1e-7 );
	cvec[ n-1 ] = 1.0;

	NagError fail;
	INIT_FAIL(fail);
	Nag_Comm* comm = NAGCOMM_NULL;

	double objf;
	nag_opt_lp( 
		n,
		nclin,
		a.data(),
		tda,
		&bl[0],
		&bu[0],
		&cvec[0],
		&x[0],
		&objf,
		options.get(),
		comm,
		&fail );

	if( fail.code != NE_NOERROR )
	{
		if( fail.code == NW_SOLN_NOT_UNIQUE )
		{
			std::cerr << "Warning: NAG solution not unique." << std::endl;
		}
		else
		{
			throw SolverException("NAG linear program solver failed: "
				+ std::string(fail.message) );
		}
	}

	double error = x[ n-1 ];
	std::vector<double> result( x.begin(), x.begin() + n-1 );

	{
		std::vector<double> err( b );
		prodAdd( A, scaled(result, -1.0), err );
		double norminf_err = norminf( err );
		assert( fabs(error - norminf_err) < 1e-4 );
	}

	return std::make_pair( error, result );
}
#endif

template <typename T>
T integerPower( const T& value, size_t exponent )
{
	T result = boost::numeric_cast<T>(1);
	for( size_t i = 0; i < exponent; ++i )
		result *= value;
	return result;
}

std::pair<float, int> fitSpline( const std::vector<float>& points,
	float predicted, float epsilon, float quadConstant
#ifdef INSTRUMENT_COMPRESSION
		, WaypointTracker& tracker
#endif
 		)
{
	// won't use startFrame as a point
	size_t numPoints = points.size();
	assert( numPoints > 0 );

	size_t order = 2;

	fortran_matrix A( numPoints, order-1 );
	
	std::vector<double> b( numPoints );

	{
#ifdef INSTRUMENT_COMPRESSION
		boost::shared_ptr<WaypointTracker::Span> span(
			tracker.newSpan("Set up matrix") );
#endif

		for( size_t i = 0; i < numPoints; ++i )
		{
			for( size_t j = 0; j < (order-1); ++j )
				A(i, j) = integerPower( i+1, j+1 );
	
			b.at(i) = points.at(i) - predicted;
			b.at(i) -= quadConstant*boost::numeric_cast<double>((i+1)*(i+1));
		}
	}

	std::vector<double> res;

	{
#ifdef INSTRUMENT_COMPRESSION
		boost::shared_ptr<WaypointTracker::Span> span(
			tracker.newSpan("Solve matrix") );
#endif
		res = leastSquares_QR( A, &b[0] );
	}

	{
#ifdef INSTRUMENT_COMPRESSION
		boost::shared_ptr<WaypointTracker::Span> span(
			tracker.newSpan("Quantize coeffs") );
#endif

		std::vector<int> intCoeffs;
		intCoeffs.reserve( res.size() );
		for( size_t i = 0; i < res.size(); ++i )
		{
			double temp = res[i] / epsilon;
			double intpart;
			double fracpart = modf( temp, &intpart );
			if( fabs(fracpart) > 0.5 )
				intpart += signum(fracpart);
	
			intCoeffs.push_back( boost::numeric_cast<int>(intpart) );
			float reconstructedCoeff = epsilon 
				* boost::numeric_cast<float>(intCoeffs.back());
	
			float error = fabs( reconstructedCoeff - res[i] );
			assert( error <= epsilon );
		}

		float maxErr = 0.0f;
		for( size_t i = 0; i < A.nrows(); ++i )
		{
			float sum = predicted;
	
			for( size_t j = 0; j < intCoeffs.size(); ++j )
			{
				float reconstructedCoeff = epsilon 
					* boost::numeric_cast<float>(intCoeffs[j]);
				sum += integerPower( i+1, j+1 ) * reconstructedCoeff;
			}
	
			sum += quadConstant*boost::numeric_cast<float>((i+1)*(i+1));
	
			float actualPos = points.at(i);
			maxErr = std::max<float>( maxErr, fabs(sum - actualPos) );
			//assert( maxErr <= res.first + (epsilon*boost::numeric_cast<double>(i+1)) );
		}

//	assert( fabs(maxErr - res.first) < 1e-4 );
		assert( intCoeffs.size() == 1 );
		return std::make_pair( maxErr, intCoeffs.back() );
	}
}

boost::tuple< float, TinyVec<int, 3>, TinyVec<unsigned char, 3> >
	fitSpline( const std::vector<vl::Vec3f>& points,
		size_t startFrame, size_t endFrame, const vl::Vec3f& predicted, float epsilon, float quadConstant
#ifdef INSTRUMENT_COMPRESSION
		, WaypointTracker& tracker
#endif
 		)
{
	TinyVec<int, 3> linearCoeffs;
	TinyVec<unsigned char, 3> quadraticCoeffs(0, 0, 0);
	float maxError = 0.0f;

	for( size_t iCoord = 0; iCoord < 3; ++iCoord )
	{
		std::vector<float> coords;
		coords.reserve( endFrame - startFrame );
		for( size_t i = startFrame + 1; i <= endFrame; ++i )
		{
			coords.push_back( (points.at(i))[iCoord] );
		}

		std::pair<float, int> res = fitSpline( 
			coords, predicted[iCoord], epsilon, 0.0f
#ifdef INSTRUMENT_COMPRESSION
			, tracker
#endif
			 );

		if( iCoord == 1 )
		{
			std::pair<float, int> newRes = fitSpline( 
				coords, predicted[iCoord], epsilon, quadConstant 
#ifdef INSTRUMENT_COMPRESSION
				, tracker
#endif
				);
			if( newRes.first < res.first )
			{
				res = newRes;
				quadraticCoeffs[iCoord] = 1;
			}
			else
			{
				quadraticCoeffs[iCoord] = 0;
			}
		}
		else
			quadraticCoeffs[iCoord] = 0;

		linearCoeffs[iCoord] = res.second;
		maxError = std::max( maxError, res.first );
	}

	return boost::make_tuple( maxError, linearCoeffs, quadraticCoeffs );
}


vl::Vec3d toScaledAxis( const Quaternion<double>& quat )
{
	vl::Vec3d axis;
	double angle;
	quat.toAxisAngle( axis, angle );
	return axis*angle;
}

/*
NxQuat toQuat( const vl::Vec3f& scaledAxis )
{
	float epsilon = 1e-5;
	NxQuat result( NxVec3(0.0f), 1.0f );
	if( vl::sqrlen(scaledAxis) > epsilon )
	{
		NxVec3 axis = toNxVec( vl::norm(scaledAxis) );
		float angle = vl::len( scaledAxis );
		result.fromAngleAxis( toDegrees(angle), axis );
	}

	return result;
}
*/


// computes the rotation R such that left*R = right
vl::Vec3d differentialRotation( const Quaternion<double>& left, const Quaternion<double>& right )
{
	assert( !isBad(left) );
	assert( !isBad(right) );

	Quaternion<double> resultRot = right * inverse(left);
	assert( !isBad(resultRot) );
    
	vl::Vec3d result = toScaledAxis( resultRot );

	{
		Quaternion<double> test = Quaternion<double>(result) * left;
		Quaternion<double> error = inverse(test) * right;
		double scalarError = arg( error );
		assert( scalarError < 1e-5 );
	}

	return result;
}

std::pair< float, TinyVec<int, 3> > findQuaternionInterp( const std::vector< Quaternion<double> >& rotations,
	size_t startFrame, size_t endFrame, const Quaternion<double>& predicted, float epsilon )
{
	const Quaternion<double>& current = rotations[endFrame];
	vl::Vec3d diff = differentialRotation( predicted, current );
	diff /= boost::numeric_cast<double>( endFrame - startFrame );

	TinyVec<int, 3> intDiff;

	// round off diff:
	for( vl::Int i = 0; i < 3; ++i )
	{
		double temp = diff[i] / epsilon;
		double intPart;
		double fracpart = modf( temp, &intPart );

		if( fabs(fracpart) > 0.5 )
			intPart += signum(fracpart);

		intDiff[i] = boost::numeric_cast<int>(intPart);
		float reconstructedCoeff = epsilon
			* boost::numeric_cast<float>(intDiff[i]);

		float error = fabs( reconstructedCoeff - diff[i] );
		assert( error <= epsilon );
		diff[i] = reconstructedCoeff;
	}

	double maxError = 0.0;
	for( size_t iFrame = startFrame+1; iFrame <= endFrame; ++iFrame )
	{
		vl::Vec3d currentDiffRot = 
			diff*boost::numeric_cast<double>( iFrame - startFrame );

		Quaternion<double> predictedRot = Quaternion<double>(currentDiffRot) * predicted;

		Quaternion<double> errorQuat = inverse(predictedRot) * rotations[iFrame];

		maxError = std::max( maxError, fabs( arg(errorQuat) ) );
	}

	return std::make_pair( maxError, intDiff );
}

template <typename T>
T log2( T value )
{
	return log(value) / log(2.0);
}

template <bool B>
struct abs_signed_selector;

template <>
struct abs_signed_selector<true>
{
	template <typename T>
	T operator()( const T& value ) const
	{
		if( value < static_cast<T>(0) )
			return -value;
		else
			return value;
	}
};

template <>
struct abs_signed_selector<false>
{
	template <typename T>
	T operator()( const T& value ) const
	{
		return value;
	}
};

template <typename T>
T generic_abs( const T& value )
{
	abs_signed_selector<std::numeric_limits<T>::is_signed> selector;
	return selector( value );
}

template <typename T>
size_t find_n( const std::vector<T>& values )
{
	size_t total = 0;
	for( typename std::vector<T>::const_iterator iter = values.begin();
		iter != values.end(); ++iter )
		total += generic_abs(*iter);
		
	double mean = boost::numeric_cast<double>( total ) / 
		boost::numeric_cast<double>( values.size() );

	size_t n = boost::numeric_cast<size_t>( std::max<double>(log2(log(2.0) * mean), 0.0) );
	return n;
}

/*
template <typename T>
std::vector<char> encode( const std::vector<T>& differences )
{
	boost::dynamic_bitset<char> result;

	const size_t n = find_n( differences );
	assert( n <= boost::numeric::bounds<char>::highest() );

	for( typename std::vector<T>::const_iterator iter = values.begin();
		iter != values.end(); ++iter )
	{
		size_t startBit = result.size();
		bool sign = (*iter < 0);
		size_t absVal = boost::numeric_cast<size_t>( generic_abs(*iter) );

		for( size_t i = 0; i < n; ++i )
		{
			result.push_back( absVal & 1 );
			absVal = absVal >> 1;
		}

		result.push_back( 1 );
		//assert( absVal < 100 );
		for( size_t i = 0; i < absVal; ++i )
			result.push_back( 0 );

		if( std::numeric_limits<T>::is_signed )
			result.push_back( sign ? 1 : 0 );

		// now need to reverse in place
		size_t endBit = result.size() - 1;
		for( size_t i = startBit; i < (result.size() - startBit) / 2; ++i )
			std::swap( result[i], result[endBit - i] );
	}

	std::vector<char> charResult;
	charResult.reserve( result.num_blocks() + 1 );
	charResult.push_back( n );
	to_block_range( result, std::back_inserter(charResult) );

	std::vector<char> asChars( result.num_blocks() );
	to_block_range( result, asChars.begin() );
	return asChars;
}
*/
template <typename T>
std::deque<char> encodeAsChars( const std::vector<T>& values, size_t n )
{
	std::deque<char> result;

	for( typename std::vector<T>::const_iterator iter = values.begin();
		iter != values.end(); ++iter )
	{
		std::deque<char> current;

		bool sign = (*iter < 0);
		size_t absVal = boost::numeric_cast<size_t>( generic_abs(*iter) );

		for( size_t i = 0; i < n; ++i )
		{
			current.push_front( absVal & 1 );
			absVal = absVal >> 1;
		}

		current.push_front( 1 );
		//assert( absVal < 100 );
		for( size_t i = 0; i < absVal; ++i )
			current.push_front( 0 );

		if( std::numeric_limits<T>::is_signed )
			current.push_front( sign ? 1 : 0 );

		std::copy( current.begin(), current.end(), 
			std::back_inserter(result) );
	}

	return result;
}

template <typename T>
std::vector<char> encode( const std::vector<T>& differences )
{
	size_t n = find_n( differences );
	std::deque<char> charArray = encodeAsChars( differences, n );

	assert( n <= boost::numeric::bounds<unsigned char>::highest() );
	std::vector<char> result;
	result.reserve( charArray.size() / 8 + 2 );
	result.push_back( n );

    for( size_t start = 0; start < charArray.size(); start += 8 )
	{
		char current = 0;
		for( size_t i = 0; (i < 8) && (start + i) < charArray.size(); ++i )
			current = current | (charArray[start+i] << (7-i));

		result.push_back( current );
	}

	return result;
}


class BitStream
{
public:
	BitStream( const std::vector<char>& values, std::vector<char>::size_type initialPos = 0 )
		: values_(values), pos_(initialPos), bit_(7) {}

	unsigned char next()
	{
		// should do something more intelligent here:
		if( pos_ >= values_.size() )
			throw IOException( "Corrupted file" );

		unsigned char result = (values_[pos_] >> bit_) & 1;
		if( bit_ == 0 )
		{
			bit_ = 7;
			++pos_;
		}
		else
		{
			--bit_;
		}

		return result;
	}

	bool hasNext() const
	{
		return (pos_ < values_.size());
	}

private:
	const std::vector<char>& values_;
	std::vector<char>::size_type pos_;
	unsigned char bit_;
};

std::vector<int> decode( const std::vector<char>& encoded, std::vector<char>::size_type pos, bool isSigned )
{
	size_t n = encoded[pos];
	pos += 1;

	std::vector<int> result;

	BitStream bitStream( encoded, pos );
	while( bitStream.hasNext() )
	{
		unsigned char sign = 0;
		if( isSigned )
			sign = bitStream.next();

		size_t leading = 0;
		if( !bitStream.hasNext() )
			return result;

		while( bitStream.next() != 1 )
		{
			if( !bitStream.hasNext() )
				return result;
			++leading;
		}

		int num = 0;
		for( size_t i = 1; i <= n; ++i )
			num = num | (boost::numeric_cast<int>(bitStream.next()) << (n-i));

		num = num | (leading << n);

		if( sign )
			result.push_back( -num );
		else
			result.push_back( num );

	}

	return result;
}

template <typename T>
std::vector<char> encode( 
	const std::vector<T>& values, 
	float initialPos, 
    float endingPos, 
	float epsilon )
{
	std::vector<char> encoded = encode( values );

	std::vector<char> result;
	result.reserve( sizeof(size_t) + sizeof(float) + sizeof(float) + sizeof(char)*encoded.size() );

	const char* initialPosPtr = reinterpret_cast<const char*>( &initialPos );
	std::copy( initialPosPtr, initialPosPtr + sizeof(float), std::back_inserter(result) );

	const char* endingPosPtr = reinterpret_cast<const char*>( &endingPos );
	std::copy( endingPosPtr, endingPosPtr + sizeof(float), std::back_inserter(result) );

	const char* epsilonPtr = reinterpret_cast<const char*>( &epsilon );
	std::copy( epsilonPtr, epsilonPtr + sizeof(float), std::back_inserter(result) );

	std::copy( encoded.begin(), encoded.end(), 
		std::back_inserter(result) );

	return result;
}

std::vector<float> decode( const std::vector<char>& values, float& initialPos, float& finalPos, bool incremental, bool isSigned = true )
{
	std::vector<char>::size_type pos = 0;

	const float* initialPosPtr = reinterpret_cast<const float*>( &values[pos] );
	initialPos = *initialPosPtr;
	pos += sizeof(float);

	const float* finalPosPtr = reinterpret_cast<const float*>( &values[pos] );
	finalPos = *finalPosPtr;
	pos += sizeof(float);

	const float* epsilonPtr = reinterpret_cast<const float*>( &values[pos] );
	float epsilon = *epsilonPtr;
	pos += sizeof(float);

	std::vector<int> intValues = decode( values, pos, isSigned );
    
	if( incremental )
		std::partial_sum( intValues.begin(), intValues.end(), intValues.begin() );

	std::vector<float> result;
	result.reserve( intValues.size() );
	for( size_t i = 0; i < intValues.size(); ++i )
	{
		result.push_back( boost::numeric_cast<double>(intValues[i]) * epsilon );
	}

	return result;
}

PiecewisePath::PiecewisePath()
{
}

PiecewisePath::PiecewisePath( const CompressedPath& path, const PiecewisePath& other, bool backwards )
{
	if( !path.positionTimes.empty() )
	{
		std::vector<int> positionTimes = 
			decode( path.positionTimes, 0, false );
		std::partial_sum( positionTimes.begin(), positionTimes.end(),
			positionTimes.begin() );

		vl::Vec3f startingPos;
  		vl::Vec3f endingPos;
		boost::array< std::vector<float>, 3 > positionLinearCoeffs;
		for( size_t i = 0; i < 3; ++i )
		{
			positionLinearCoeffs[i] = decode( path.positionLinearCoeffs[i], startingPos[i], endingPos[i], true, true );
			assert( positionLinearCoeffs[i].size() == positionTimes.size() );
		}

		boost::array< std::vector<float>, 3 > positionQuadraticCoeffs;
		for( size_t i = 0; i < 3; ++i )
		{
			float dummy1, dummy2;
			if( !path.positionQuadraticCoeffs[i].empty() )
				positionQuadraticCoeffs[i] = decode( path.positionQuadraticCoeffs[i], dummy1, dummy2, false, false );
			else
				positionQuadraticCoeffs[i] = std::vector<float>( positionTimes.size(), 0.0f );

			assert( positionQuadraticCoeffs[i].size() == positionTimes.size() );
		}

		vl::Vec3f current = startingPos;
        size_t prevTime = path.startFrame;
		boost::array<CoeffList, 3> allCoeffs = {{ vl::vl_0, vl::vl_0, vl::vl_0 }};
		positionCoeffs_.reserve( positionTimes.size() );
		for( size_t i = 0; i < positionTimes.size(); ++i )
		{
			size_t currentTime = path.startFrame + positionTimes[i];
			for( size_t iCoord = 0; iCoord < 3; ++iCoord )
			{
				CoeffList& coeffs = allCoeffs[iCoord];
				double timeDiff = boost::numeric_cast<double>(currentTime - prevTime);
				current[iCoord] += coeffs[1]*timeDiff + coeffs[2]*timeDiff*timeDiff;
			}

			for( size_t iCoord = 0; iCoord < 3; ++iCoord )
			{
				CoeffList& coeffs = allCoeffs[iCoord];
				coeffs[0] = current[iCoord];
				coeffs[1] = (positionLinearCoeffs[iCoord])[i];
				coeffs[2] = (positionQuadraticCoeffs[iCoord])[i];
			}

			positionCoeffs_.push_back( std::make_pair(currentTime, allCoeffs) );
			prevTime = currentTime;
		}

        for( int i = 0; i < 3; ++i )
            (positionCoeffs_.back().second[i])[0] = endingPos[i];
	}

	if( !path.rotationTimes.empty() )
	{
		std::vector<int> rotationTimes = 
			decode( path.rotationTimes, 0, false );
		std::partial_sum( rotationTimes.begin(), rotationTimes.end(),
			rotationTimes.begin() );

		boost::array< std::vector<float>, 3 > rotationCoeffs;
		vl::Vec3f startingRot;
        vl::Vec3f endingRot;
		for( size_t i = 0; i < 3; ++i )
			rotationCoeffs[i] = decode( path.rotationCoeffs[i], startingRot[i], endingRot[i], true, true );

		Quaternion<double> current( toVec3d(startingRot) );
		Quaternion<double> endRot( toVec3d(endingRot) );
		size_t prevTime = path.startFrame;
		vl::Vec3f differentialRotation( vl::vl_0 );
		for( size_t i = 0; i < rotationTimes.size(); ++i )
		{
			size_t currentTime = path.startFrame + rotationTimes[i];
			double timeDiff = boost::numeric_cast<double>( currentTime - prevTime );
			current = Quaternion<double>( timeDiff * toVec3d(differentialRotation) ) * current;

			for( size_t iCoord = 0; iCoord < 3; ++iCoord )
				differentialRotation[iCoord] = (rotationCoeffs[iCoord])[i];

			rotationCoeffs_.push_back( 
				boost::make_tuple(
					currentTime,
					current, 
					differentialRotation ) );

			prevTime = currentTime;
		}

        rotationCoeffs_.back().get<1>() = endRot;
	}

	if( !other.positionCoeffs_.empty() )
	{
		if( positionCoeffs_.empty() )
		{
			this->positionCoeffs_ = other.positionCoeffs_;
		}
		else if( backwards )
        {
            // backwards sim, put other stuff after my stuff
            typedef std::vector<TimedCoeffList>::const_iterator Iter;
            Iter iter = 
		        std::upper_bound( other.positionCoeffs_.begin(), other.positionCoeffs_.end(), 
			        positionCoeffs_.back(), Compare1st<TimedCoeffList>() );

            if( iter != other.positionCoeffs_.end() )
            {
                assert( iter != other.positionCoeffs_.begin() );
                this->positionCoeffs_.reserve( positionCoeffs_.size() + 
                    std::distance( iter, other.positionCoeffs_.end() ) );

                // our last frame is kind of a placeholder, and we don't want to use its velocity
                // to extrapolate.  As a result, we will stick an extra frame on the end which 
                // does the extrapolation correctly using the original motion.
                Iter prev = iter - 1;
                int frame = positionCoeffs_.back().first + 1;
                vl::Vec3f eval = other.position( frame );
                boost::array<CoeffList, 3> coeffs = prev->second;
                for( size_t i = 0; i < 3; ++i )
                    (coeffs[i])[0] = eval[i];
                positionCoeffs_.push_back( TimedCoeffList( frame, coeffs ) );

                std::copy( iter, other.positionCoeffs_.end(), 
                    std::back_inserter( this->positionCoeffs_ ) );
            }
        }
        else
        {
            // regular sim, my stuff goes at the tail end of other stuff
		    std::vector<TimedCoeffList>::const_iterator iter = 
			    std::lower_bound( other.positionCoeffs_.begin(), other.positionCoeffs_.end(), 
				    positionCoeffs_.front(), Compare1st<TimedCoeffList>() );

		    std::vector<TimedCoeffList> temp( other.positionCoeffs_.begin(), iter );
		    std::copy( this->positionCoeffs_.begin(), this->positionCoeffs_.end(), 
			    std::back_inserter(temp) );
		    temp.swap( this->positionCoeffs_ );
        }
	}

	if( !other.rotationCoeffs_.empty() )
	{
		if( rotationCoeffs_.empty() )
		{
			this->rotationCoeffs_ = other.rotationCoeffs_;
		}
		else if( backwards )
        {
            // backwards sim, put other stuff after my stuff
            typedef std::vector<TimedDifferentialRotation>::const_iterator Iter;
	        Iter iter = 
		        std::upper_bound( other.rotationCoeffs_.begin(), other.rotationCoeffs_.end(), 
			        rotationCoeffs_.back(), CompareTimedRot() );
            if( iter != other.rotationCoeffs_.end() )
            {
                // if the other coeffs start well after mine end, then something is seriously wrong
                assert( iter != other.rotationCoeffs_.begin() );
                this->rotationCoeffs_.reserve( rotationCoeffs_.size() + 
                    std::distance( iter, other.rotationCoeffs_.end() ) + 1 );

                // our last frame is kind of a placeholder, and we don't want to use its velocity
                // to extrapolate.  As a result, we will stick an extra frame on the end which 
                // does the extrapolation correctly using the original motion.
                Iter prev = iter - 1;
                int frame = rotationCoeffs_.back().get<0>() + 1;
                Quaternion<double> eval = other.rotation( frame );
                rotationCoeffs_.push_back( TimedDifferentialRotation( frame, eval, prev->get<2>() ) );

                std::copy( iter, other.rotationCoeffs_.end(), 
                    std::back_inserter( this->rotationCoeffs_ ) );
            }
        }
        else
        {
            // regular sim, my stuff goes at the tail end of other stuff
		    std::vector<TimedDifferentialRotation>::const_iterator iter = 
			    std::lower_bound( other.rotationCoeffs_.begin(), other.rotationCoeffs_.end(), 
				    rotationCoeffs_.front(), CompareTimedRot() );

		    std::vector<TimedDifferentialRotation> temp( other.rotationCoeffs_.begin(), iter );
		    std::copy( this->rotationCoeffs_.begin(), this->rotationCoeffs_.end(), 
			    std::back_inserter(temp) );
		    temp.swap( this->rotationCoeffs_ );
        }
	}

	assert( !rotationCoeffs_.empty() );
	assert( !positionCoeffs_.empty() );
}

void PiecewisePath::dumpToArrays(
		std::vector<int>& positionTimes,
		std::vector<float>& positionCoeffs,
		std::vector<int>& rotationTimes,
		std::vector<float>& rotations,
		std::vector<float>& rotationCoeffs ) const
{
	positionTimes.resize( positionCoeffs_.size() );
	positionCoeffs.resize( 9*positionCoeffs_.size() );
	for( size_t i = 0; i < positionCoeffs_.size(); ++i )
	{
		positionTimes[i] = positionCoeffs_[i].first;
		for( size_t j = 0; j < 3; ++j )
			for( size_t k = 0; k < 3; ++k )
				positionCoeffs[9*i + 3*j + k] = (positionCoeffs_[i].second[j])[k];
	}

	rotationTimes.resize( rotationCoeffs_.size() );
	rotations.resize( 4*rotationCoeffs_.size() );
	rotationCoeffs.resize( 3*rotationCoeffs_.size() );
	for( size_t i = 0; i < rotationCoeffs_.size(); ++i )
	{
		rotationTimes[i] = rotationCoeffs_[i].get<0>();

		boost::array<double, 4> q = rotationCoeffs_[i].get<1>().get( 
			Quaternion<double>::ORDER_WXYZ );
		for( size_t j = 0; j < 4; ++j )
			rotations[4*i+j] = q[j];

		for( size_t j = 0; j < 3; ++j )
			rotationCoeffs[3*i + j] = (rotationCoeffs_[i].get<2>())[j];
	}
}

PiecewisePath::PiecewisePath(
		const std::vector<int>& positionTimes,
		const std::vector<float>& positionCoeffs,
		const std::vector<int>& rotationTimes,
		const std::vector<float>& rotations,
		const std::vector<float>& rotationCoeffs )
{
	this->positionCoeffs_.resize( positionTimes.size() );
	assert( positionCoeffs.size() == 9*positionTimes.size() );
	for( size_t i = 0; i < positionTimes.size(); ++i )
	{
		this->positionCoeffs_[i].first = positionTimes[i];
		for( size_t j = 0; j < 3; ++j )
			for( size_t k = 0; k < 3; ++k )
				(positionCoeffs_[i].second[j])[k] = positionCoeffs[9*i + 3*j + k];
	}

	this->rotationCoeffs_.resize( rotationTimes.size() );
	assert( rotations.size() == 4*rotationTimes.size() );
	assert( rotationCoeffs.size() == 3*rotationTimes.size() );
	for( size_t i = 0; i < rotationCoeffs_.size(); ++i )
	{
		this->rotationCoeffs_[i].get<0>() = rotationTimes[i];

		boost::array<float, 4> q;
		std::copy( rotations.begin() + 4*i, rotations.begin() + 4*(i+1), &q[0] );
		rotationCoeffs_[i].get<1>() = Quaternion<double>( &q[0],
			Quaternion<double>::ORDER_WXYZ );

		for( size_t j = 0; j < 3; ++j )
			(rotationCoeffs_[i].get<2>())[j] = rotationCoeffs[3*i + j];
	}
}

void dumpAttribute( 
	std::ostream& os,
	const char* attributeName, 
	const std::vector<double>& values,
	size_t nValues )
{
	assert( (values.size() % nValues) == 0 );
	os << "\tsetAttr -s " << nValues << " \"" << attributeName << 
		"[0:" << (nValues-1) << "]\"";
	size_t i = 0;
	for( std::vector<double>::const_iterator itr = values.begin();
		itr != values.end(); ++itr )
	{
		if( ((i+1) % 6) == 0 )
			os << "\n\t\t";
		else
			os << " ";
		os << *itr;
		++i;
	}
	os << ";\n";
}

void dumpAttribute( 
	std::ostream& os,
	const char* attributeName, 
	const char* value,
	size_t nValues )
{
	os << "\tsetAttr -s " << nValues << " \"" << attributeName << 
		"[0:" << (nValues-1) << "]\"";
	for( size_t i = 0; i < nValues; ++i )
	{
		if( ((i+1) % 6) == 0 )
			os << "\n\t\t";
		else
			os << " ";
		os << value;
	}
	os << ";\n";
}

double toDegrees( double radians )
{
	return radians * (180.0 / M_PI);
}

float PiecewisePath::angularEnergy() const
{
	float result = 0.0;
	for( std::vector<TimedDifferentialRotation>::const_iterator coeffItr = rotationCoeffs_.begin();
		coeffItr != rotationCoeffs_.end(); ++coeffItr )
	{
		vl::Vec3f velocity = coeffItr->get<2>();
		if( (coeffItr+1) != rotationCoeffs_.end() )
		{
			float timeDiff = boost::numeric_cast<float>( (coeffItr+1)->get<0>() - coeffItr->get<0>() );
			result += timeDiff * vl::dot( velocity, velocity );
		}
	}

	return result;
}

void PiecewisePath::toMayaAsciiFile( std::ostream& ofs, const std::string& name, size_t frameRate, PiecewisePath::FrameType startFrame, PiecewisePath::FrameType endFrame ) const
{
    startFrame = std::max( this->startFrame(), startFrame );
	endFrame = std::min( this->endFrame(), endFrame );
	assert( startFrame < endFrame );

    size_t frameSkip;
    if( frameRate >= 30 )
        frameSkip = frameRate / 30;
    else
        frameSkip = 1;

	for( size_t iTranslate = 0; iTranslate < 3; ++iTranslate )
	{
		std::string translateAttrName = name + "_translate";
		translateAttrName.push_back( 'X' + iTranslate );
		ofs << "createNode animCurveTL -n \"" << translateAttrName << "\";\n";

		std::vector<double> points;        points.reserve( (this->endFrame() - this->startFrame())*2 );

		for( FrameType frame = startFrame; frame < endFrame; frame += frameSkip )
		{
			double time = 30.0*boost::numeric_cast<double>(frame - startFrame) / 
				boost::numeric_cast<double>(frameRate);

			vl::Vec3f pos = this->position( frame );
			points.push_back( time );
			points.push_back( pos[iTranslate] );
		}

		dumpAttribute( ofs, ".ktv", points, points.size()/2 );
		ofs << "connectAttr \"" << translateAttrName << ".o\" \"" << name << ".t" << std::string(1, 'x' + iTranslate) << "\";\n";
	}
	/*
	for( size_t iTranslate = 0; iTranslate < 3; ++iTranslate )
	{
		std::string translateAttrName = name + "_translate";
		translateAttrName.push_back( 'X' + iTranslate );
		ofs << "createNode animCurveTL -n \"" << translateAttrName << "\";\n";

		std::vector<double> points;        points.reserve( positionCoeffs_.size() * 2 );
		std::vector<double> tangentInX;    tangentInX.reserve( positionCoeffs_.size() );
		std::vector<double> tangentInY;    tangentInY.reserve( positionCoeffs_.size() );
		std::vector<double> tangentOutX;   tangentOutX.reserve( positionCoeffs_.size() );
		std::vector<double> tangentOutY;   tangentOutY.reserve( positionCoeffs_.size() );

		tangentInY.push_back( 0.0 );
		tangentInX.push_back( 1.0 );
		for( std::vector<TimedCoeffList>::const_iterator coeffItr = positionCoeffs_.begin();
			coeffItr != positionCoeffs_.end(); ++coeffItr )
		{
			vl::Vec3f coeffs = coeffItr->second[iTranslate];
			double time = 30.0*boost::numeric_cast<double>(coeffItr->first) / 
				boost::numeric_cast<double>(frameRate);

			// position
			points.push_back( time );
			points.push_back( coeffs[0] );

			// now tangent
			if( (coeffItr+1) != positionCoeffs_.end() )
			{
				size_t nextTime = 30.0*boost::numeric_cast<double>((coeffItr+1)->first) / 
					boost::numeric_cast<double>(frameRate);

				double diff = nextTime - time;
				vl::Vec2d p0( boost::numeric_cast<double>(time), coeffs[0] );
				vl::Vec2d p1( boost::numeric_cast<double>(nextTime), 
					coeffs[0] + diff*coeffs[1] + diff*diff*coeffs[2] );
				vl::Vec2d c( time + 0.5*diff, coeffs[0] + 0.5*diff*coeffs[1] );

				vl::Vec2d c1 = (2.0/3.0)*c + (1.0/3.0)*p0;
				vl::Vec2d c2 = (2.0/3.0)*c + (1.0/3.0)*p1;

				vl::Vec2d outTangent = c1 - p0;
				vl::Vec2d inTangent = p0 - c2;
				tangentOutX.push_back( outTangent[0] );
				tangentOutY.push_back( outTangent[1] );
				tangentInX.push_back( inTangent[0] );
				tangentInY.push_back( inTangent[1] );
			}
		}
		tangentOutY.push_back( 0.0 );
		tangentOutX.push_back( 1.0 );

		dumpAttribute( ofs, ".ktv", points, positionCoeffs_.size() );

		dumpAttribute( ofs, ".ktl", "no", positionCoeffs_.size() );
		dumpAttribute( ofs, ".kwl", "no", positionCoeffs_.size() );

		dumpAttribute( ofs, ".kix", tangentInX, positionCoeffs_.size() );
		dumpAttribute( ofs, ".kiy", tangentInY, positionCoeffs_.size() );
		dumpAttribute( ofs, ".kox", tangentOutX, positionCoeffs_.size() );
		dumpAttribute( ofs, ".koy", tangentOutY, positionCoeffs_.size() );
		dumpAttribute( ofs, ".ktv", points, positionCoeffs_.size() );

		ofs << "connectAttr \"" << translateAttrName << ".o\" \"" << name << ".t" << std::string(1, 'x' + iTranslate) << "\";\n";
	}
	*/

	std::vector< std::vector<double> > coeffs( 3 );
	for( size_t i = 0; i < 3; ++i )
		coeffs.at(i).reserve( endFrame - startFrame );

	for( FrameType frame = startFrame; frame < endFrame; frame += frameSkip )
	{
		double time = 30.0*boost::numeric_cast<double>(frame - startFrame) / 
			boost::numeric_cast<double>(frameRate);

		vl::Vec3d ang = toMayaRotation( this->rotation( frame ) );
		for( size_t i = 0; i < 3; ++i )
		{
			if( !coeffs[i].empty() )
			{
				while( fabs((ang[i] - 360.0) - coeffs[i].back()) < fabs(ang[i] - coeffs[i].back()) )
					ang[i] -= 360.0;

				while( fabs((ang[i] + 360.0) - coeffs[i].back()) < fabs(ang[i] - coeffs[i].back()) )
					ang[i] += 360.0;
			}

			coeffs[i].push_back( time );
			coeffs[i].push_back( ang[i] );
		}
	}

	for( size_t i = 0; i < 3; ++i )
	{
		std::string rotateAttrName = name + "_Rotate";
		rotateAttrName.push_back( 'X' + i );
		ofs << "createNode animCurveTA -n \"" << rotateAttrName << "\";\n";
		dumpAttribute( ofs, ".ktv", coeffs[i], coeffs[i].size()/2 );

		ofs << "connectAttr \"" << rotateAttrName << ".o\" \"" << name << ".r" << std::string(1, 'x' + i) << "\";\n";
	}

	/*
	std::vector< std::vector<double> > coeffs( 4 );
	for( size_t i = 0; i < 4; ++i )
		coeffs.at(i).reserve( rotationCoeffs_.size() );

	for( std::vector<TimedDifferentialRotation>::const_iterator coeffItr = rotationCoeffs_.begin();
		coeffItr != rotationCoeffs_.end(); ++coeffItr )
	{
		double time = 30.0*boost::numeric_cast<double>(coeffItr->get<0>()) / 
			boost::numeric_cast<double>(frameRate);

		Quaternion<double> quat = coeffItr->get<1>();

		for( size_t i = 0; i < 4; ++i )
			coeffs[i].push_back( time );

		coeffs[0].push_back( quat.w() );
		coeffs[1].push_back( quat.x() );
		coeffs[2].push_back( quat.y() );
		coeffs[3].push_back( quat.z() );
	}

	for( size_t i = 0; i < 4; ++i )
	{
		std::string rotateAttrName = name + "_Rotate";
		rotateAttrName.push_back( 'W' + i );
		ofs << "createNode animCurveTU -n \"" << rotateAttrName << "\";\n";
		ofs << "\tsetAttr \".roti\" 3;\n";
		dumpAttribute( ofs, ".ktv", coeffs[i], rotationCoeffs_.size() );

		ofs << "connectAttr \"" << rotateAttrName << ".o\" \"" << name << ".rq" << std::string(1, 'w' + i) << "\";\n";
	}
	*/
}

void PiecewisePath::pathToMayaAsciiFile( std::ostream& ofs, const std::string& name, const vl::Vec3d& pointInLocalFrame ) const
{
	const double epsilon = 1e-3;
	std::deque<vl::Vec3d> points;
	for( FrameType iFrame = this->startFrame(); iFrame != this->endFrame(); ++iFrame )
	{
		vl::Vec3f pos = this->position( iFrame );
		Quaternion<double> q = this->rotation( iFrame );
		vl::Mat3d rot = q.toRotMatd();

		vl::Vec3d localPos = rot*pointInLocalFrame + toVec3d(pos);
		if( points.empty() || vl::len( localPos - points.back() ) > epsilon )
			points.push_back( localPos );
	}

	ofs << "createNode transform -n \"" << name << "\";\n";
	ofs << "createNode nurbsCurve -n \"" << name << "_shape\" -p \"" << name << "\";\n";
	ofs << "setAttr -k off \".v\";\n";
	ofs << "setAttr \".cc\" -type \"nurbsCurve\"\n";

	// First line: degree, number of spans, form (0=open, 1=closed, 2=periodic), rational (yes/no), dimension
	ofs << "\t1 " << (points.size()-1) << " 0 no 3\n";

	// Second line: number of knots, list of knot values
	ofs << "\t" << points.size();
	for( size_t i = 0; i < points.size(); ++i )
		ofs << " " << i;
	ofs << "\n";

	// Third line: number of CVs
	ofs << points.size();

	// Fourth and later lines: CV positions in x,y,z (and w if rational)
	for( size_t i = 0; i < points.size(); ++i )
	{
		const vl::Vec3d& localPos = points.at(i);
		ofs << "\t" << localPos[0] << " " << localPos[1] << " " << localPos[2] << "\n";
	}
	ofs << "\t;\n";
}

vl::Vec3f PiecewisePath::linearVelocity( const float time ) const
{
	TimedCoeffList tempValue;
	tempValue.first = boost::numeric_cast<FrameType>(time);
	std::vector<TimedCoeffList>::const_iterator iter = 
		std::upper_bound( positionCoeffs_.begin(), positionCoeffs_.end(), tempValue, Compare1st<TimedCoeffList>() );

    if( iter != positionCoeffs_.begin() )
		--iter;

	float diff = time - boost::numeric_cast<float>( iter->first );
    diff = std::max( diff, 0.0f );
	vl::Vec3f result;
	for( size_t iCoord = 0; iCoord < 3; ++iCoord )
	{
		const CoeffList& coeffs = iter->second[iCoord];
		result[iCoord] = coeffs[1] + 2.0*coeffs[2]*diff;
	}
	return result;
}

vl::Vec3f PiecewisePath::position( const float time ) const
{
	TimedCoeffList tempValue;
	tempValue.first = boost::numeric_cast<FrameType>(time);
	std::vector<TimedCoeffList>::const_iterator iter = 
		std::upper_bound( positionCoeffs_.begin(), positionCoeffs_.end(), tempValue, Compare1st<TimedCoeffList>() );

	// the first coeff must have time 0, so there shouldn't be anything before it:
	if( iter != positionCoeffs_.begin() )
	    --iter;

	float diff = time - boost::numeric_cast<float>( iter->first );
    diff = std::max( diff, 0.0f );
	vl::Vec3f result;
	for( size_t iCoord = 0; iCoord < 3; ++iCoord )
	{
		const CoeffList& coeffs = iter->second[iCoord];
		result[iCoord] = coeffs[0] + coeffs[1]*diff + coeffs[2]*diff*diff;
	}
	return result;
}

vl::Vec3f PiecewisePath::angularVelocity( const float time ) const
{
	TimedDifferentialRotation tempValue;
	tempValue.get<0>() = boost::numeric_cast<FrameType>(time);
	std::vector<TimedDifferentialRotation>::const_iterator iter = 
		std::upper_bound( rotationCoeffs_.begin(), rotationCoeffs_.end(), tempValue, CompareTimedRot() );

	if( iter != rotationCoeffs_.begin() )
		--iter;

	return iter->get<2>();
}

Quaternion<double> PiecewisePath::rotation( const float time ) const
{
	TimedDifferentialRotation tempValue;
	tempValue.get<0>() = boost::numeric_cast<FrameType>(time);
	std::vector<TimedDifferentialRotation>::const_iterator iter = 
		std::upper_bound( rotationCoeffs_.begin(), rotationCoeffs_.end(), tempValue, CompareTimedRot() );

	if( iter != rotationCoeffs_.begin() )
	    --iter;

	double diff = time - boost::numeric_cast<double>( iter->get<0>() );
    diff = std::max( diff, 0.0 );
	Quaternion<double> result = Quaternion<double>(diff*toVec3d(iter->get<2>())) * iter->get<1>();
	return result;
}


PiecewisePath::FrameType PiecewisePath::startFrame() const
{
	if( this->positionCoeffs_.empty() )
		return boost::numeric::bounds<size_t>::highest();

	return std::min(this->positionCoeffs_.front().first, 
		this->rotationCoeffs_.front().get<0>());
}

PiecewisePath::FrameType PiecewisePath::endFrame() const
{
	if( this->positionCoeffs_.empty() )
		return boost::numeric::bounds<size_t>::lowest();

	return std::max(this->positionCoeffs_.back().first,
		this->rotationCoeffs_.back().get<0>());
}

template <typename T>
void testEncodeDecode( const std::vector<T>& testValues )
{
	std::vector<char> encoded = encode( testValues );
	std::vector<int> decoded = decode( encoded, 0, std::numeric_limits<T>::is_signed );

/*
	if( !testValues.empty() )
	{
		MATFile values( "//volley/cdtwigg/values.mat", "blah" );
		values.add( "A", std::vector<double>(testValues.begin(), testValues.end()) );
		values.add( "B", std::vector<double>(decoded.begin(), decoded.end()) );
	}
*/


	assert( decoded.size() == testValues.size() );
	for( size_t i = 0; i < decoded.size(); ++i )
		assert( decoded[i] == testValues[i] );
}

std::vector<CompressedPath> compress( const std::deque<Simulation::State>& states, float quadConstant, size_t frameRate, const float maxPositionError, const float maxRotationError )
{
	std::vector<CompressedPath> result;
	result.reserve( states.front().stateCount() );

	size_t objectCount = states.front().stateCount();
	PercentageUpdater updater( objectCount, "Compressing" );
	assert( !states.empty() );
	for( size_t iObject = 0; iObject < objectCount; ++iObject )
	{
		updater.update( iObject );

		std::deque<RigidDynamicState> objectStates;
		for( std::deque<Simulation::State>::const_iterator itr = states.begin();
			itr != states.end(); ++itr )
		{
			objectStates.push_back( itr->state(iObject) );
		}

		result.push_back( compress(objectStates, quadConstant, 0, frameRate, maxPositionError, maxRotationError ) );
	}

	return result;
}

//#define CHECK_GAPS

CompressedPath compress( 
	const std::deque<RigidDynamicState>& states, 
	float quadConstant,
	size_t startOffset,
    size_t frameRate,
	const float maxPositionError,
	const float maxRotationError )
{
#ifdef INSTRUMENT_COMPRESSION
	static WaypointTracker::CumulativeTimes cumulativeTimes;
	WaypointTracker tracker( cumulativeTimes );
#endif

	assert( !states.empty() );

	// @todo factor compression into a separate function
	// we will encode each object separately
	const float epsilon = maxPositionError;
	const float coeffsEpsilon = epsilon / 50.0;

	const double rotEpsilon = maxRotationError;
	const double rotCoeffEpsilon = rotEpsilon / 50.0;

	// @todo remove this requirement:
	assert( states.size() >= 2 );

#ifdef CHECK_GAPS
	std::deque<size_t> gaps;
#endif

	CompressedPath path;
    path.startFrame = startOffset;

	{
		std::vector<size_t> startTimes;
		boost::array< std::vector<int>, 3 > linearCoeffs;
		boost::array< std::vector<unsigned char>, 3 > quadCoeffs;

		std::vector< vl::Vec3f > points;
		points.reserve( states.size() );

		for( size_t iState = 0; iState < states.size(); ++iState )
			points.push_back( toVec3f(states[iState].position()) );

		size_t startFrame = 0;
		vl::Vec3f predicted = points.front();

		while( startFrame < (states.size()-1) )
		{
			{
				size_t endFrame = startFrame + 1;
				assert( endFrame < states.size() );

				float error;
				TinyVec<int, 3> linearCoeff;
				TinyVec<unsigned char, 3> quadraticCoeff;

				while( endFrame < states.size() )
				{
					TinyVec<int, 3> linearCoeffNew = linearCoeff;
					TinyVec<unsigned char, 3> quadraticCoeffNew = quadraticCoeff;

					// we will fit a spline to the first n=5 points, and after that
					//   we will simply extrapolate; this will significantly reduce
					//   the runtime cost
					if( endFrame - startFrame < 10 )
					{
						boost::tie( error, linearCoeffNew, quadraticCoeffNew ) =
							fitSpline( points, startFrame, endFrame, predicted, coeffsEpsilon, quadConstant
#ifdef INSTRUMENT_COMPRESSION
									, tracker
#endif
									);
						if( error > epsilon )
							break;
					}
					else
					{
						vl::Vec3f origPoint = points.at( endFrame );

						// extrapolate value using previously computed coeffs
						size_t iFrame = endFrame - startFrame;
						vl::Vec3f sum = predicted;

						for( size_t iCoord = 0; iCoord < 3; ++iCoord )
						{
							float reconstructedCoeff = 
								boost::numeric_cast<float>( linearCoeffNew[iCoord] ) * coeffsEpsilon; 
							sum[iCoord] += reconstructedCoeff * iFrame;
							sum[iCoord] += boost::numeric_cast<float>( quadraticCoeffNew[iCoord] )
								* quadConstant * boost::numeric_cast<float>(iFrame*iFrame);
						}

						if( vl::len(sum - origPoint) > epsilon )
						{
							break;
						}
					}

					++endFrame;
					linearCoeff = linearCoeffNew;
					quadraticCoeff = quadraticCoeffNew;
				}

//				std::cout << "endFrame: " << endFrame << ", startFrame: " << startFrame << ", states.size(): " << states.size() << std::endl;

				// endFrame needs to refer to an actual frame:
				assert( endFrame > 0 );
				--endFrame;
				assert( endFrame > startFrame );

				for( size_t iCoord = 0; iCoord < 3; ++iCoord )
				{
					linearCoeffs[iCoord].push_back( linearCoeff[iCoord] );
					quadCoeffs[iCoord].push_back( quadraticCoeff[iCoord] );
				}

#ifdef CHECK_GAPS
				gaps.push_back( endFrame - startFrame );
#endif

				assert( startTimes.empty() || startFrame > startTimes.back() );
				startTimes.push_back( startFrame );
	
				size_t iFrame = endFrame - startFrame;
				vl::Vec3f sum = predicted;

				for( size_t iCoord = 0; iCoord < 3; ++iCoord )
				{
					float reconstructedCoeff = 
						boost::numeric_cast<float>( linearCoeffs[iCoord].back() ) * coeffsEpsilon; 
					sum[iCoord] += reconstructedCoeff * iFrame;
					sum[iCoord] += boost::numeric_cast<float>(quadCoeffs[iCoord].back())
						* quadConstant * boost::numeric_cast<float>(iFrame*iFrame);
				}

				predicted = sum;

				vl::Vec3f actualPos = points.at(startFrame + iFrame);
				float err = vl::len( predicted - actualPos );
				assert( fabs(err) <= sqrt(3.0)*epsilon );

				startFrame += iFrame;
			}
		}

        // It is fairly important that we get the ending velocity correct
		startTimes.push_back( states.size() );
        vl::Vec3f endingVelocity = states.back().linearVelocity();
		for( size_t iCoord = 0; iCoord < 3; ++iCoord )
		{
            float v = endingVelocity[iCoord] / static_cast<float>(frameRate);

            // want to round toward 0:
            int rounded = boost::numeric_cast<int>( v < 0 ? v + 1.0 : v );
			linearCoeffs[iCoord].push_back( rounded );
			quadCoeffs[iCoord].push_back( 0 );
		}

        /*
		// need to have a marker for the stop:
		startTimes.push_back( states.size() );
		for( size_t iCoord = 0; iCoord < 3; ++iCoord )
		{
			linearCoeffs[iCoord].push_back( 0 );
			quadCoeffs[iCoord].push_back( 0 );
		}
        */

		{
#ifdef INSTRUMENT_COMPRESSION
			boost::shared_ptr<WaypointTracker::Span> span(
				tracker.newSpan("Encoding") );
#endif

			for( size_t iCoord = 0; iCoord < 3; ++iCoord )
			{
				std::vector<int> differences( linearCoeffs[iCoord].size() );
				std::adjacent_difference( linearCoeffs[iCoord].begin(), linearCoeffs[iCoord].end(), 
					differences.begin() );
				assert( linearCoeffs[iCoord].size() == startTimes.size() );

				path.positionLinearCoeffs[iCoord] =
					encode( differences, (points.front())[iCoord], (points.back())[iCoord], coeffsEpsilon );
#ifdef _DEBUG
				std::vector<int> decoded = decode( path.positionLinearCoeffs[iCoord], 3*sizeof(float), true );
				assert( decoded.size() == differences.size() );
				for( size_t i = 0; i < differences.size(); ++i )
					assert( decoded[i] == differences[i] );
#endif
			}

			{
				path.positionQuadraticCoeffs[ 1 ] =
					encode( quadCoeffs[1], 0.0f, 0.0f, quadConstant );
#ifdef _DEBUG
				std::vector<int> decoded = decode( path.positionQuadraticCoeffs[ 1 ], 3*sizeof(float), false );
				assert( decoded.size() == quadCoeffs[1].size() );
				for( size_t i = 0; i < quadCoeffs[1].size(); ++i )
					assert( decoded[i] == (quadCoeffs[1])[i] );
#endif
			}

			{
				std::vector<size_t> differences( startTimes.size() );
				std::adjacent_difference( startTimes.begin(), startTimes.end(), 
					differences.begin() );

				path.positionTimes = encode( differences );
#ifdef _DEBUG
				std::vector<int> decoded = decode( path.positionTimes, 0, false );
				assert( decoded.size() == differences.size() );
				for( size_t i = 0; i < differences.size(); ++i )
					assert( decoded[i] == differences[i] );
#endif
			}

			TinyVec<size_t, 3> linearCoeffsSize;
			for( size_t iCoord = 0; iCoord < 3; ++iCoord )
				linearCoeffsSize[iCoord] = path.positionLinearCoeffs[ iCoord ].size();
			size_t quadraticCoeffsSize = path.positionQuadraticCoeffs[1].size();
			size_t timesSize = path.positionTimes.size();

			size_t bytes = linearCoeffsSize[0] + linearCoeffsSize[1] + linearCoeffsSize[2]
				+ quadraticCoeffsSize + timesSize;
			size_t original = 3*4*states.size();

			double ratio = boost::numeric_cast<double>( original ) /
				boost::numeric_cast<double>( bytes );
		}
	}

	{
#ifdef INSTRUMENT_COMPRESSION
		boost::shared_ptr<WaypointTracker::Span> span(
			tracker.newSpan("Rotation") );
#endif

		// now for angular component
		std::deque< TinyVec<int, 3> > coeffs;
		std::deque<size_t> startTimes;

		std::vector< Quaternion<double> > rotations;
		rotations.reserve( states.size() );
		for( size_t iState = 0; iState < states.size(); ++iState )
			rotations.push_back( 
				Quaternion<double>(states[iState].orientation()) );

		/*
#ifdef _DEBUG
		boost::shared_ptr<MATFile> test( 
			new MATFile("//volley/cdtwigg/compTest.mat", "test compression") );
		{
			fortran_matrix m( rotations.size(), 4 );
			for( size_t i = 0; i < rotations.size(); ++i )
			{
				Quaternion<double> q = rotations[i];
				m(i, 0) = q.x();
				m(i, 1) = q.y();
				m(i, 2) = q.z();
				m(i, 3) = q.w();
			}

			test->add( "R", m );
		}
#endif
		*/


		size_t startFrame = 0;
		vl::Vec3f initialRot;
        vl::Vec3f finalRot;

		{
			vl::Vec3d axis;
			double angle;
			rotations.front().toAxisAngle(axis, angle);
			initialRot = toVec3f( angle*axis );
			assert( !isBad(initialRot) );

            rotations.back().toAxisAngle(axis, angle);
			finalRot = toVec3f( angle*axis );
			assert( !isBad(finalRot) );
		}

		Quaternion<double> predicted( toVec3d(initialRot) );

		while( startFrame < (states.size()-1) )
		{
			size_t high = startFrame + 1;
			size_t low = startFrame;

			while( high < states.size() )
			{
				float error = 
					findQuaternionInterp( rotations, startFrame, high, predicted, rotCoeffEpsilon ).first;

				if( error > rotEpsilon )
					break;

				// double it
				low = high;
				high += (high - startFrame);
			}

			// binary search implementation borrowed from 
			//   http://www.tbray.org/ongoing/When/200x/2003/03/22/Binary
			high = std::min( high, states.size() );

			while( high - low > 1 )
			{
				size_t probe = (high + low) / 2;

				float error = 
					findQuaternionInterp( rotations, startFrame, probe, predicted, rotCoeffEpsilon ).first;

				if( error > epsilon )
					high = probe;
				else
					low = probe;
			}

			size_t endFrame = std::max(low, startFrame+1);

			float error;
			TinyVec<int, 3> currentCoeff;
			boost::tie( error, currentCoeff ) =
				findQuaternionInterp( rotations, startFrame, endFrame, predicted, rotCoeffEpsilon );

			coeffs.push_back( currentCoeff );
			startTimes.push_back( startFrame );

			size_t iFrame = endFrame - startFrame;

			vl::Vec3d reconstructedDiff;
			for( vl::Int i = 0; i < 3; ++i )
				reconstructedDiff[i] = boost::numeric_cast<double>(currentCoeff[i]) * rotCoeffEpsilon;

			vl::Vec3d currentDiffRot = 
					reconstructedDiff*boost::numeric_cast<double>( iFrame );
			predicted = Quaternion<double>(currentDiffRot) * predicted;

			Quaternion<double> errorQuat = inverse(predicted) * rotations[endFrame];
			assert( fabs(arg(errorQuat)) <= rotEpsilon );

			startFrame += iFrame;
		}

        vl::Vec3f endingVelocity = states.back().angularVelocity();
        TinyVec<int, 3> roundedVelocity;
		for( size_t iCoord = 0; iCoord < 3; ++iCoord )
		{
            float v = endingVelocity[iCoord] / static_cast<float>(frameRate);
            int rounded = boost::numeric_cast<int>( v < 0 ? v - 0.5 : v + 0.5 );
            roundedVelocity[ iCoord ] = rounded;
		}

		startTimes.push_back( states.size() );
		coeffs.push_back( roundedVelocity );

		/*
#ifdef _DEBUG
		{
			fortran_matrix tmp( coeffs.size(), 3 );
			for( size_t i = 0; i < coeffs.size(); ++i )
				for( size_t j = 0; j < 3; ++j )
					tmp(i, j) = (coeffs[i])[j];

			test->add( "coeffs", tmp );
			test->add( "startTimes", std::vector<int>(startTimes.begin(), startTimes.end()) );
		}
#endif
		*/

		for( size_t iCoord = 0; iCoord < 3; ++iCoord )
		{
			std::vector<int> currentCoeffs;
			currentCoeffs.reserve( coeffs.size() );
			for( std::deque< TinyVec<int, 3> >::const_iterator iter = coeffs.begin();
				iter != coeffs.end(); ++iter )
			{
				currentCoeffs.push_back( (*iter)[iCoord] );
			}

			std::vector<int> differences( currentCoeffs.size() );
			std::adjacent_difference( currentCoeffs.begin(), currentCoeffs.end(), 
				differences.begin() );

			path.rotationCoeffs[iCoord] =
				encode( differences, initialRot[iCoord], finalRot[iCoord], rotCoeffEpsilon );
		}

		{
			std::vector<size_t> differences( startTimes.size() );
			std::adjacent_difference( startTimes.begin(), startTimes.end(), 
				differences.begin() );

			path.rotationTimes = encode( differences );
		}
	}

	TinyVec<size_t, 3> coeffsSize;
	for( size_t iCoord = 0; iCoord < 3; ++iCoord )
		coeffsSize[iCoord] = path.rotationCoeffs[ iCoord ].size();
	size_t timesSize = path.rotationTimes.size();

	size_t bytes = coeffsSize[0] + coeffsSize[1] + coeffsSize[2]
		+ timesSize;
	size_t original = 3*4*states.size();

	double ratio = boost::numeric_cast<double>( original ) /
		boost::numeric_cast<double>( bytes );

#ifdef _DEBUG
	{
		PiecewisePath reconstructedPath( path );
	}
#endif

	/*
#ifdef _DEBUG
	PiecewisePath reconstructedPath( path );

	for( size_t iState = 0; iState < states.size(); ++iState )
	{
		vl::Vec3f position = reconstructedPath.position( iState+startOffset );
		Quaternion<double> rotation = reconstructedPath.rotation( iState+startOffset );

		vl::Vec3f originalPos = toVec3f(states[iState].position());
		Quaternion<double> originalRot( states[iState].orientation() );

		vl::Vec3f errorPos = position - originalPos;
		assert( vl::len(errorPos) <= sqrt(3.0)*epsilon );

		Quaternion<double> errorRot = inverse( originalRot ) * rotation;
		assert( arg(errorRot) <= rotEpsilon );
	}
#endif
	*/

#ifdef INSTRUMENT_COMPRESSION
	tracker.dump(std::cout);
#endif

#ifdef CHECK_GAPS
	std::deque<size_t>::iterator mid = gaps.begin() + (gaps.size()/2);
	std::nth_element( gaps.begin(), mid, gaps.end() );
	size_t medianGap = *mid;
	std::cerr << "Median gap: " << medianGap << std::endl;
#endif

	return path;

	/*
	for( size_t iCoord = 0; iCoord < 3; ++iCoord )
	{
		size_t order = 2;
		if( iCoord == 1 )
			order = 3;

		--order;

		// for time
		++order;
		size_t numRows = intRepresentation[iCoord].size() / order;
		fortran_matrix A( numRows, order );
		for( size_t i = 0; i < numRows; ++i )
			for( size_t j = 0; j < order; ++j )
				A(i, j) = (intRepresentation[iCoord])[ i*order + j ];

		fortran_matrix orig( states.size(), this->dynamicObjects_.size() );
		for( size_t i = 0; i <  states.size(); ++i )
			for( size_t j = 0; j < orig.ncols(); ++j )
				orig(i, j) = (states.at(i).state(j).position())[ iCoord ];

		char which = 'x' + iCoord;
		std::string matFilename( "//volley/cdtwigg/coeffs_" );
		matFilename.push_back( which );
		matFilename.append( ".mat" );
		MATFile matFile( matFilename, "coord" );
		matFile.add( "A", A );
		matFile.add( "orig", orig );

		size_t origFloats = this->dynamicObjects_.size() * states.size();
		size_t newFloats = A.nrows() * (order - 1);
		double ratio = boost::numeric_cast<double>(origFloats) / boost::numeric_cast<double>(newFloats);
		std::cerr << "Coord " << iCoord << ": ratio: " << ratio << std::endl;
	}
	*/

	/*
	typedef boost::numeric::converter<
		int
		, double
		, boost::numeric::conversion_traits<int, float>
		, boost::numeric::silent_overflow_handler
		, boost::numeric::Trunc<float> > Rounder;

	std::deque<vl::Vec3f> differenceSamples;

	// @todo factor compression into a separate function
	// we will encode each object separately
	const float epsilon = 1e-4;
	for( size_t iObject = 0; iObject < this->dynamicObjects_.size(); ++iObject )
	{
		std::deque<vl::Vec3f> currentFrontier( 3, 
			toPosition(states.front().extendedState(iObject)) );

		for( size_t iFrame = 0; iFrame < states.size(); ++iFrame )
		{
			const char* newState = states[iFrame].extendedState( iObject );

			// we use a second-order predictor here
			vl::Vec3f predictedPos = 
				3*currentFrontier[0]
					- 3*currentFrontier[1]
					+ currentFrontier[2];
			vl::Vec3f actualPos = 
				toPosition(newState);

			vl::Vec3f difference = actualPos - predictedPos;
			vl::Vec3f storedVec;
			for( vl::Int j = 0; j < 3; ++j )
			{
				float dividedDiff = difference[j] / epsilon;
				int diffInt = Rounder::convert( dividedDiff );
				difference[j] = boost::numeric_cast<float>(diffInt) * epsilon;
				storedVec[j] = predictedPos[j] + difference[j];
			}
			differenceSamples.push_back( difference );

			currentFrontier.push_front( storedVec );
			currentFrontier.pop_back();
		}
	}
	*/

	/*
	fortran_matrix mat( differenceSamples.size(), 3 );
	for( size_t i = 0; i < differenceSamples.size(); ++i )
		for( vl::Int j = 0; j < 3; ++j )
			mat( i, j ) = (differenceSamples[i])[j];

	wtf, mate?

	try
	{
		boost::mutex::scoped_lock lock( hdf5Mutex );
		HDF5File outFile( "samples.hdf5" );
		outFile.add( "S", mat );
	}
	catch( Exception& e )
	{
		std::string msg = e.message();
		std::cerr << msg << std::endl;
	}
	*/

	/*
	std::ofstream ofs( "samples.dat", std::ios::app );
	for( std::deque<vl::Vec3f>::const_iterator iter = differenceSamples.begin();
		iter != differenceSamples.end(); ++iter )
	{
		for( vl::Int i = 0; i < 3; ++i )
			ofs << (*iter)[i] << " ";

		ofs << "\n";
	}
	*/
}

void testCompression()
{
	RandomStream random;

	for( size_t i = 0; i < 3; ++i )
	{
		testEncodeDecode( std::vector<int>(i, 0) );
		testEncodeDecode( std::vector<int>(i, 1) );
		testEncodeDecode( std::vector<size_t>(i, 0) );
		testEncodeDecode( std::vector<size_t>(i, 1) );

		std::vector<size_t> test(i, 0);
		test.push_back( 1 );
		testEncodeDecode(test);
	}

	for( size_t iTest = 0; iTest < 10; ++iTest )
	{
		std::vector<int> testValues;
		size_t num = random.uniform( 10, 100 );
		for( size_t i = 0; i < num; ++i )
			testValues.push_back( boost::numeric_cast<int>(random.normal( 0.0, 60.0 )) );

		testEncodeDecode( testValues );
	}

	for( size_t iTest = 0; iTest < 10; ++iTest )
	{
		std::vector<size_t> testValues;
		size_t num = random.uniform( 10, 100 );
		for( size_t i = 0; i < num; ++i )
			testValues.push_back( random.uniform( 0, 50 ) );

		testEncodeDecode( testValues );
	}

#ifdef NAG
	// test least infinity norm stuff
	{
		size_t numElts = 5;
		size_t order = 3;
		fortran_matrix A( numElts, order );
		for( size_t i = 0; i < numElts; ++i )
		{
			for( size_t j = 0; j < order; ++j )
			{
				A(i, j) = std::pow(
					boost::numeric_cast<double>(i), 
					boost::numeric_cast<double>(j) );
			}
		}

		std::vector<double> b( 5 );
		for( size_t i = 0; i < b.size(); ++i )
			b.at(i) = random.normal();

		leastInfinityNorm( A, b );
	}
#endif
}

void computeRateDistortionCurve( const std::deque<Simulation::State>& states,
	const float quadConstant,
    size_t frameRate,
	const std::string& matFilename )
{
	assert( !states.empty() );

	std::deque< std::vector<double> > samples;

	for( double k = -1.0; k > -4.0; k -= 0.01 )
	{
		// compress
		double maxErr = pow( 10.0, k );
		std::vector<CompressedPath> compressed = compress( states, quadConstant, frameRate, maxErr, maxErr );
		assert( compressed.size() == 1 );

		// decompress
		std::vector<PiecewisePath> piecewisePaths( compressed.begin(), compressed.end() );

		// compute error
		std::vector<double> l2Error( piecewisePaths.size(), 0.0 );
		std::vector<double> lInfError( piecewisePaths.size(), 0.0 );
		for( size_t iState = 0; iState < states.size(); ++iState )
		{
			for( size_t iObject = 0; iObject < piecewisePaths.size(); ++iObject )
			{
				vl::Vec3d compressedPos = toVec3d( piecewisePaths.at(iObject).position( iState ) );
				vl::Vec3d originalPos = toVec3d( states.at(iState).state(iObject).position() );
				vl::Vec3d diff = compressedPos - originalPos;
				l2Error.at(iObject) += dot( diff, diff );
				lInfError.at(iObject) = std::max( vl::len(diff), lInfError.at(iObject) );
			}
		}

		size_t numBits = 0;
		for( size_t iObject = 0; iObject < compressed.size(); ++iObject )
		{
			numBits += 8*compressed.at(iObject).positionTimes.size();
			for( size_t j = 0; j < 3; ++j )
				numBits += 8*compressed.at(iObject).positionLinearCoeffs[j].size() + 32;

			for( size_t j = 0; j < 3; ++j )
				numBits += 8*compressed.at(iObject).positionQuadraticCoeffs[j].size() + 32;
		}

		// since we need to keep track of sizes of the arrays to decompress
		numBits += 5*4;

		double seconds = states.back().time() - states.front().time();

		samples.push_back( std::vector<double>(5) );
		samples.back().at(0) = k;
		samples.back().at(1) = maxErr;
		samples.back().at(2) = std::accumulate( l2Error.begin(), l2Error.end(), 0.0 );
		samples.back().at(3) = norminf( lInfError );
		samples.back().at(4) = boost::numeric_cast<double>(numBits) / 
			boost::numeric_cast<double>(states.size());
	}

	fortran_matrix mat( samples.size(), 5 );
	for( size_t iSample = 0; iSample < samples.size(); ++iSample )
		for( size_t j = 0; j < samples.at(iSample).size(); ++j )
			mat(iSample, j) = samples.at(iSample).at(j);

	MATFile matFile( matFilename, "rate/distortion curves" );
	matFile.add( "rate", mat );
}

} // namespace planning

