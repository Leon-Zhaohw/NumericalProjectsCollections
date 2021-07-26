#ifndef __UTIL_H__

#define __UTIL_H__

// Global utility functions


#include <gsl/gsl_math.h>

#include <vl/VL.h>
#include <vl/VLf.h>
#include <vl/VLd.h>

#include <boost/array.hpp>
#include <boost/static_assert.hpp>

#include <algorithm>

template <typename T>
T signum(T value)
{
	if (value >= 0.0) return static_cast<T>(1.0); else return static_cast<T>(-1.0);
}

template <typename T>
bool isBad( const T& );

template <>
bool isBad(const double& number);

template <>
bool isBad(const float& number);

template <>
bool isBad(const vl::Vec4d& vec);

template <>
bool isBad(const vl::Vec4f& vec);

template <>
bool isBad(const vl::Vec3d& vec);

template <>
bool isBad(const vl::Vec3f& vec);

template <>
bool isBad(const vl::Vec2d& vec);

template <>
bool isBad(const vl::Vec2f& vec);

template <>
bool isBad(const vl::Mat3d& mat);

template <>
bool isBad(const vl::Mat3f& mat);

template <>
bool isBad(const vl::Mat4d& mat);

template <>
bool isBad(const vl::Mat4f& mat);

template <typename const_iterator>
bool isBad(const_iterator start, const_iterator end)
{
	return std::find_if( start, end, &isBad ) != end;
}

#include <vector>
template <typename T>
bool isBad(const std::vector<T>& vec)
{
	typedef bool (*BadFcn)(const T& t);
	BadFcn fcn = &isBad<T>;
	return std::find_if( vec.begin(), vec.end(), fcn ) != vec.end();
}

template <int v>
struct Int2Type
{
	enum { value = v };
};

#include <string>
class Named
{
public:
	explicit Named( const std::string& name )
		: name_( name )							{}
	std::string name() const					{ return name_; }
	void setName( const std::string& name )		{ name_ = name; }
	virtual ~Named()							{}

private:
	std::string name_;
};

template <typename NamedClass>
struct NamedComparator
{
	bool operator()( const NamedClass& first, 
		const NamedClass& second )
	{
		return std::less<std::string>()( first.name(), second.name() );
	}
};

struct NamedPtrComparator
{
	template <typename FirstPtr, typename SecondPtr>
	bool operator()( const FirstPtr first, 
		const SecondPtr second )
	{
		return std::less<std::string>()( first->name(), second->name() );
	}
};

template <typename T, unsigned int N>
class TinyVec
{
public:
	explicit TinyVec()
	{
		for( unsigned int i = 0; i < N; ++i )
			vec_[i] = 0;
	}

	explicit TinyVec( T i )
	{
		BOOST_STATIC_ASSERT( N == 1 );
		vec_[0] = i;
	}


	explicit TinyVec( T i, T j )
	{
		BOOST_STATIC_ASSERT( N == 2 );
		vec_[0] = i;
		vec_[1] = j;
	}

	explicit TinyVec( T i, T j, T k )
	{
		BOOST_STATIC_ASSERT( N == 3 );
		vec_[0] = i;
		vec_[1] = j;
		vec_[2] = k;
	}

	explicit TinyVec( T i, T j, T k, T l )
	{
		BOOST_STATIC_ASSERT( N == 4 );
		vec_[0] = i;
		vec_[1] = j;
		vec_[2] = k;
		vec_[3] = l;
	}

	T operator[]( unsigned int index ) const
	{
		assert( index < N );
		return vec_[index];
	}

	T& operator[]( unsigned int index )
	{
		assert( index < N );
		return vec_[index];
	}

	inline TinyVec<T, N>& operator*=( const T& rhs )
	{
		for( unsigned int i = 0; i < N; ++i )
			vec_[i] *= rhs;
		return *this;
	}

	inline TinyVec<T, N>& operator/=( const T& rhs )
	{
		for( unsigned int i = 0; i < N; ++i )
			vec_[i] /= rhs;
		return *this;
	}

	inline TinyVec<T, N>& operator+=( const TinyVec<T, N>& rhs )
	{
		for( unsigned int i = 0; i < N; ++i )
			vec_[i] += rhs[i];
		return *this;
	}

	inline TinyVec<T, N>& operator-=( const TinyVec<T, N>& rhs )
	{
		for( unsigned int i = 0; i < N; ++i )
			vec_[i] -= rhs[i];
		return *this;
	}

	const boost::array<T, N>& asArray() const
	{
		return vec_;
	}

private:
	boost::array<T, N> vec_;
};

template <typename T, unsigned int N>
TinyVec<T, N> operator+( const TinyVec<T, N>& lhs, const TinyVec<T, N>& rhs )
{
	TinyVec<T, N> result(lhs);
	result += rhs;
	return result;
}

template <typename T, unsigned int N>
TinyVec<T, N> operator-( const TinyVec<T, N>& lhs, const TinyVec<T, N>& rhs )
{
	TinyVec<T, N> result(lhs);
	result -= rhs;
	return result;
}

template <typename T, unsigned int N>
TinyVec<T, N> operator/( const TinyVec<T, N>& lhs, const T& rhs )
{
	TinyVec<T, N> result( lhs );
	result /= rhs;
	return result;
}

template <typename T, unsigned int N>
TinyVec<T, N> operator*( const TinyVec<T, N>& lhs, const T& rhs )
{
	TinyVec<T, N> result( lhs );
	result *= rhs;
	return result;
}

template <typename _Pair>
class IgnoreSecondComparator
{
public:
	bool operator()( const _Pair& left, const _Pair& right )
	{
		return std::less<typename _Pair::first_type>()( left.first, right.first );
	}
};


#endif
