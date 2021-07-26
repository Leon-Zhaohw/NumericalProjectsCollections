#ifndef __MECHMODEL_H__
#define __MECHMODEL_H__

#include "novodexWrappers.h"
#include "physicsFwd.h"

#include <deque>
#include <iosfwd>
#include <sstream>
#include <string>
#include <iomanip>

namespace planning {

class MechModelObject
{
public:
	virtual void dump( std::ostream& os, const std::string& prefix ) const = 0;
};

template <typename T>
class StringForMechModelValue
{
public:
	std::string operator()( const T& value )
	{
		std::ostringstream oss;
		oss << value;
		return oss.str();
	}
};

template <>
class StringForMechModelValue<std::string>
{
public:
	std::string operator()( const std::string& value )
	{
		std::string result;
		result.append( "\"" );
		result.append( value );
		result.append( "\"" );
		return result;
	}
};

template <>
class StringForMechModelValue< std::pair<std::string, std::string> >
{
public:
	std::string operator()( const std::pair<std::string, std::string>& value )
	{
		StringForMechModelValue<std::string> conv;
		std::string result;
		result.append( conv( value.first ) );
		result.append( "." );
		result.append( conv(value.second) );
		return result;
	}
};

template <>
class StringForMechModelValue<double>
{
public:
	std::string operator()( const double& value )
	{
		std::ostringstream oss;
		oss << std::setprecision(10) << value;
		return oss.str();
	}
};

template <>
class StringForMechModelValue<unsigned int>
{
public:
	std::string operator()( const unsigned int& value )
	{
		std::ostringstream oss;
		oss << value;
		return oss.str();
	}
};

template <>
class StringForMechModelValue<int>
{
public:
	std::string operator()( const int& value )
	{
		std::ostringstream oss;
		oss << value;
		return oss.str();
	}
};

template <>
class StringForMechModelValue<vl::Vec2d>
{
public:
	std::string operator()( const vl::Vec2d& value )
	{
		std::ostringstream oss;
		StringForMechModelValue<double> conv;
		oss <<  "(" << conv( value[0] ) 
			<< ", " << conv( value[1] ) << ")";
		return oss.str();
	}
};

template <>
class StringForMechModelValue<vl::Vec3d>
{
public:
	std::string operator()( const vl::Vec3d& value )
	{
		std::ostringstream oss;
		StringForMechModelValue<double> conv;
		oss <<  "(" << conv( value[0] ) 
			<< ", " << conv( value[1] ) 
			<< ", " << conv( value[2] ) << ")";
		return oss.str();
	}
};

template <>
class StringForMechModelValue<vl::Vec4d>
{
public:
	std::string operator()( const vl::Vec4d& value )
	{
		std::ostringstream oss;
		StringForMechModelValue<double> conv;
		oss <<  "(" << conv( value[0] ) 
			<< ", " << conv( value[1] ) 
			<< ", " << conv( value[2] ) 
			<< ", " << conv( value[3] ) << ")";
		return oss.str();
	}
};

template <>
class StringForMechModelValue<vl::Vec2f>
{
public:
	std::string operator()( const vl::Vec2f& value )
	{
		StringForMechModelValue<vl::Vec2d> conv;
		return conv( toVec2d(value) );
	}
};

template <>
class StringForMechModelValue<vl::Vec3f>
{
public:
	std::string operator()( const vl::Vec3f& value )
	{
		StringForMechModelValue<vl::Vec3d> conv;
		return conv( toVec3d(value) );
	}
};

template <>
class StringForMechModelValue<vl::Vec4f>
{
public:
	std::string operator()( const vl::Vec4f& value )
	{
		StringForMechModelValue<vl::Vec4d> conv;
		return conv( toVec4d(value) );
	}
};

template <>
class StringForMechModelValue< std::vector<vl::Vec3f> >
{
public:
	std::string operator()( const std::vector<vl::Vec3f>& value )
	{
		std::ostringstream oss;
		oss << "( ";
		for( std::vector<vl::Vec3f>::const_iterator itr = value.begin();
			itr != value.end(); ++itr )
		{
			if( itr != value.begin() )
				oss << ", ";
			oss << "(" << (*itr)[0] << ", " << (*itr)[1] << ", " << (*itr)[2] << ")";
		}
		oss << " )";
		return oss.str();
	}
};

template <>
class StringForMechModelValue<bool>
{
public:
	std::string operator()( const bool& value )
	{
		if( value )
			return "true";
		else
			return "false";
	}
};

template <typename T, unsigned int N>
class StringForMechModelValue< TinyVec<T, N> >
{
public:
	std::string operator()( const TinyVec<T, N>& value )
	{
		std::ostringstream oss;
		oss << "(";
		StringForMechModelValue<T> conv;
		for( unsigned int i = 0; i < (N-1); ++i )
			oss << conv( value[i] ) << ", ";
		oss << conv( value[N - 1] ) << ")";
		return oss.str();
	}

};

class MechModelValuePair
	: public MechModelObject
{
public:
	virtual ~MechModelValuePair() {}

	template <typename ValueType>
	MechModelValuePair( const std::string& name, const ValueType& value )
		: name_(name), value_( StringForMechModelValue<ValueType>()( value ) )
	{
	}

	void dump( std::ostream& os, const std::string& prefix ) const;

private:
	const std::string name_;
	const std::string value_;
};

class MechModelValueList
	: public MechModelObject
{
public:
	virtual ~MechModelValueList() {}

	MechModelValueList( const std::string& name )
		: name_(name)
	{
	}

	template <typename T>
	void add( const T& value )
	{
		StringForMechModelValue<T> converter;
		values_.push_back( converter( value ) );
	}

	void dump( std::ostream& os, const std::string& prefix ) const;

private:
	const std::string name_;
	std::vector<std::string> values_;
};

template <typename T>
boost::shared_ptr<MechModelValuePair> mechModelValuePair( const std::string& name, const T& value )
{
	return boost::shared_ptr<MechModelValuePair>(
		new MechModelValuePair(name, value) );
}

class MechModelList
	: public MechModelObject
{
public:
	virtual ~MechModelList() {}

	typedef boost::shared_ptr<const MechModelObject> MechModelObjectPtr;

	void add( MechModelObjectPtr child )
	{
		children_.push_back( child );
	}

	void dump( std::ostream& os, const std::string& prefix ) const;

private:
	typedef std::deque<MechModelObjectPtr> MechModelObjectPtrList;
	MechModelObjectPtrList children_;
};

class MechModelScope
	: public MechModelList
{
public:
	virtual ~MechModelScope() {}
	typedef boost::shared_ptr<const MechModelObject> MechModelObjectPtr;
	void dump( std::ostream& os, const std::string& prefix ) const;
};

class MechModelDeclaredObject
	: public MechModelObject
{
public:
	MechModelDeclaredObject( const std::string& type, const std::string& name, bool compressible = true );
	virtual ~MechModelDeclaredObject();

	void add( MechModelObjectPtr object );
	void dump( std::ostream& os, const std::string& prefix ) const;
	void setName( const std::string& name );

private:
	const std::string type_;
	std::string name_;

	typedef std::deque<MechModelObjectPtr> MechModelObjectPtrList;
	MechModelObjectPtrList children_;
	bool compressible_;
};

} // namespace planning

#endif

