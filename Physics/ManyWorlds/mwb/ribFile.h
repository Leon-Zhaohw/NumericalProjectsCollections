#ifdef WIN32
#pragma once
#endif

#include "physicsFwd.h"
#include <iosfwd>

namespace planning {

class RibObject
{
public:
	virtual void dump( std::ostream& os, const std::string& prefix ) const = 0;
	virtual ~RibObject();
};

template <typename T>
class StringForRibValue
{
public:
	std::string operator()( const T& value )
	{
		return boost::lexical_cast<std::string>(value);
	}
};

template <>
class StringForRibValue< vl::Mat4f >
{
public:
	std::string operator()( const vl::Mat4f& mat )
	{
		std::ostringstream oss;
		oss << "[";
		for( vl::Int jCol = 0; jCol < 4; ++jCol )
		{
			for( vl::Int iRow = 0; iRow < 4; ++iRow )
			{
				oss << mat[iRow][jCol];
				if( iRow != 3 )
					oss << " ";
			}

			if( jCol != 3 )
				oss << "   ";
		}
		oss << "]";
		return oss.str();
	}
};

template <>
class StringForRibValue< vl::Mat4d >
{
public:
	std::string operator()( const vl::Mat4d& mat )
	{
		std::ostringstream oss;
		oss << "[";
		for( vl::Int jCol = 0; jCol < 4; ++jCol )
		{
			for( vl::Int iRow = 0; iRow < 4; ++iRow )
			{
				oss << mat[iRow][jCol];
				if( iRow != 3 )
					oss << " ";
			}

			if( jCol != 3 )
				oss << "   ";
		}
		oss << "]";
		return oss.str();
	}
};


template <>
class StringForRibValue< vl::Vec3f >
{
public:
	std::string operator()( const vl::Vec3f& vec )
	{
		std::ostringstream oss;
		oss << "[" << vec[0];
		for( vl::Int i = 1; i < 3; ++i )
			oss << " " << vec[i];
		oss << "]";
		return oss.str();
	}
};

template <>
class StringForRibValue< std::vector<float> >
{
public:
	std::string operator()( const std::vector<float>& values )
	{
		std::ostringstream oss;
		oss << "[";
		for( size_t iValue = 0; iValue < values.size(); ++iValue )
			oss << values[iValue] << " ";
		oss << "]";
		return oss.str();
	}
};

template <>
class StringForRibValue<const char*>
{
public:
	std::string operator()( const char* value )
	{
		std::ostringstream oss;
		oss << "\"";
		oss << value;
		oss << "\"";
		return oss.str();
	}
};

template <>
class StringForRibValue<std::string>
{
public:
	std::string operator()( const std::string& value )
	{
		std::ostringstream oss;
		oss << "\"";
		oss << value;
		oss << "\"";
		return oss.str();
	}
};


class RibObjectWithAttributes
	: public RibObject
{
public:
	RibObjectWithAttributes( const std::string& name );

	template <typename ValueType>
	void addAttribute( const ValueType& attribute )
	{
		StringForRibValue<ValueType> converter;
		attributes_.push_back( converter(attribute) );
	}

	void dump( std::ostream& os, const std::string& prefix ) const;

private:
	std::deque<std::string> attributes_;
	std::string name_;
};

class RibMetaComment
	: public RibObject
{
public:
	RibMetaComment( const std::string& text )
		: text_(text) {}
	void dump( std::ostream& os, const std::string& prefix ) const;

private:
	std::string text_;
};

class RibComment
	: public RibObject
{
public:
	RibComment( const std::string& text )
		: text_(text) {}
	void dump( std::ostream& os, const std::string& prefix ) const;

private:
	std::string text_;
};

class RibNamedBlock
	: public RibObject
{
public:
	RibNamedBlock( const std::string& name );
	void addChild( RibObjectPtr child );
	void dump( std::ostream& os, const std::string& prefix ) const;

	template <typename ValueType>
	void addAttribute( const ValueType& attribute )
	{
		StringForRibValue<ValueType> converter;
		attributes_.push_back( converter(attribute) );
	}

private:
	std::string name_;
	
	std::deque<std::string> attributes_;
	typedef std::deque<RibObjectPtr> RibObjectList;
	RibObjectList children_;
};

class RibBlock
	: public RibObject
{
public:
	void addChild( RibObjectPtr child );
	void dump( std::ostream& os, const std::string& prefix ) const;

private:
	typedef std::deque<RibObjectPtr> RibObjectList;
	RibObjectList children_;
};

template <typename T>
RibObjectPtr ribObjectWithAttribute( const std::string& name, const T& attributeValue )
{
	boost::shared_ptr<RibObjectWithAttributes> result( new RibObjectWithAttributes(name) );
	result->addAttribute( attributeValue );
	return result;
}

} // namespace planning
