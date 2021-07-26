#include "stdafx.h"


#include "mechModel.h"

#include <iostream>

namespace planning {

std::string quoted( const std::string& in )
{
	std::string result;
	result.reserve( in.size() + 2 );
	result.push_back( '\"' );
	for( std::string::const_iterator itr = in.begin();
		itr != in.end(); ++itr )
	{
		if( *itr == '\"' )
			result.push_back( '\\' );
		result.push_back( *itr );
	}

	result.push_back( '\"' );
	return result;
}

MechModelDeclaredObject::MechModelDeclaredObject( const std::string& type, const std::string& name, bool compressible )
	: type_(type), name_(name), compressible_(compressible)
{
}

MechModelDeclaredObject::~MechModelDeclaredObject()
{
}

void MechModelDeclaredObject::setName( const std::string& name )
{
	this->name_ = name;
}

void MechModelDeclaredObject::add( MechModelObjectPtr object )
{
	children_.push_back( object );
}

void MechModelDeclaredObject::dump( std::ostream& os, const std::string& prefix ) const
{
	os << prefix << type_;
	if( !name_.empty() )
		os << " " << quoted(name_);

	if( compressible_ && children_.size() == 1 )
	{
		MechModelObjectPtr item = children_[0];
		if( dynamic_cast<const MechModelValuePair*>(item.get()) )
		{
			int lengthSoFar = prefix.size() + type_.size() + name_.size();
			const int lineLength = 20;
			os << std::string(std::max(0, lineLength - lengthSoFar), ' ') << "{ ";
			item->dump( os, std::string() );
			os << " }";
			return;
		}
	}

	os << "\n";
	os << prefix << "{\n";
	for( MechModelObjectPtrList::const_iterator itr = children_.begin();
		itr != children_.end();
		++itr )
	{
		(*itr)->dump( os, prefix + std::string("  ") );
		os << "\n";
	}
	os << prefix << "}  // end " << type_ << "\n";
}


void MechModelScope::dump( std::ostream& os, const std::string& prefix ) const
{
	os << prefix << "{\n";
	MechModelList::dump( os, prefix + "  " );
	os << prefix << "}\n\n";
}

void MechModelList::dump( std::ostream& os, const std::string& prefix ) const
{
	for( MechModelObjectPtrList::const_iterator itr = children_.begin();
		itr != children_.end();
		++itr )
	{
		(*itr)->dump( os, prefix );
		os << "\n";
	}
	os << "\n";
}

void MechModelValuePair::dump( std::ostream& os, const std::string& prefix ) const
{
	os << prefix << name_ << " " << value_ << ";";
}

void MechModelValueList::dump( std::ostream& os, const std::string& prefix ) const
{
	os << prefix << name_;
	for( size_t iValue = 0; iValue < values_.size(); ++iValue )
		os << " " << values_[iValue];
	os << ";";
}


} // namespace planning
