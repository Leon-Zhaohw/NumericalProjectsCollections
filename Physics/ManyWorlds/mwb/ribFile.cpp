#include "stdafx.h"
#include "ribFile.h"

namespace planning {


RibObject::~RibObject()
{
}

RibObjectWithAttributes::RibObjectWithAttributes(const std::string& name)
	: name_(name)
{
}

void RibObjectWithAttributes::dump( std::ostream& os, const std::string& prefix ) const
{
	os << prefix;
	os << name_;
	for( std::deque<std::string>::const_iterator attributeItr = this->attributes_.begin();
		attributeItr != this->attributes_.end(); ++attributeItr )
	{
		os << " " << *attributeItr;
	}
	os << "\n";
}

void RibMetaComment::dump( std::ostream& os, const std::string& prefix ) const
{
	os << prefix << "##" << this->text_ << "\n";
}

void RibComment::dump( std::ostream& os, const std::string& prefix ) const
{
	os << prefix << "# " << this->text_ << "\n";
}

RibNamedBlock::RibNamedBlock( const std::string& name )
	: name_(name)
{
}

void RibNamedBlock::addChild( RibObjectPtr child )
{
	this->children_.push_back( child );
}

void RibNamedBlock::dump( std::ostream& os, const std::string& prefix ) const
{
	os << prefix << name_ << "Begin";
	for( std::deque<std::string>::const_iterator attributeItr = this->attributes_.begin();
		attributeItr != this->attributes_.end(); ++attributeItr )
	{
		os << " " << *attributeItr;
	}
	os << "\n";

	std::for_each( this->children_.begin(), this->children_.end(),
		boost::bind( &RibObject::dump, _1, boost::ref(os), std::string("  ") + prefix ) );
	os << prefix << name_ << "End\n\n";
}

void RibBlock::addChild( RibObjectPtr child )
{
	this->children_.push_back( child );
}

void RibBlock::dump( std::ostream& os, const std::string& prefix ) const
{
	std::for_each( this->children_.begin(), this->children_.end(),
		boost::bind( &RibObject::dump, _1, boost::ref(os), std::string("  ") + prefix ) );
}

} // namespace planning


