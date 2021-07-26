#include "stdafx.h"

#include "twigg/keyable.h"
#include "twigg/interpolation.h"

#include <boost/spirit/core.hpp>
#include <boost/spirit/attribute.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/utility/chset.hpp>
#include <boost/spirit/utility/escape_char.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/utility/loops.hpp>
#include <boost/spirit/iterator.hpp>
#include <boost/spirit/error_handling/exceptions.hpp>

#include <iostream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

AttributeWithKeys::AttributeWithKeys( KeyablePtr keyable )
	: Named( keyable->name() ), keyable_(keyable)
{
}

AttributeWithKeys::~AttributeWithKeys()
{
}

void AttributeWithKeys::addKey( const KeyframeType& keyframe )
{
	typedef KeyframeList::iterator KeyframeListIter;
	typedef std::pair<KeyframeListIter, KeyframeListIter> KeyframeListIterPr;
	KeyframeListIterPr pr = std::equal_range( keys_.begin(), keys_.end(), keyframe );
	KeyframeListIter newPos = keys_.erase( pr.first, pr.second );
	keys_.insert( newPos, keyframe );

	interpolator_.reset();
}

void AttributeWithKeys::setKey( TimeType time )
{
	KeyframeType keyframe( time, keyable_->get() );
	addKey( keyframe );
}

void AttributeWithKeys::deleteKey( TimeType time )
{
	typedef KeyframeList::iterator KeyframeListIter;
	typedef std::pair<KeyframeListIter, KeyframeListIter> KeyframeListIterPr;
	KeyframeType keyframe( time, keyable_->get() );
	KeyframeListIterPr pr = std::equal_range( keys_.begin(), keys_.end(), keyframe );
	keys_.erase( pr.first, pr.second );

	interpolator_.reset();
}

AttributeWithKeys::TimeType AttributeWithKeys::nextKey( TimeType time ) const
{
	typedef KeyframeList::const_iterator KeyframeListIter;
	KeyframeListIter iter = std::upper_bound( keys_.begin(), keys_.end(),
		KeyframeType( time, 0.0 ) );
	if( iter == keys_.end() )
		return lastKey();
	return iter->time();
}

AttributeWithKeys::TimeType AttributeWithKeys::prevKey( TimeType time ) const
{
	typedef KeyframeList::const_iterator KeyframeListIter;
	KeyframeListIter iter = std::lower_bound( keys_.begin(), keys_.end(),
		KeyframeType( time, 0.0 ) );
	if( iter == keys_.begin() )
		return firstKey();
	return (iter-1)->time();
}

AttributeWithKeys::TimeType AttributeWithKeys::firstKey() const
{
	if( keys_.empty() )
		return TimeType(0);

	return keys_.front().time();
}

AttributeWithKeys::TimeType AttributeWithKeys::lastKey() const
{
	if( keys_.empty() )
		return TimeType(0);

	return keys_.back().time();
}

void AttributeWithKeys::setTime( TimeType time )
{
	if( keys_.empty() )
		return;

	if( time >= keys_.back().time() )
	{
		keyable_->set( keys_.back().value() );
		return;
	}

	if( time <= keys_.front().time() )
	{
		keyable_->set( keys_.front().value() );
		return;
	}

	if( !interpolator_ )
	{
		std::vector<double> x; x.reserve( keys_.size() );
		std::vector<double> y; y.reserve( keys_.size() );

		for( KeyframeList::const_iterator keyItr = keys_.begin();
			keyItr != keys_.end(); ++keyItr )
		{
			x.push_back( keyItr->time() );
			y.push_back( keyItr->value() );
		}

		assert( x.size() == y.size() );
#ifdef USE_GSL
		const gsl_interp_type *interpType;
		if( x.size() < 3 )
			interpType = gsl_interp_linear;
		else if( x.size() < 5 )
			interpType = gsl_interp_cspline;
		else
			interpType = gsl_interp_akima;
#endif
		interpolator_.reset( new Interpolator(&x[0], &y[0], x.size()
#ifdef USE_GSL
			, interpType
#endif
			) );
	}

	keyable_->set( (*interpolator_)( time ) );
}

void AttributeWithKeys::write( std::ostream& out )
{
	/*
	out << "keys " << keys_.size() << "\n";
	for( KeyframeList::const_iterator keyItr = keys_.begin(); 
		keyItr != keys_.end(); ++keyItr )
	{
		out << "key " << keyItr->time() << " " << keyItr->value() << "\n";
	}
	*/

	// use binary instead
	int numKeys = keys_.size();
	out.write( reinterpret_cast<const char*>(&numKeys), sizeof(int) );

	for( KeyframeList::const_iterator keyItr = keys_.begin(); 
		keyItr != keys_.end(); ++keyItr )
	{
		int time = keyItr->time();
		float value = keyItr->value();
		out.write( reinterpret_cast<const char*>(&time), sizeof(int) );
		out.write( reinterpret_cast<const char*>(&value), sizeof(float) );
	}
}

void AttributeWithKeys::read( std::istream& ifs )
{
	this->keys_.clear();

	int numKeys;
	ifs.read( reinterpret_cast<char*>( &numKeys ), sizeof( int ) );

	for( int i = 0; i < numKeys; ++i )
	{
		int time;
		float value;
		ifs.read( reinterpret_cast<char*>(&time), sizeof(int) );
		ifs.read( reinterpret_cast<char*>(&value), sizeof(float) );

		this->keys_.push_back( KeyframeType(time, value) );
	}

	/*
	using namespace boost::spirit;

	std::streamsize maxLine = 100000;
	std::vector<char> line( maxLine, 0 );
	while( strlen(&line[0]) == 0 )
		ifs.getline( &line[0], maxLine );
	int numKeys = 0;
	parse_info<char const*> info = parse( &line[0], &line[strlen(&line[0])],
		str_p("keys") >> int_p[ assign(numKeys) ], space_p );
	if( !info.hit )
		throw CreationException( "Expected: 'keys'" );

	for( unsigned int i = 0; i < numKeys; ++i )
	{
		ifs.getline( &line[0], maxLine );
		const char* startLine = &line[0];
		while( (startLine < &line[maxLine]) && (isspace(*startLine)) && (*startLine != 0) )
			++startLine;
		const char* endLine = startLine;
		while( (endLine < &line[maxLine]) && (*endLine != 0) )
			++endLine;

		if( startLine == endLine )
			continue;

		int time;
		float value;
		info = parse( startLine, endLine,
			str_p("key") >> int_p[ assign(time) ] >> real_p[ assign(value) ], space_p );
		if( !info.hit )
			throw CreationException( "Expected: 'key'" );

		this->keys_.push_back( KeyframeType(time, value) );
	}
	*/
}

AttributeWithKeys::keyframe_iterator AttributeWithKeys::begin_keys() const
{
	return keys_.begin();
}

AttributeWithKeys::keyframe_iterator AttributeWithKeys::end_keys() const
{
	return keys_.end();
}

const KeyableAttribute<float>& AttributeWithKeys::keyable() const
{
	return *this->keyable_;
}
