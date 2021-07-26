#ifdef WIN32
#pragma once
#endif

#ifndef __KEYABLE_H__
#define __KEYABLE_H__

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <deque>
#include <vector>
#include <iosfwd>


template <typename TimeType, typename ValueType>
class Keyframe
{
public:
	Keyframe( TimeType time, ValueType value )
		: time_(time), value_(value) {}

	bool operator<( const Keyframe<TimeType, ValueType>& other ) const
	{
		return std::less<TimeType>()( time_, other.time_ );
	}

	TimeType time() const		{ return time_; }
	ValueType value() const		{ return value_; }

private:
	TimeType time_;
	ValueType value_;
};

template <typename T>
class KeyableAttribute
	: public Named
{
public:
	KeyableAttribute( const std::string& name )
		: Named( name ) {}
	virtual ~KeyableAttribute() {}

	virtual void set( const T& value ) = 0;
	virtual T get() const = 0;
};

class Interpolator;
// This is going to be hardcoded for floats for the moment
//   due to the GSL interpolation code.

class AttributeWithKeys
	: public Named
{
public:
	typedef float TimeType;
	typedef Keyframe<TimeType, float> KeyframeType;

	typedef boost::shared_ptr< KeyableAttribute<float> > KeyablePtr;
	AttributeWithKeys( KeyablePtr keyable );
	~AttributeWithKeys();

	// sets the key using the current values of the attribute
	void setKey( TimeType time );
	void addKey( const KeyframeType& key );
	void deleteKey( TimeType time );

	// uses the keyframes to interpolate a value for the given time.
	void setTime( TimeType time );

	// Informational stuff:
	// gets the first key after a given time
	TimeType nextKey( TimeType time ) const;

	// gets the first key previous to a given time
	TimeType prevKey( TimeType time ) const;

	TimeType firstKey() const;
	TimeType lastKey() const;

	void write( std::ostream& out );
	void read( std::istream& file );

	typedef std::deque< KeyframeType >::const_iterator keyframe_iterator;
	keyframe_iterator begin_keys() const;
	keyframe_iterator end_keys() const;

	const KeyableAttribute<float>& keyable() const;

	typedef std::deque< KeyframeType > KeyframeList;
	KeyframeList keys()                      { return this->keys_; }
	void setKeys( const KeyframeList& keys ) { this->keys_ = keys; }

private:

	KeyablePtr keyable_;
	KeyframeList keys_;

	mutable boost::scoped_ptr<Interpolator> interpolator_;
};

template <typename T>
class HasKeyable
{
public:
	typedef boost::shared_ptr< KeyableAttribute<T> > KeyableAttributePtr;
	typedef std::deque<KeyableAttributePtr> AttributeList;

	virtual ~HasKeyable() {}
	virtual AttributeList keyable() = 0;
};

template <typename T>
class SimpleKeyable
	: public KeyableAttribute<T>
{
public:
	SimpleKeyable( const std::string& name, T& value )
		: KeyableAttribute<T>(name), value_(value) {}

	void set( const T& value )
	{
		value_ = value;
	}

	T get() const
	{
		return value_;
	}

private:
	T& value_;
};

#endif 

