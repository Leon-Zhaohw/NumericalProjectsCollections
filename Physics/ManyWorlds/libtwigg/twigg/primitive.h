#ifndef __PRIMITIVE_H__
#define __PRIMITIVE_H__

#include <boost/array.hpp>
#include "twigg/util.h"

// Abstract for triangle or line (or point, I suppose)
template <unsigned int N>
class GeometricPrimitive
{
public:
	typedef GeometricPrimitive<N> BaseType;
	typedef TinyVec<unsigned int, N> TinyVecType;

	GeometricPrimitive()
	{
	}

	GeometricPrimitive( const TinyVecType& v )
	{
		for( unsigned int iVertex = 0; iVertex < N; ++iVertex )
			vertices_[iVertex] = v[iVertex];
	}

	unsigned int operator[]( unsigned int iVertex ) const
	{
		return vertices_[iVertex];
	}

	unsigned int vertexCount() const	{ return N; }

	bool operator ==( const GeometricPrimitive<N>& other ) const
	{
		for( unsigned int iVertex = 0; iVertex < N; ++iVertex )
		{
			if( std::find( vertices_.begin(), vertices_.end(), other[iVertex] ) == vertices_.end() )
				return false;
		}

		return true;
	}

	bool operator<( const GeometricPrimitive<N>& other ) const
	{
		return std::lexicographical_compare( 
			this->vertices_.begin(), this->vertices_.end(), 
			other.vertices_.begin(), other.vertices_.end() );
	}

private:
	boost::array<unsigned int, N> vertices_;

protected:
	void setVertex( unsigned int iVertex, unsigned int value )
	{
		vertices_[iVertex] = value;
	}
};

class Triangle : public GeometricPrimitive<3>
{
public:
	Triangle( unsigned int v0, unsigned int v1, unsigned int v2 )
		: GeometricPrimitive<3>( TinyVec<unsigned int, 3>(v0, v1, v2) ) {}

		Triangle( const BaseType::TinyVecType& v )
		: GeometricPrimitive<3>(v) {}

	Triangle() {}
};

class Edge : public GeometricPrimitive<2>
{
public:
	Edge(unsigned int v0, unsigned int v1)
		: GeometricPrimitive<2>( TinyVec<unsigned int, 2>( 
			std::min(v0, v1),
			std::max(v0, v1) ) ) { assert( v0 != v1 ); }

	Edge( const BaseType::TinyVecType& v )
		: GeometricPrimitive<2>(v) {}
};

class Point : public GeometricPrimitive<1>
{
public:
	Point(unsigned int v0)
		: GeometricPrimitive<1>( TinyVec<unsigned int, 1>(v0) ) {}

	Point(const BaseType::TinyVecType& v)
		: GeometricPrimitive<1>(v) {}
};

#endif


