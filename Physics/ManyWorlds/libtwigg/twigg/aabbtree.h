#ifndef __AABBTREE_H__

#define __AABBTREE_H__

#include <vector>
#include <queue>
#include <stack>
#include "twigg/boundingbox.h"
#include "twigg/util.h"
#include "twigg/screenSpace.h"

class CollisionObj;
class CollisionHandler;
class Cloth;

const unsigned int splitAxes = 2;
const unsigned int degree = 1 << splitAxes;

template< typename Vec3Array >
class PositionWrapper
{
public:
	PositionWrapper( const Vec3Array& positions )
		: positions_(positions) {}

	template <typename Primitive>
	BoundingBox3d operator()( const Primitive& primitive ) const
	{
		// hopefully uses return value optimization
		BoundingBox3d bounds;
		for( unsigned int iVertex = 0; iVertex < primitive.vertexCount(); ++iVertex )
			bounds.expand( positions_[ primitive[iVertex] ] );

		return bounds;
	}

private:
	const Vec3Array& positions_;
};

template <typename Primitive, typename BoundingBox = BoundingBox3d>
class AABBTree
{
public:
	typedef Primitive PrimitiveType;

private:
	class AABBNode
	{
	public:
		AABBNode(const Primitive& primitive)
			: isLeaf_(true), primitive_(primitive)
		{
		}

		void addChild( unsigned int childLocation )
		{
			if( isLeaf_ )
			{
				isLeaf_ = false;
				childrenStart_ = childLocation;
			}

			childrenEnd_ = childLocation + 1;
		}

		~AABBNode()
		{
		}

		BoundingBox bounds() const				{ return bounds_; }

		bool leaf() const						{ return isLeaf_; }

		unsigned int children_begin() const		{ assert( !leaf() ); return childrenStart_; }
		unsigned int children_end() const		{ assert( !leaf() ); return childrenEnd_; }

		Primitive primitive() const				{ return primitive_; }

		void expand( const AABBNode& node )		{ bounds_.expand( node.bounds_ ); }

		template <typename T>
		void expand( const T& amount )			{ bounds_.expand( amount ); }
		void reset()							{ bounds_.nullify(); }

	private:
		BoundingBox bounds_;
		bool isLeaf_;

		// Due to the way we order the nodes in the array, all the nodes in
		// [childrenStart, childrenEnd) are children
		unsigned int childrenStart_;
		unsigned int childrenEnd_;

		Primitive primitive_;
	};

	// We'll put our tree in an array, just as in
	// [van den Bergen, 1998]).
public:
	// @todo do something to fix this!
	typedef std::deque<AABBNode> NodeArray;
	NodeArray nodes_;

private:
	class SeparatingAxis
	{
	public:
		int axis;
		double length;

		SeparatingAxis(int ax, const typename BoundingBox::point_type& diff)
			: axis( ax ), length( diff[ax] )
		{
		}

		bool operator< ( const SeparatingAxis& other ) const
		{
			return length < other.length;
		}
	};

	template <typename BoundsForPrimitive>
	class PrimitiveComparator
	{
	public:
		PrimitiveComparator( const SeparatingAxis& axis, const BoundsForPrimitive& bounds )
			: axis_( axis ), bounds_(bounds)
		{
		}

		// projects each primitive onto the current separating axis (using the bbox midpoint)
		//   and returns the lesser according to that measure
		bool operator()( const Primitive& lhs, const Primitive& rhs ) const
		{
			return projectionMidpoint(lhs) <
				projectionMidpoint(rhs);
		}

	private:
		const SeparatingAxis axis_;
		const BoundsForPrimitive& bounds_;

		double projectionMidpoint(const Primitive& primitive) const
		{
			BoundingBox bbox = bounds_(primitive);
			return 0.5*(bbox.minimum()[axis_.axis] + bbox.maximum()[axis_.axis]);
		}
	};


public:
	template <typename PrimitiveList, typename PrimitiveBounds>
	AABBTree( const PrimitiveList& primitives, 
		const PrimitiveBounds& boundsForPrimitive )
	{
		// the rest of this function will assume that the set of objects is nonempty
		if( primitives.empty() )
			return;

		// first element of the pair tells us what the node's parent is; 
		//   second element tells us all the primitives that are children of
		//   the node.
		std::queue< std::pair<unsigned int, PrimitiveList > > nodeQueue;

		// the top node in the tree will have all primitives as its children
		nodeQueue.push( std::make_pair(0, primitives) );

		while( !nodeQueue.empty() )
		{
			// create a node:
			typename NodeArray::iterator nodeItr =
				nodes_.insert( nodes_.end(), 
					AABBNode( *nodeQueue.front().second.begin() ) );

			// assign its parent
			if( nodeItr != nodes_.begin() )
			{
				int parent = nodeQueue.front().first;
				nodes_[parent].addChild( std::distance( nodes_.begin(), nodeItr ) );
			}

			PrimitiveList primitives = nodeQueue.front().second;
			nodeQueue.pop();

			// if there is only one child, then the node is already a leaf node
			if( primitives.size() == 1 )
				continue;

			if( primitives.size() <= degree )
			{
				// Base case; add a single node for each
				// primitive and return.
				for( typename PrimitiveList::const_iterator primitiveItr = primitives.begin();
					primitiveItr != primitives.end();
					++primitiveItr )
				{
					PrimitiveList oneTriangleVec(1, *primitiveItr);
					nodeQueue.push( std::make_pair( std::distance(nodes_.begin(), nodeItr),
						oneTriangleVec ) );
				}
				continue;
			}

			// find the bounds for the entire set of primitives
			BoundingBox bounds;
			for( typename PrimitiveList::const_iterator primitiveItr = primitives.begin();
				primitiveItr != primitives.end();
				++primitiveItr )
			{
				bounds.expand( boundsForPrimitive( *primitiveItr ) );
			}

			// Now, subdivide the remaining primitives according to the 
			//   largest axes (we allow up to splitAxes splits, resulting
			//   in 2^splitAxes children)
			std::vector<SeparatingAxis> axes;
			for( size_t iAxis = 0; iAxis < 3; ++iAxis )
				axes.push_back(
					SeparatingAxis(iAxis,
						bounds.maximum() - bounds.minimum() ) );

			std::sort( axes.begin(), axes.end() );

			// initialize the face subsets to all faces.  We will progressively
			//   split these into smaller and smaller subsets
			std::vector< PrimitiveList > faceSubsets(1, primitives);

			for( typename std::vector<SeparatingAxis>::reverse_iterator axesItr = axes.rbegin();
				axesItr != axes.rend();
				++axesItr )
			{
				if( static_cast<unsigned int>(std::distance(axes.rbegin(), axesItr)) + 1 > splitAxes )
					break;

				// for each of the current subsets
				std::vector< PrimitiveList > newSubsets;
				newSubsets.reserve( 2*faceSubsets.size() );
				for( typename std::vector< PrimitiveList >::iterator subsetItr = faceSubsets.begin();
					subsetItr != faceSubsets.end();
					++subsetItr )
				{
					PrimitiveList& subset = *subsetItr;
					typename PrimitiveList::size_type m = subset.size()/2;
					typename PrimitiveList::iterator half = subset.begin() + m;
					
					assert( !subset.empty() && m < subset.size() );
					assert( half != subset.end() ); // should never be true as subset is nonempty

					// nth_element using this comparator will find the median using the heuristic
					//   specified in PrimitiveComparator (midpoint of projection onto given axis)
					PrimitiveComparator<PrimitiveBounds> comparator(*axesItr, boundsForPrimitive);
#ifdef WIN32
					std::nth_element( subset.begin(), half, subset.end(), comparator );
#else
					// weird bug in gcc nth_element causes crashes here, if someone can
					//   tell me why I'd be greatly appreciative.  partial_sort has the
					//   same effective result at significantly greater cost
					std::partial_sort( subset.begin(), half+1, subset.end(), comparator );
#endif

					// split subset in two, stick new subsets into new subset list
					newSubsets.push_back(
						PrimitiveList(subset.begin(), half));
					newSubsets.push_back(
						PrimitiveList(half, subset.end()));
				}

				faceSubsets.swap( newSubsets );
			}

			// Now we have our subsets.  We can stick the new children into the node queue:
			for( typename std::vector< PrimitiveList >::iterator subsetItr = faceSubsets.begin();
				subsetItr != faceSubsets.end();
				++subsetItr )
			{
				if( subsetItr->empty() )
					continue;

				nodeQueue.push( std::make_pair( std::distance(nodes_.begin(), nodeItr),
					*subsetItr ) );
			}
		}

#ifdef CHECK_TOUCHED
		// this is just a sanity check to make sure we've placed every primitive into some leaf node
		std::deque<bool> touched( primitives.size(), false );
#endif

		std::stack<unsigned int> nodeStack;
		nodeStack.push(0);
		while( !nodeStack.empty() )
		{
			unsigned int current = nodeStack.top();
			nodeStack.pop();

			if( nodes_[current].leaf() )
			{
#ifdef CHECK_TOUCHED
				int iPrimitive = std::distance( primitives.begin(),
					std::find( primitives.begin(), primitives.end(), 
						nodes_[current].primitive() ) );

				touched[ iPrimitive ] = true;
#endif
			}
			else
			{
				for( unsigned int iChild = nodes_[current].children_begin();
					iChild < nodes_[current].children_end();
					++iChild )
				{
					nodeStack.push( iChild );
				}
			}
		}

#ifdef CHECK_TOUCHED
		for( unsigned int iPrimitive = 0; iPrimitive < touched.size(); ++iPrimitive )
		{
			if( !touched[iPrimitive] )
				std::cout << "Not touched! " << iPrimitive << std::endl;
			//else
			//	std::cout << "ok." << iPrimitive << std::endl;
		}
#endif

		update( boundsForPrimitive, 0.0 );
	}

	BoundingBox bounds() const
	{
		return nodes_.begin()->bounds();
	}

	template <typename BoundsForPrimitive>
	void update( const BoundsForPrimitive& boundsForPrimitive, double padding )
	{
		for( typename NodeArray::reverse_iterator nodeItr = nodes_.rbegin();
			nodeItr != nodes_.rend();
			++nodeItr )
		{
			nodeItr->reset();

			if( nodeItr->leaf() )
			{
				nodeItr->expand( boundsForPrimitive( nodeItr->primitive() ) );
				nodeItr->expand( padding );
			}
			else
			{
				for( unsigned int iChild = nodeItr->children_begin();
					iChild < nodeItr->children_end();
					++iChild )
				{
					nodeItr->expand( nodes_[iChild] );
				}
			}
		}
	}

	template <typename CollisionHandler, typename RightAABBTree>
	bool intersectTree( const RightAABBTree& rhs, CollisionHandler handler ) const
	{
		if( nodes_.empty() )
			return false;

		std::stack< std::pair< unsigned int, unsigned int > > toHandle;
		toHandle.push( std::make_pair(0,0) );
		bool result = false;
		while( ! toHandle.empty() )
		{
			std::pair<unsigned int, unsigned int> current = toHandle.top();
			toHandle.pop();
			if( intersects( this->nodes_[ current.first ].bounds(), 
						rhs.nodes_[ current.second ].bounds(), 0.0 ) )
			{
				if( nodes_[ current.first ].leaf() &&
					rhs.nodes_[ current.second ].leaf() )
				{
					result = handler( nodes_[ current.first ].primitive(),
						rhs.nodes_[ current.second ].primitive() ) || result;
				}
				else if( nodes_[ current.first ].leaf() )
				{
					for( unsigned int iChild = rhs.nodes_[ current.second ].children_begin();
						iChild != rhs.nodes_[ current.second ].children_end();
						++iChild )
					{
						toHandle.push( std::make_pair( current.first, iChild ) );
					}
				}
				else if( rhs.nodes_[ current.second ].leaf() )
				{
					for( unsigned int iChild = nodes_[ current.first ].children_begin();
						iChild != nodes_[ current.first ].children_end();
						++iChild )
					{
						toHandle.push( std::make_pair( iChild, current.second ) );
					}
				}
				else
				{
					for( unsigned int iChild = nodes_[ current.first ].children_begin();
						iChild != nodes_[ current.first ].children_end();
						++iChild )
					{
						for( unsigned int jChild = rhs.nodes_[ current.second ].children_begin();
							jChild != rhs.nodes_[ current.second ].children_end();
							++jChild )
						{
							toHandle.push( std::make_pair( iChild, jChild ) );
						}
					}
				}
			}
		}

		return result;
	}

	template <typename Handler, typename Ray>
	bool intersectRay( const Ray& ray, Handler& handler ) const
	{
		if( nodes_.empty() )
			return false;

		std::stack< unsigned int > toHandle;
		toHandle.push( 0 );
		bool result = false;
		while( !toHandle.empty() )
		{
			unsigned int current = toHandle.top();
			toHandle.pop();
			double t0, t1;
			if( this->nodes_[current].bounds().intersect( ray, t0, t1 ) )
			{
				if( nodes_[current].leaf() )
				{
					if( handler( nodes_[ current ].primitive(), ray ) )
						return true;
				}
				else
				{
					for( unsigned int iChild = nodes_[ current ].children_begin();
						iChild != nodes_[ current ].children_end();
						++iChild )
						toHandle.push( iChild );
				}
			}
		}

		return false;
	}

	template <typename Handler, typename BoundingBox2>
	bool intersect( const BoundingBox2& box, const ScreenSpaceConverter& converter, Handler& handler ) const
	{
		if( nodes_.empty() )
			return false;

		std::stack< unsigned int > toHandle;
		toHandle.push( 0 );
		bool result = false;
		while( !toHandle.empty() )
		{
			unsigned int current = toHandle.top();
			toHandle.pop();
			BoundingBox2 projectedBox;
			boost::array<typename BoundingBox::point_type, 2> params = 
				{{	this->nodes_[current].bounds().minimum(), 
					this->nodes_[current].bounds().maximum() }};
			for( unsigned int i = 0; i <= 1; ++i )
				for( unsigned int j = 0; j <= 1; ++j )
					for( unsigned int k = 0; k <= 1; ++k )
					{
						projectedBox.expand( 
							strip( converter.toScreenSpace(
							vl::Vec3d( (params[i])[0], (params[j])[1], (params[k])[2] ) ) ) ); 
					}

			if( intersects( projectedBox, box ) )
			{
				if( nodes_[current].leaf() )
				{
					if( handler(nodes_[ current ].primitive()) )
						return true;
				}
				else
				{
					for( unsigned int iChild = nodes_[ current ].children_begin();
						iChild != nodes_[ current ].children_end();
						++iChild )
						toHandle.push( iChild );
				}
			}
		}

		return false;
	}
};

#endif

