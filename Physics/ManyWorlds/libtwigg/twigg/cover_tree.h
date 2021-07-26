
#include <boost/pool/object_pool.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/utility.hpp>

#include <cassert>
#include <cmath>
#include <vector>
#include <deque>

// utility classes
template <typename Pair>
struct compare1st
{
	bool operator()( const Pair& left, const Pair& right ) const
	{
		std::less<typename Pair::first_type> comparator;
		return comparator(left.first, right.first);
	}
};

template <typename Pair>
struct less_than_predicate
{
public:
	less_than_predicate( const typename Pair::first_type& value )
		: value_( value ) {}

	bool operator()( const Pair& left ) const
	{
		std::less<typename Pair::first_type> comparator;
		return comparator(left.first, value_);
	}

private:
	typename Pair::first_type value_;
};


template <typename Point, typename Distance>
class cover_tree
	: boost::noncopyable
{
public:
	typedef Point point_type;
	typedef Distance distance_type;
	typedef typename Point::value_type real_type;
	typedef short int level_type;


private:
	class node_type
	{
		template <class NodeType>
		class child_iterator_base
			:	public boost::iterator_adaptor<
				child_iterator_base<NodeType>,
				NodeType*,
				boost::use_default,
				boost::forward_traversal_tag>
		{
		private:
			struct enabler {};

		public:
			child_iterator_base()
				: child_iterator_base::iterator_adaptor_(0) {}

			explicit child_iterator_base(NodeType* p)
				: child_iterator_base::iterator_adaptor_(p) {}

			template <class OtherNodeType>
			child_iterator_base( child_iterator_base<OtherNodeType> const& other,
				typename boost::enable_if<
					boost::is_convertible<OtherNodeType*, NodeType*>,
					enabler >::type = enabler() )
				: child_iterator_base::iterator_adaptor_(other.base()) {}

		private:
			friend class boost::iterator_core_access;
			void increment() { this->base_reference() = this->base()->nextSibling_; }
		};

	public:
		node_type( const point_type& p, const level_type level )
			: point_( p ), level_( level ), firstChild_(0), nextSibling_(0)    {}

		level_type level() const              { return level_; }
		const point_type& point() const       { return point_; }

		typedef child_iterator_base<node_type> child_iterator;
		typedef child_iterator_base<const node_type> const_child_iterator;

		child_iterator children_begin()
		{
			return child_iterator( this->firstChild_ );
		}

		child_iterator children_end()
		{
			return child_iterator();
		}

		const_child_iterator children_begin() const
		{
			return const_child_iterator( this->firstChild_ );
		}

		const_child_iterator children_end() const
		{
			return const_child_iterator();
		}

		void add_child(node_type* child)
		{
			if( firstChild_ == 0 )
				firstChild_ = child;
			else
			{
				child->nextSibling_ = this->firstChild_;
				this->firstChild_ = child;
			}
		}


		private:
			point_type point_;
			level_type level_;

			node_type* firstChild_;
			node_type* nextSibling_;
		};


public:
	template <typename Sequence>
	cover_tree( const Sequence& points )
		: root_(0), numPoints_(0)
	{
		for( typename Sequence::const_iterator iter = points.begin();
			iter != points.end(); ++iter )
		{
			this->insert( *iter );
		}
	}

	cover_tree()
		: root_(0), numPoints_(0)
	{
	}

	void insert( const point_type& p );

	void check_invariants() const;

	size_t size() const
	{
		return numPoints_;
	}

	std::vector<point_type> k_nearest_neighbors( const point_type& point, size_t n = 1 ) const;
	std::vector<point_type> epsilon_nearest_neighbors( const point_type& point, real_type epsilon ) const;

private:
	void check_child_invariants(const node_type* node) const;

	typedef std::pair< real_type, const node_type* > distance_node_type;

	class epsilon_bound
	{
		real_type epsilon_;

	public:
		epsilon_bound( const real_type epsilon )
			: epsilon_(epsilon) {}

		real_type operator()( std::deque< distance_node_type >& Q ) const
		{
			return epsilon_;
		}
	};

	class epsilon_select
	{
		real_type epsilon_;

	public:
		epsilon_select( const real_type epsilon )
			: epsilon_(epsilon) {}

		typename std::vector<point_type> operator()( 
			typename std::deque< distance_node_type >& allCoverSets ) const
		{
			less_than_predicate<distance_node_type> pred( epsilon_ );
			typename std::deque< distance_node_type >::iterator mid = 
				std::partition( allCoverSets.begin(), allCoverSets.end(), pred );
			std::sort( allCoverSets.begin(), mid, compare1st<distance_node_type>() );

			typename std::vector<point_type> result;
			result.reserve( std::distance( allCoverSets.begin(), mid ) );
			for( typename std::deque< distance_node_type >::const_iterator iter = allCoverSets.begin();
				iter != mid; ++iter )
			{
				result.push_back( iter->second->point() );
			}

			return result;
		}
	};

	class knn_bound
	{
		size_t k_;

	public:
		knn_bound( const size_t k )
			: k_(k) {}

		real_type operator()( std::deque< distance_node_type >& Q ) const
		{
			assert( !Q.empty() );
			typename std::deque< distance_node_type >::iterator kth = 
				Q.begin() + (std::min)(k_ - 1, Q.size() - 1);
			std::nth_element( Q.begin(), 
				kth, 
				Q.end(), 
				compare1st<distance_node_type>() );
			return kth->first;
		}
	};

	class knn_select
	{
		size_t k_;

	public:
		knn_select( const size_t k )
			: k_(k) {}

		std::vector<point_type> operator()( std::deque< distance_node_type >& allCoverSets ) const
		{
			typename std::deque< distance_node_type >::iterator kth = 
				allCoverSets.begin() + (std::min)(k_, allCoverSets.size());
			std::partial_sort( allCoverSets.begin(), kth, allCoverSets.end(),
				compare1st<distance_node_type>() );

			std::vector<point_type> result;
			result.reserve( k_ );
			for( typename std::deque< distance_node_type >::const_iterator iter = allCoverSets.begin();
				iter != kth; ++iter )
			{
				result.push_back( iter->second->point() );
			}

			return result;
		}
	};

	template <typename Bound, typename Select>
	std::vector<point_type> generic_neighbor_search( 
		const point_type& point, 
		const Bound& bound, 
		const Select& select ) const;


	node_type* insert( const point_type& p, const std::deque<node_type*>& Q_i, level_type level );
	// insert all children of this node.
	//   useful for rebuilding tree
	void insert( node_type* node );

//	real_type distance( const point_type& p, const std::deque<node_type*>& Q_i );
	real_type distance( const point_type& p, const point_type& q ) const
	{
		distance_type dist;
		return dist(p, q);
	}

	real_type scale( level_type level ) const;

	boost::object_pool<node_type> nodePool_;
	node_type* root_;

	size_t numPoints_;
};


template <typename Point, typename Distance>
void cover_tree<Point, Distance>::insert( const Point& p )
{
	if( root_ == 0 )
	{
		// start out with level = -\infty
		// when we get our next point it will be too far away
		//   and so we will automagically build the tree with the 
		//   correct scale.
		root_ = nodePool_.construct( 
			p, 
			boost::numeric::bounds<level_type>::lowest() );

		numPoints_ = 1;
	}
	else
	{
		// if the new point cannot be a child of the root,
		// then we need to induce a parent/child relation between 
		// root_ and p.
		real_type d = distance(root_->point(), p);
		if( d < std::numeric_limits<real_type>::epsilon() )
			return;

		if( d > scale(root_->level()) )
		{
			level_type level = boost::numeric_cast<level_type>( 
#ifdef _WIN32
				_logb( d ) 
#else
				logb( d )
#endif
					) + 1;
			node_type* oldRoot = nodePool_.construct( p, level );
			std::swap( root_, oldRoot );
			numPoints_ = 1;

			// now, dfs and reinsert all the points into the tree.
			insert( oldRoot );
		}
		else
		{
			std::deque<node_type*> Q_i( 1, root_ );
			node_type* parent = this->insert( p, Q_i, root_->level() );
			assert( parent != 0 );

			++numPoints_;
		}
	}
}


template <typename Point, typename Distance>
std::vector<Point> cover_tree<Point, Distance>::epsilon_nearest_neighbors( 
	const Point& p, real_type epsilon ) const
{
	assert( epsilon >= 0.0 );

	epsilon_bound bound( epsilon );
	epsilon_select select( epsilon );

	return generic_neighbor_search( p, bound, select );
}

template <typename Point, typename Distance>
std::vector<Point> cover_tree<Point, Distance>::k_nearest_neighbors( const Point& p, size_t k ) const
{
	assert( k > 0 );

	knn_bound bound( k );
	knn_select select( k );

	std::vector<Point> result = 
		generic_neighbor_search( p, bound, select );
	assert( result.size() == (std::min)(k, this->size()) );
	return result;
}

template <typename Point, typename Distance>
template <typename Bound, typename Select>
std::vector<Point> cover_tree<Point, Distance>::generic_neighbor_search( 
	const Point& p, 
	const Bound& bound, 
	const Select& select ) const
{
	if( root_ == 0 )
		return std::vector<Point>();

	std::deque< distance_node_type > allCoverSets;
	std::deque< distance_node_type > Q_i;

	Q_i.push_back( std::make_pair( distance(root_->point(), p), root_ ) );
	while( !Q_i.empty() )
	{
		std::deque< distance_node_type > Q;
		level_type level = Q_i.front().second->level();

		// consider the set of children of Q_i
		// Q = { Children(q) : q \in Q_i }
		for( typename std::deque< distance_node_type >::const_iterator Q_i_iter = Q_i.begin();
			Q_i_iter != Q_i.end(); ++Q_i_iter )
		{
			const node_type* node = Q_i_iter->second;
			assert( node->level() == level );

			for( typename node_type::const_child_iterator childItr = node->children_begin();
				childItr != node->children_end(); ++childItr )
			{
				Q.push_back( distance_node_type(distance(childItr->point(), p), boost::addressof(*childItr)) );
			}
		}

		if( !Q.empty() )
		{
			// form the next cover set:
			//   Q_{i-1} = { q \ in Q : Bound(Q) < d(p, Q) + 2^i }
			real_type minDist = bound(Q) + scale(level);
			less_than_predicate<distance_node_type> pred( minDist );
			Q.erase( std::partition( Q.begin(), Q.end(), pred ), Q.end() );
		}

		std::copy( Q_i.begin(), Q_i.end(), 
			std::back_inserter( allCoverSets ) );
		std::swap( Q_i, Q );
	}

	return select( allCoverSets );
}

template <typename Point, typename Distance>
typename cover_tree<Point, Distance>::real_type 
	cover_tree<Point, Distance>::scale( typename cover_tree<Point, Distance>::level_type level ) const
{
	return boost::numeric_cast<real_type>( ldexp( 1.0, level ) );
}

template <typename Point, typename Distance>
typename cover_tree<Point, Distance>::node_type* 
	cover_tree<Point, Distance>::insert( 
		const Point& p, 
		const std::deque<node_type*>& Q_i, 
		typename cover_tree<Point, Distance>::level_type level )
{
	real_type scale = this->scale( level );

	bool parentFound = false;
	for( typename std::deque<node_type*>::const_iterator Q_i_iter = Q_i.begin();
		Q_i_iter != Q_i.end(); ++Q_i_iter )
	{
		if( distance( (*Q_i_iter)->point(), p ) <= scale )
		{
			parentFound = true;
			break;
		}
	}

	if( !parentFound )
		return 0;

	std::deque< node_type* > Q_i_minus_1;
	// set         Q = {Children(Q) : q \in Q_i}
	//       Q_{i-1} = {q \in Q : d(p, q) < 2^i}
	for( typename std::deque<node_type*>::const_iterator Q_i_iter = Q_i.begin();
		Q_i_iter != Q_i.end(); ++Q_i_iter )
	{
		for( typename node_type::child_iterator Q_iter = (*Q_i_iter)->children_begin();
			Q_iter != (*Q_i_iter)->children_end(); ++Q_iter )
		{
			if( distance(Q_iter->point(), p) < scale )
				Q_i_minus_1.push_back( boost::addressof(*Q_iter) );
		}
	}

	// if Insert(p, Q_{i-1}, i-1) = "no parent found"
	node_type* result = insert( p, Q_i_minus_1, level - 1 );
	if( result == 0 )
	{
		for( typename std::deque<node_type*>::const_iterator Q_i_iter = Q_i.begin();
			Q_i_iter != Q_i.end(); ++Q_i_iter )
		{
			// if d(p, Q_i) <= 2^i
			// pick q \in Q_i satisfying d(p, q) < 2^i
			if( distance(p, (*Q_i_iter)->point()) <= scale )
			{
				node_type* child = nodePool_.construct(
					p, level-1 );

				// insert p into Children(q)
				(*Q_i_iter)->add_child( child );

				// return "parent found"
				return *Q_i_iter;
			}
		}
	}

	return result;       // else return "no parent found"
}

template <typename Point, typename Distance>
void cover_tree<Point, Distance>::insert( node_type* node )
{
	insert( node->point() );

	for( typename node_type::child_iterator iter = node->children_begin();
		iter != node->children_end(); ++iter )
	{
		insert( boost::addressof(*iter) );
	}
}

template <typename Point, typename Distance>
void cover_tree<Point, Distance>::check_child_invariants( 
	const typename cover_tree<Point, Distance>::node_type* node ) const
{
	// (cover tree invariant)
	// (2) For every p \in C_{i-1}, there exists a q \in C_i satisfying d(p, q) <= 2^i, and exactly one
	//        such q is a parent of p
	for( typename node_type::const_child_iterator childIter = node->children_begin();
		childIter != node->children_end(); ++childIter )
	{
		assert( childIter->level() == (node->level() - 1) );
		assert( distance( childIter->point(), node->point() ) <= scale( node->level() ) );
	}
}

template <typename Point, typename Distance>
void cover_tree<Point, Distance>::check_invariants() const
{
	check_child_invariants( root_ );

	// now, do a bfs and check that the separation invariant holds
	typedef std::deque<const node_type*> NodeQueue;
	NodeQueue nodes;
	nodes.push_back( root_ );

	while( !nodes.empty() )
	{
		const node_type* currentNode = nodes.front();
		nodes.pop_front();

		// every time we see a node, we compare to all the other nodes
		//   in its level, this will in effect perform the all pairs check
		typedef typename NodeQueue::iterator NodeItr;
		for( NodeItr iter = nodes.begin(); 
			iter != nodes.end() && (*iter)->level() == currentNode->level();
			++iter )
		{
			const node_type* otherNode = *iter;
			assert( distance( currentNode->point(), otherNode->point() ) > scale( currentNode->level() ) );
		}

		for( typename node_type::const_child_iterator childIter = currentNode->children_begin();
			childIter != currentNode->children_end(); ++childIter )
		{
			nodes.push_back( boost::addressof(*childIter) );
		}
	}
}

