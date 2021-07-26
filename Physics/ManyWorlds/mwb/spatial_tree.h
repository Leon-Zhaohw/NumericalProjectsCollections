#ifdef WIN32
#pragma once
#endif

namespace planning {

// a spatial_tree represents a quadtree/octree/higher-d version
// for now, we will assume that the contained data can be stored in 
//   hash sets.  later, if necessary, we will relax this constraint
template <typename real, size_t num_dimensions, typename data_type>
class spatial_tree
{
public:
	typedef STLEXT hash_set<data_type> data_set;

	spatial_tree( size_t maxDepth )
		: maxDepth_( maxDepth ), nodePool_( new boost::object_pool<spatial_tree_node>() )
	{
		root_ = nodePool_->construct();
	}

	template <typename sequence>
	void insert( const sequence& position, const data_type& data )
	{
		// first, verify that the position is in the set
		for( size_t i = 0; i < num_dimensions; ++i )
		{
			if( position[i] < 0.0 || position[i] > 1.0 )
				return;
		}

		root_->insert( position, data, maxDepth_, *nodePool_ );
	}

	// includePartial specifies whether to include objects that are in boxes that are only
	//   partially included in the set
	template <typename sequence>
	data_set query( const sequence& lowerPosition, const sequence& higherPosition, bool includePartial ) const
	{
		data_set result;
		root_->query( lowerPosition, higherPosition, includePartial, result, this->max_depth() );
		return result;
	}

	void clear()
	{
		nodePool_.reset( new boost::object_pool<spatial_tree_node>() );
		root_ = nodePool_->construct();
	}

	size_t max_depth() const
	{
		return this->maxDepth_;
	}

private:
	typedef boost::array<real, num_dimensions> position_type;

	class spatial_tree_node
	{
		BOOST_STATIC_CONSTANT( size_t, numChildren = (2 << num_dimensions) );

	public:
		spatial_tree_node()
		{
			for( size_t i = 0; i < numChildren; ++i )
				children_[i] = 0;

			//std::fill( children_.begin(), children_.end(), 0 );
		}

		template <typename sequence>
		void query( const sequence& lowerPosition, 
			const sequence& higherPosition, 
			bool includePartial, 
			data_set& result, 
			size_t remainingDepth )
		{
			bool contained = true;
			bool intersects = true;
			for( size_t i = 0; i < num_dimensions; ++i )
			{
				if( lowerPosition[i] > 0.0 || higherPosition[i] < 1.0 )
					contained = false;
				if( lowerPosition[i] > 1.0 || higherPosition[i] < 0.0 )
					intersects = false;
			}

			if( !intersects )
				return;

			if( contained )
			{
				result.insert( data_.begin(), data_.end() );
				return;
			}

			// recurse on children
			bool hasChildren = false;
			for( size_t i = 0; i < children_array::static_size; ++i )
			{
				if( children_[i] == 0 )
					continue;

				hasChildren = true;
				boost::array<real, num_dimensions> childLower;
				boost::array<real, num_dimensions> childUpper;
				for( size_t jDim = 0; jDim < num_dimensions; ++jDim )
				{
					if( (i >> (num_dimensions - jDim - 1)) & 1 )
					{
						childLower[jDim] = 2.0*(lowerPosition[jDim] - 0.5);
						childUpper[jDim] = 2.0*(higherPosition[jDim] - 0.5);
					}
					else
					{
						childLower[jDim] = 2.0*(lowerPosition[jDim]);
						childUpper[jDim] = 2.0*(higherPosition[jDim]);
					}
				}

				children_[i]->query( childLower, childUpper, includePartial, result, remainingDepth-1 );
			}

			if( !hasChildren && includePartial )
			{
				assert( remainingDepth == 0 || data_.empty() );
				result.insert( data_.begin(), data_.end() );
			}
		}

		template <typename sequence>
		void insert( const sequence& position, const data_type& data, size_t remainingDepth, 
			boost::object_pool<spatial_tree_node>& nodePool )
		{
			data_.insert( data );
			if( remainingDepth == 0 )
				return;

			typename children_array::size_type whichChild = 0;
			boost::array<real, num_dimensions> childPos;
			for( size_t i = 0; i < num_dimensions; ++i )
			{
				whichChild = whichChild << 1;

				if( position[i] > 0.5 )
				{
					whichChild |= 1;
					childPos[i] = 2.0*(position[i] - 0.5);
				}
				else
				{
					whichChild |= 0;
					childPos[i] = 2.0*position[i];
				}
			}

			if( children_[whichChild] == 0 )
				children_[whichChild] = nodePool.construct();

			children_[whichChild]->insert( 
				childPos,
				data, 
				remainingDepth - 1,
				nodePool );
		}

	private:
		typedef boost::array<spatial_tree_node*, numChildren> children_array;
		children_array children_;
		data_set data_;
	};

	size_t maxDepth_;
	spatial_tree_node* root_;
	boost::scoped_ptr< boost::object_pool<spatial_tree_node> > nodePool_;
};

void test_spatial_tree();

} // namespace planning

