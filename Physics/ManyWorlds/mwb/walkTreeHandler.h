#include "planningUI.h"

#ifndef __WALKTREEHANDLER_H__
#define __WALKTREEHANDLER_H__

#include <algorithm>

namespace planning
{


struct ElementToHandle
{
	ElementToHandle( SceneGraphElementPtr elt,
		const vl::Mat4f& trans, UsefulBits& bits )
		: object(elt), transform(trans), bits(bits)
	{
	}

	SceneGraphElementPtr object;
	vl::Mat4f transform;
	UsefulBits bits;
};


template <typename Handler>
void MainFrame::walkTree( Handler& handler, const Simulation::State& state ) const
{
	std::vector<ConstPhysicsObjectPtr> dynamicObjects;
	if( this->simulation_ )
		dynamicObjects = simulation_->dynamicObjects();
	else if( this->currentTree_ )
		dynamicObjects = currentTree_->dynamicObjects();

	typedef STLEXT hash_map< std::string, size_t > NameToIdMap;
	NameToIdMap dynamicsObjectMap;
	for( size_t i = 0; i < dynamicObjects.size(); ++i )
	{
		std::pair<NameToIdMap::iterator, bool> result = 
			dynamicsObjectMap.insert( 
				NameToIdMap::value_type( dynamicObjects.at(i)->name(), i ) );
		assert( result.second );
	}

	typedef STLEXT hash_set< const SceneGraphElement* > SelectedSet;
	SelectedSet selectedSet;
	SelectedSet involvedInJointSet;
	{
		std::deque<SceneGraphElementPtr> selected = this->getSelected();
		for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
			selectedItr != selected.end(); ++selectedItr )
		{
			selectedSet.insert( selectedItr->get() );
		}

		for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
			selectedItr != selected.end(); ++selectedItr )
		{
			JointPtr joint = boost::dynamic_pointer_cast<const Joint>( *selectedItr );
			if( joint )
			{
				std::vector<ConstPhysicsObjectPtr> objects = joint->objects();
				for( std::vector<ConstPhysicsObjectPtr>::const_iterator objectItr = objects.begin();
					objectItr != objects.end(); ++objectItr )
				{
					involvedInJointSet.insert( objectItr->get() );
				}
			}
		}
	}

	// we will use this for joints:
	typedef STLEXT hash_map<const PhysicsObject*, const char*> ObjectToStateMap;
	ObjectToStateMap objectToState;

	std::deque<ElementToHandle> toHandle( 1, 
		ElementToHandle(
			this->scene_->root(),
			this->scene_->root()->transform(),
			UsefulBits() ) );
	while( !toHandle.empty() )
	{
		ElementToHandle current = toHandle.back();
		toHandle.pop_back();

		boost::shared_ptr<const Joint> joint = 
			boost::dynamic_pointer_cast<const Joint>( current.object );
		ConstPhysicsObjectPtr physicsObject = 
			boost::dynamic_pointer_cast<const PhysicsObject>( current.object );

		UsefulBits currentBits = current.bits;
		if( !currentBits.changedFromTree() )
		{
			currentBits.setChangedFromTree( isModifiedFromTree( current.object ) );
		}

		if( !currentBits.selected() )
		{
			SelectedSet::const_iterator selectedItr = 
				selectedSet.find( current.object.get() );
			currentBits.setSelected(selectedItr != selectedSet.end());
			currentBits.setCurrentSelected(selectedItr != selectedSet.end());
		}
		else
		{
			currentBits.setCurrentSelected(false);
		}

		if( !currentBits.involvedInConstraint() && this->selectedConstraint_ )
		{
			for( Constraint::object_iterator objectItr = selectedConstraint_->objects_begin();
				objectItr != selectedConstraint_->objects_end(); ++objectItr )
			{
				if( *objectItr < dynamicObjects.size() && 
					current.object->name() == dynamicObjects.at( *objectItr )->name() )
					currentBits.setInvolvedInConstraint();
			}
		}

		vl::Mat4f currentTransform;
		if( joint && !currentBits.changedFromTree() )
		{
			// joints get treated specially
			std::vector<ConstPhysicsObjectPtr> objects = joint->objects();
			std::vector<const char*> stateVec( objects.size(), 0 );

			typedef const char* StatePtr;
			bool hasState = false;
			for( size_t iObject = 0; iObject < objects.size(); ++iObject )
			{
				ObjectToStateMap::const_iterator itr = 
					objectToState.find( objects.at(iObject).get() );
				if( itr != objectToState.end() )
				{
					stateVec[iObject] = itr->second;
					hasState = true;
				}
			}

			if( hasState )
			{
				const char* statePtr = reinterpret_cast<const char*>(&stateVec[0]);
				handler( current.object, vl::vl_1, currentBits, statePtr );
			}
			else
			{
                currentTransform = 
					current.transform * joint->localTransform( 0 );
				handler( current.object, currentTransform, currentBits, 0 );
			}
		}
		else if( physicsObject && !currentBits.changedFromTree() )
		{
			if( !currentBits.involvedInJoint() )
			{
				SelectedSet::const_iterator involvedInJointItr = 
					involvedInJointSet.find( current.object.get() );
				currentBits.setInvolvedInJoint(
					involvedInJointItr != involvedInJointSet.end());
			}

			NameToIdMap::const_iterator mapItr = 
				dynamicsObjectMap.find( physicsObject->name() );
			if( mapItr != dynamicsObjectMap.end() && mapItr->second < state.stateCount() )
			{
				const char* statePtr = state.extendedState( mapItr->second );
				objectToState.insert( ObjectToStateMap::value_type(physicsObject.get(), statePtr) );
				currentTransform = physicsObject->localTransform( statePtr );
				handler( current.object, currentTransform, currentBits, statePtr );
			}
			else
			{
				currentTransform = 
					current.transform * physicsObject->localTransform( 0 );
				handler( current.object, currentTransform, currentBits, 0 );
			}
		}
		else
		{
			currentTransform = 
				current.transform * current.object->localTransform( 0 );
			handler( current.object, currentTransform, currentBits, 0 );
		}
	

		size_t stackSize = toHandle.size();
		for( SceneGraphElement::child_iterator itr = current.object->children_begin();
			itr != current.object->children_end(); ++itr )
		{
			toHandle.push_back( ElementToHandle( 
				*itr, currentTransform, currentBits ) );
		}
		// want to evaluate in the original order
		std::reverse( toHandle.begin() + stackSize, toHandle.end() );
	}
}

} // namespace planning

#endif
