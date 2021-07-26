#include "stdafx.h"

#include "planningUI.h"
#include "hierarchyView.h"
#include "physicsObject.h"
#include "scene.h"

#include "bitmaps/jointImage.xpm"
#include "bitmaps/cube.xpm"
#include "bitmaps/fusedCube.xpm"
#include "bitmaps/group.xpm"

namespace planning
{

BEGIN_EVENT_TABLE(HierarchyFrame, wxFrame)
	EVT_CLOSE(HierarchyFrame::OnClose)
	EVT_TREE_SEL_CHANGED(HierarchyFrame::TREE, HierarchyFrame::OnChangeSelection)
	EVT_KEY_DOWN(HierarchyFrame::OnKey)
	EVT_TREE_BEGIN_RDRAG(HierarchyFrame::TREE, HierarchyFrame::OnBeginDrag)
	EVT_TREE_END_DRAG(HierarchyFrame::TREE, HierarchyFrame::OnEndDrag)
	EVT_TREE_BEGIN_LABEL_EDIT(HierarchyFrame::TREE, HierarchyFrame::OnStartLabelEdit)
	EVT_TREE_END_LABEL_EDIT(HierarchyFrame::TREE, HierarchyFrame::OnEndLabelEdit)
END_EVENT_TABLE()

HierarchyFrame::HierarchyFrame(wxFrame *frame, const wxString& title, const wxPoint& pos, const wxSize& size,
        long style)
	: wxFrame(frame, -1, title, pos, size, style), imageList_(20, 20)
{
	tree_ = new HierarchyControl( this, TREE, wxDefaultPosition, wxDefaultSize, 
		wxTR_MULTIPLE | wxTR_EXTENDED | wxTR_HAS_BUTTONS | 
		wxTR_ROW_LINES | wxTR_HIDE_ROOT | 
		wxTR_LINES_AT_ROOT | wxTR_EDIT_LABELS );
	tree_->AddRoot( wxT("root") );

	this->childIsSelectedBackgroundColor_ = wxColour( 200, 200, 200 );
	this->defaultBackgroundColor_ = tree_->GetItemBackgroundColour( tree_->GetRootItem() );
	const wxColour childIsSelectedColour( 200, 200, 200 );

	wxBitmap groupImage(group_xpm);
	imageList_.Add( groupImage );

	wxBitmap cubeImage(cube_xpm);
	imageList_.Add( cubeImage );

	wxBitmap fusedCubeImage(fusedCube_xpm);
	imageList_.Add( fusedCubeImage );

	wxBitmap jointImage(jointImage_xpm);
	imageList_.Add( jointImage );

	tree_->SetImageList( &imageList_ );

}

HierarchyFrame::~HierarchyFrame()
{
}

HierarchyFrame::SceneGraphElementWrapper::SceneGraphElementWrapper(SceneGraphElementPtr elt)
	: element(elt)
{
}

void HierarchyFrame::OnKey(wxKeyEvent& event)
{
	tree_->OnKey(event);
}

HierarchyFrame::SceneGraphElementWrapper::SceneGraphElementWrapper()
{
}

HierarchyFrame::SceneGraphElementWrapper::~SceneGraphElementWrapper()
{
}

void HierarchyFrame::OnClose(wxCloseEvent& event)
{
	if( event.CanVeto() )
	{
		this->Hide();
		event.Veto();
	}
	else
	{
		this->Destroy();
	}
}

void HierarchyFrame::buildHierarchy( wxTreeItemId treeId, SceneGraphElementPtr elt )
{
	eltToItem_[elt.get()] = treeId;

	wxTreeItemIdValue cookie;
	wxTreeItemId currentTreeChild = tree_->GetFirstChild( treeId, cookie );

	for( SceneGraphElement::child_iterator currentSceneChild = elt->children_begin();
		currentSceneChild != elt->children_end(); ++currentSceneChild )
	{
		wxTreeItemId myId;
		wxTreeItemId tempTreeChild = currentTreeChild;
		for( size_t numChecked = 0; numChecked < 3 && tempTreeChild.IsOk(); 
			++numChecked, tempTreeChild = tree_->GetNextSibling(tempTreeChild) )
		{
			wxTreeItemData* itemData = tree_->GetItemData(tempTreeChild);
			const SceneGraphElementWrapper* eltWrapper = 
				dynamic_cast<const SceneGraphElementWrapper*>( itemData );
			assert( eltWrapper != 0 );
			if( eltWrapper->element == *currentSceneChild )
			{
				// check if name matches
				std::string currentItemText( tree_->GetItemText(tempTreeChild).c_str() );
				if( eltWrapper->element->name() != currentItemText )
					tree_->SetItemText( tempTreeChild, wxT(eltWrapper->element->name().c_str()) );

				// clear out the dead items
				for( size_t i = 0; i < numChecked; ++i )
					tree_->Delete( tree_->GetPrevSibling(tempTreeChild) );

				currentTreeChild = tree_->GetNextSibling(tempTreeChild);
				myId = tempTreeChild;
				break;
			}
		}

		if( !myId.IsOk() )
		{
			int image = 0;
			if( dynamic_cast<const Joint*>( currentSceneChild->get() ) )
				image = 3;
			else if( dynamic_cast<const CombinedObject*>( currentSceneChild->get() ) )
				image = 2;
			else if( dynamic_cast<const PhysicsObject*>( currentSceneChild->get() ) )
				image = 1;

			// doesn't exist, need to create new entry
			// insert before current spot:
			// this is gross because wxWidgets doesn't include an "insert before current"
			//   operator
			if( !currentTreeChild.IsOk() )
			{
				myId = tree_->AppendItem( 
						treeId, (*currentSceneChild)->name(), image, -1, 
						new SceneGraphElementWrapper(*currentSceneChild) );
			}
			else
			{
				wxTreeItemId prev = tree_->GetPrevSibling( currentTreeChild );
				if( prev.IsOk() )
				{
					myId = tree_->InsertItem( 
						treeId, prev, (*currentSceneChild)->name(), image, -1, 
						new SceneGraphElementWrapper(*currentSceneChild) );
				}
				else
				{
					myId = tree_->InsertItem( 
						treeId, (size_t) 0, (*currentSceneChild)->name(), image, -1, 
						new SceneGraphElementWrapper(*currentSceneChild) );
				}
			}
		}

		this->buildHierarchy( myId, *currentSceneChild );
	}

	// clear out any other items that might be there
	if( currentTreeChild.IsOk() )
	{
		while( tree_->GetNextSibling( currentTreeChild ).IsOk() )
			tree_->Delete( tree_->GetNextSibling( currentTreeChild ) );

		tree_->Delete( currentTreeChild );
	}
}

void HierarchyFrame::setScene( boost::shared_ptr<Scene> scene, 
	const std::deque<SceneGraphElementPtr>& selected )
{
	this->buildHierarchy( tree_->GetRootItem(), scene->root() );
	this->setSelected( selected );
}

class Freezer
{
public:
	Freezer( wxWindow* window )
		: window_(window)
	{
		window_->Freeze();
	}

	~Freezer()
	{
		window_->Thaw();
	}

private:
	wxWindow* window_;
};

// sometimes we will need to force setting the selection because wxWidgets insists on changing the
//   selection without notifying us via a CHANGE_SELECTION message
void HierarchyFrame::setSelected( const std::deque<SceneGraphElementPtr>& selected, bool force )
{
	Freezer freezer( this->tree_ );

	if( currentSelected_ == selected && !force )
		return;

	tree_->UnselectAll();

	for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
		selectedItr != selected.end(); ++selectedItr )
	{
		SceneGraphElementToTreeItemMap::const_iterator mapItr = 
			eltToItem_.find( selectedItr->get() );
		// this will happen when we e.g. add an object to the scene, and
		//   the selection gets changed before we are notified about
		//   the scene changing.
		//   The easiest way is to simply ignore it because it will get 
		//   fixed once the action is complete
		if( mapItr == eltToItem_.end() )
			continue;

		tree_->SelectItem( mapItr->second );
	}

	this->updateGray();
	currentSelected_ = selected;
}

void HierarchyFrame::updateGray()
{
	updateGray( tree_->GetRootItem() );
}

bool equalColors( const wxColour& left, const wxColour& right )
{
	return left.Blue() == right.Blue()
		&& left.Red() == right.Red()
		&& left.Green() == right.Green();
}

bool HierarchyFrame::updateGray( wxTreeItemId id )
{
	bool childSelected = false;
	wxTreeItemIdValue cookie;
	for( wxTreeItemId childItem = tree_->GetFirstChild(id, cookie);
		childItem.IsOk(); childItem = tree_->GetNextChild(id, cookie) )
	{
		childSelected = this->updateGray( childItem ) || childSelected;
	}

	if( id == tree_->GetRootItem() )
		return childSelected;

	if( childSelected && !equalColors(tree_->GetItemBackgroundColour(id), childIsSelectedBackgroundColor_) )
		tree_->SetItemBackgroundColour(id, childIsSelectedBackgroundColor_);
	else if( !childSelected && !equalColors(tree_->GetItemBackgroundColour(id), defaultBackgroundColor_) )
		tree_->SetItemBackgroundColour(id, defaultBackgroundColor_);

	return childSelected || tree_->IsSelected( id );
}

void HierarchyFrame::OnChangeSelection(wxTreeEvent& event)
{
	// hack hack hack
	// for some reason, it insists on sending us TWO events every
	//  time the user clicks something: one to disable the previous
	//  selection, and one to enable the new selection.  We will 
	//  alsway ignore the first and only respond to the second.
	if( !event.GetItem().IsOk() )
	{
		event.Skip();
		return;
	}

	wxTreeItemId current = tree_->GetItemParent( event.GetItem() );
	while( current.IsOk() && current != tree_->GetRootItem() )
	{
		if( tree_->IsSelected( current ) )
			tree_->UnselectItem( current );
		current = tree_->GetItemParent( current );
	}

	Freezer freezer( this->tree_ );

	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	assert( parent );

	std::deque<SceneGraphElementPtr> newSelection;

	wxArrayTreeItemIds selection;
	size_t numSel = tree_->GetSelections(selection);
	for( size_t i = 0; i < selection.GetCount(); ++i )
	{
		wxTreeItemId id = selection.Item( i );
		if( id == tree_->GetRootItem() )
			continue;

		wxTreeItemData* itemData = tree_->GetItemData(id);
		const SceneGraphElementWrapper* eltWrapper = 
			dynamic_cast<const SceneGraphElementWrapper*>( itemData );
		assert( eltWrapper != 0 );
        newSelection.push_back( eltWrapper->element );
	}

	if( newSelection == this->currentSelected_ )
		return;

	this->updateGray();
	this->Update();

	this->currentSelected_ = newSelection;
	ActionPtr changeSelectionAction( 
		new ChangeSelectionAction(
			parent,
			parent->getSelected(),
			cleanSelected(newSelection),
			ConstraintPtr(), ConstraintPtr() ) );
	parent->performAction( changeSelectionAction );
}

void HierarchyFrame::OnBeginDrag(wxTreeEvent& event)
{
	event.Allow();	
}

void HierarchyFrame::OnEndDrag(wxTreeEvent& event)
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	wxPoint point = event.GetPoint();
	wxTreeItemId id = event.GetItem();//tree_->HitTest( point, flags );
	SceneGraphElementPtr newParent;
	if( id.IsOk() )
	{
		wxTreeItemData* itemData = tree_->GetItemData(id);
		const SceneGraphElementWrapper* eltWrapper = 
			dynamic_cast<const SceneGraphElementWrapper*>( itemData );
		assert( eltWrapper != 0 );
		newParent = eltWrapper->element;
	}
	else
	{
		newParent = parent->scene()->root();
	}

	// need to check if it's okay to drop here:
	// not okay to drop a parent on its own child
	std::deque<SceneGraphElementPtr> selected = parent->getSelected();
	for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
		selectedItr != selected.end(); ++selectedItr )
	{
		SceneGraphElementPtr current = newParent;
		while( current )
		{
			if( current == *selectedItr )
			{
				event.Veto();
				return;
			}

			current = current->parent().lock();
		}
	}

	try
	{
		std::deque<ActionPtr> actions = parent->reparent(
			selected,
			std::deque<SceneGraphElementPtr>( selected.size(), newParent ),
			true );
		boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
		parent->performAction( compoundAction );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error reparenting: " << e.message();
		errMsg( oss.str(), this );	
	}
}

void HierarchyFrame::OnStartLabelEdit(wxTreeEvent& event)
{
	event.Allow();
}

void HierarchyFrame::OnEndLabelEdit(wxTreeEvent& event)
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	std::string newName(event.GetLabel().c_str());

	ScenePtr scene = parent->scene();

	wxTreeItemId id = event.GetItem();
	if( id == tree_->GetRootItem() )
	{
		event.Veto();
		return;
	}

	wxTreeItemData* itemData = tree_->GetItemData(id);
	const SceneGraphElementWrapper* eltWrapper = 
		dynamic_cast<const SceneGraphElementWrapper*>( itemData );
	assert( eltWrapper != 0 );

	if( newName.empty() || eltWrapper->element->name() == newName )
	{
		// just return to old name.
		event.Veto();
		return;
	}
	else if( scene->object( newName ) )
	{
		event.Veto();
		return;
	}
	else
	{
		// new name is okay
		boost::shared_ptr<Action> nameChangeAction(
			new NameChangeAction( eltWrapper->element, scene, eltWrapper->element->name(), newName ) );
		parent->performAction( nameChangeAction );
		event.Allow();
		return;
	}
}

HierarchyControl::HierarchyControl( wxWindow* parent, wxWindowID id, 
	const wxPoint& pos, const wxSize& size, 
	long style, const wxValidator& validator, 
	const wxString& name )
	: wxTreeCtrl( parent, id, pos, size, style, validator, name )
{
}

HierarchyControl::~HierarchyControl()
{
}

BEGIN_EVENT_TABLE(HierarchyControl, wxTreeCtrl)
	EVT_KEY_DOWN(HierarchyControl::OnKey)
END_EVENT_TABLE()

void HierarchyControl::OnKey( wxKeyEvent& event )
{
	MainFrame* frame = dynamic_cast<MainFrame*>(this->GetParent()->GetParent());
	HierarchyFrame* parent = dynamic_cast<HierarchyFrame*>(this->GetParent());
		if( event.ControlDown() )
	{
		switch( event.GetKeyCode() )
		{
		case 'a':
		case 'A':
			frame->showProperties( !frame->arePropertiesShown() );
			break;
		case 'd':
		case 'D':
			frame->duplicate( frame->getSelected() );
			break;
		case 'C':
		case 'c':
			frame->copy( frame->getSelected() );
			break;
		case 'v':
		case 'V':
			frame->paste();
			break;
		case 'z':
		case 'Z':
			frame->undo();
			break;
		case 'y':
		case 'Y':
			frame->redo();
			break;
		default:
			event.Skip();
			return;
		}
	}
	else if( event.AltDown() )
	{
		event.Skip();
		return;
	}
	else
	{
		switch( event.GetKeyCode() )
		{
		case WXK_UP:
			frame->upHierarchy();
			parent->setSelected( frame->getSelected(), true );
			break;
		case WXK_DOWN:
			frame->downHierarchy();
			parent->setSelected( frame->getSelected(), true );
			break;
		case WXK_DELETE:
		{
			std::deque<SceneGraphElementPtr> selected = frame->getSelected();
			if( !selected.empty() )
				frame->deleteObjects( selected );
			break;
		}
		default:
			event.Skip();
			return;
		}
	}
}

} // namespace planning
