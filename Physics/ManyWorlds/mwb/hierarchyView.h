#ifndef __HIERARCHYVIEW_H__
#define __HIERARCHYVIEW_H__

#include "physicsFwd.h"

#include <wx/treectrl.h>

#include <hash_map>

namespace planning
{


class HierarchyControl;

class HierarchyFrame
	: public wxFrame
{
public:
	HierarchyFrame(wxFrame *frame, const wxString& title, const wxPoint& pos, const wxSize& size,
        long style = wxDEFAULT_FRAME_STYLE);
	~HierarchyFrame();

	void setScene( boost::shared_ptr<Scene> scene, const std::deque<SceneGraphElementPtr>& selected );
	void setSelected( const std::deque<SceneGraphElementPtr>& selected, bool force = true );

	void OnClose(wxCloseEvent& event);
	void OnChangeSelection(wxTreeEvent& event);
	void OnKey(wxKeyEvent& event);

	void OnBeginDrag(wxTreeEvent& event);
	void OnEndDrag(wxTreeEvent& event);

	void OnStartLabelEdit(wxTreeEvent& event);
	void OnEndLabelEdit(wxTreeEvent& event);

private:
	void buildHierarchy( wxTreeItemId treeId, SceneGraphElementPtr elt );

	void updateGray();
	bool updateGray( wxTreeItemId id );

	struct SceneGraphElementWrapper
		: public wxTreeItemData
	{
		SceneGraphElementWrapper(SceneGraphElementPtr elt);
		SceneGraphElementWrapper();
		~SceneGraphElementWrapper();

        SceneGraphElementPtr element;
	};

	enum 
	{
		TREE,
	};

	HierarchyControl* tree_;
	boost::shared_ptr<Scene> scene_;

	typedef STLEXT hash_map< const SceneGraphElement*, wxTreeItemId, 
		hash_ptr<const SceneGraphElement*> > 
		SceneGraphElementToTreeItemMap;
	SceneGraphElementToTreeItemMap eltToItem_;

	wxColour defaultBackgroundColor_;
	wxColour childIsSelectedBackgroundColor_;

	std::deque<SceneGraphElementPtr> currentSelected_;

	wxImageList imageList_;

	DECLARE_EVENT_TABLE()
};

class HierarchyControl
	: public wxTreeCtrl
{
public:
	HierarchyControl( wxWindow* parent, wxWindowID id, 
		const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, 
		long style = wxTR_HAS_BUTTONS, const wxValidator& validator = wxDefaultValidator, 
		const wxString& name = "listCtrl" );
	~HierarchyControl();

	void OnKey( wxKeyEvent& event );

private:
	DECLARE_EVENT_TABLE()
};

} // namespace planning

#endif
