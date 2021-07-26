#include "stdafx.h"

#include "GLView.h"
#include "planningUI.h"
#include "walkTreeHandler.h"
#include "constraints.h"
#include "sphere_tree.h"

#include "twigg/renderUtils.h"

#include "Wm4Ray3.h"
#include "Wm4DistVector2Line2.h"
#include "Wm4DistLine3Line3.h"
#include "Wm4DistLine3Circle3.h"
#include "Wm4DistVector3Line3.h"

// used by the inertia tensor stuff
#include <mkl_lapack.h>

#include <wx/image.h>

#include <iostream>
#include <iomanip>

namespace planning {

BEGIN_EVENT_TABLE(GLView, wxGLCanvas)
    EVT_SIZE(GLView::OnSize)
    EVT_PAINT(GLView::OnPaint)
    EVT_ERASE_BACKGROUND(GLView::OnEraseBackground)
    EVT_MOUSE_EVENTS(GLView::OnMouse)
	EVT_KEY_DOWN(GLView::OnKey)
	EVT_KEY_UP(GLView::OnKeyUp)
END_EVENT_TABLE()

GLView::GLView(wxWindow *parent, wxWindowID id,
    const wxPoint& pos, const wxSize& size, long style, const wxString& name)
	:	wxGLCanvas(parent, id, pos, size, style, name),
		frameTime_(0.0),
		headsUp_(true),
		selecting_(false)
{
	spacePressed_ = false;

	viewportWidth_ = size.GetWidth();
	viewportHeight_ = size.GetHeight();
}

GLView::GLView(wxWindow *parent, wxGLContext* sharedContext, wxWindowID id,
    const wxPoint& pos, const wxSize& size, long style, const wxString& name)
	:	wxGLCanvas(parent, sharedContext, id, pos, size, style, name),
		frameTime_(0.0),
		headsUp_(true),
		selecting_(false)
{
	spacePressed_ = false;

	viewportWidth_ = size.GetWidth();
	viewportHeight_ = size.GetHeight();
}


GLView::~GLView()
{
}

template <typename T>
struct identity
{
	T operator()( const T& t ) const
	{
		return t;
	}
};

/*
void GLView::splatNodes(SimulationTreePtr tree, const std::deque<SimulationTree::Path>& nodes )
{
	// take out the lock before we grab the matrix so we know it's valid
	boost::try_mutex::scoped_lock lock( quadTreeMutex_ );
	vl::Mat4f matrix = this->camera_->camera()->projectiveTransform() 
		* this->camera_->camera()->modelViewTransform();

	std::vector<size_t> dynamicObjectIds = 
		objectListToDynamicObjectIds( selected_, tree->dynamicObjects() );

	for( std::deque<SimulationTree::Path>::const_iterator nodeIter = nodes.begin();
		nodeIter != nodes.end(); ++nodeIter )
	{
		const SimulationNode* origNode = nodeIter->tail;
		const SimulationNode* node = nodeIter->tail;

		while( node != 0 )
		{
			for( size_t iObject = 0; iObject < dynamicObjectIds.size(); ++iObject )
			{
				if( this->quadTreeUpdateStop_ )
					return;

				size_t objectId = dynamicObjectIds[iObject];
				const vl::Vec3f position = toVec3f(node->state().state( iObject ).position());
				const vl::Vec3f transformedPosition = 
					vl::xform( matrix, position );
				boost::array<float, 2> normalizedPos = 
					{{ 0.5*( transformedPosition[0] + 1.0 ), 0.5*( transformedPosition[1] + 1.0 ) }};

				this->quadTree_.insert( normalizedPos, origNode );
			}

			node = node->parent();
		}
	}
}
*/

void GLView::setCamera( CameraWrapperPtr camera )
{
	this->camera_ = camera;
	this->Refresh(FALSE);
}

CameraWrapperPtr GLView::cameraWrapper()
{
	return this->camera_;
}

Camera& GLView::camera()
{
	assert( camera_ );
	return *this->camera_->camera();
}

const Camera& GLView::camera() const
{
	assert( camera_ );
	return *this->camera_->camera();
}

void GLView::errMsg(const std::string& message)
{
	planning::errMsg( message, this );
}

void GLView::clearCache()
{
	this->Refresh( FALSE );
}

void GLView::OnKeyUp(wxKeyEvent& event)
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	time_t keyReleaseTime = clock();

	if( event.GetKeyCode() == WXK_SPACE )
	{
		if( (keyReleaseTime - keyPressTime_) < 300 )
		{
			parent->setCurrentView( this, true );
		}

		this->spacePressed_ = false;
		this->Refresh(FALSE);
		return;
	}

	event.Skip();
}

void GLView::OnKey(wxKeyEvent& event)
{
	this->SetFocus();

	keyPressTime_ = clock();
	if( event.GetKeyCode() == WXK_SPACE )
	{
		this->spacePressed_ = true;
		this->Refresh(FALSE);
		return;
	}

	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	if( event.ControlDown() )
	{
		switch( event.GetKeyCode() )
		{
		case 'a':
		case 'A':
			parent->showProperties( !parent->arePropertiesShown() );
			break;
		case 'd':
		case 'D':
			parent->duplicate( parent->getSelected() );
			break;
		case 'C':
		case 'c':
			parent->copy( parent->getSelected() );
			break;
		case 'v':
		case 'V':
			parent->paste();
			break;
		case 'z':
		case 'Z':
			parent->undo();
			break;
		case 'y':
		case 'Y':
			parent->redo();
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
		case WXK_DELETE:
		{
			std::deque<SceneGraphElementPtr> selected = parent->getSelected();
			if( !selected.empty() )
				parent->deleteObjects( selected );
			else if( parent->selectedConstraint() )
				parent->deleteConstraint( parent->selectedConstraint() );
			break;
		}
		case WXK_UP:
			parent->upHierarchy();
			break;
		case WXK_DOWN:
			parent->downHierarchy();
			break;
		case 'q':
		case 'Q':
			parent->setObjectMode( MainFrame::MODE_OBJECT_SELECT );
			break;
		case 'w':
		case 'W':
			parent->setObjectMode( MainFrame::MODE_OBJECT_TRANSLATE );
			break;
		case 'e':
		case 'E':
			parent->setObjectMode( MainFrame::MODE_OBJECT_ROTATE );
			break;
		case 'r':
		case 'R':
			parent->setObjectMode( MainFrame::MODE_OBJECT_SCALE );
			break;
		case WXK_INSERT:
			parent->setObjectMode( MainFrame::MODE_OBJECT_CHANGEPIVOT );
			break;
		case 's':
		case 'S':
			setCameraKey();
			break;
		case 'd':
		case 'D':
			deleteCameraKey();
			break;
		case ',':
		case '<':
			prevCameraKey();
			break;
		case '.':
		case '>':
			nextCameraKey();
			break;
		case 'f':
		case 'F':
			{
				BoundingBox3f bounds;
				std::deque<SceneGraphElementPtr> selected = parent->getSelected();
				if( selected.empty() )
					bounds = parent->scene()->bounds();
				else
				{
					for( std::deque<SceneGraphElementPtr>::const_iterator itr = selected.begin();
						itr != selected.end(); ++itr )
					{
						bounds.expand( (*itr)->bounds() );
					}
				}

				const Path* path = parent->getPath();
				if( path )
				{
					std::vector<size_t> dynamicIds = parent->objectListToDynamicObjectIds( 
						parent->getSelected(), parent->getTree()->dynamicObjects() );

					for( std::vector<size_t>::const_iterator itr = dynamicIds.begin();
						itr != dynamicIds.end(); ++itr )
					{
						std::vector<vl::Vec3f> points = path->obbTree(*itr).points();
						for( std::vector<vl::Vec3f>::const_iterator pointItr = points.begin();
							pointItr != points.end(); ++pointItr )
						{
							bounds.expand( *pointItr );
						}
					}
				}

				HasKeyable<float>::AttributeList keyable =
					this->camera().keyable();
				boost::ptr_deque<AttributeWithKeys> atts;
				for( HasKeyable<float>::AttributeList::const_iterator itr = keyable.begin();
					itr != keyable.end(); ++itr )
				{
					atts.push_back( new AttributeWithKeys(*itr) );
				}

				std::for_each( atts.begin(), atts.end(), 
					boost::bind( &AttributeWithKeys::setKey, _1, 0.0 ) );
				this->camera().frame( bounds );
				std::for_each( atts.begin(), atts.end(), 
					boost::bind( &AttributeWithKeys::setKey, _1, 1.0 ) );

				size_t numSteps = 10;
				for( size_t i = 0; i <= numSteps; ++i )
				{
					std::for_each( atts.begin(), atts.end(), 
						boost::bind( &AttributeWithKeys::setTime, _1, 
							boost::numeric_cast<float>(i) / boost::numeric_cast<float>(numSteps) ) );
					this->Refresh( FALSE );
					this->Update();
				}

				break;
			}
#ifdef WINDOWS_MEDIA
		case WXK_NUMPAD_INSERT:
		case WXK_NUMPAD0:
			{
				parent->recordMovie();
				break;
			}
#endif
#ifdef SOUND
		case WXK_NUMPAD_DELETE:
		case WXK_NUMPAD_DECIMAL:
			{
				parent->recordSound();
				break;
			}		
#endif
		default:
			event.Skip();
			return;
		}
	}

	parent->updateViews();
}

void GLView::nextCameraKey()
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	float time = this->camera_->nextKey( parent->getTime() );
	parent->setTime( time, true );
}

void GLView::prevCameraKey()
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	float time = this->camera_->prevKey( parent->getTime() );
	parent->setTime( time, true );
}

void GLView::setCameraKey()
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	this->camera_->setKey( parent->getTime() );
}

void GLView::deleteCameraKey()
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	this->camera_->deleteKey( parent->getTime() );
	camera_->setTime( parent->getTime() );
}

void initGlew()
{
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		const GLubyte* e = glewGetErrorString(err);
		std::ostringstream oss;
		oss << e;
		errMsg( oss.str(), 0 );
	}
}


boost::once_flag glewOnce = BOOST_ONCE_INIT;

void InitGL(const Camera& camera)
{
	boost::call_once(&initGlew, glewOnce);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glDisable(GL_AUTO_NORMAL);
	glDisable(GL_NORMALIZE);
	glDisable(GL_TEXTURE_2D);
	glShadeModel(GL_SMOOTH);

	{
		// set up lighting
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);

		// set up lighting
		vl::Vec4f lightAmbient0( 0.2, 0.2, 0.2, 1.0 );
		vl::Vec4f lightPosition0( 1, -2, 1, 0 );
		vl::Vec4f lightDiffuse0( 0.5, 0.5, 0.5, 1.0 );
		vl::Vec4f lightPosition1( -2, -1, 5, 0 );
		vl::Vec4f lightDiffuse1( 0.5, 0.5, 0.5, 1 );
		vl::Vec4f lightPosition2 = vl::Vec4f(
			vl::norm(camera.getPosition() - camera.getLookAt()), 0);
		vl::Vec4f lightDiffuse2( 0.3, 0.3, 0.3, 1 );

		glLightfv( GL_LIGHT0, GL_POSITION, norm( lightPosition0 ).Ref() );
		glLightfv( GL_LIGHT0, GL_DIFFUSE, lightDiffuse0.Ref() );
		glLightfv( GL_LIGHT0, GL_AMBIENT, lightAmbient0.Ref() );
		glLightfv( GL_LIGHT1, GL_POSITION, norm( lightPosition1 ).Ref() );
		glLightfv( GL_LIGHT1, GL_DIFFUSE, lightDiffuse1.Ref() );
		glLightfv( GL_LIGHT2, GL_POSITION, lightPosition2.Ref() );
		glLightfv( GL_LIGHT2, GL_DIFFUSE, lightDiffuse2.Ref() );

		glEnable( GL_LIGHT0 );
		glEnable( GL_LIGHT1 );
		glEnable( GL_LIGHT2 );
	}
}


void grid()
{
	// draw direction arrows
	vl::Vec3 origin( -1, 0, -1 );
	for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
	{
		vl::Vec4f color( 0.0, 0.0, 0.0, 1.0 );
		color[iAxis] = 1.0;

		vl::Vec3 orientation(vl::vl_0);
		orientation[iAxis] = 1.0;
		drawShadedArrow( origin, orientation, 0.05, color );
	}


	glDisable(GL_LIGHTING);

	// Draw the grid lines.
	unsigned int x = 0;
	unsigned int y = 2;

	for( int i = 0; i < 2; ++i )
	{
		vl::Vec3 dir_x = vl::vl_0;
		dir_x[x] = 1.0;

		vl::Vec3 dir_y = vl::vl_0;
		dir_y[y] = 1.0;

		{
			glColor3ub( 0, 0, 0 );
			glLineWidth( 2.0 );

			GLActionHandler lines( GL_LINES );
			glVertex3dv( (-dir_y).Ref() );
			glVertex3dv( (dir_y).Ref() );
		}

		{
			glColor3ub( 128, 128, 128 );
			glLineWidth( 1.0 );

			GLActionHandler lines( GL_LINES );
			for( int i = -10; i <= 10; ++i )
			{
				if( i == 0 )
					continue;

				vl::Vec3 begin = 0.1 * boost::numeric_cast<double>(i) * dir_x - dir_y;
				vl::Vec3 end   = 0.1 * boost::numeric_cast<double>(i) * dir_x + dir_y;

				glVertex3dv( begin.Ref() );
				glVertex3dv( end.Ref() );
			}
		}

		std::swap(x, y);
	}
}

struct RenderOnlyDynamicObjectHandler
{
	void operator()( ConstSceneGraphElementPtr object, const vl::Mat4f& transform, const UsefulBits& bits, const char* statePtr )
	{
		boost::shared_ptr<const PhysicsObject> po = boost::dynamic_pointer_cast<const PhysicsObject>( object );
		if( !po )
			return;

		bool isDynamicObject = !po->isStatic();

		ConstSceneGraphElementPtr current = object;
		while( current )
		{
			boost::shared_ptr<const CombinedObject> fused = boost::dynamic_pointer_cast<const CombinedObject>( object );

			if( fused )
				isDynamicObject = !fused->isStatic();
			current = current->parent().lock();
		}

		if( !isDynamicObject )
			return;

		SceneGraphElement::RenderState renderState;
		renderState.alpha = 1.0;
		renderState.polygonMode = PhysicsObject::RenderState::FILL;

		object->render(renderState, transform, statePtr);

	}
};

struct RenderObjectHandler
{
	const MainFrame* parent;

	void operator()( ConstSceneGraphElementPtr object, const vl::Mat4f& transform, const UsefulBits& bits, const char* statePtr )
	{
		const PhysicsObject* physicsObject = dynamic_cast<const PhysicsObject*>( object.get() );
		if( (physicsObject != 0) && (parent->visualizeInertiaTensor()) 
			&& !physicsObject->isStatic() && bits.currentSelected() )
		{
			const int workSize = 512;
			dMass mass = physicsObject->getOdeMassProps();
			vl::Vec3f centerOfMass = vl::Vec3f( mass.c[0], mass.c[1], mass.c[2] );
			dMassTranslate( &mass, -mass.c[0], -mass.c[1], -mass.c[2] );

			fortran_matrix a( 3, 3 );
			for( int i = 0; i < 3; ++i )
				for( int j = 0; j < 3; ++j )
					a(i, j) = mass.I[ i*4 + j ];

			std::vector<double> w;
			boost::tie(a, w) = symmetricJacobi( a );

			double trace = w[0] + w[1] + w[2];
			vl::Vec3d scales;
			for( vl::Int i = 0; i < 3; ++i )
				scales[i] = sqrt( 5.0*(trace - 2*w[i])/2.0*mass.mass );

			ConstMaterialPtr mat = physicsObject->material();
			double density = mat->density();
			double s = pow( 3.0*mass.mass / (4.0*M_PI*mat->density()*scales[0]*scales[1]*scales[2]), 1.0/3.0 );
			for( vl::Int i = 0; i < 3; ++i )
				scales[i] *= s;


			vl::Mat3f localRot( vl::vl_1 );
			for( vl::Int i = 0; i < 3; ++i )
				for( vl::Int j = 0; j < 3; ++j )
					localRot[i][j] = a(i, j) * scales[j];
			{
				vl::Mat3f tmp = vl::trans(localRot);
				vl::Vec3f z = vl::cross( tmp[0], tmp[1] );
				if( vl::dot( tmp[2], z ) < 0.0f )
				{
					std::swap( tmp[0], tmp[1] );
					localRot = vl::trans(tmp);
				}
			}

			vl::Mat4f localTransform( vl::vl_1 );
			for( vl::Int i = 0; i < 3; ++i )
				for( vl::Int j = 0; j < 3; ++j )
					localTransform[i][j] = localRot[i][j];

			for( vl::Int i = 0; i < 3; ++i )
				localTransform[i][3] = centerOfMass[i];

			vl::Mat4f transform = vl::trans(physicsObject->rigidTransform() * localTransform);

			vl::Vec3f color( 1.0f, 1.0f, 0.0f );
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
			glEnable( GL_LIGHTING );
			glEnable( GL_NORMALIZE );
			const GLfloat no_mat[] = { 0.0f, 0.0f, 0.0f, 0.0f };
			glMaterialfv(GL_FRONT, GL_AMBIENT, color.Ref());
			glMaterialfv(GL_FRONT, GL_DIFFUSE, color.Ref());
			glMaterialfv(GL_FRONT, GL_SPECULAR, no_mat);
			glMaterialfv(GL_FRONT, GL_SHININESS, no_mat);
			glMaterialfv(GL_FRONT, GL_EMISSION, no_mat );

			GLMatrixStackHandler pushMatrix;
			glMultMatrixf( transform.Ref() );
			GLUQuadricWrapper fillQuadric;
			gluQuadricDrawStyle( fillQuadric, GLU_FILL );

			const int sphereSubdivisions = 20;
			gluSphere( fillQuadric, 1.0, sphereSubdivisions, sphereSubdivisions );
		}

		SceneGraphElement::RenderState renderState;
		renderState.alpha = 1.0;

		MainFrame::RenderStyle renderStyle = parent->renderStyle();
		switch( renderStyle )
		{
		case MainFrame::STYLE_POLYS:
		case MainFrame::STYLE_BOXES:
			renderState.polygonMode = PhysicsObject::RenderState::FILL;
			if( parent->wireframeOnShaded() )
				renderState.wireframeColor = vl::Vec4f( 40, 43, 112, 100 ) / 255.0;
			break;
		case MainFrame::STYLE_WIREFRAME:
			renderState.polygonMode = PhysicsObject::RenderState::LINES;
			renderState.wireframeColor = vl::Vec4f( 40, 43, 112, 100 ) / 255.0;
			break;
		case MainFrame::STYLE_POINTS:
			renderState.polygonMode = PhysicsObject::RenderState::POINTS;
			break;
		}

		if( parent->visualizeConvexHulls() )
			renderState.renderHulls = true;

		if( bits.selected() && parent->wireframeOnSelected() )
			renderState.wireframeColor = vl::Vec4f( 91, 231, 162, 255 ) / 255.0;
		else if( bits.involvedInJoint() || bits.involvedInConstraint() )
			renderState.wireframeColor = vl::Vec4f( 224, 40, 213, 255 ) / 255.0;

		if( bits.changedFromTree() )
		{
			renderState.overrideColor = true;
			renderState.newColor = vl::Vec3f( 227.0, 12.0, 12.0 ) / 255.0;
		}


		object->render(renderState, transform, statePtr);
	}
};

void GLView::OnPaint( wxPaintEvent& event )
{
	// We will keep a running average of the last 30 frames.
	frameTimes_.push_front( clock() );
	while( frameTimes_.size() > 30 )
		frameTimes_.pop_back();

	// must always be here
    wxPaintDC dc(this);

#ifndef __WXMOTIF__
    if (!GetContext()) return;
#endif

    SetCurrent();
   
	time_t start = clock();

	wxWindow* p = this->GetParent();

	const MainFrame* parent = dynamic_cast<const MainFrame*>(this->GetParent());

	if( !statesToPaint_.empty() )
	{
		glClearColor( 0.0, 0.0, 0.0, 0.0 );
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);
		assert( tileRenderer_ );

		camera().applyViewingTransform();
		InitGL(camera());

		RenderOnlyDynamicObjectHandler handler;
		for( size_t i = 0; i < statesToPaint_.size(); ++i )
			parent->walkTree( handler, statesToPaint_.at(i) );

	    glFlush();

		return;
	}

	vl::Vec3f backgroundColor = parent->backgroundColor();
	glClearColor( backgroundColor[0], backgroundColor[1], backgroundColor[2], 0.0 );
	
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	//glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	if( !tileRenderer_ )
	{
		setupViewport();

		camera().applyProjectiveTransform( 0 );
	}

	camera().applyViewingTransform();

	InitGL(camera());

	/*
	{
		glDisable(GL_LIGHTING);
		glColor4f( 1.0, 0.0, 0.0, 1.0 );
		GLActionHandler lines( GL_LINES );
		glVertex3dv( this->nearTestPt_.Ref() );
		glVertex3dv( this->farTestPt_.Ref() );
	}
	{
		GLActionHandler point( GL_POINTS );
		vl::Vec3d mid = 0.5*(this->nearTestPt_+this->farTestPt_);
		glVertex3dv( mid.Ref() );
	}
	*/

	if( parent->renderGrid() )
		grid();

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glEnable(GL_LIGHTING);

	/*
	{
		PhysicsObject::RenderState renderState;
		renderState.alpha = 1.0;
		renderState.polygonMode = PhysicsObject::RenderState::FILL;
		if( parent->wireframeOnShaded() )
			renderState.wireframeColor = vl::Vec4f( 40, 43, 112, 100 ) / 255.0;
		simulation_->render(renderState);
	}
	*/

	Simulation::State state = parent->state();


	ConstSimulationTreePtr currentTree = parent->getTree();
	if( parent->getMode() == MainFrame::MODE_SOLVE && !parent->playing() && !this->spacePressed_ && currentTree ) 
	{
		// display currently selected simulations
		glDisable( GL_LIGHTING );

        /*
		std::vector<size_t> dynamicObjectIds = 
			parent->objectListToDynamicObjectIds( parent->getSelected(), currentTree->dynamicObjects() );
            */

		std::vector<size_t> dynamicObjectIds(1, 0);

        /*
		if( parent->getPath() )
		{
			// paint selected path:
			glColor4f( 1.0, 1.0, 0.0, 1.0 );
			glLineWidth(2.0);

			GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
			ReaderWriterLock::ScopedReaderLock lock( parent->selectedLineCacheMutex() );

			const std::deque<MainFrame::LineCache>& selectedLineCache = parent->selectedLineCache();
			if( !selectedLineCache.empty() )
			{
				for( std::vector<size_t>::const_iterator obItr = dynamicObjectIds.begin();
					obItr != dynamicObjectIds.end(); ++obItr )
				{
					const MainFrame::LineCache& lineCache = selectedLineCache.at( *obItr );
					if( lineCache.empty() )
						continue;

					glVertexPointer(3, GL_FLOAT, 0, &lineCache[0]);
					glDrawArrays( GL_LINE_STRIP, 0, lineCache.size() );
				}
			}

			glLineWidth(1.0);
		}
        */

		/*
		if( parent->backgroundStyle() == MainFrame::BACKGROUND_WHITE )
			glColor4f( 0.0, 0.0, 0.0, 1.0 );
		else
//			glColor4f( 0.6, 0.6, 0.6, 1.0 );
			glColor4ub( 0, 108, 255, 255 );
			*/

		glLineWidth( 2.0 );


		{
			for( std::vector<size_t>::const_iterator obItr = dynamicObjectIds.begin();
				obItr != dynamicObjectIds.end(); ++obItr )
			{
				size_t iObject = *obItr;
				ReaderWriterLock::ScopedReaderLock lock( parent->lineCacheMutex() );
				if( iObject >= parent->lineCache().size() )
					continue;

				const std::vector<vl::Vec3f>& lineCache = parent->lineCache().at(iObject);
				const std::vector<vl::Vec3f>& lineCacheColors = parent->lineCacheColors().at(iObject);
				const std::vector<GLint>& lineCacheIndices = parent->lineCacheStarts().at(iObject);
				const std::vector<GLsizei>& lineCacheCounts = parent->lineCacheCounts().at(iObject);

				if( !lineCacheIndices.empty() )
				{
					GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
					glVertexPointer(3, GL_FLOAT, 0, &lineCache[0]);

					GLEnableClientStateHandler colorArray(GL_COLOR_ARRAY);
					glColorPointer(3, GL_FLOAT, 0, &lineCacheColors[0]);

					if( glMultiDrawArrays )
					{
						glMultiDrawArrays( GL_LINE_STRIP, 
							const_cast<GLint*>(&lineCacheIndices[0]), 
							const_cast<GLsizei*>(&lineCacheCounts[0]), 
							lineCacheIndices.size() );
					}
					else
					{
						for( size_t i = 0; i < lineCacheIndices.size(); ++i )
						{
							GLint first = lineCacheIndices[i];
							GLint size = lineCacheCounts[i];
							vl::Vec3f tmp1 = lineCache[first];
							vl::Vec3f tmp2 = lineCache[first + size - 1];

							glDrawArrays( GL_LINE_STRIP, first, size );

						}
					}
				}
			}
		}

		/*
		{
			glDisable( GL_LIGHTING );
			glColor4f( 1.0f, 1.0f, 0.0f, 1.0f );

			// Visualize velocities
			for( std::vector<size_t>::const_iterator obItr = dynamicObjectIds.begin();
				obItr != dynamicObjectIds.end(); ++obItr )
			{
				size_t iObject = *obItr;
				if( iObject >= state.stateCount() )
					continue;

				vl::Vec3f position = state.state(iObject).position();
				vl::Vec3f velocity = state.state(iObject).linearVelocity();

				GLActionHandler lines( GL_LINES );
				glVertex3fv( (position + 10.0f*velocity).Ref() );
				glVertex3fv( (position - 10.0f*velocity).Ref() );
			}
		}
		*/
	}

	RenderObjectHandler handler;
	handler.parent = parent;
	parent->walkTree( handler, state );

	std::pair<SceneGraphElementPtr, std::deque<unsigned int> > selectedPoints =
		parent->getSelectedPoints();
	if( selectedPoints.first && parent->getMode() == MainFrame::MODE_EDIT &&
		parent->getObjectMode() == MainFrame::MODE_OBJECT_SELECTPOINTS )
	{
		const char* statePtr = 0;

		/*
		ConstPhysicsObjectPtr physObj = 
			boost::dynamic_pointer_cast<const PhysicsObject>( selectedPoints.first );
		if( physObj )
		{
			std::vector<ConstPhysicsObjectPtr> dynamicObjects;
			if( this->simulation_ )
				dynamicObjects = simulation_->dynamicObjects();
			else if( this->currentTree_ )
				dynamicObjects = currentTree_->dynamicObjects();

			std::vector<size_t> dynamicObjectIds = 
				parent->objectListToDynamicObjectIds(  
					std::vector<ConstPhysicsObjectPtr>(1, physObj),
					dynamicObject );
			if( !dynamicObjectIds.empty() )
				statePtr = state.state( dynamicObjectIds.front() );
		}
		*/

		vl::Mat4f transform = selectedPoints.first->transform( statePtr );
		selectedPoints.first->renderPoints( selectedPoints.second, transform, statePtr );
	}

	{
		// display all constraints
		ConstSimulationTreePtr tree = parent->getTree();
		if( tree )
		{
			glDisable( GL_LIGHTING );
			std::deque<ConstraintPtr> constraints = tree->constraints();
			ConstraintPtr selectedConstraint = parent->selectedConstraint();
			for( std::deque<ConstraintPtr>::const_iterator constraintItr = constraints.begin();
				constraintItr != constraints.end(); ++constraintItr )
			{
				if( *constraintItr == selectedConstraint )
					(*constraintItr)->render(true);
				else
					(*constraintItr)->render(false);
			}
		}
	}

	boost::shared_ptr<const MousingAction> mouseAction = parent->mouseAction();
	if( !parent->isSimulating() && mouseAction )
	{
		CameraTransformConverter converter( this->camera(), this->viewportWidth_, this->viewportHeight_ );
		mouseAction->render(converter);
	}

	if( parent->visualizeContactPoints() )
	{
		GLUQuadricWrapper fillQuadric;
		gluQuadricDrawStyle( fillQuadric, GLU_FILL );

		glPolygonMode( GL_FRONT, GL_FILL );
		std::vector<Simulation::ContactPoint> contactPoints = parent->contactPoints();
		for( std::vector<Simulation::ContactPoint>::const_iterator iter = contactPoints.begin();
			iter != contactPoints.end(); ++iter )
		{
			{
				GLMatrixStackHandler pushMatrix;
				glTranslatef( iter->position[0], iter->position[1], iter->position[2] );
				gluSphere( fillQuadric, 0.1, 10, 10 );
			}

			drawShadedArrow( toVec3d(iter->position), 
				toVec3d(iter->normal), 
				1.0, 
				vl::Vec4f(1.0, 0.0, 0.0, 1.0), 
				false );
		}
	}

	time_t end = clock();
	{
		time_t elapsed = end - start;
		double newFrameTime = boost::numeric_cast<double>( elapsed ) 
			/ boost::numeric_cast<double>( CLOCKS_PER_SEC );
		
		const double alpha = 0.1;
		frameTime_ = alpha*newFrameTime + (1.0-alpha)*frameTime_;
	}

	if( headsUp_ )
	{
		glDisable(GL_CULL_FACE);
		glDisable(GL_LIGHTING);
		GLDisableHandler depthTest( GL_DEPTH_TEST );

		if( parent->GetCanvas() == this )
		{
			glMatrixMode (GL_PROJECTION);
			glLoadIdentity();
			gluOrtho2D(0, 1, 0, 1);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			glLineWidth( 4 );
			glColor4f( 0.2, 0.2, 1.0, 1.0 );
			GLActionHandler lineLoop( GL_LINE_LOOP );
			glVertex2f( 0, 0 );
			glVertex2f( 0, 1 );
			glVertex2f( 1, 1 );
			glVertex2f( 1, 0 );
		}

		GLOverlayHandler overlay;

		glColor4ub(255, 255, 255, 255);

		typedef std::pair<std::string, std::string> Value;
		std::deque<Value> values;

/*
		{
			std::ostringstream oss;
			oss << setiosflags(std::ios::fixed) << std::setprecision(1) << fps();
			values.push_back( Value( "fps", oss.str() ) );
		}
		{
			std::ostringstream oss;
			oss << setiosflags(std::ios::fixed) << std::setprecision(1) << (1000.0*frameTime_) << "ms";
			values.push_back( Value( "time/frame", oss.str() ) );
		}
		{
			std::ostringstream oss;
			oss << simulationTree_->nodesExpanded();
			values.push_back( Value( "nodes expanded", oss.str() ) );
		}
*/
		/*
		debugging stuff

		if( !this->frames_.empty() )
		{
			{
				std::ostringstream oss;
				NxVec3 pos = this->frames_[frameNo_].state( 0 ).position();
				oss << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")";
				values.push_back( Value( "position", oss.str() ) );
			}
			{
				std::ostringstream oss;
				NxVec3 mo = this->frames_[frameNo_].state( 0 ).linearMomentum();
				oss << "(" << mo[0] << ", " << mo[1] << ", " << mo[2] << ")";
				values.push_back( Value( "linear momentum", oss.str() ) );
			}
		}
		*/

//		printStrings( vl::Vec3(-90.0, 90.0, 0.0), values );

		std::string name = this->camera_->name();
		printString( name.c_str(), vl::Vec3(-90.0, 90.0, 0.0), ALIGN_HORIZ_LEFT, ALIGN_VERT_TOP, true );
	}

	if( selecting_ )
	{
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, viewportWidth_, 0, viewportHeight_, -1.0, 1.0);

		GLDisableHandler disableLighting(GL_LIGHTING);
		GLDisableHandler depthTest(GL_DEPTH_TEST);

        if( parent->mouseMode() == MainFrame::MODE_MOUSE_SKETCH_ACTION ||
            parent->mouseMode() == MainFrame::MODE_MOUSE_SKETCH_VELOCITY )
        {
            if( !sketchPoints_.empty() )
            {
                glColor4ub(255, 253, 118, 255);
                GLEnableClientStateHandler vertexArray( GL_VERTEX_ARRAY );
                glVertexPointer( 2, GL_FLOAT, 0, &sketchPoints_[0] );
                glDrawArrays( GL_LINE_STRIP, 0, sketchPoints_.size() );
            }
        }
        else
        {
		    boost::scoped_ptr<GLEnableHandler> lineStipple;
		    if( parent->mouseMode() == MainFrame::MODE_MOUSE_SELECT )
		    {
			    glColor3ub(255,255,255);

			    lineStipple.reset( new GLEnableHandler( GL_LINE_STIPPLE ) );
			    glLineWidth( 1.0 );
			    glLineStipple(1, 0xF0F0);
		    }
		    else if( parent->mouseMode() == MainFrame::MODE_MOUSE_ADDITIVE_CONSTRAINT )
		    {
			    glColor4ub(40, 255, 40, 100);
			    glLineWidth( 3.0 );
		    }
		    else if( parent->mouseMode() == MainFrame::MODE_MOUSE_SUBTRACTIVE_CONSTRAINT )
		    {
			    glColor4ub(255, 40, 40, 100);
			    glLineWidth( 3.0 );
		    }
		    else if( parent->mouseMode() == MainFrame::MODE_MOUSE_REFINE )
		    {
			    glColor4ub(255, 255, 40, 100);
			    glLineWidth(3.0);
		    }

		    if( parent->mouseMode() == MainFrame::MODE_MOUSE_ADDITIVE_CONSTRAINT
			    || parent->mouseMode() == MainFrame::MODE_MOUSE_SUBTRACTIVE_CONSTRAINT
			    || parent->mouseMode() == MainFrame::MODE_MOUSE_REFINE )
		    {
			    GLEnableHandler blend( GL_BLEND );
			    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			    BoundingBox2d bounds;
			    bounds.expand( this->selectingAnchor_ );
			    bounds.expand( this->selectingFree_ );
			    vl::Vec2d lower = bounds.minimum();
			    vl::Vec2d upper = bounds.maximum();
			    {
				    GLActionHandler quads(GL_QUADS);
				    glVertex2d( lower[0], lower[1] );
				    glVertex2d( lower[0], upper[1] );
				    glVertex2d( upper[0], upper[1] );
				    glVertex2d( upper[0], lower[1] );
			    }
		    }

		    {
			    GLActionHandler lineLoop( GL_LINE_LOOP );
			    glVertex2d( selectingAnchor_[0], selectingAnchor_[1] );
			    glVertex2d( selectingAnchor_[0], selectingFree_[1] );
			    glVertex2d( selectingFree_[0], selectingFree_[1] );
			    glVertex2d( selectingFree_[0], selectingAnchor_[1] );
		    }
        }
	}

    glFlush();

	if( !this->tileRenderer_ )
		SwapBuffers();
}

void GLView::OnSize(wxSizeEvent& event)
{
    // this is also necessary to update the context on some platforms
    wxGLCanvas::OnSize(event);
    
    // set GL viewport (not called by wxGLCanvas::OnSize on all platforms...)
    GetClientSize(&viewportWidth_, &viewportHeight_);
#ifndef __WXMOTIF__
    if (GetContext())
#endif
    {
//        SetCurrent();
		this->Refresh(FALSE);
    }
}

void GLView::OnEraseBackground(wxEraseEvent& event)
{
    // Do nothing, to avoid flashing on MSW
}

class SelectClosestPath
{
public:
	SelectClosestPath( 
		const vl::Vec2f& clickedPoint,
		const vl::Mat4f& mat,
		const SimulationTree& tree, 
		const boost::dynamic_bitset<>& paths,
		const std::vector<size_t>& dynamicObjectIds,
		std::pair<float, const Path*>& result,
		size_t start,
		size_t end,
		float maxZ )
		:	clickedPoint_(clickedPoint),
			matrix_(mat),
			tree_(tree),
			paths_(paths),
			dynamicObjectIds_(dynamicObjectIds),
			start_(start),
			end_(end),
			result_(result),
			maxZ_(maxZ) {}

	void operator()()
	{
		for( size_t i = start_; i < end_; ++i )
		{
			if( !paths_[i] )
				continue;

			const Path* path = tree_.path( i );
			for( std::vector<size_t>::const_iterator obIter = dynamicObjectIds_.begin();
				obIter != dynamicObjectIds_.end(); ++obIter )
			{
				float obbTreeDist = path->obbTree(*obIter).dist( clickedPoint_, 
					matrix_, this->result_.first, this->maxZ_ );
				if( obbTreeDist < this->result_.first )
				{
                    this->result_.second = path;
					this->result_.first = obbTreeDist;
				}
			}
		}
	}

private:
	vl::Vec2f clickedPoint_;
	vl::Mat4f matrix_;
	const boost::dynamic_bitset<>& paths_;
	const SimulationTree& tree_;
	const std::vector<size_t>& dynamicObjectIds_;
	size_t start_;
	size_t end_;
	std::pair<float, const Path*>& result_;
	float maxZ_;
};

struct SelectBoxHandler
{
	SelectBoxHandler( const BoundingBox2d& b, const ScreenSpaceConverter& c )
		: box( b ), converter( c ) {}

	std::deque<SceneGraphElementPtr> newSelected;
	BoundingBox2d box;
	const ScreenSpaceConverter& converter;

	void operator()( SceneGraphElementPtr element, const vl::Mat4f& transform, const UsefulBits& bits, const char* statePtr )
	{
		if( element->intersect( box, converter, transform, statePtr ) )
			newSelected.push_back( element );
	}
};

struct SelectRayHandler
{
	Ray3d ray;
	double t;
	SceneGraphElementPtr intersected;

	SelectRayHandler( const Ray3d& r )
		: ray(r), t( boost::numeric::bounds<double>::highest() ) {}

	void operator()( SceneGraphElementPtr object, const vl::Mat4f& transform, UsefulBits bits, const char* statePtr )
	{
		double tTemp;
		if( object->intersect( ray, tTemp, transform, statePtr ) )
		{
			if( tTemp < t )
			{
				intersected = object;
				t = tTemp;
			}
		}
	}
};

void GLView::OnMouse( wxMouseEvent& event )
{
 	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	wxSize sz(GetClientSize());

	camera_->camera()->setRatio( boost::numeric_cast<double>(viewportWidth_) 
		/ boost::numeric_cast<double>(viewportHeight_) );
	CameraTransformConverter converter( this->camera(), this->viewportWidth_, this->viewportHeight_ );

	boost::shared_ptr<MousingAction> mouseAction = parent->mouseAction();

	vl::Mat4f matrix = this->camera_->camera()->projectiveTransform() 
		* this->camera_->camera()->modelViewTransform();

	if( event.ButtonDown() )
	{
		this->SetFocus();
		parent->setCurrentView( this, false );
		mouseCapture_.reset( new MouseCapture(this) );

		if( event.AltDown() )
		{
			if( event.LeftIsDown() )
				camera().clickMouse(kActionRotate, event.GetX(), event.GetY() );
			else if( event.MiddleIsDown() )
				camera().clickMouse(kActionTranslate, event.GetX(), event.GetY() );
			else if( event.RightIsDown() )
				camera().clickMouse(kActionZoom, event.GetX(), event.GetY() );
		}
		else
		{
			if( mouseAction && 
				mouseAction->pressMouse( vl::Vec2d(event.GetX(), viewportHeight_ - event.GetY()), converter ) )
			{
				parent->updateViews();
				return;
			}
            else if( parent->mouseMode() == MainFrame::MODE_MOUSE_SKETCH_ACTION ||
                parent->mouseMode() == MainFrame::MODE_MOUSE_SKETCH_VELOCITY )
            {
                this->sketchPoints_.clear();
                this->times_.clear();

                this->sketchPoints_.push_back( 
                    vl::Vec2f( event.GetX(), viewportHeight_ - event.GetY() ) );
                this->times_.push_back( LARGE_INTEGER() );
                BOOL ret = QueryPerformanceCounter( &this->times_.back() );
                assert( ret != 0 );

                selecting_ = true;
				selectingAnchor_ = vl::Vec2d(event.GetX(), viewportHeight_ - event.GetY());
				selectingFree_ = selectingAnchor_;

                Refresh( FALSE );
                return;
            }
			else
			{
				selecting_ = true;
				selectingAnchor_ = vl::Vec2d(event.GetX(), viewportHeight_ - event.GetY());
				selectingFree_ = selectingAnchor_;

				Refresh( FALSE );
				return;
			}
		}
	}
	else if( event.Dragging() )
	{
		if( camera().moving() )
		{
			camera().dragMouse(event.GetX(), event.GetY());
			Refresh(FALSE);
			return;
		}

        if( parent->mouseMode() == MainFrame::MODE_MOUSE_SKETCH_ACTION ||
            parent->mouseMode() == MainFrame::MODE_MOUSE_SKETCH_VELOCITY )
        {
            this->sketchPoints_.push_back( 
                vl::Vec2f( event.GetX(), viewportHeight_ - event.GetY() ) );
            this->times_.push_back( LARGE_INTEGER() );
            BOOL ret = QueryPerformanceCounter( &this->times_.back() );
            assert( ret != 0 );

            Refresh( FALSE );
            return;
        }
		else if( selecting_ )
		{
			vl::Vec2d oldSelectingFree = selectingFree_;
			selectingFree_ = vl::Vec2d(event.GetX(), viewportHeight_ - event.GetY());

			if( parent->getMode() == MainFrame::MODE_SOLVE 
				&& (parent->mouseMode() == MainFrame::MODE_MOUSE_ADDITIVE_CONSTRAINT ||
					parent->mouseMode() == MainFrame::MODE_MOUSE_SUBTRACTIVE_CONSTRAINT) )
			{
				BoundingBox2f box;
				box.expand( vl::Vec2f( 
					selectingFree_[0] / boost::numeric_cast<float>(this->viewportWidth_),
					selectingFree_[1] / boost::numeric_cast<float>(this->viewportHeight_) ) );
				box.expand( vl::Vec2f( 
					selectingAnchor_[0] / boost::numeric_cast<float>(this->viewportWidth_),
					selectingAnchor_[1] / boost::numeric_cast<float>(this->viewportHeight_) ) );
			}

			Refresh(FALSE);
			return;
		}

		if( mouseAction )
		{
			mouseAction->moveMouse( vl::Vec2d(event.GetX(), viewportHeight_ - event.GetY()), converter );
			parent->updateViews();
			return;
		}
	}
	else if( event.ButtonUp() )
	{
		mouseCapture_.reset();

		if( camera().moving() )
		{
			camera().releaseMouse( event.GetX(), event.GetY() );
			parent->updateViews();
//			updateQuadTree();
			return;
		}

		if( mouseAction && mouseAction->moving() )
		{
			mouseAction->releaseMouse( vl::Vec2d(event.GetX(), viewportHeight_ - event.GetY()), converter );
			ActionPtr action = mouseAction->action();
			parent->performAction(action);

			parent->setSelected( parent->getSelected() );
			parent->updateViews();
			return;
		}

		/*
		{
			vl::Vec2d clickedPoint( event.GetX(), viewportHeight_ - event.GetY() );
			vl::Vec2f normalized( 
				(clickedPoint[0] - 0.0f) * (2.0f / static_cast<float>(this->viewportWidth_)) - 1.0f,
				(clickedPoint[1] - 0.0f) * (2.0f / static_cast<float>(this->viewportHeight_)) - 1.0f );
			std::ofstream ofs( "clicks.bin", std::ios::binary | std::ios::app );
			vl::Mat4f modelViewTransform = this->camera_->camera()->modelViewTransform();
			vl::Mat4f projectiveTransform = this->camera_->camera()->projectiveTransform();
			ofs.write( reinterpret_cast<const char*>( modelViewTransform.Ref() ), sizeof( vl::Mat4f ) );
			ofs.write( reinterpret_cast<const char*>( projectiveTransform.Ref() ), sizeof( vl::Mat4f ) );
			ofs.write( reinterpret_cast<const char*>( normalized.Ref() ), sizeof( vl::Vec2f ) );
		}
		*/

		if( !selecting_ )
			return;

        SimulationTreePtr tree = parent->getTree();
        if( tree &&
            (parent->mouseMode() == MainFrame::MODE_MOUSE_SKETCH_ACTION ||
             parent->mouseMode() == MainFrame::MODE_MOUSE_SKETCH_VELOCITY) &&
             sketchPoints_.size() > 3 &&
             vl::len( sketchPoints_.front() - sketchPoints_.back() ) > 3.0 )
        {
            // need to decide if user indicated flat or bounce
            fortran_matrix A( this->sketchPoints_.size(), 3 );
            fortran_matrix B( this->sketchPoints_.size(), 3 );
            std::vector<double> x( A.nrows() );
            std::vector<double> y( A.nrows() );

            LARGE_INTEGER frequency;
            int result = QueryPerformanceFrequency( &frequency );
            assert( result != 0 );

            double maxTime = 0;
            for( size_t i = 0; i < A.nrows(); ++i )
            {
                double time = boost::numeric_cast<double>( this->times_[i].QuadPart - this->times_.front().QuadPart ) /
                    boost::numeric_cast<double>( frequency.QuadPart );
                maxTime = std::max( time, maxTime );
                A(i, 0) = 1;
                A(i, 1) = time;
                A(i, 2) = time*time;

                vl::Vec2f point = this->sketchPoints_.at( i );
                B(i, 0) = 1;
                B(i, 1) = point[0];
                B(i, 2) = point[0]*point[0];

                x.at(i) = point[0];
                y.at(i) = point[1];
            }

            {
                MATFile matFile( "curves.mat", "curves" );
                matFile.add( "A", A );
                matFile.add( "B", B );
                matFile.add( "x", x );
                matFile.add( "y", y );
            }

            fortran_matrix A1 = A;

            // obviously would be faster to do these at the same time:
            std::vector<double> xCoeffs = leastSquares_QR( A, &x[0] );
            std::vector<double> yCoeffs = leastSquares_QR( A1, &y[0] );
            std::vector<double> x_vs_yCoeffs = leastSquares_QR( B, &y[0] );
            assert( xCoeffs.size() == 3 );
            assert( yCoeffs.size() == 3 );
            assert( x_vs_yCoeffs.size() == 3 );

			std::vector<size_t> dynamicObjectIds = 
				parent->objectListToDynamicObjectIds( parent->getSelected(), tree->dynamicObjects() );

            // now, what is our threshold on the quadratic term?
            double val = fabs( x_vs_yCoeffs.at(2) );
            if( val < 0.001 )
            {
                // need to figure out the direction
                vl::Vec2f minPos( xCoeffs[0], yCoeffs[0] );
                vl::Vec2f maxPos( 
                    xCoeffs[0] + xCoeffs[1]*maxTime + xCoeffs[2]*maxTime*maxTime,
                    yCoeffs[0] + yCoeffs[1]*maxTime + yCoeffs[2]*maxTime*maxTime );

                std::vector<ConstPhysicsObjectPtr> dynOb = tree->dynamicObjects();
                Simulation::State state = parent->state();

                // Mark selected objects as "slide"
                for( std::vector<size_t>::const_iterator itr = dynamicObjectIds.begin();
                    itr != dynamicObjectIds.end(); ++itr )
                {
                    SamplingProperties props;
                    props.slide = true;

                    // now, we project the direction onto the plane containing the object
                    // for the moment we'll assume the center of mass corresponds to the 
                    // object's origin, which really _should_ be the case if we want decent
                    // compression, but is not always the case...
                    // @todo fix this assumption
                    vl::Vec3f objPos = vl::xform( dynOb.at( *itr )->transform(), vl::Vec3f( vl::vl_0 ) );
                    if( *itr < state.stateCount() )
                        objPos = state.state( *itr ).position();

                    vl::Vec2f screenSpacePositions[] = { minPos, maxPos };
                    vl::Vec3d worldSpacePositions[2];
                    for( int i = 0; i < 2; ++i )
                    {
                        vl::Vec2f screenPos = screenSpacePositions[i];
                        vl::Vec3d start = converter.toWorldSpace( vl::Vec3d(screenPos[0], screenPos[1], 0.0) );
                        vl::Vec3d end = converter.toWorldSpace( vl::Vec3d(screenPos[0], screenPos[1], 1.0) );
                        vl::Vec3d p = start;
                        vl::Vec3d d = end - start;
                        double t = (objPos[1] - p[1]) / d[1];
                        worldSpacePositions[i] = p + t*d;
                    }

                    vl::Vec3d dir = vl::norm( worldSpacePositions[1] - worldSpacePositions[0] );

                    props.frictionDirection = toVec3f( dir );
                    tree->setSamplingProperties( *itr, props );
                }
            }
            else
            {
                // Mark selected objects as "bounce"
                for( std::vector<size_t>::const_iterator itr = dynamicObjectIds.begin();
                    itr != dynamicObjectIds.end(); ++itr )
                {
                    SamplingProperties props;
                    props.bounce = true;
                    tree->setSamplingProperties( *itr, props );
                }
            }
        }

		std::deque<SceneGraphElementPtr> newSelected;

		selecting_ = false;
		selectingFree_ = vl::Vec2( event.GetX(), viewportHeight_ - event.GetY() );

        /*
		setupViewport();
		camera().applyProjectiveTransform( 0 );
		camera().applyViewingTransform();
        */

		bool isClick = (len( selectingAnchor_ - selectingFree_ ) < 3.0);
		if( parent->mouseMode() == MainFrame::MODE_MOUSE_SELECT || isClick )
		{
			// pick an object
			if( isClick )
			{
				vl::Vec2d clickedPoint( event.GetX(), viewportHeight_ - event.GetY() );
				vl::Vec3d near_v = converter.toWorldSpace( vl::Vec3d( clickedPoint, 0.0 ) );
				vl::Vec3d far_v = converter.toWorldSpace( vl::Vec3d( clickedPoint, 1.0 ) );
				vl::Vec3d dir = vl::norm( far_v - near_v );

				/*
				this->nearTestPt_ = near_v;
				this->farTestPt_ = far_v;
				this->Refresh(FALSE);
				*/

				Ray3d ray( near_v, dir );
				SelectRayHandler handler( ray );
				parent->walkTree( handler, parent->state()  );

				float intersectedZ = 1.0;

				ConstraintPtr newSelectedConstraint;
				ConstSimulationTreePtr tree = parent->getTree();
				if( tree )
				{
					double min_t = handler.t;
					std::deque<ConstraintPtr> constraints = tree->constraints();
					for( size_t iConstraint = 0; iConstraint < constraints.size(); ++iConstraint )
					{
						double tmp_t;
						if( constraints[iConstraint]->intersect( ray, tmp_t ) 
							&& (tmp_t < min_t) )
						{
							newSelectedConstraint = constraints[iConstraint];
							min_t = tmp_t;
						}
					}
				}

				if( handler.intersected || newSelectedConstraint )
				{
					vl::Vec3d intersectedPoint = ray.at(handler.t);
					vl::Vec3d intersectedPoint_clipSpace( converter.toScreenSpace(intersectedPoint) );

					intersectedZ = intersectedPoint_clipSpace[2];
					newSelected.push_back( handler.intersected );
				}

				// must be within 5 pixels to select:
				float minPathDist = 10.0f;

				if( parent->getMode() == MainFrame::MODE_SOLVE
					&& !parent->playing() 
					&& !this->spacePressed_)
				{
					vl::Vec2f clickedPoint2f( clickedPoint[0] - boost::numeric_cast<float>(viewportWidth_) / 2.0f, 
						clickedPoint[1] - boost::numeric_cast<float>(viewportHeight_) / 2.0f );

					vl::Mat4f viewportMat = matrix;
					for( vl::Int i = 0; i < 4; ++i )
					{
						viewportMat[0][i] *= boost::numeric_cast<float>(viewportWidth_) / 2.0f;
						viewportMat[1][i] *= boost::numeric_cast<float>(viewportHeight_) / 2.0f;
					}

					SimulationTreePtr currentTree = parent->getTree();
					boost::dynamic_bitset<> paths = 
						currentTree->activePaths();
					std::vector<size_t> dynamicObjectIds = 
						parent->objectListToDynamicObjectIds( parent->getSelected(), currentTree->dynamicObjects() );

					size_t numThreads = 1;
					if( dynamicObjectIds.size() * paths.size() > 1000 )
						numThreads = numberOfProcessors;

					std::vector< std::pair<float, const Path*> > closestPaths( numThreads,
						std::make_pair( minPathDist, (const Path*) 0 ) );
					boost::ptr_vector<SelectClosestPath> selectPathsThreadFunctions;
					selectPathsThreadFunctions.reserve( numThreads );

					for( size_t i = 0; i < numThreads; ++i )
					{
						selectPathsThreadFunctions.push_back( new SelectClosestPath( 
								clickedPoint2f,
								viewportMat,
								*currentTree,
								paths,
								dynamicObjectIds,
								closestPaths[i],
								i*paths.size() / numThreads,
								(i+1)*paths.size() / numThreads, intersectedZ ) );
					}

					if( selectPathsThreadFunctions.size() > 1 )
					{
						boost::thread_group threads;
						std::for_each( selectPathsThreadFunctions.begin(), selectPathsThreadFunctions.end(),
							boost::bind( &boost::thread_group::create_thread, boost::ref(threads), _1 ) );
						threads.join_all();
					}
					else
					{
						(selectPathsThreadFunctions.front())();
					}

					double closest = minPathDist;
					const Path* closestPath = 0;
					for( size_t i = 0; i < closestPaths.size(); ++i )
					{
						if( closestPaths[i].second && 
							closestPaths[i].first < closest )
						{
							closestPath = closestPaths[i].second;
							closest = closestPaths[i].first;
						}
					}

					if( closestPath )
					{
						boost::shared_ptr<Action> setNodeAction( 
							new SetPathSelectionAction(
								std::make_pair(currentTree, closestPath),
								std::make_pair(currentTree, parent->getPath()),
								parent ) );
						parent->performAction( setNodeAction );
						return;
					}
				}
				else if( parent->getObjectMode() == MainFrame::MODE_OBJECT_SELECTPOINTS )
				{
					// if selecting points, we'll see if we can select any in the currently
					// selected object
					MainFrame::SelectedPoints currentSelected = parent->getSelectedPoints();
					if( currentSelected.first )
					{
						// todo this does not use the correct transform always
						std::pair<bool, unsigned int> newSelectedPoint = 
							currentSelected.first->selectPoint( clickedPoint, converter, currentSelected.first->transform() );
						if( newSelectedPoint.first )
						{
							std::deque<unsigned int> newSelectedPoints( 1, newSelectedPoint.second );
							if( event.ShiftDown() )
							{
								std::deque<unsigned int> output;
								std::set_symmetric_difference( currentSelected.second.begin(),
									currentSelected.second.end(),
									newSelectedPoints.begin(),
									newSelectedPoints.end(),
									std::back_inserter(output) );
								output.swap( newSelectedPoints );
							}

							ActionPtr changeSelectionAction( 
								new ChangeSelectionAction(
									parent,
									parent->getSelected(), parent->getSelected(),
									parent->selectedConstraint(), parent->selectedConstraint(),
									currentSelected, MainFrame::SelectedPoints(currentSelected.first, newSelectedPoints) ) );
							parent->performAction( changeSelectionAction );
							return;
						}
					}
				}

				if( newSelectedConstraint )
				{
					ActionPtr changeSelectionAction( 
						new ChangeSelectionAction(
							parent,
							parent->getSelected(), std::deque<SceneGraphElementPtr>(),
							parent->selectedConstraint(), newSelectedConstraint ) );
					parent->performAction( changeSelectionAction );
					return;
				}
			}
			else
			{
				BoundingBox2d box;
				box.expand( selectingAnchor_ );
				box.expand( selectingFree_ );

				SelectBoxHandler handler( box, converter );
				parent->walkTree( handler, parent->state()  );
				newSelected = handler.newSelected;
			}

			newSelected = parent->convertSelected( newSelected );

			std::deque<SceneGraphElementPtr> oldSelected = parent->getSelected();
			newSelected = parent->convertSelected( newSelected );
			if( event.ShiftDown() )
			{
				std::deque<SceneGraphElementPtr> result;

				// need to perform a symmetric difference between the two sets,
				//   but retain ordering
				// we first transfer all the members from the current selected set
				//   into the result set provided they don't appear in the second
				//   set.  We will use a hash_set on all elements of new_selected
				//   to figure this out
				typedef STLEXT hash_set<const SceneGraphElement*, 
					hash_ptr<const SceneGraphElement*> > SceneGraphElementHash;
				{
					SceneGraphElementHash newSelSet( 
						boost::make_transform_iterator(newSelected.begin(), 
							boost::bind(&SceneGraphElementPtr::get, _1) ),
						boost::make_transform_iterator(newSelected.end(), 
							boost::bind(&SceneGraphElementPtr::get, _1) ) );
					for( std::deque<SceneGraphElementPtr>::const_iterator oldItr = oldSelected.begin();
						oldItr != oldSelected.end(); ++oldItr )
					{
						SceneGraphElementHash::iterator hashItr = newSelSet.find( oldItr->get() );
						if( hashItr == newSelSet.end() )
							result.push_back(*oldItr);
					}
				}

				// now, we transfer any elements of the second set that aren't in the
				//   first set
				{
					SceneGraphElementHash oldSelSet( 
						boost::make_transform_iterator(oldSelected.begin(), 
							boost::bind(&SceneGraphElementPtr::get, _1) ),
						boost::make_transform_iterator(oldSelected.end(), 
							boost::bind(&SceneGraphElementPtr::get, _1) ) );
					for( std::deque<SceneGraphElementPtr>::const_iterator newItr = newSelected.begin();
						newItr != newSelected.end(); ++newItr )
					{
						SceneGraphElementHash::iterator hashItr = oldSelSet.find( newItr->get() );
						if( hashItr == oldSelSet.end() )
							result.push_back(*newItr);
					}
				}

				newSelected.swap( result );
			}

			std::deque<SceneGraphElementPtr> cleaned = 
				cleanSelected(newSelected);
			if( cleaned == oldSelected )
			{
				this->Refresh(FALSE);
				return;
			}

			assert( cleaned.size() != newSelected.size() 
				|| newSelected == cleaned );

			ActionPtr changeSelectionAction( 
				new ChangeSelectionAction(
					parent,
					oldSelected,
					cleaned,
					parent->selectedConstraint(), ConstraintPtr() ) );
			parent->performAction( changeSelectionAction );
		}

		BoundingBox3f box;
		box.expand( vl::Vec3f( 
			2.0*selectingFree_[0] / this->viewportWidth_ - 1.0,
			2.0*selectingFree_[1] / this->viewportHeight_ - 1.0,
			-1.0 ) );
		box.expand( vl::Vec3f( 
			2.0*selectingAnchor_[0] / this->viewportWidth_ - 1.0,
			2.0*selectingAnchor_[1] / this->viewportHeight_ - 1.0,
			1.0 ) );

		if( parent->getMode() == MainFrame::MODE_SOLVE 
			&& (!parent->playing())
			&& (parent->mouseMode() == MainFrame::MODE_MOUSE_REFINE) )
		{
			float earliestTime = boost::numeric::bounds<float>::highest();
			float latestTime = boost::numeric::bounds<float>::lowest();

			// assume compressed path:
			SimulationTreePtr currentTree = parent->getTree();
			const Path* path = parent->getPath();
			std::vector<PiecewisePath> piecewisePaths = parent->piecewisePaths();
			assert( !piecewisePaths.empty() || path == 0 );

			boost::dynamic_bitset<> headInPath( piecewisePaths.size(), 0 );
			boost::dynamic_bitset<> tailInPath( piecewisePaths.size(), 0 );

			std::vector<size_t> dynamicObjectIds = 
				parent->objectListToDynamicObjectIds( parent->getSelected(), currentTree->dynamicObjects() );
			if( !piecewisePaths.empty() && !dynamicObjectIds.empty() )
			{
				for( std::vector<size_t>::const_iterator dynamicObjectItr = dynamicObjectIds.begin();
					dynamicObjectItr != dynamicObjectIds.end(); ++dynamicObjectItr )
				{
					const PiecewisePath& path = piecewisePaths.at( *dynamicObjectItr );
					for( long iFrame = path.startFrame(); iFrame < path.endFrame(); ++iFrame )
					{
						vl::Vec3f pos = path.position( iFrame );
						vl::Vec3f projected = vl::xform( matrix, pos );
						if( contains(box, projected) )
						{
							earliestTime = std::min<float>( iFrame, earliestTime );
							latestTime = std::max<float>( iFrame, latestTime );

							if( iFrame == path.startFrame() )
								headInPath[ *dynamicObjectItr ] = true;
							if( iFrame == (path.endFrame() - 1) )
								tailInPath[ *dynamicObjectItr ] = true;
                        }
					}
				}
			}

			// @todo fix this
			// @todo why?
			std::deque<SceneGraphElementPtr> selected = parent->getSelected();
			std::vector<ConstPhysicsObjectPtr> dynamicObjects = currentTree->dynamicObjects();
			std::deque<ConstPhysicsObjectPtr> refinedObjects;

			for( std::vector<ConstPhysicsObjectPtr>::const_iterator itr = dynamicObjects.begin();
				itr != dynamicObjects.end(); ++itr )
			{
				ConstSceneGraphElementPtr current = *itr;
				while( current )
				{
					if( std::find_if( selected.begin(), selected.end(), 
						NameEquals(current->name()) ) != selected.end() )
					{
						refinedObjects.push_back( *itr );
						break;
					}

					current = current->parent().lock();
				}
			}

			// Our heuristic is that if the user selects the head of a path
			// he wants to refine backwards, unless he selects the whole path
			if( headInPath.any() && !tailInPath.any() )
			{
				parent->startSubtree( currentTree, 
					piecewisePaths, 
					parent->piecewisePathsFrameRate(),
					latestTime / boost::numeric_cast<float>(parent->piecewisePathsFrameRate()),
					refinedObjects,
					true );
				
			}
			else if( earliestTime < boost::numeric::bounds<float>::highest() )
			{
				parent->startSubtree( currentTree, 
					piecewisePaths, 
					parent->piecewisePathsFrameRate(),
					earliestTime / boost::numeric_cast<float>(parent->piecewisePathsFrameRate()),
					refinedObjects,
					false );
			}
		}

		if( parent->getMode() == MainFrame::MODE_SOLVE 
			&& (!parent->playing())
			&& (parent->mouseMode() == MainFrame::MODE_MOUSE_ADDITIVE_CONSTRAINT ||
				parent->mouseMode() == MainFrame::MODE_MOUSE_SUBTRACTIVE_CONSTRAINT) &&
				len( selectingAnchor_ - selectingFree_ ) > 15.0 )
		{
			SimulationTreePtr currentTree = parent->getTree();
			/*
			if( querying_ && updateQuadTreeThread_ )
				updateQuadTreeThread_->join();
				*/

			vl::Mat4f matrix = this->camera_->camera()->projectiveTransform() 
				* this->camera_->camera()->modelViewTransform();

			std::vector<size_t> dynamicObjectIds = 
				parent->objectListToDynamicObjectIds( parent->getSelected(), currentTree->dynamicObjects() );
			// don't know what to do with a constraint on 0 objects:
			if( dynamicObjectIds.empty() )
				return;			
			ConstraintPtr constraint;
			std::string constraintType;
			if( parent->mouseMode() == MainFrame::MODE_MOUSE_ADDITIVE_CONSTRAINT )
			{
				constraint.reset( new PositiveProjectiveConstraint(dynamicObjectIds, matrix, 
					BoundingBox2f( strip(box.minimum()), strip(box.maximum()) ) ) );
				constraintType = "positive";
			}
			else
			{
				constraint.reset( new NegativeProjectiveConstraint(dynamicObjectIds, matrix, 
					BoundingBox2f( strip(box.minimum()), strip(box.maximum()) ) ) );
				constraintType = "negative";
			}

			{
				BOOL result;
				LARGE_INTEGER before;
				result = QueryPerformanceCounter( &before );
				assert( result );

				SimulationTreePtr tree = parent->getTree();
				boost::shared_ptr<Action> addConstraintAction(
					new AddConstraintAction(constraint, tree) );
				parent->performAction( addConstraintAction );

				LARGE_INTEGER after;
				result = QueryPerformanceCounter( &after );
				assert( result );

				LARGE_INTEGER frequency;
				result = QueryPerformanceFrequency( &frequency );
				assert( result );

				double time = boost::numeric_cast<double>( after.QuadPart - before.QuadPart ) 
					/ boost::numeric_cast<double>( frequency.QuadPart );
				std::ostringstream oss;
				oss << "New " << constraintType << " constraint; evaluated in " << (time*1000.0) << "ms";
				parent->log( oss.str() );
			}
		}

		parent->updateViews();
	}
}

void GLView::setupViewport() const
{
	camera_->camera()->setRatio( boost::numeric_cast<double>(viewportWidth_) 
		/ boost::numeric_cast<double>(viewportHeight_) );
	glViewport(0, 0, viewportWidth_, viewportHeight_);
}


float GLView::fps() const
{
	if( this->frameTimes_.size() < 2 )
		return std::numeric_limits<float>::quiet_NaN();

	double elapsed = boost::numeric_cast<double>(frameTimes_.front() - frameTimes_.back()) / 
						boost::numeric_cast<double>(CLOCKS_PER_SEC);
	return boost::numeric_cast<double>(frameTimes_.size()) / elapsed;
}


template <typename T>
class SwapBack
{
public:
	SwapBack( T& first, T& second )
		: first_( first ), second_( second )
	{
		std::swap( first_, second_ );
	}

	~SwapBack()
	{
		std::swap( first_, second_ );
	}

private:
	T& first_;
	T& second_;
};

wxImage GLView::dumpFrame(
	int width,
	int height,
	const std::vector< Simulation::State >& states )
{
	if( states.empty() )
		return wxImage( width, height, true );

	tileRenderer_.reset( new TileRendererWrapper() );
	trTileSize(tileRenderer_->get(), this->viewportWidth_, this->viewportHeight_, 0);
	trImageSize( tileRenderer_->get(), width, height );

	wxImage image( width, height );
	trImageBuffer( tileRenderer_->get(), GL_RGB, GL_UNSIGNED_BYTE, image.GetData() );

	camera().setRatio( boost::numeric_cast<double>(width)
		/ boost::numeric_cast<double>(height) );
	camera().applyProjectiveTransform( tileRenderer_->get() );

	statesToPaint_ = states;

	int more = 1;
	while( more )
	{
		trBeginTile( tileRenderer_->get() );

		Refresh();
		Update();

		more = trEndTile( tileRenderer_->get() );
	}

	tileRenderer_.reset();
	statesToPaint_.clear();

	return image.Mirror(false);
}

void GLView::dumpFrame( 
	int width, 
	int height, 
	const std::string& filename )
{
	//glViewport(0, 0, width, height);

	{
		bool headsUp = false;
		SwapBack<bool> heads( headsUp, this->headsUp_ );

		tileRenderer_.reset( new TileRendererWrapper() );
		trTileSize(tileRenderer_->get(), this->viewportWidth_, this->viewportHeight_, 0);
		trImageSize( tileRenderer_->get(), width, height );

		wxImage image( width, height );
		trImageBuffer( tileRenderer_->get(), GL_RGB, GL_UNSIGNED_BYTE, image.GetData() );

		camera().setRatio( boost::numeric_cast<double>(width)
			/ boost::numeric_cast<double>(height) );
		camera().applyProjectiveTransform( tileRenderer_->get() );

		int more = 1;
		while( more )
		{
			trBeginTile( tileRenderer_->get() );

			Refresh();
			Update();

			more = trEndTile( tileRenderer_->get() );
		}

		image = image.Mirror(false);
		if( !image.SaveFile( wxT(filename.c_str()) ) )
		{
			throw IOException( "Unable to write to image '" + filename + "'." );
		}

		tileRenderer_.reset();
	}
}

void GLView::dumpFrame( 
	int width, 
	int height, 
	BYTE* buffer )
{
	{
		bool headsUp = false;
		SwapBack<bool> heads( headsUp, this->headsUp_ );

		tileRenderer_.reset( new TileRendererWrapper() );
		trTileSize(tileRenderer_->get(), this->viewportWidth_, this->viewportHeight_, 0);
		trImageSize( tileRenderer_->get(), width, height );

		trImageBuffer( tileRenderer_->get(), GL_BGR, GL_UNSIGNED_BYTE, buffer );

		camera().setRatio( boost::numeric_cast<double>(width)
			/ boost::numeric_cast<double>(height) );
		camera().applyProjectiveTransform( tileRenderer_->get() );

		int more = 1;
		while( more )
		{
			trBeginTile( tileRenderer_->get() );

			Refresh();
			Update();

			more = trEndTile( tileRenderer_->get() );
		}

		tileRenderer_.reset();
		SwapBuffers();
	}
}


/*
void GLView::writeRibFile( const std::string& filename )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	std::string prefix( filename );
	std::string::size_type dotPos = filename.rfind(".");
	if( dotPos != std::string::npos )
		prefix = filename.substr(0, dotPos);

	std::string linesFilename = prefix + ".lines";
	std::ofstream ofs( linesFilename.c_str() );
	if( !ofs )
		throw IOException( "Unable to open file '" + linesFilename + "' for writing." );

	if( parent->getTree() )
	{
		std::vector<size_t> dynamicObjectIds = 
			objectListToDynamicObjectIds( selected_, parent->getTree()->dynamicObjects() );
		for( std::vector<size_t>::const_iterator obIter = dynamicObjectIds.begin();
			obIter != dynamicObjectIds.end(); ++obIter )
		{
			std::vector<vl::Vec3f> points;
			assert( lineCache_[*obIter].size() % 2 == 0 );
			for( size_t i = 0; i < this->lineCache_[*obIter].size(); i += 2 )
			{
				vl::Vec3f firstPoint = this->lineCache_[*obIter].at( i );
				vl::Vec3f secondPoint = this->lineCache_[*obIter].at( i+1 );
				if( !points.empty() && firstPoint != points.back() )
				{
					ofs << "\n";
					ofs << firstPoint[0] << " " << firstPoint[1] << " " << firstPoint[2] << " ";
					points.clear();
				}

				ofs << secondPoint[0] << " " << secondPoint[1] << " " << secondPoint[2] << " ";
			}
			ofs << std::endl;
		}
	}
}
*/

CameraTransformConverter::CameraTransformConverter( 
	const Camera& camera, int viewportWidth, int viewportHeight )
	:	projectiveTransform_(toMat4d(camera.projectiveTransform())),
		modelViewTransform_(toMat4d(camera.modelViewTransform())),
		viewportWidth_(viewportWidth), 
		viewportHeight_(viewportHeight)
{
	this->invProjectiveTransform_ = vl::inv( projectiveTransform_ );
	this->invModelViewTransform_ = vl::inv( modelViewTransform_ );
}

vl::Vec3d CameraTransformConverter::toScreenSpace( const vl::Vec3d& worldSpace ) const
{
	vl::Vec3d eyeCoord = vl::xform( modelViewTransform_, worldSpace );
	vl::Vec3d ndCoord = vl::xform( projectiveTransform_, eyeCoord );
	vl::Vec3d result( 
		(ndCoord[0] + 1.0f)*(0.5f*this->viewportWidth_),
		(ndCoord[1] + 1.0f)*(0.5f*this->viewportHeight_), 
		(ndCoord[2] + 1.0f)*(0.5f) );
	return result;
}

vl::Vec3d CameraTransformConverter::toWorldSpace( const vl::Vec3d& screenSpace ) const
{
	vl::Vec3d ndCoord(
		screenSpace[0]/(0.5f*this->viewportWidth_) - 1.0f,
		screenSpace[1]/(0.5f*this->viewportHeight_) - 1.0f,
		2.0f*screenSpace[2] - 1.0f );
	vl::Vec3d eyeCoord = vl::xform( invProjectiveTransform_, ndCoord );
	vl::Vec3d result = vl::xform( invModelViewTransform_, eyeCoord );
	return result;
}

} // namespace planning
