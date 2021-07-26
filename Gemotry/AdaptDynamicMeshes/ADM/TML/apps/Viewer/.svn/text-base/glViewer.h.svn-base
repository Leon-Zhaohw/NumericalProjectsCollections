#ifndef _GL_VIEWER_
#define _GL_VIEWER_ 1

#ifndef GL_EXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

#include <FL/Fl.h>
#include <FL/Fl_Gl_Window.h>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

//#include <Magick++.h>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>

class glViewer : public Fl_Gl_Window
{
 protected:
  GLuint displaylist_, mode_, nmodes_;
  
  float cx_, cy_, cz_;
  GLfloat rot_[16];
  GLfloat scale_, zoom_, angle_;
  bool xrot_;
  
  bool first_, saving_;
  int  animating_, counter_;
  
 public:
  glViewer(int width, int height, const char* label)
    : Fl_Gl_Window(0,0,width,height,label),
    mode_(0), nmodes_(1), 
    cx_(0.), cy_(0.), cz_(0.), scale_(1.), 
    zoom_(-4.), angle_(0.), xrot_(false), 
    first_(true), saving_(false), animating_(0), counter_(0)
    {
      for (int i=0; i<16; i++) rot_[i] = (i%5==0)? 1 : 0;
      resizable(this);
    }
  
  virtual ~glViewer() { glDeleteLists(displaylist_, nmodes_); }
  
  /*************************************************************/
  
  virtual void DrawModel()
    { glCallList(displaylist_+mode_); }

  virtual void draw()
    {
      if (!valid())
	{
	  GLfloat aspect = (GLfloat)w()/(GLfloat)h();
	  glViewport(0,0,w(),h());
	  glMatrixMode(GL_PROJECTION);
	  glLoadIdentity();
	  gluPerspective(45.0,aspect,1.0,10.0);
	}
  
      if (first_)
	{
	  SetPosition();
	  InitializeGL();
	  CreateDisplayLists();
	  first_ = false;
	}
      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

      glPushMatrix();
      
      RigidTransform();
      DrawModel();
      
      UpdateRot();
      glPopMatrix();

      glFlush();
      
      if (animating_>0 && saving_) saveASpng();
    }
  
  virtual int handle(int e)
    {
      switch (e)
	{
	case FL_KEYBOARD:
	  switch (Fl::event_key())
	    {
	    case 'h':
	      KeyboardOptions();
	      break;
	      //SAVING
	    case 'i':
	      saving_ = !saving_;
	      break;
	    case 'o':
	      saveASpng();
	      break;
	    case 'p':
	      saveASoff();
	      break;
	      //RESET
	    case 'u':
	      for (int i=0; i<16; i++) rot_[i] = (i%5==0)? 1 : 0;
	      redraw();
	      break;
	      //ANIMATION BUTTONS
	    case ',':
	      animating_ = 0;
	      break;
	    case '.':
	      animating_ = 1;
	      break;	 
	      //MOVE SCENE
	    case FL_Down:
	      angle_ = 2;
	      xrot_ = true;
	      redraw(); 
	      break;
	    case FL_Up:
	      angle_ = -2;
	      xrot_ = true;
	      redraw(); 
	      break;
	    case FL_Left:
	      angle_ = -2;
	      xrot_ = false;
	      redraw(); 
	      break;
	    case FL_Right:
	      angle_ = 2;
	      xrot_ = false;
	      redraw(); 
	      break;
	    case FL_Page_Up:   
	      zoom_ += 0.05; 
	      redraw(); 
	      break;
	    case FL_Page_Down: 
	      zoom_ -= 0.05; 
	      redraw(); 
	      break;
	      
	    default: return Fl_Gl_Window::handle(e);
	    };
	  break;
	case FL_RELEASE: return 1;
	case FL_PUSH:    return 1;
	case FL_DRAG:    return 1;
	default: return Fl_Gl_Window::handle(e);
	};
      return 1;
    }
  
  virtual void saveASoff() = 0;
  
  void saveASpng()
    {
      /*
      char filename[256];
      sprintf(filename,"%s-m%.1d-n%.3d.jpg",label(),mode_,counter_);
      counter_++;
      
      float* pixels = new float[3*w()*h()];
      glReadPixels(0,0,w(),h(), GL_RGB,GL_FLOAT, pixels);
      
      Magick::Image img(w(),h(), "RGB",Magick::FloatPixel, pixels);
      img.flip();
      
      img.write(filename);
      
      std::cout << "Saving .. " << filename << std::endl;
      
      delete[] pixels;
      */
    }

  virtual void KeyboardOptions()
    {
      std::cout << "\t h - help"                      << std::endl;
      std::cout << "\t , - stop animation"            << std::endl;
      std::cout << "\t . - play animation"            << std::endl;
      std::cout << "\t i - saving animation"          << std::endl;
      std::cout << "\t o - save screenshot"           << std::endl;
      std::cout << "\t p - save model as off"         << std::endl;
      std::cout << "\t u - reset rotation"            << std::endl;
      std::cout << "     @ Arrow to rotate and page_up(down) to z-translation @" << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    }
  
 protected:
  virtual void InitializeGL()
    {
      glClearColor(1.,1.,1.,0.);
            
      glShadeModel(GL_SMOOTH);
      
      glEnable(GL_DEPTH_TEST);

      glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
      
      glLightf(GL_LIGHT0,GL_LINEAR_ATTENUATION,1.0);
    }
  
  virtual void UpdateRot()
    {
      angle_ = 0;
      glGetFloatv(GL_MODELVIEW_MATRIX,rot_);
      rot_[12] = rot_[13] = rot_[14] = 0;
    }
  
  virtual void RigidTransform()
    {
      glTranslatef(0,0,zoom_);
      if (angle_!=0.)
	if (xrot_) glRotatef(angle_,1.0,0.0,0.0);
	else       glRotatef(angle_,0.0,1.0,0.0);
      glMultMatrixf(rot_);
    }

  virtual void Centering()
    {
      glScalef(scale_,scale_,scale_);
      glTranslatef(-cx_,-cy_,-cz_);
    }
  
  virtual void CreateDisplayLists() = 0;
  
  virtual void SetPosition()        = 0;
  
  virtual void Color(float t, float min, float zero, float max)
    {
      float r, g, b, alpha;
      if (std::fabs(t-zero)<=1e-12)
	{ r = 0.; g = 1.; b = 0.; }
      else 
	if (t>zero)
	  {
	    // green -> yellow -> red (bezier)
	    alpha = (t - zero) / (max - zero);
	    r = alpha*(2 - alpha);
	    g = 1. - alpha*alpha;
	    b = 0.;
	  } else {
	    // blue -> cyan -> green
	    alpha = (t - min) / (zero - min);
	    r = 0.;
	    g = alpha*(2 - alpha);
	    b = 1. - alpha*alpha;
	  }
      
      glColor3f(r,g,b);
    }
};

#endif
