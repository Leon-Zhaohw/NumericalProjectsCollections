#ifndef VISUALIZER_H
#define VISUALIZER_H

// OpenGL
//
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <Eigen/Core>
#include <string>

typedef Eigen::Vector3d Point3d;
typedef Eigen::Vector2d Point2d;



class Visualizer
{

public:
    Visualizer (int w, int h);

    // GLUT event handlers
	//
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mousePress(int button, int state, int x, int y);
    void mouseMove(int x, int y);
	void keyPress(unsigned char k);

    // Camera Control
	//
    
	void resetView(); // show the entire mesh
	
    void updateGL();
    void initMesh();

	Point3d project(const Point3d& p);
	Point3d unproject(double sx, double sy, double sz);
	
    int width;  // window width
    int height; // window height

	void renderText(Point3d p, std::string str, bool _small); // render Text using glut
	std::string convertInt(int number);
	std::string convertDouble(double number);

    // Temporary variables for selection
	//
    bool isSelecting;
    bool isDown;
    Point2d initialPos;
    
    Eigen::Vector2d mousePos;
    
};

#endif

