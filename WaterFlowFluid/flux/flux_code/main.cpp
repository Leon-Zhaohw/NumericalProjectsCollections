/*
**	main.cpp
**
**	Created by Anonymouse Authors
**
**	Permission is hereby granted, free of charge, to any person obtaining a copy of
**	this software and associated documentation files (the "Software"), to deal in
**	the Software without restriction, including without limitation the rights to use,
**	copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
**	Software, and to permit persons to whom the Software is furnished to do so,
**	subject to the following conditions:
**
**	The above copyright notice and this permission notice shall be included in all copies
**	or substantial portions of the Software.
**
**	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
**	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
**	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
**	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
**	CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
**	OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//
// How to compile:
//
// macOS: g++ main.cpp -std=c++14 -framework OpenGL -lglfw -Wno-deprecated-declarations
// Linux: g++ main.cpp -std=c++14 -lGL -lglfw -Wno-deprecated-declarations
//
#include "facecutter.h"
#include <GLFW/glfw3.h>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <unistd.h>
//
std::vector<std::vector<point2> > g_polygons;
//
static void process() {
	//
	float area_sum (0.0f);
	auto output_func = [&](const vertex2 *polygon, int count) {
		//
		std::vector<point2> tmp(count);
		for( int i=0; i<count; ++i ) {
			tmp[i] = polygon[i].xi;
		}
		//
		g_polygons.push_back(std::vector<point2>(tmp));
	};
	//
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<float> distr(-0.5,0.5);
	//
	const float scale (1.5);
	point2 offsets[4] = {
		0.5*point2(-scale,-scale),
		0.5*point2(scale,-scale),
		0.5*point2(scale,scale),
		0.5*point2(-scale,scale),
	};
	const float translate = distr(eng);
	for( int i=0; i<4; ++i ) {
		offsets[i] = offsets[i] + point2(distr(eng)+translate,distr(eng)+translate);
	}
	//
	const vertex2 polygon[4] = {
		vertex2(point2(-0.9,-0.9),offsets[0]),
		vertex2(point2(0.9,-0.9),offsets[1]),
		vertex2(point2(0.9,0.9),offsets[2]),
		vertex2(point2(-0.9,0.9),offsets[3])
	};
	//
	g_polygons.clear();
	subdivide(polygon,4,output_func);
}
//
static void draw() {
	//
	glColor3f(0.25,0.25,0.4);
	for( auto polygon : g_polygons ) {
		glBegin(GL_TRIANGLE_FAN);
		for( auto v : polygon ) glVertex2f(v.x[0],v.x[1]);
		glEnd();
	}
	//
	glLineWidth(3.0);
	glColor3f(1.0,1.0,1.0);
	for( auto polygon : g_polygons ) {
		glBegin(GL_LINE_LOOP);
		for( auto v : polygon ) glVertex2f(v.x[0],v.x[1]);
		glEnd();
	}
};
//
int main( int argc, const char *argv[] ) {
	//
	assert(glfwInit());
	GLFWwindow* window = glfwCreateWindow(500,500,"Face Subdividion Example",nullptr,nullptr);
	assert(window);
	glfwMakeContextCurrent(window);
	glClearColor(0.0,0.0,0.0,1.0);
	//
	while (!glfwWindowShouldClose(window)) {
		//
		glClear(GL_COLOR_BUFFER_BIT);
		process();
		draw();
		//
		glfwSwapBuffers(window);
		glfwPollEvents();
		usleep(300000);
	}
	//
	glfwTerminate();
	return 0;
}
//