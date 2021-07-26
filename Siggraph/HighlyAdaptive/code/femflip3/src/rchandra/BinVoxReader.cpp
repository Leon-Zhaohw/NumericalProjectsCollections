#include "BinVoxReader.h"
#include <stdlib.h>
//
// This example program reads a .binvox file and writes
// an ASCII version of the same file called "voxels.txt"
//
// 0 = empty voxel
// 1 = filled voxel
// A newline is output after every "dim" voxels (depth = height = width = dim)
//
// Note that this ASCII version is not supported by "viewvox" and "thinvox"
//
// The x-axis is the most significant axis, then the z-axis, then the y-axis.
//

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

typedef unsigned char byte;

static int version;
static int depth, height, width;
static int size;
static byte *voxels = 0;
static float tx, ty, tz;
static float scale;

int read_binvox(string filespec)
{
	ifstream *input = new ifstream(filespec.c_str(), ios::in | ios::binary);
	//
	// read header
	//
	string line;
	*input >> line;  // #binvox
	if (line.compare("#binvox") != 0) {
		cout << "Error: first line reads [" << line << "] instead of [#binvox]" << endl;
		delete input;
		return 0;
	}
	*input >> version;
	cout << "reading binvox version " << version << endl;
	
	depth = -1;
	int done = 0;
	while(input->good() && !done) {
		*input >> line;
		if (line.compare("data") == 0) done = 1;
		else if (line.compare("dim") == 0) {
			*input >> depth >> height >> width;
		}
		else if (line.compare("translate") == 0) {
			*input >> tx >> ty >> tz;
		}
		else if (line.compare("scale") == 0) {
			*input >> scale;
		}
		else {
			cout << "  unrecognized keyword [" << line << "], skipping" << endl;
			char c;
			do {  // skip until end of line
				c = input->get();
			} while(input->good() && (c != '\n'));
			
		}
	}
	if (!done) {
		cout << "  error reading header" << endl;
		return 0;
	}
	if (depth == -1) {
		cout << "  missing dimensions in header" << endl;
		return 0;
	}
	
	size = width * height * depth;
	voxels = new byte[size];
	if (!voxels) {
		cout << "  error allocating memory" << endl;
		return 0;
	}
	
	//
	// read voxel data
	//
	byte value;
	byte count;
	int index = 0;
	int end_index = 0;
	int nr_voxels = 0;
	
	input->unsetf(ios::skipws);  // need to read every byte now (!)
	*input >> value;  // read the linefeed char
	
	while((end_index < size) && input->good()) {
		*input >> value >> count;
		
		if (input->good()) {
			end_index = index + count;
			if (end_index > size) return 0;
			for(int i=index; i < end_index; i++) voxels[i] = value;
			
			if (value) nr_voxels += count;
			index = end_index;
		}  // if file still ok
		
	}  // while
	
	input->close();
	cout << "  read " << nr_voxels << " voxels" << endl;
	
	return 1;
	
}

unsigned char* loadBinVox( const char *path, FLOAT64 translate[3], FLOAT64 &sc )
{	
	if (!read_binvox(path)) {
		cout << "Error reading [" << path << "]" << endl << endl;
		exit(1);
	}
	
#if 0
	//
	// now write the data to as ASCII
	//
	ofstream *out = new ofstream("voxels.txt");
	if(!out->good()) {
		cout << "Error opening [voxels.txt]" << endl << endl;
		exit(1);
	}
	
	cout << "Writing voxel data to ASCII file..." << endl;
	
	*out << "#binvox ASCII data" << endl;
	*out << "dim " << depth << " " << height << " " << width << endl;
	*out << "translate " << tx << " " << ty << " " << tz << endl;
	*out << "scale " << scale << endl;
	*out << "data" << endl;
	
	for(int i=0; i < size; i++) {
		*out << (char) (voxels[i] + '0') << " ";
		if (((i + 1) % width) == 0) *out << endl;
	}
	
	out->close();
	
	cout << "done" << endl << endl;
#else
	translate[0] = tx;
	translate[1] = ty;
	translate[2] = tz;
	sc = scale;
	return voxels;
#endif
}
