//////////////////////////////////////////////////////////////////////
// This file is part of Wavelet Turbulence.
// 
// Wavelet Turbulence is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Wavelet Turbulence is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Wavelet Turbulence.  If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2008 Theodore Kim and Nils Thuerey
// 
//////////////////////////////////////////////////////////////////////
//
#ifndef IMAGE_H
#define IMAGE_H

#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <zlib.h>

//////////////////////////////////////////////////////////////////////
// NT helper functions
//////////////////////////////////////////////////////////////////////
template < class T > inline T ABS( T a ) { 
	return (0 < a) ? a : -a ; 
}

template < class T > inline void SWAP_POINTERS( T &a, T &b ) { 
	T temp = a;
	a = b;
	b = temp;
}

template < class T > inline void CLAMP( T &a, T b=0., T c=1.) { 
	if(a<b) { a=b; return; }
	if(a>c) { a=c; return; }
}

template < class T > inline T MIN( T a, T b) {
	return (a < b) ? a : b;
}

template < class T > inline T MAX( T a, T b) {
	return (a > b) ? a : b;
}

template < class T > inline T MAX3( T a, T b, T c) {
	T max = (a > b) ? a : b;
	max = (max > c) ? max : c;
	return max;
}

template < class T > inline float MAX3V( T vec) {
	float max = (vec[0] > vec[1]) ? vec[0] : vec[1];
	max = (max > vec[2]) ? max : vec[2];
	return max;
}

template < class T > inline float MIN3V( T vec) {
	float min = (vec[0] < vec[1]) ? vec[0] : vec[1];
	min = (min < vec[2]) ? min : vec[2];
	return min;
}

//////////////////////////////////////////////////////////////////////
// PNG, POV-Ray, and PBRT output functions
//////////////////////////////////////////////////////////////////////
#include <png.h>

namespace IMAGE {
  static int writePng(const char *fileName, unsigned char **rowsp, int w, int h, bool normalize)
  {
    // defaults 
    const int colortype = PNG_COLOR_TYPE_RGBA;
    const int bitdepth = 8;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_bytep *rows = rowsp;

    FILE *fp = NULL;
    std::string doing = "open for writing";
    if (!(fp = fopen(fileName, "wb"))) goto fail;

    if(!png_ptr) {
      doing = "create png write struct";
      if (!(png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL))) goto fail;
    }
    if(!info_ptr) {
      doing = "create png info struct";
      if (!(info_ptr = png_create_info_struct(png_ptr))) goto fail;
    }

    if (setjmp(png_jmpbuf(png_ptr))) goto fail;
    doing = "init IO";
    png_init_io(png_ptr, fp);
    doing = "write header";
    png_set_IHDR(png_ptr, info_ptr, w, h, bitdepth, colortype, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    doing = "write info";
    png_write_info(png_ptr, info_ptr);
    doing = "write image";
    png_write_image(png_ptr, rows);
    doing = "write end";
    png_write_end(png_ptr, NULL);
    doing = "write destroy structs";
    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose( fp );
    return 0;

  fail:	
    std::cerr << "writePng: could not "<<doing<<" !\n";
    if(fp) fclose( fp );
    if(png_ptr || info_ptr) png_destroy_write_struct(&png_ptr, &info_ptr);
    return -1;
  }

  /////////////////////////////////////////////////////////////////////////////////
  // write a numbered PNG file out, padded with zeros up to three zeros
  /////////////////////////////////////////////////////////////////////////////////
  static void dumpNumberedPNG(int counter, std::string prefix, float* field, int xRes, int yRes)
  {
    char buffer[256];
    sprintf(buffer,"%04i", counter);
    std::string number = std::string(buffer);

    unsigned char pngbuf[xRes*yRes*4];
    unsigned char *rows[yRes];
    float *pfield = field;
    for (int j=0; j<yRes; j++) {
      for (int i=0; i<xRes; i++) {
        float val = *pfield;
        if(val>1.) val=1.;
        if(val<0.) val=0.;
        pngbuf[(j*xRes+i)*4+0] = (unsigned char)(val*255.); 
        pngbuf[(j*xRes+i)*4+1] = (unsigned char)(val*255.); 
        pngbuf[(j*xRes+i)*4+2] = (unsigned char)(val*255.); 
        pfield++;
        pngbuf[(j*xRes+i)*4+3] = 255;
      }
      rows[j] = &pngbuf[(yRes-j-1)*xRes*4];
    }
    std::string filenamePNG = prefix + number + std::string(".png");
    writePng(filenamePNG.c_str(), rows, xRes, yRes, false);
    printf("Writing %s\n", filenamePNG.c_str());
  }
};


#endif
