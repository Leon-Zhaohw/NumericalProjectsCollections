/*
 *  write_bmp.cpp
 */

#include "write_bmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct header {
	char bfType[3];
	int bfSize;
	int bfReserved;
	int bfOffBits;
	
	int biSize;
	int biWidth;
	int biHeight;
	short biPlanes;
	short biBitCount;
	int biCompression;
	int biSizeImage;
	int biXPelsPerMeter;
	int biYPelsPerMeter;
	int biClrUsed;
	
	int biClrImportant;
};

template <class T>
static void my_fwrite( T *data, size_t size, size_t n, FILE * fp ) {
	int x=1; // 0x00000001
	if (*(char*)&x)
	{
		// Little Endian
		fwrite( data, size, n, fp );
	}
	else
	{
		//Big Endian
		T value = *data;
		for( int i=0; i<size; i++ )
		{
			int v = (value >> (8*i)) & 0x000000ff;
			putc(v, fp);
		}
	}
}

void write_bmp( const char *filename, unsigned char *image, int width, int height, bool invertY ) {
	
	int bytesPerLine;
	unsigned char *lines;
	
	FILE *fp;
	header bmph;
	
	bytesPerLine = (3 * (width + 1) / 4) * 4;
	
	strcpy(bmph.bfType, "BM");
	bmph.bfOffBits = 54;
	bmph.bfSize = bmph.bfOffBits + bytesPerLine * height;
	bmph.bfReserved = 0;
	bmph.biSize = 40;
	bmph.biWidth = width;
	bmph.biHeight = height;
	bmph.biPlanes = 1;
	bmph.biBitCount = 24;
	bmph.biCompression = 0;
	bmph.biSizeImage = bytesPerLine * height;
	bmph.biXPelsPerMeter = 0;
	bmph.biYPelsPerMeter = 0;
	bmph.biClrUsed = 0;
	bmph.biClrImportant = 0;
	
	fp = fopen (filename, "wb");
	if (! fp ) return;
	
	fwrite(bmph.bfType, 2, 1, fp );
	my_fwrite<int>(&bmph.bfSize, 4, 1, fp );
	my_fwrite<int>(&bmph.bfReserved, 4, 1, fp );
	my_fwrite<int>(&bmph.bfOffBits, 4, 1, fp );
	my_fwrite<int>(&bmph.biSize, 4, 1, fp ) ;
	my_fwrite<int>(&bmph.biWidth, 4, 1, fp );
	my_fwrite<int>(&bmph.biHeight, 4, 1, fp );
	my_fwrite<short>(&bmph.biPlanes, 2, 1, fp );
	my_fwrite<short>(&bmph.biBitCount, 2, 1, fp );
	my_fwrite<int>(&bmph.biCompression, 4, 1, fp );
	my_fwrite<int>(&bmph.biSizeImage, 4, 1, fp );
	my_fwrite<int>(&bmph.biXPelsPerMeter, 4, 1, fp );
	my_fwrite<int>(&bmph.biYPelsPerMeter, 4, 1, fp );
	my_fwrite<int>(&bmph.biClrUsed, 4, 1, fp );
	my_fwrite<int>(&bmph.biClrImportant, 4, 1, fp );
	
	lines = (unsigned char *)malloc(bytesPerLine * height);
	for ( int v=0; v < height; v++ ) {
		int i;
		if( invertY ) i = v;
		else i = height - v - 1;
		
		for ( int j = 0; j < width; j++)
		{
			int pos = 4 * (width * i + j);
			lines[3*j+i*bytesPerLine] = image[pos + 2];
			lines[3*j+i*bytesPerLine+1] = image[pos + 1];
			lines[3*j+i*bytesPerLine+2] = image[pos];
		}
	}
	fwrite(lines, bytesPerLine * height, 1, fp );
	
	free(lines);
	fclose(fp);
	return;
}
