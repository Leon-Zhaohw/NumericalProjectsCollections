/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
 */
///////////////////////////////////////////////////////////////////////
// This is just a nice interface on Andrew Selle's code:
// http://physbam.stanford.edu/~aselle/code/make_quicktime.cpp
//
// To use, call "addLuminanceFrame" for each frame of your movie
// (it will set dimensions based on the first frame passed in)
// and "writeMovie" to write the MOV out when you're done.
//
// Or, to grab a frame from GL, call "addFrameGL"
///////////////////////////////////////////////////////////////////////

#ifndef QUICKTIME_MOVIE_H
#define QUICKTIME_MOVIE_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cassert>
#include <vector>
#include <jpeglib.h>

// enables OpenGL screengrabs
#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

typedef unsigned int uint;
typedef unsigned short ushort;

class QT_ATOM
{
    FILE *fp;
    long start_offset;
    const char* type;
public:
    QT_ATOM(FILE* fp,const char* type)
        :fp(fp),type(type)
    {
        start_offset=ftell(fp);
        uint dummy;
        fwrite(&dummy,4,1,fp);
        fputs(type,fp);
    }

    ~QT_ATOM()
    {
        uint atom_size=ftell(fp)-start_offset;
        uint atom_size_endian=htonl(atom_size);
        fseek(fp,start_offset,SEEK_SET);
        fwrite(&atom_size_endian,4,1,fp);
        fseek(fp,0,SEEK_END);
    }

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
class QUICKTIME_MOVIE {

public:  
  QUICKTIME_MOVIE() {
    _width = -1;
    _height = -1;
    _totalFrames = 0;
  };

  ~QUICKTIME_MOVIE() {
    for (unsigned int x = 0; x < _frameRows.size(); x++)
      delete[] _frameRows[x];
  };

  ////////////////////////////////////////////////////////////////////////
  // add a new frame to the movie, assuming it is luminance, [0,1]
  ////////////////////////////////////////////////////////////////////////
  void addLuminanceFrame(const float* image, const int& width, const int& height) 
  {
    if (_width == -1 && _height == -1)
    {
      _width = width;
      _height = height;
    }
    assert(width == _width);
    assert(height == _height);

    for (int y = 0; y < _height; y++)
    {
      JSAMPLE* row = new JSAMPLE[3 * _width];

      for (int x = 0; x < _width; x++)
      {
        float sample = image[x + y * _width];
        sample = (sample > 1.0) ? 1.0 : sample;
        sample = (sample < 0.0) ? 0.0 : sample;
       
        unsigned char scaled = (unsigned char)(sample * 255); 
        row[3 * x] = scaled;
        row[3 * x + 1] = scaled;
        row[3 * x + 2] = scaled;
      }

      _frameRows.push_back(row);
    }
    _totalFrames++;
  };

  ////////////////////////////////////////////////////////////////////////
  // Grab the current OpenGL frame and add it to the movie
  //
  // The original screengrab code is from GLVU:
  // http://www.cs.unc.edu/~walk/software/glvu/
  ////////////////////////////////////////////////////////////////////////
  void addFrameGL()
  {
    GLint OldReadBuffer;
    glGetIntegerv(GL_READ_BUFFER,&OldReadBuffer);
    glReadBuffer(GL_FRONT);

    GLint OldPackAlignment;
    glGetIntegerv(GL_PACK_ALIGNMENT,&OldPackAlignment); 
    glPixelStorei(GL_PACK_ALIGNMENT,1);

    // get the screen pixels
    int width = glutGet((GLenum)GLUT_WINDOW_WIDTH);
    int height = glutGet((GLenum)GLUT_WINDOW_HEIGHT);
    int NumPixels = width * height;
    GLubyte* Pixels = new GLubyte[NumPixels*3];
    if (Pixels==NULL) { printf("UNABLE TO ALLOC PIXEL READ ARRAY!\n"); return; }
    glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,Pixels);

    // store the screen pixels
    if (_width == -1 && _height == -1)
    {
      _width = width;
      _height = height;
    }
    assert(width == _width);
    assert(height == _height);
    for (int y = 0; y < _height; y++)
    {
      JSAMPLE* row = new JSAMPLE[3 * _width];

      // invert y, because of pixel ordering
      for (int x = 0; x < _width * 3; x++)
        row[x] = Pixels[x + (height - 1 - y) * (3 * _width)];

      _frameRows.push_back(row);
    }

    // clean up
    delete[] Pixels;

    glPixelStorei(GL_PACK_ALIGNMENT,OldPackAlignment);
    glReadBuffer((GLenum)OldReadBuffer);

    _totalFrames++;
  }

  ////////////////////////////////////////////////////////////////////////
  // dump out the final movie
  ////////////////////////////////////////////////////////////////////////
  void writeMovie(const char* filename)
  {
    std::cout << " Writing movie " << filename << "..."; flush(std::cout);
    unsigned char EndianTest[2]={0,1};
    big_endian=*(short*)EndianTest==1;
    //std::cout<<"big endian is "<<big_endian<<std::endl;

    FILE *fp=fopen(filename,"w");

    const int frames_per_second=30;

    const int frames = _totalFrames;
    const int width = _width;
    const int height = _height;

    std::vector<int> samplesizes;
    std::vector<int> offsets;

    // Write the samples (i.e. the mdat part in the quicktime)
    {long mdat_begin=ftell(fp);
    QT_ATOM mdat(fp,"mdat");

        struct jpeg_compress_struct cinfo;
        struct jpeg_error_mgr jerr;
     
        // this is where image data is set
        int currentFrameRow = 0;
        for(int i = 0; i < frames; i++){
          long initial_pos=ftell(fp);
          offsets.push_back(initial_pos-mdat_begin);
          cinfo.err=jpeg_std_error(&jerr);

          jpeg_create_compress(&cinfo);
          jpeg_stdio_dest(&cinfo,fp);
          
          cinfo.image_width=width;
          cinfo.image_height=height;
          cinfo.input_components=3;
          cinfo.in_color_space=JCS_RGB;
          jpeg_set_defaults(&cinfo);

          jpeg_set_quality(&cinfo,95,TRUE);
          jpeg_start_compress(&cinfo,TRUE);
          
          //int row_stride = cinfo.image_width * 3;
          //JSAMPLE* row = new JSAMPLE[row_stride];
          //unsigned char rowi=i;
          while(cinfo.next_scanline < cinfo.image_height)
          {
            JSAMPROW row_pointer[]={_frameRows[currentFrameRow]};
            /*
            int j;
            for(j = 0; j < row_stride; j++) 
              row[j] = rowi;
            rowi = (rowi+1) % 256;
            */
            jpeg_write_scanlines(&cinfo,row_pointer,1);
            currentFrameRow++;
          }
          jpeg_finish_compress(&cinfo);
          jpeg_destroy_compress(&cinfo);
          samplesizes.push_back(ftell(fp)-initial_pos);
        }
    }

    // Write the header
    {QT_ATOM a(fp,"moov");
        {QT_ATOM a(fp,"mvhd");
            char c=0;
            Write(fp,c); // version
            Write(fp,c);Write(fp,c);Write(fp,c); // reserved
            Write(fp,(uint)0); // creation time
            Write(fp,(uint)0); // modification time
            Write(fp,(uint)frames_per_second); // time rate 1/30th second
            Write(fp,(uint)frames); // duration
            Write(fp,(uint)0x10000); // preferred rate 16bit fixed pt 1.0
            Write(fp,(ushort)0x100); // full volume
            Write(fp,(uint)0);Write(fp,(uint)0);Write(fp,(ushort)0x0); // 10 bytes padded
            Write_Identity_Matrix(fp);
            Write(fp,(uint)0); // preview time
            Write(fp,(uint)0); // preview duration
            Write(fp,(uint)0); // poster time 
            Write(fp,(uint)0); // selection time
            Write(fp,(uint)frames); // selection duration
            Write(fp,(uint)0); // current time
            Write(fp,(uint)2);} // next track
        {QT_ATOM a(fp,"trak");
            {QT_ATOM a(fp,"tkhd");
                Write(fp,(uint)0xf); // flag visibble
                Write(fp,(uint)0); // creation time
                Write(fp,(uint)0); // modification time
                Write(fp,(uint)1); // track id
                Write(fp,(uint)0); // reserved
                Write(fp,(uint)frames); // duration
                Write(fp,(uint)0);Write(fp,(uint)0); // reserved
                Write(fp,(ushort)0); // layer
                Write(fp,(ushort)0); // alternative group
                Write(fp,(ushort)0x100); // volume
                Write(fp,(ushort)0); // reserved
                Write_Identity_Matrix(fp);
                Write(fp,(uint)width<<16); // width
                Write(fp,(uint)height<<16);} // height
            {QT_ATOM a(fp,"edts");
              {QT_ATOM a(fp,"elst");
                Write(fp,(uint)0); // version flags
                Write(fp,(uint)0);}} // 1 entry

            {QT_ATOM a(fp,"mdia");
                {QT_ATOM a(fp,"mdhd");
                    Write(fp,(uint)0x0); // version/flag visibble
                    Write(fp,(uint)0); // creation time
                    Write(fp,(uint)0); // modified time
                    Write(fp,(uint)frames_per_second); // time scale
                    Write(fp,(uint)frames); // duration
                    Write(fp,(ushort)0); // english language
                    Write(fp,(ushort)0xffff);} // quality
                {QT_ATOM a(fp,"hdlr");
                    Write(fp,(uint)0x0); // version/flags
                    fputs("mhlrvide",fp);
                    Write(fp,(uint)0); // component manufacture
                    Write(fp,(uint)0); // component flags
                    Write(fp,(uint)0); // component flags mask
                    Write(fp,(char)0); // component name
                    fputs("Linux Video Media Handler",fp);}
                {QT_ATOM a(fp,"minf");
                  {QT_ATOM a(fp,"vmhd");
                    Write(fp,(uint)0x0001); // version/flags set 1 for compatibility
                    Write(fp,(ushort)0x40); // graphics mode copy
                    Write(fp,(ushort)0x8000); // unused graphics mode opcolor
                    Write(fp,(ushort)0x8000); // unused graphics mode opcolor
                    Write(fp,(ushort)0x8000);} // unused graphics mode opcolor
                  {QT_ATOM a(fp,"hdlr");
                    Write(fp,(uint)0x0); // version/flags
                    fputs("dhlralis",fp);
                    Write(fp,(uint)0); // component manufacture
                    Write(fp,(uint)0); // component flags
                    Write(fp,(uint)0); // component flags mask
                    Write(fp,(char)0); // component name
                    fputs("Linux Alias Data Handler",fp);}
                  {QT_ATOM a(fp,"dinf");
                    {QT_ATOM a(fp,"dref");
                      Write(fp,(uint)0x0); // vvvf version flags
                      Write(fp,(uint)0x1); // 1 entry
                      {QT_ATOM a(fp,"alis");
                        Write(fp,(uint)1);}}}
                  {QT_ATOM a(fp,"stbl");
                    {QT_ATOM a(fp,"stsd");
                      Write(fp,(uint)0); // version and flags
                      Write(fp,(uint)1); // 1 entry
                      {QT_ATOM a(fp,"jpeg");
                        Write(fp,(uint)0); //reserved
                        Write(fp,(ushort)0); //reserved
                        Write(fp,(ushort)1); // data reference index
                        // write video specific data
                        Write(fp,(ushort)0); // version
                        Write(fp,(ushort)0); // revision level
                        fputs("lnux",fp); // vendor
                        Write(fp,(uint)100); //temporal quality (max)
                        Write(fp,(uint)258); //spatial quality (max)
                        Write(fp,(ushort)width); // width of image
                        Write(fp,(ushort)height); // height of image
                        Write(fp,(uint)0x00480000); // height of image (72dpi)
                        Write(fp,(uint)0x00480000); // height of image (72dpi)
                        Write(fp,(uint)0); // data size (must be zero)
                        Write(fp,(ushort)1); // frames per sample (usually 1)
                        const char* descript="Quicktime for Linux";
                        Write(fp,(char)strlen(descript)); 
                        fputs(descript,fp); // compressor
                        for(unsigned int i=0;i<32-strlen(descript)-1;i++) Write(fp,(char)0);
                        Write(fp,(ushort)24); // color depth
                        Write(fp,(ushort)(short)-1); // use default color table id
                      }
                    }
                  
                    {QT_ATOM a(fp,"stts");
                      Write(fp,(uint)0); // version and flags
                      Write(fp,(uint)1); // 1 entry
                      Write(fp,(uint)frames); // all frames have same duration
                      Write(fp,(uint)1); // duration is one time unit
                    }
                    {QT_ATOM a(fp,"stsc");
                      Write(fp,(uint)0); // version and flags
                      Write(fp,(uint)1); // 1 entry
                      Write(fp,(uint)1); // first sample to use in chunk
                      Write(fp,(uint)1); // number of samples per chunk
                      Write(fp,(uint)1); // index of descriptor (points to stsd)
                    }
                    {QT_ATOM a(fp,"stsz");
                      Write(fp,(uint)0); // version and flags
                      Write(fp,(uint)0); // sample size (non-uniform so zero and table follows)
                      Write(fp,(uint)frames); // one entry per frame
                      for(unsigned int i=0;i<samplesizes.size();i++)  Write(fp,(uint)samplesizes[i]);
                    }
                    {QT_ATOM a(fp,"stco");
                      Write(fp,(uint)0); // version and flags
                      Write(fp,(uint)frames); // one entry per frame
                      for(unsigned int i=0;i<samplesizes.size();i++)  Write(fp,(uint)offsets[i]); // offset from begin of file
                    }
                  }
                }
            }
        }
    }
    fclose(fp);
    std::cout << " done." << std::endl;
  }

private:
  // video dimensions
  int _width;
  int _height;
  std::vector<JSAMPLE*> _frameRows;
  int _totalFrames;

  bool big_endian;

  template<class T>
  inline void Swap_Endianity(T& x)
  { 
    assert(sizeof(T)<=8);
    if (big_endian) return;
    if (sizeof(T) > 1) {
      T old = x;
      for (unsigned int k = 1; k <= sizeof(T); k++) 
        ((char*)&x)[k-1]=((char*)&old)[sizeof(T)-k];
    }
  }

  void Write(FILE* fp, char num)
  {
    fwrite(&num,sizeof(num),1,fp);
  }

  void Write(FILE* fp, ushort num)
  {
    Swap_Endianity(num);
    fwrite(&num,sizeof(num),1,fp);
  }

  void Write(FILE* fp, uint num)
  {
    Swap_Endianity(num);
    fwrite(&num,sizeof(num),1,fp);
  }

  void Write_Identity_Matrix(FILE* fp)
  {
    Write(fp,(uint)0x10000);Write(fp,(uint)0x00000);Write(fp,(uint)0); // 16.16 fixed pt
    Write(fp,(uint)0x00000);Write(fp,(uint)0x10000);Write(fp,(uint)0); // 16.16 fixed pt
    Write(fp,(uint)0x00000);Write(fp,(uint)0x00000);Write(fp,(uint)0x40000000); // 2.30 fixed pt
  }

};

#endif
