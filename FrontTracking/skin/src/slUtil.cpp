// Copyright (c) 2011, Regents of the University of Utah
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "slUtil.H"
#include <iostream>
#include <fstream>
#include <float.h>
#include "slUtil.H"

using namespace std;

double pow2(int p) {
  switch (p) {
	case -20: return 9.53674e-07;
	case -19: return 1.90735e-06;
	case -18: return 3.8147e-06;
	case -17: return 7.62939e-06;
	case -16: return 1.52588e-05;
	case -15: return 3.05176e-05;
	case -14: return 6.10352e-05;
	case -13: return 0.0001220703125;
	case -12: return 0.000244140625;
	case -11: return 0.00048828125;
	case -10: return 0.0009765625;
	case -9: return 0.001953125;
	case -8: return 0.00390625;
	case -7: return 0.0078125;
	case -6: return 0.015625;
	case -5: return 0.03125;
	case -4: return 0.0625;
	case -3: return 0.125;
	case -2: return 0.25;
	case -1: return 0.5;
	case 0: return 1;
	case 1: return 2;
	case 2: return 4;
	case 3: return 8;
	case 4: return 16;
	case 5: return 32;
	case 6: return 64;
	case 7: return 128;
	case 8: return 256;
	case 9: return 512;
	case 10: return 1024;
	case 11: return 2048;
	case 12: return 4096;
	case 13: return 8192;
	case 14: return 16384;
	case 15: return 32768;
	case 16: return 65536;
	case 17: return 131072;
	case 18: return 262144;
	case 19: return 524288;
	case 20: return 1048576;
  default:
    double ret = 1;
    if (abs(p) == p)
      for (int i=0; i<abs(p); i++)
				ret *= 2.0;
    else
      for (int i=0; i<abs(p); i++)
				ret /= 2.0;
    return ret;
  }
}


istream &operator>>(std::istream &strm, std::vector<int> &l) {
  unsigned int j,s = 0;
	std::vector<int>::iterator i;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> s;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.resize(s);
    for(i=l.begin(),j=0;i!=l.end();i++,j++) {
      strm >> (*i);
      if (j != (s-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const std::vector<int> &l) {
  strm << "[";
  strm << l.size();
  strm << ":";
  unsigned int j;
	std::vector<int>::const_iterator i;
  for(i=l.begin(),j=0;i!=l.end();i++,j++) {
    strm << *i;
    if (j != (l.size()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

void sl::read(ifstream &in, std::vector<int> &l) {
	unsigned int s;
	std::vector<int>::iterator i;
	in.read((char *)&s, sizeof(int));
	l.resize(s);
	for (i=l.begin(); i!=l.end(); i++) {
		in.read((char*)(&(*i)), sizeof(int));
	}
}

void sl::write(ofstream &out, const std::vector<int> &l) {
	unsigned int s=l.size();
	std::vector<int>::const_iterator i;
	out.write((char*)&s, sizeof(int));
	for (i=l.begin(); i!=l.end(); i++) {
		out.write((char*) (&(*i)), sizeof(int));
	}
}

void randomize(std::vector<int> &list) {
	int i=list.size(), j, tmp;

	for (i=0; i>0;) {
		j = rand()%(i--);
		tmp = list[i];
		list[i] = list[j];
		list[j] = tmp;
	}
		
}

//////////////////////////////////////////////////////////////////
// doubleList
//////////////////////////////////////////////////////////////////
istream &operator>>(std::istream &strm, std::vector<double> &l) {
	unsigned int j,s=0;
	std::vector<double>::iterator i;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);
	
  eatStr("[",strm);
  strm >> s;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.resize(s);
    for(i=l.begin(),j=0;i!=l.end();i++,j++) {
      strm >> *i;
      if (j != (s-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const std::vector<double> &l) {
  unsigned int j;
	std::vector<double>::const_iterator i;
  strm << "[";
  strm << l.size();
  strm << ":";
  for(i=l.begin(),j=0;i!=l.end();i++,j++) {
    strm << *i;
    if (j != (l.size()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

void sl::read(ifstream &in, std::vector<double> &l) {
	unsigned int s;
	std::vector<double>::iterator i;
	in.read((char *)&s, sizeof(int));
	l.resize(s);
	for (i=l.begin(); i!=l.end(); i++) {
		in.read((char*)(&(*i)), sizeof(double));
	}
}

void sl::write(ofstream &out, const std::vector<double> &l) {
	unsigned int s=l.size();
	std::vector<double>::const_iterator i;
	out.write((char *)&s, sizeof(int));
	for (i=l.begin(); i!=l.end(); i++) {
		out.write((char*) (&(*i)), sizeof(double));
	}
}


//////////////////////////////////////////////////////////////////
// SlTri
//////////////////////////////////////////////////////////////////
std::istream &operator>>(std::istream &strm, SlTri &t) {
  ios::fmtflags orgFlags = strm.setf(ios::skipws);
  eatStr("[",strm);
  strm >> t[0];
  eatChar(',',strm);
  strm >> t[1];
  eatChar(',',strm);
  strm >> t[2];
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

std::ostream &operator<<(std::ostream &strm, const SlTri &t) {
  strm << "[";
  strm << t[0];
  strm << ",";
  strm << t[1];
  strm << ",";
  strm << t[2];
  strm << "]";
  return strm;
}
void sl::read(ifstream &in, SlTri &t) {
	in.read((char*) &(t[0]), sizeof(int));
	in.read((char*) &(t[1]), sizeof(int));
	in.read((char*) &(t[2]), sizeof(int));
}
void sl::write(ofstream &out, const SlTri &_t) {
	SlTri t(_t);
	out.write((char*) &(t[0]), sizeof(int));
	out.write((char*) &(t[1]), sizeof(int));
	out.write((char*) &(t[2]), sizeof(int));
}

//////////////////////////////////////////////////////////////////
// SlTriList
//////////////////////////////////////////////////////////////////
istream &operator>>(std::istream &strm, std::vector<SlTri> &l) {
  unsigned int j,s = 0;
	std::vector<SlTri>::iterator i;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> s;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.resize(s);
    for(i=l.begin(),j=0;i!=l.end();i++,j++) {
      strm >> *i;
      if (j != (s-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm, const std::vector<SlTri> &l) {
  unsigned int j;
	std::vector<SlTri>::const_iterator i;
  strm << "[";
  strm << l.size();
  strm << ":";
  for(i=l.begin(),j=0;i!=l.end();i++,j++) {
    strm << *i;
    if (j != (l.size()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

void sl::read(ifstream &in, std::vector<SlTri> &l) {
	unsigned int s;
	std::vector<SlTri>::iterator i;
	in.read((char *)&s, sizeof(int));
	l.resize(s);
	for (i=l.begin(); i!=l.end(); i++) {
		sl::read(in, *i);
	}
}

void sl::write(ofstream &out, const std::vector<SlTri> &l) {
	unsigned int s=l.size();
	std::vector<SlTri>::const_iterator i;
	out.write((char*)&s, sizeof(int));
	for (i=l.begin(); i!=l.end(); i++) {
		sl::write(out, *i);
	}
}


//////////////////////////////////////////////////////////////////
// SlVecList
//////////////////////////////////////////////////////////////////
istream &operator>>(std::istream &strm, std::vector<SlVector3> &l) {
  unsigned int j,s = 0;
	std::vector<SlVector3>::iterator i;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> s;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.resize(s);
    for(i=l.begin(),j=0;i!=l.end();i++,j++) {
      strm >> *i;
      if (j != (s-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const std::vector<SlVector3> &l) {
  unsigned int j;
	std::vector<SlVector3>::const_iterator i;
  strm << "[";
  strm << l.size();
  strm << ":";
  for(i=l.begin(),j=0;i!=l.end();i++,j++) {
    strm << *i;
    if (j != (l.size()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

void sl::read(ifstream &in, std::vector<SlVector3> &l) {
	int s=0;
	std::vector<SlVector3>::iterator i;
	in.read((char*)&s, sizeof(int));
	l.resize(s);
	for (i=l.begin(); i!=l.end(); i++) {
		in.read((char*)&((*i)[0]), sizeof(double));
		in.read((char*)&((*i)[1]), sizeof(double));
		in.read((char*)&((*i)[2]), sizeof(double));
	}
}

void sl::write(ofstream &out, const std::vector<SlVector3> &l) {
	int s=l.size();
	std::vector<SlVector3>::const_iterator i;
	out.write((char*)&s, sizeof(int));
	for (i=l.begin(); i!=l.end(); i++) {
		SlVector3 x(*i);
		out.write((char*)&(x[0]), sizeof(double));
		out.write((char*)&(x[1]), sizeof(double));
		out.write((char*)&(x[2]), sizeof(double));
	}
}

//////////////////////////////////////////////////////////////////
// SlInt3List
//////////////////////////////////////////////////////////////////
istream &operator>>(std::istream &strm, SlInt3 &v) {
	ios::fmtflags orgFlags = strm.setf(ios::skipws);
  eatChar('[',strm);
  strm >> v[0];
  eatChar(',',strm);
  strm >> v[1];
  eatChar(',',strm);
  strm >> v[2];
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}
  
ostream &operator<<(ostream &strm,const SlInt3 &v) {
  strm << "[";
	strm << v[0]; strm << ",";
	strm << v[1]; strm << ",";
	strm << v[2]; strm << "]";
  return strm;
}

istream &operator>>(std::istream &strm, std::vector<SlInt3> &l) {
  unsigned int j,s = 0;
	std::vector<SlInt3>::iterator i;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> s;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.resize(s);
    for(i=l.begin(),j=0;i!=l.end();i++,j++) {
      strm >> *i;
      if (j != (s-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const std::vector<SlInt3> &l) {
  unsigned int j;
	std::vector<SlInt3>::const_iterator i;
  strm << "[";
  strm << l.size();
  strm << ":";
  for(i=l.begin(),j=0;i!=l.end();i++,j++) {
    strm << *i;
    if (j != (l.size()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

void sl::read(ifstream &in, std::vector<SlInt3> &l) {
	int s;
	std::vector<SlInt3>::iterator i;
	in.read((char*)&s, sizeof(int));
	l.resize(s);
	for (i=l.begin(); i<l.end(); i++) {
		in.read((char*)&((*i)[0]), sizeof(int));
		in.read((char*)&((*i)[1]), sizeof(int));
		in.read((char*)&((*i)[2]), sizeof(int));
	}
}

void sl::write(ofstream &out, const std::vector<SlInt3> &l) {
	int s=l.size();
	std::vector<SlInt3>::const_iterator i;
	out.write((char*)&s, sizeof(int));
	for (i=l.begin(); i<l.end(); i++) {
		SlInt3 x(*i);
		out.write((char*)&(x[0]), sizeof(int));
		out.write((char*)&(x[1]), sizeof(int));
		out.write((char*)&(x[2]), sizeof(int));
	}
}

//////////////////////////////////////////////////////////////////
// SlInt6List
//////////////////////////////////////////////////////////////////
istream &operator>>(std::istream &strm, SlInt6 &v) {
	ios::fmtflags orgFlags = strm.setf(ios::skipws);
  eatChar('[',strm);
  strm >> v[0];
  eatChar(',',strm);
  strm >> v[1];
  eatChar(',',strm);
  strm >> v[2];
  eatChar(',',strm);
  strm >> v[3];
  eatChar(',',strm);
  strm >> v[4];
  eatChar(',',strm);
  strm >> v[5];
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}
  
ostream &operator<<(ostream &strm,const SlInt6 &v) {
  strm << "[";
	strm << v[0]; strm << ",";
	strm << v[1]; strm << ",";
	strm << v[2]; strm << ",";
	strm << v[3]; strm << ",";
	strm << v[4]; strm << ",";
	strm << v[5]; strm << "]";
  return strm;
}

void SlInt6::read(ifstream &in) {
	in.read((char*)&((*this)[0]), sizeof(int));
	in.read((char*)&((*this)[1]), sizeof(int));
	in.read((char*)&((*this)[2]), sizeof(int));
	in.read((char*)&((*this)[3]), sizeof(int));
	in.read((char*)&((*this)[4]), sizeof(int));
	in.read((char*)&((*this)[5]), sizeof(int));
}

void SlInt6::write(std::ofstream &out) const {
	SlInt6 x(*this);
	out.write((char*)&(x[0]), sizeof(int));
	out.write((char*)&(x[1]), sizeof(int));
	out.write((char*)&(x[2]), sizeof(int));
	out.write((char*)&(x[3]), sizeof(int));
	out.write((char*)&(x[4]), sizeof(int));
	out.write((char*)&(x[5]), sizeof(int));
}

istream &operator>>(std::istream &strm, std::vector<SlInt6> &l) {
  unsigned int i,s = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);
	std::vector<SlInt6>::iterator c;

  eatStr("[",strm);
  strm >> s;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.resize(s);
    for(c=l.begin(),i=0;c!=l.end();c++,i++) {
      strm >> (*c);
      if (i != (s-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const std::vector<SlInt6> &l) {
	std::vector<SlInt6>::const_iterator c;
  strm << "[";
  strm << l.size();
  strm << ":";
  unsigned int i;
  for(c=l.begin(),i=0;c!=l.end();c++,i++) {
    strm << (*c);
    if (i != (l.size()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

void sl::read(ifstream &in, std::vector<SlInt6> &l) {
	int s;
	std::vector<SlInt6>::iterator c;
	in.read((char *)&s, sizeof(int));
	l.resize(s);
	for (c=l.begin(); c!=l.end(); c++) (*c).read(in);
}

void sl::write(std::ofstream &out, const std::vector<SlInt6> &l) {
	int s = l.size();
	std::vector<SlInt6>::const_iterator c;
	out.write((char *)&s, sizeof(int));
	for (c=l.begin(); c<l.end(); c++) (*c).write(out);
}

//////////////////////////////////////////////////////////////////
// Useful operations
//////////////////////////////////////////////////////////////////
void computeNormals(const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles, 
										std::vector<SlVector3> &vertexNormals, std::vector<SlVector3> &faceNormals) {
	std::vector<SlVector3>::iterator n;
	std::vector<SlTri>::const_iterator t;
	faceNormals.resize(triangles.size());
	vertexNormals.resize(meshPts.size());
	for (n=vertexNormals.begin(); n!=vertexNormals.end(); n++) {
		(*n) = 0.0;
	}
	for (n=faceNormals.begin(),t=triangles.begin(); n!=faceNormals.end(); n++,t++) {
		(*n) = cross(meshPts[(*t)[1]]-meshPts[(*t)[0]], 
								 meshPts[(*t)[2]]-meshPts[(*t)[0]]);
		vertexNormals[(*t)[0]]+=(*n);
		vertexNormals[(*t)[1]]+=(*n);
		vertexNormals[(*t)[2]]+=(*n);
	}

	for (n=vertexNormals.begin(); n!=vertexNormals.end(); n++)
		normalize(*n);
}

void computeFaceNormals(const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles, 
												std::vector<SlVector3> &faceNormals) {
	std::vector<SlVector3>::iterator n;
	std::vector<SlTri>::const_iterator t;
	faceNormals.resize(triangles.size());
	for (n=faceNormals.begin(),t=triangles.begin(); n!=faceNormals.end(); n++,t++) {
		(*n) = cross(meshPts[(*t)[1]]-meshPts[(*t)[0]], 
								 meshPts[(*t)[2]]-meshPts[(*t)[0]]);
	}
}

void bvtDump(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles) {
	ofstream fout (fname, ios::out | ios::binary);
	int version = 1;
	fout.write((char*)&version, sizeof(int));
	sl::write(fout, meshPts);
	sl::write(fout, triangles);
	fout.close();
}

void bvtRead(char *fname, std::vector<SlVector3> &meshPts, std::vector<SlTri> &triangles) {
	int version;
	ifstream fin (fname, ios::in | ios::binary);
	fin.read((char*)&version, sizeof(int));
	if (version == 1) {
		sl::read(fin, meshPts);
		sl::read(fin, triangles);
	} else {
		cerr << "Incompatible bvt file: "<<fname<<endl;
		exit(-1);
	}
	fin.close();
}

void gtDump(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles) {
	ofstream out;
	std::vector<SlVector3>::const_iterator p,n;
	std::vector<SlTri>::const_iterator t;
	std::vector<SlVector3> vertexNormals;
	std::vector<SlVector3> faceNormals;
	computeNormals(meshPts,triangles,vertexNormals,faceNormals);
	out.open(fname);
	out<<"vertices: "<<meshPts.size()<<std::endl;
	out<<"faces: "<<triangles.size()<<std::endl;

	for (p=meshPts.begin(),n=vertexNormals.begin(); p!=meshPts.end(); p++,n++) 
		out<<"v "<<(*p)[0]<<" "<<(*p)[1]<<" "<<(*p)[2]<<" "<<(*n)[0]<<" "<<(*n)[1]<<" "<<(*n)[2]<<std::endl;
	
	for (t=triangles.begin(); t!=triangles.end(); t++) 
		out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;

	out.close();
}

void bgtDump(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles) {
	ofstream out;
	std::vector<SlVector3>::const_iterator p,n;
	std::vector<SlTri>::const_iterator t;
	std::vector<SlVector3> vertexNormals;
	std::vector<SlVector3> faceNormals;
	computeNormals(meshPts,triangles,vertexNormals,faceNormals);
	out.open(fname, ios::out | ios::binary);
	int foo = meshPts.size();
	out.write((char*)&foo, sizeof(int));
	foo = triangles.size();
	out.write((char*)&foo, sizeof(int));
	double bar;
	for (p=meshPts.begin(),n=vertexNormals.begin(); p!=meshPts.end(); p++,n++) {
		bar = (*p)[0];
		out.write((char*)&bar, sizeof(double));
		bar = (*p)[1];
		out.write((char*)&bar, sizeof(double));
		bar = (*p)[2];
		out.write((char*)&bar, sizeof(double));
		bar = (*n)[0];
		out.write((char*)&bar, sizeof(double));
		bar = (*n)[1];
		out.write((char*)&bar, sizeof(double));
		bar = (*n)[2];
		out.write((char*)&bar, sizeof(double));
	}
	for (t=triangles.begin(); t!=triangles.end(); t++) {
		foo = (*t)[0];
		out.write((char*)&foo, sizeof(int));
		foo = (*t)[1];
		out.write((char*)&foo, sizeof(int));
		foo = (*t)[2];
		out.write((char*)&foo, sizeof(int));
	}
	out.close();
}

void gtRead(char *fname, std::vector<SlVector3> &meshPts, std::vector<SlTri> &triangles) {
	int i,j;
	ifstream in;
	std::vector<SlVector3>::iterator p,n;
	std::vector<SlTri>::iterator t;
	std::vector<SlVector3> vertexNormals;
	in.open(fname);
	eatStr("vertices:", in);
	in>>i;
	meshPts.resize(i);
	vertexNormals.resize(i);
	eatStr("faces:",in);
	in>>j;
	triangles.resize(j);
	for (p=meshPts.begin(),n=vertexNormals.begin(); p!=meshPts.end(); p++,n++) {
		eatChar('v', in);
		in>>(*p)[0]>>(*p)[1]>>(*p)[2]>>(*n)[0]>>(*n)[1]>>(*n)[2];
	}
	for (t=triangles.begin(); t!=triangles.end(); t++) {
		eatChar('f', in);
		in>>(*t)[0]>>(*t)[1]>>(*t)[2];
		(*t)[0]--;
		(*t)[1]--;
		(*t)[2]--;
	}
	in.close();
}

void objDump(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles) {
	ofstream out;
	std::vector<SlVector3>::const_iterator p;
	std::vector<SlTri>::const_iterator t;

	out.open(fname);

	for (p=meshPts.begin(); p!=meshPts.end(); p++) 
		out<<"v "<<(*p)[0]<<" "<<(*p)[1]<<" "<<(*p)[2]<<std::endl;
	
	for (t=triangles.begin(); t!=triangles.end(); t++) 
		out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;

	out.close();
}

bool readObjFile(char *fname, std::vector<SlVector3> &pts, std::vector<SlTri> &triangles)
{
    char c[500];
    SlVector3 lc(DBL_MAX), uc(-DBL_MAX);

    int numVertices=0, numFaces=0;
    bool normals = false, texture = false;
    int tint;
    char ch;
    int p, q, r;
    double x, y, z;
    std::vector<SlVector3>::iterator v;
    std::vector<SlTri>::iterator t;
    std::ifstream in1(fname, std::ios::in);
    if (!in1.is_open()) {
        return false;
    }
    in1.flags(in1.flags() & ~std::ios::skipws);

  while (in1>>ch) {
    if (ch == 'v') {
      in1>>ch;
      if (ch == ' ') numVertices++;
      else if (ch == 'n') normals = true;
      else if (ch == 't') texture = true;
      else std::cerr<<"error \'"<<ch<<"\'"<<std::endl;
    } else if (ch == '#') {
        while (in1 >> ch && ch != '\n'); // Read to the end of the line.
    } else if (ch == 'f') numFaces++;
  }
  in1.close();

    pts.resize(numVertices);
    triangles.resize(numFaces);
    v = pts.begin();
    t = triangles.begin();

    std::ifstream in(fname, std::ios::in);
    if (!in.is_open()) {
        return false;
    }

  while (in>>ch) {
    if (ch == '#') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'g') {
      in.getline(c,500);
      continue;
    }
    if (ch == 's') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'm') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'u') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'v') {
      ch = in.peek();
      if (ch != 't' && ch != 'n') {
				in>>x>>y>>z;
				(*v).set(x,y,z);
				if ((*v)[0] < lc[0]) lc[0] = (*v)[0];
				if ((*v)[1] < lc[1]) lc[1] = (*v)[1];
				if ((*v)[2] < lc[2]) lc[2] = (*v)[2];
				if ((*v)[0] > uc[0]) uc[0] = (*v)[0];
				if ((*v)[1] > uc[1]) uc[1] = (*v)[1];
				if ((*v)[2] > uc[2]) uc[2] = (*v)[2];
				v++;
      } else {
				in.getline(c, 500);
      }
      continue;
    }
    if (ch == 'f') {
      if (normals && texture) {
				in>>p>>ch>>tint>>ch>>tint>>q>>ch>>tint>>ch>>tint>>r>>ch>>tint>>ch>>tint;
      } else if (normals) {
				in>>p>>ch>>ch>>tint>>q>>ch>>ch>>tint>>r>>ch>>ch>>tint;
      } else if (texture) {
				in>>p>>ch>>tint>>q>>ch>>tint>>r>>ch>>tint;
      } else {
				in>>p>>q>>r;
      }
      (*t)[0] = p-1;
      (*t)[1] = q-1;
      (*t)[2] = r-1;
      t++;
      continue;
    }
  }
  in.close();
  return true;
}
