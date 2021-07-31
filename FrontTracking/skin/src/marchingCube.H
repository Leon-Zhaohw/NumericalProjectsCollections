#ifndef marchingCube_H
#define marchingCube_H

#include <map>
#include "slUtil.H"

class MarchingCell;

class MarchingCube {
public:
	MarchingCube(int nx, int ny, int nz, double h, const SlVector3 &lc);
	void buildTriangleMesh(const SlArray3D<double> &phi, std::vector<SlTri> &triangles, std::vector<SlVector3> &meshPts);
	int buildCellMesh(MarchingCell &cell,const SlArray3D<double> &phi, std::map<std::pair<int,int>,int> &edgeIntersection,
					   std::vector<SlTri> &triangles, std::vector<SlVector3> &meshPts);
	int vertexInterp(const MarchingCell &cell, int a, int b, 
									 std::map<std::pair<int,int>,int> &edgeIntersection, std::vector<SlVector3> &meshPts);
	int buildTriangles(int cubeindex, MarchingCell &cell, std::map<std::pair<int,int>,int> &edgeIntersection, 
										 std::vector<SlTri> &triangles, std::vector<SlVector3> &meshPts);

protected:
	static int edgeTable[256];
	static int triTable[256][16];
	int nx_, ny_, nz_;
	double h_;
	SlVector3 lc_, uc_;
};

class MarchingCell {
public:
	MarchingCell(int nx, int ny, int nz, double hi, const SlVector3 &lc);
	void set(int i, int j, int k);
	int pointIndex(int c) const;
	SlVector3 point(int c) const;
	SlInt3 pointCoord(int c) const;

	int nx_, ny_, nz_, nynz_;
	double h_, h2_;
	SlVector3 lc_;
	int x, y, z;
	double phiCache[8];
};

#endif
