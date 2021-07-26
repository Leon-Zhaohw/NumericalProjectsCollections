#include "stdafx.h"

#include "twigg/objfile.h"
#include "twigg/volumeIntegrals.h"
#include "twigg/boundingbox.h"

#include <boost/multi_array.hpp>
#include <vector>

/* projection integrals */
struct ProjectionIntegrals
{
	double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;
};

/* face integrals */
struct FaceIntegrals
{
	double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;
};

/* volume integrals */
struct VolumeIntegrals
{
	double T0, T1[3], T2[3], TP[3];
};

struct Indices
{
	int A;   /* alpha */
	int B;   /* beta */
	int C;   /* gamma */
};


/*
   ============================================================================
   data structures
   ============================================================================
*/

struct POLYHEDRON;

struct FACE
{
	FACE( int n )
		: verts(n)
	{
	}

	size_t numVerts() const
	{
		return verts.size();
	}

	boost::array<double, 3> norm;
	double w;
	std::vector<int> verts;
	POLYHEDRON *poly;
};

struct POLYHEDRON
{
	typedef boost::multi_array<double, 2> VertexArray;
	typedef std::vector<FACE> FaceArray;

	POLYHEDRON( size_t numVerts )
		: verts( boost::extents[numVerts][3] )
	{
	}

	FaceArray::size_type numFaces() const
	{
		return faces.size();
	}

	VertexArray::size_type numVerts() const
	{
		return (verts.shape())[0];
	}

	VertexArray verts;
	FaceArray faces;
};

template <typename T>
T SQR(T x)
{
	return (x)*(x);
}

template <typename T>
T CUBE(T x)
{
	return (x)*(x)*(x);
}

/* compute various integrations over projection of face */
ProjectionIntegrals compProjectionIntegrals(FACE *f, Indices indices)
{
  ProjectionIntegrals result;

  double a0, a1, da;
  double b0, b1, db;
  double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
  double a1_2, a1_3, b1_2, b1_3;
  double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
  double Cab, Kab, Caab, Kaab, Cabb, Kabb;
  int i;

  result.P1 = result.Pa = result.Pb = result.Paa = 
	  result.Pab = result.Pbb = result.Paaa = 
	  result.Paab = result.Pabb = result.Pbbb = 0.0;

  for (i = 0; i < f->numVerts(); i++) {
    a0 = f->poly->verts[f->verts[i]][indices.A];
    b0 = f->poly->verts[f->verts[i]][indices.B];
    a1 = f->poly->verts[f->verts[(i+1) % f->numVerts()]][indices.A];
    b1 = f->poly->verts[f->verts[(i+1) % f->numVerts()]][indices.B];
    da = a1 - a0;
    db = b1 - b0;
    a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
    b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
    a1_2 = a1 * a1; a1_3 = a1_2 * a1; 
    b1_2 = b1 * b1; b1_3 = b1_2 * b1;

    C1 = a1 + a0;
    Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
    Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
    Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
    Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
    Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
    Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

    result.P1 += db*C1;
    result.Pa += db*Ca;
    result.Paa += db*Caa;
    result.Paaa += db*Caaa;
    result.Pb += da*Cb;
    result.Pbb += da*Cbb;
    result.Pbbb += da*Cbbb;
    result.Pab += db*(b1*Cab + b0*Kab);
    result.Paab += db*(b1*Caab + b0*Kaab);
    result.Pabb += da*(a1*Cabb + a0*Kabb);
  }

  result.P1 /= 2.0;
  result.Pa /= 6.0;
  result.Paa /= 12.0;
  result.Paaa /= 20.0;
  result.Pb /= -6.0;
  result.Pbb /= -12.0;
  result.Pbbb /= -20.0;
  result.Pab /= 24.0;
  result.Paab /= 60.0;
  result.Pabb /= -60.0;

  return result;
}

FaceIntegrals compFaceIntegrals(FACE *f, Indices indices)
{
	FaceIntegrals result;
  double *n, w;
  double k1, k2, k3, k4;

  ProjectionIntegrals p = compProjectionIntegrals(f, indices);

  w = f->w;
  n = &f->norm[0];
  k1 = 1 / n[indices.C]; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

  result.Fa = k1 * p.Pa;
  result.Fb = k1 * p.Pb;
  result.Fc = -k2 * (n[indices.A]*p.Pa + n[indices.B]*p.Pb + w*p.P1);

  result.Faa = k1 * p.Paa;
  result.Fbb = k1 * p.Pbb;
  result.Fcc = k3 * (SQR(n[indices.A])*p.Paa + 2*n[indices.A]*n[indices.B]*p.Pab 
	 + SQR(n[indices.B])*p.Pbb
	 + w*(2*(n[indices.A]*p.Pa + n[indices.B]*p.Pb) + w*p.P1));

  result.Faaa = k1 * p.Paaa;
  result.Fbbb = k1 * p.Pbbb;
  result.Fccc = -k4 * (CUBE(n[indices.A])*p.Paaa + 3*SQR(n[indices.A])*n[indices.B]*p.Paab 
	   + 3*n[indices.A]*SQR(n[indices.B])*p.Pabb + CUBE(n[indices.B])*p.Pbbb
	   + 3*w*(SQR(n[indices.A])*p.Paa + 2*n[indices.A]*n[indices.B]*p.Pab 
	   + SQR(n[indices.B])*p.Pbb)
	   + w*w*(3*(n[indices.A]*p.Pa + n[indices.B]*p.Pb) + w*p.P1));

  result.Faab = k1 * p.Paab;
  result.Fbbc = -k2 * (n[indices.A]*p.Pabb + n[indices.B]*p.Pbbb + w*p.Pbb);
  result.Fcca = k3 * (SQR(n[indices.A])*p.Paaa 
	 + 2*n[indices.A]*n[indices.B]*p.Paab + SQR(n[indices.B])*p.Pabb
	 + w*(2*(n[indices.A]*p.Paa + n[indices.B]*p.Pab) + w*p.Pa));

  return result;
}

VolumeIntegrals compVolumeIntegrals(POLYHEDRON *p)
{
  const size_t X = 0;
  const size_t Y = 1;
  const size_t Z = 2;

  VolumeIntegrals result;
  FACE *f;
  double nx, ny, nz;
  int i;

  result.T0 = result.T1[X] = result.T1[Y] = result.T1[Z] 
     = result.T2[X] = result.T2[Y] = result.T2[Z] 
     = result.TP[X] = result.TP[Y] = result.TP[Z] = 0;

  for (i = 0; i < p->numFaces(); i++) {
    Indices indices;
    f = &p->faces[i];

    nx = fabs(f->norm[X]);
    ny = fabs(f->norm[Y]);
    nz = fabs(f->norm[Z]);
    if (nx > ny && nx > nz) indices.C = X;
    else indices.C = (ny > nz) ? Y : Z;
    indices.A = (indices.C + 1) % 3;
    indices.B = (indices.A + 1) % 3;

    FaceIntegrals faceIntegrals = compFaceIntegrals(f, indices);

    result.T0 += f->norm[X] * 
		((indices.A == X) ? faceIntegrals.Fa : ((indices.B == X) ? faceIntegrals.Fb : faceIntegrals.Fc));

    result.T1[indices.A] += f->norm[indices.A] * faceIntegrals.Faa;
    result.T1[indices.B] += f->norm[indices.B] * faceIntegrals.Fbb;
    result.T1[indices.C] += f->norm[indices.C] * faceIntegrals.Fcc;
    result.T2[indices.A] += f->norm[indices.A] * faceIntegrals.Faaa;
    result.T2[indices.B] += f->norm[indices.B] * faceIntegrals.Fbbb;
    result.T2[indices.C] += f->norm[indices.C] * faceIntegrals.Fccc;
    result.TP[indices.A] += f->norm[indices.A] * faceIntegrals.Faab;
    result.TP[indices.B] += f->norm[indices.B] * faceIntegrals.Fbbc;
    result.TP[indices.C] += f->norm[indices.C] * faceIntegrals.Fcca;
  }

  result.T1[X] /= 2; result.T1[Y] /= 2; result.T1[Z] /= 2;
  result.T2[X] /= 3; result.T2[Y] /= 3; result.T2[Z] /= 3;
  result.TP[X] /= 2; result.TP[Y] /= 2; result.TP[Z] /= 2;

  return result;
}

MassProperties computeMassProperties( const ObjFile& objFile, const double density, const vl::Vec3d& scale )
{
	const size_t X = 0;
	const size_t Y = 1;
	const size_t Z = 2;

	BoundingBox3d bounds;

	// first vertices
	POLYHEDRON p( objFile.vertexCount() );
	std::vector<vl::Vec3d> vertexPositions = objFile.vertexPositions();
	for( POLYHEDRON::VertexArray::index iVertex = 0; iVertex < objFile.vertexCount(); ++iVertex )
	{
		for( size_t iCoord = 0; iCoord < 3; ++iCoord )
			p.verts[iVertex][iCoord] = (vertexPositions[iVertex])[iCoord] * scale[iCoord];

		bounds.expand( scale*vertexPositions[iVertex] );
	}

	// now faces:
	for( size_t iGroup = 0; iGroup < objFile.numGroups(); ++iGroup )
	{
		ObjFile::Group group = objFile.group( iGroup );
		for( unsigned int iFace = 0; iFace < group.faceCount(); ++iFace )
		{
			ObjFile::Face<unsigned int> face = group.face( iFace );
			FACE f( face.vertexCount() );
			f.poly = &p;
			for( size_t iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
				f.verts[iVertex] = face.vertex(iVertex).position() - 1;

			/* compute face normal and offset w from first 3 vertices */
			double dx1 = p.verts[f.verts[1]][X] - p.verts[f.verts[0]][X];
			double dy1 = p.verts[f.verts[1]][Y] - p.verts[f.verts[0]][Y];
			double dz1 = p.verts[f.verts[1]][Z] - p.verts[f.verts[0]][Z];
			double dx2 = p.verts[f.verts[2]][X] - p.verts[f.verts[1]][X];
			double dy2 = p.verts[f.verts[2]][Y] - p.verts[f.verts[1]][Y];
			double dz2 = p.verts[f.verts[2]][Z] - p.verts[f.verts[1]][Z];
			double nx = dy1 * dz2 - dy2 * dz1;
			double ny = dz1 * dx2 - dz2 * dx1;
			double nz = dx1 * dy2 - dx2 * dy1;
			double len = sqrt(nx * nx + ny * ny + nz * nz);
			f.norm[X] = nx / len;
			f.norm[Y] = ny / len;
			f.norm[Z] = nz / len;
			f.w = - f.norm[X] * p.verts[f.verts[0]][X]
				- f.norm[Y] * p.verts[f.verts[0]][Y]
				- f.norm[Z] * p.verts[f.verts[0]][Z];

			p.faces.push_back( f );
		}
	}

	VolumeIntegrals v = compVolumeIntegrals(&p);

	MassProperties result;

	result.mass = density * v.T0;
	assert( result.mass > 0.0 );

	// compute center of mass
	result.centerOfMass[X] = v.T1[X] / v.T0;
	result.centerOfMass[Y] = v.T1[Y] / v.T0;
	result.centerOfMass[Z] = v.T1[Z] / v.T0;
	assert( contains( bounds, result.centerOfMass ) );

	// compute inertia tensor
	result.inertiaTensor[X][X] = density * (v.T2[Y] + v.T2[Z]);
	result.inertiaTensor[Y][Y] = density * (v.T2[Z] + v.T2[X]);
	result.inertiaTensor[Z][Z] = density * (v.T2[X] + v.T2[Y]);
	result.inertiaTensor[X][Y] = result.inertiaTensor[Y][X] = density * v.TP[X];
	result.inertiaTensor[Y][Z] = result.inertiaTensor[Z][Y] = density * v.TP[Y];
	result.inertiaTensor[Z][X] = result.inertiaTensor[X][Z] = density * v.TP[Z];

	// translate inertia tensor to center of mass
	const vl::Vec3d& r = result.centerOfMass;
	result.inertiaTensor[X][X] -= result.mass * (r[Y]*r[Y] + r[Z]*r[Z]);
	result.inertiaTensor[Y][Y] -= result.mass * (r[Z]*r[Z] + r[X]*r[X]);
	result.inertiaTensor[Z][Z] -= result.mass * (r[X]*r[X] + r[Y]*r[Y]);
	result.inertiaTensor[X][Y] = result.inertiaTensor[Y][X] -= result.mass * r[X] * r[Y]; 
	result.inertiaTensor[Y][Z] = result.inertiaTensor[Z][Y] -= result.mass * r[Y] * r[Z]; 
	result.inertiaTensor[Z][X] = result.inertiaTensor[X][Z] -= result.mass * r[Z] * r[X]; 

	return result;
}

