#ifdef WIN32
#pragma once
#endif

struct MassProperties
{
	double mass;
	vl::Vec3d centerOfMass;
	vl::Mat3d inertiaTensor;
};

class ObjFile;


// WARNING: inertia tensor is computed in the center of mass frame
//   rather than world frame.  Use e.g. dMassTranslate to correct this
MassProperties computeMassProperties( const ObjFile& objFile, 
									 const double density, 
									 const vl::Vec3d& scale = vl::Vec3d(vl::vl_1) );

