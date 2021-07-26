#ifndef __ODEWRAPPERS_H__

#define __ODEWRAPPERS_H__

#include <ode/ode.h>

#include "twigg/primitive.h"

#if (ODE_VERSION > 5)
#define dCreateCCylinder dCreateCapsule
#define dMassSetCappedCylinder dMassSetCapsule
#endif

namespace planning {

class ODEWorldWrapper
{
public:
	ODEWorldWrapper();
	~ODEWorldWrapper();

	dWorldID get() const;
	
private:
	dWorldID id_;
};

typedef boost::shared_ptr<ODEWorldWrapper> ODEWorldWrapperPtr;

class ODEBodyWrapper
{
public:
	ODEBodyWrapper(ODEWorldWrapperPtr world);
	~ODEBodyWrapper();

	dBodyID get() const;

private:
	dBodyID id_;
	ODEWorldWrapperPtr world_;
};

#if dTRIMESH_ENABLED
class ODETriMeshDataWrapper
{
public:
	ODETriMeshDataWrapper( const std::vector<vl::Vec3f>& positions, const std::vector<Triangle>& triangles);

private:
	std::vector<vl::Vec3f> positions_;
	std::vector<vl::Vec3f> faceNormals_;
	std::vector<int> triangles_;

	friend class ODETriMeshWrapper;
};

class ODETriMeshWrapper
{
public:
	ODETriMeshWrapper( boost::shared_ptr<ODETriMeshDataWrapper> data );
	~ODETriMeshWrapper();

	dTriMeshDataID get() const;

private:
	dTriMeshDataID id_;
	boost::shared_ptr<ODETriMeshDataWrapper> data_;
};

typedef boost::shared_ptr<ODETriMeshDataWrapper> ODETriMeshDataWrapperPtr;
typedef boost::shared_ptr<ODETriMeshWrapper> ODETriMeshWrapperPtr;
#endif

class ODEInitializer
{
public:
	ODEInitializer();
	~ODEInitializer();
};

} // namespace planning

#endif

