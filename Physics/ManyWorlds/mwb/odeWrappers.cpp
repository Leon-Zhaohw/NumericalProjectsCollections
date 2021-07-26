#include "stdafx.h"
#include "odeWrappers.h"

namespace planning {

ODEWorldWrapper::ODEWorldWrapper()
	: id_(dWorldCreate())
{
}

ODEWorldWrapper::~ODEWorldWrapper()
{
	dWorldDestroy( id_ );
}

dWorldID ODEWorldWrapper::get() const
{
	return id_;
}

ODEBodyWrapper::ODEBodyWrapper( boost::shared_ptr<ODEWorldWrapper> world )
	: id_( dBodyCreate(world->get()) )
{
}

ODEBodyWrapper::~ODEBodyWrapper()
{
	dBodyDestroy( id_ );
}

dBodyID ODEBodyWrapper::get() const
{
	return id_;
}

#if dTRIMESH_ENABLED
ODETriMeshDataWrapper::ODETriMeshDataWrapper( 
	const std::vector<vl::Vec3f>& positions, 
	const std::vector<Triangle>& triangles)
{
	positions_ = positions;
	faceNormals_.reserve( triangles.size() );

	triangles_.reserve( 3*triangles.size() );
	for( std::vector<Triangle>::const_iterator triItr = triangles.begin();
		triItr != triangles.end(); ++triItr )
	{
		Triangle tri = *triItr;
		for( unsigned int i = 0; i < 3; ++i )
			triangles_.push_back( tri[i] );

		faceNormals_.push_back( vl::norm( 
			vl::cross(positions[ tri[1] ] - positions[ tri[0] ],
				positions[ tri[2] ] - positions[ tri[0] ]) ) );
	}
}

ODETriMeshWrapper::ODETriMeshWrapper( boost::shared_ptr<ODETriMeshDataWrapper> data )
	: data_( data )
{
	id_ = dGeomTriMeshDataCreate();

	dGeomTriMeshDataBuildSingle(id_,
		&data_->positions_[0], sizeof( vl::Vec3f ), data_->positions_.size(),
		&data_->triangles_[0], data_->triangles_.size(), 3*sizeof(int)
		/* ,&data_->faceNormals_[0] */ 
		);
}

ODETriMeshWrapper::~ODETriMeshWrapper()
{
	dGeomTriMeshDataDestroy(id_);
}

dTriMeshDataID ODETriMeshWrapper::get() const
{
	return id_;
}
#endif

boost::once_flag initODEOnce = BOOST_ONCE_INIT;
boost::once_flag closeODEOnce = BOOST_ONCE_INIT;

ODEInitializer::ODEInitializer()
{
#if (ODE_VERSION > 7)
	boost::call_once(&dInitODE, initODEOnce);
#endif
}

ODEInitializer::~ODEInitializer()
{
#if (ODE_VERSION > 7)
	boost::call_once(&dCloseODE, initODEOnce);
#endif
}

} // namespace planning


