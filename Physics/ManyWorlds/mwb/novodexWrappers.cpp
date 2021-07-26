#include "stdafx.h"
#include "novodexWrappers.h"

#ifdef USE_NOVODEX
bool isBad( const NxVec3& vec )
{
	return (isBad(vec[0]) || isBad(vec[1]) || isBad(vec[2]));
}

bool isBad( const NxQuat& quat )
{
	return (isBad(quat.x) || isBad(quat.y) || isBad(quat.z) || isBad(quat.w));
}

namespace planning {

boost::recursive_mutex novodexCreateDeleteMutex;

NxVec3 toNxVec( const vl::Vec3d& v )
{
	return NxVec3( v[0], v[1], v[2] );
}

NxVec3 toNxVec( const vl::Vec3f& v )
{
	return NxVec3( v.Ref() );
}

NxQuat toNxQuat( const vl::Vec4d& v )
{
	NxQuat result;
	result.setXYZW( v[0], v[1], v[2], v[3] );
	return result;
}

NxQuat toNxQuat( const vl::Vec4f& v )
{
	NxQuat result;
	result.setXYZW( v.Ref() );
	return result;
}

float defaultSkinWidth()
{
	return 0.01;
}

NxPhysicsSDKPtr instantiateSDK( NxUserAllocator* myAllocator, NxUserOutputStream* myOutputStream )
{
	boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );

	NxPhysicsSDK* sdk = 
		NxCreatePhysicsSDK(
			NX_PHYSICS_SDK_VERSION,
			myAllocator,
			myOutputStream);


	if( sdk == 0 )
	{
		std::ostringstream oss;
		oss << "Unable to instantiate Novodex physics SDK; wrong DLL version?\n";
		oss << "Expected version: " << NX_PHYSICS_SDK_VERSION;
		throw CreationException( oss.str() );
	}

	NxPhysicsSDKPtr result( sdk, 
		boost::mem_fn(&NxPhysicsSDK::release) );
//	result->setParameter( NX_BOUNCE_THRESHOLD, -0.01f );
	result->setParameter( NX_MAX_ANGULAR_VELOCITY, 100.0f );
	result->setParameter( NX_COLL_VETO_JOINTED, 0.0f );
	result->setParameter( NX_SKIN_WIDTH, defaultSkinWidth() );

	return result;
}

NxSceneDeleter::NxSceneDeleter( NxPhysicsSDKPtr ptr )
	: parent_(ptr)
{
}

void NxSceneDeleter::operator()( NxScene* scene )
{
	if( scene == 0 )
		return;

	boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
	parent_->releaseScene( *scene );
}

NxActorDeleter::NxActorDeleter( NxScenePtr ptr )
	: parent_(ptr)
{
}

void NxActorDeleter::operator()( NxActor* actor )
{
	if( actor == 0 )
		return;

	boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
	parent_->releaseActor( *actor );
}

NxJointDeleter::NxJointDeleter( NxScenePtr ptr )
	: parent_(ptr)
{
}

void NxJointDeleter::operator()( NxJoint* joint )
{
	if( joint == 0 )
		return;

	boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
	parent_->releaseJoint( *joint );
}

NxMaterialDeleter::NxMaterialDeleter( NxScenePtr parent )
	: parent_(parent)
{
}

void NxMaterialDeleter::operator()( NxMaterial* mat )
{
	if( mat == 0 )
		return;

	boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
	parent_->releaseMaterial( *mat );
}

NxTriangleMeshDeleter::NxTriangleMeshDeleter( NxPhysicsSDKPtr ptr )
	: parent_(ptr)
{
}

void NxTriangleMeshDeleter::operator()( NxTriangleMesh* mesh )
{
	if( mesh == 0 )
		return;

	boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
	parent_->releaseTriangleMesh( *mesh );
}


NxConvexMeshDeleter::NxConvexMeshDeleter( NxPhysicsSDKPtr ptr )
	: parent_(ptr)
{
}

void NxConvexMeshDeleter::operator()( NxConvexMesh* mesh )
{
	if( mesh == 0 )
		return;

	boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
	parent_->releaseConvexMesh( *mesh );
}

MemoryWriteBuffer::MemoryWriteBuffer() : currentSize(0), maxSize(0), data(NULL)
	{
	}

MemoryWriteBuffer::~MemoryWriteBuffer()
	{
	NX_DELETE_ARRAY(data);
	}

void MemoryWriteBuffer::clear()
	{
	currentSize = 0;
	}

NxStream& MemoryWriteBuffer::storeByte(NxU8 b)
	{
	storeBuffer(&b, sizeof(NxU8));
	return *this;
	}
NxStream& MemoryWriteBuffer::storeWord(NxU16 w)
	{
	storeBuffer(&w, sizeof(NxU16));
	return *this;
	}
NxStream& MemoryWriteBuffer::storeDword(NxU32 d)
	{
	storeBuffer(&d, sizeof(NxU32));
	return *this;
	}
NxStream& MemoryWriteBuffer::storeFloat(NxReal f)
	{
	storeBuffer(&f, sizeof(NxReal));
	return *this;
	}
NxStream& MemoryWriteBuffer::storeDouble(NxF64 f)
	{
	storeBuffer(&f, sizeof(NxF64));
	return *this;
	}
NxStream& MemoryWriteBuffer::storeBuffer(const void* buffer, NxU32 size)
	{
	NxU32 expectedSize = currentSize + size;
	if(expectedSize > maxSize)
		{
		maxSize = expectedSize + 4096;

		NxU8* newData = new NxU8[maxSize];
		NX_ASSERT(newData!=NULL);

		if(data)
			{
			memcpy(newData, data, currentSize);
			delete[] data;
			}
		data = newData;
		}
	memcpy(data+currentSize, buffer, size);
	currentSize += size;
	return *this;
	}


MemoryReadBuffer::MemoryReadBuffer(const NxU8* data) : buffer(data)
	{
	}

MemoryReadBuffer::~MemoryReadBuffer()
	{
	// We don't own the data => no delete
	}

NxU8 MemoryReadBuffer::readByte() const
	{
	NxU8 b;
	memcpy(&b, buffer, sizeof(NxU8));
	buffer += sizeof(NxU8);
	return b;
	}

NxU16 MemoryReadBuffer::readWord() const
	{
	NxU16 w;
	memcpy(&w, buffer, sizeof(NxU16));
	buffer += sizeof(NxU16);
	return w;
	}

NxU32 MemoryReadBuffer::readDword() const
	{
	NxU32 d;
	memcpy(&d, buffer, sizeof(NxU32));
	buffer += sizeof(NxU32);
	return d;
	}

float MemoryReadBuffer::readFloat() const
	{
	float f;
	memcpy(&f, buffer, sizeof(float));
	buffer += sizeof(float);
	return f;
	}

double MemoryReadBuffer::readDouble() const
	{
	double f;
	memcpy(&f, buffer, sizeof(double));
	buffer += sizeof(double);
	return f;
	}

void MemoryReadBuffer::readBuffer(void* dest, NxU32 size) const
	{
	memcpy(dest, buffer, size);
	buffer += size;
	}

} // namespace planning

#endif

