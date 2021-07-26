#ifdef WIN32
#pragma once
#endif

#ifndef __NOVODEXWRAPPERS_H__
#define __NOVODEXWRAPPERS_H__

#ifdef USE_NOVODEX
#include "NxMat33.h"
#include "NxVec3.h"
#include "NxQuat.h"

#include "NxPhysics.h"
#include "NxStream.h"

bool isBad( const NxVec3& vec );
bool isBad( const NxQuat& quat );

namespace vl
{

inline vl::Vec3f toVec3f( const NxVec3& v )
{
	return vl::Vec3f( v[0], v[1], v[2] );
}

inline vl::Vec3d toVec3d( const NxVec3& v )
{
	return vl::Vec3d( v[0], v[1], v[2] );
}

}

namespace planning {

// Novodex is not thread-safe, especially in its memory management code.
//   So wherever we create or delete something, we need to take out a
//   lock:
extern boost::recursive_mutex novodexCreateDeleteMutex;

NxVec3 toNxVec( const vl::Vec3d& v );
NxVec3 toNxVec( const vl::Vec3f& v );
NxQuat toNxQuat( const vl::Vec4d& v );
NxQuat toNxQuat( const vl::Vec4f& v );

template <typename Real>
NxVec3 toNxVec3( const Real* v )
{
	return NxVec3( v[0], v[1], v[2] );
}

typedef boost::shared_ptr<NxPhysicsSDK> NxPhysicsSDKPtr;
typedef boost::shared_ptr<NxScene> NxScenePtr;
typedef boost::shared_ptr<NxActor> NxActorPtr;
typedef boost::shared_ptr<NxMaterial> NxMaterialPtr;
typedef boost::shared_ptr<NxJoint> NxJointPtr;

NxPhysicsSDKPtr instantiateSDK( NxUserAllocator* myAllocator, NxUserOutputStream* myOutputStream );

float defaultSkinWidth();

class NxSceneDeleter
{
public:
	NxSceneDeleter( NxPhysicsSDKPtr ptr );
	void operator()( NxScene* scene );

private:
	NxPhysicsSDKPtr parent_;
};

class NxActorDeleter
{
public:
	NxActorDeleter( NxScenePtr scene );
	void operator()( NxActor* actor );

private:
	NxScenePtr parent_;
};

class NxJointDeleter
{
public:
	NxJointDeleter( NxScenePtr scene );
	void operator()( NxJoint* actor );

private:
	NxScenePtr parent_;
};

class NxTriangleMeshDeleter
{
public:
	NxTriangleMeshDeleter( NxPhysicsSDKPtr ptr );
	void operator()( NxTriangleMesh* mesh );

private:
	NxPhysicsSDKPtr parent_;
};

class NxConvexMeshDeleter
{
public:
	NxConvexMeshDeleter( NxPhysicsSDKPtr ptr );
	void operator()( NxConvexMesh* mesh );

private:
	NxPhysicsSDKPtr parent_;
};

class NxMaterialDeleter
{
public:
	NxMaterialDeleter( NxScenePtr parent );
	void operator()( NxMaterial* mat );

private:
	NxScenePtr parent_;
};

class MemoryWriteBuffer : public NxStream
	{
	public:
								MemoryWriteBuffer();
	virtual						~MemoryWriteBuffer();
				void			clear();

	virtual		NxU8			readByte()								const	{ NX_ASSERT(0);	return 0;	}
	virtual		NxU16			readWord()								const	{ NX_ASSERT(0);	return 0;	}
	virtual		NxU32			readDword()								const	{ NX_ASSERT(0);	return 0;	}
	virtual		float			readFloat()								const	{ NX_ASSERT(0);	return 0.0f;}
	virtual		double			readDouble()							const	{ NX_ASSERT(0);	return 0.0;	}
	virtual		void			readBuffer(void* buffer, NxU32 size)	const	{ NX_ASSERT(0);				}

	virtual		NxStream&		storeByte(NxU8 b);
	virtual		NxStream&		storeWord(NxU16 w);
	virtual		NxStream&		storeDword(NxU32 d);
	virtual		NxStream&		storeFloat(NxReal f);
	virtual		NxStream&		storeDouble(NxF64 f);
	virtual		NxStream&		storeBuffer(const void* buffer, NxU32 size);

				NxU32			currentSize;
				NxU32			maxSize;
				NxU8*			data;
	};

class MemoryReadBuffer : public NxStream
	{
	public:
								MemoryReadBuffer(const NxU8* data);
	virtual						~MemoryReadBuffer();

	virtual		NxU8			readByte()								const;
	virtual		NxU16			readWord()								const;
	virtual		NxU32			readDword()								const;
	virtual		float			readFloat()								const;
	virtual		double			readDouble()							const;
	virtual		void			readBuffer(void* buffer, NxU32 size)	const;

	virtual		NxStream&		storeByte(NxU8 b)							{ NX_ASSERT(0);	return *this;	}
	virtual		NxStream&		storeWord(NxU16 w)							{ NX_ASSERT(0);	return *this;	}
	virtual		NxStream&		storeDword(NxU32 d)							{ NX_ASSERT(0);	return *this;	}
	virtual		NxStream&		storeFloat(NxReal f)						{ NX_ASSERT(0);	return *this;	}
	virtual		NxStream&		storeDouble(NxF64 f)						{ NX_ASSERT(0);	return *this;	}
	virtual		NxStream&		storeBuffer(const void* buffer, NxU32 size)	{ NX_ASSERT(0);	return *this;	}

	mutable		const NxU8*		buffer;
	};

} // namespace planning
#endif


#endif

