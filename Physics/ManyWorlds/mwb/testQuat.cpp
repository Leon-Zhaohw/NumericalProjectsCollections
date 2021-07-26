#include "stdafx.h"

#include "quat.h"

#include "twigg/random.h"

Quaternion<float> randomQuat(RandomStream& stream)
{
	vl::Vec3f vec;
	for( size_t i = 0; i < 3; ++i )
		vec[i] = stream.normal();

	Quaternion<float> result( vl::norm(vec), stream.uniform(0.0, M_PI) );
	result.normalize();
	return result;
}

void testQuat()
{
	RandomStream stream;
	for( size_t i = 0; i < 30; ++i )
	{
		Quaternion<float> left = randomQuat(stream);
		Quaternion<float> right = randomQuat(stream);

		Quaternion<float> test = left*right;

		Quaternion<float> nxTest = 
			Quaternion<float>(left * right);

		Quaternion<float> error = inverse( test ) * nxTest;
		float scalarError = arg( error );
		assert( scalarError < 1e-4 );
	}

    // Test slerp code
    {
        Quaternion<float> q1 = randomQuat(stream);
        Quaternion<float> q2 = randomQuat(stream);
        
        Quaternion<float> result = slerp(q1, q2, 0.0f);
        Quaternion<float> diff = inverse(result) * q1;
        assert( arg(diff) < 1e-4 );

        result = slerp(q1, q2, 1.0f);
        diff = inverse(result) * q2;
        assert( arg(diff) < 1e-4 );

        q1 = Quaternion<float>( vl::Vec3d(0.0, 0.0, 0.0) );
        q2 = Quaternion<float>( vl::Vec3d(0.0, M_PI, 0.0) );
        result = slerp( q1, q2, 0.5f );
        Quaternion<float> q3( vl::Vec3d(0.0, 0.5*M_PI, 0.0) );
        diff = inverse(result) * q3;
        assert( arg(diff) < 1e-4 );
    }
}
