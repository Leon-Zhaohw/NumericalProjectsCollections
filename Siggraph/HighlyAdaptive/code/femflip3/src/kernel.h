/*
 *	kernel.h
 *	
 *	Created by Ryoichi Ando on 1/8/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include <math.h>

#ifndef _KERNEL_H_
#define _KERNEL_H_

namespace kernel {
    static FLOAT64 smooth_kernel( FLOAT64 r2, FLOAT64 h ) {
		r2 = r2/(h*h);
		return powf(fmax(0.0,1.0-r2),3.0);
	}
    static FLOAT64 sharp_kernel( FLOAT64 r2, FLOAT64 h ) {
		return r2 ? fmax( h*h/r2 - 1.0, 0.0 ) : 1e20;
	}
}

#endif
