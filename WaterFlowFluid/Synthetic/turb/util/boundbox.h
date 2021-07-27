/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Bounding box class
 *
 *****************************************************************************/

#ifndef BOUNDBOX_H
#define BOUNDBOX_H

#include "vectorbase.h"


namespace DDF { 

//! bbox for a vector
template<class VecClass>
class BboxVec {
	public:
		BboxVec() : mStart(0), mEnd(0) { };
		BboxVec(VecClass s, VecClass e) :
			mStart(s), mEnd(e) { };
		~BboxVec() {};

		bool contains(VecClass p) {
			for(int i=0; i<3; i++) if(p[i]<mStart[i]) return false;
			for(int i=0; i<3; i++) if(p[i]>mEnd[i]) return false;
// print stats of pressure min/max
			return true;
		};

		inline VecClass& getStart() { return mStart; }
		inline VecClass& getEnd() { return mEnd; }

		std::string toString() {
			std::ostringstream out;
			out <<"bbox s:"<<mStart
			        <<" e:"<<mEnd <<" ";
			return out.str();
		};

	protected:
		// region
		VecClass mStart;
		VecClass mEnd;
		
}; // BboxVec

typedef BboxVec<Vec3>  BboxVecr; 
typedef BboxVec<nVec3i>  BboxVeci; 

} // namespace DDF 

#endif


