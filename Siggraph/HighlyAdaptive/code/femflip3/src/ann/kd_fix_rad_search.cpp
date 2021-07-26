//----------------------------------------------------------------------
// File:			kd_fix_rad_search.cpp
// Programmer:		Sunil Arya and David Mount
// Description:		Standard kd-tree fixed-radius kNN search
// Last modified:	05/03/05 (Version 1.1)
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------
// History:
//	Revision 1.1  05/03/05
//		Initial release
//----------------------------------------------------------------------

#include "kd_fix_rad_search.h"			// kd fixed-radius search decls

//----------------------------------------------------------------------
//	Approximate fixed-radius k nearest neighbor search
//		The squared radius is provided, and this procedure finds the
//		k nearest neighbors within the radius, and returns the total
//		number of points lying within the radius.
//
//		The method used for searching the kd-tree is a variation of the
//		nearest neighbor search used in kd_search.cpp, except that the
//		radius of the search ball is known.  We refer the reader to that
//		file for the explanation of the recursive search procedure.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	annkFRSearch - fixed radius search for k nearest neighbors
//----------------------------------------------------------------------

std::vector<ANNidx>
ANNkd_tree::annkFRSearch(
	ANNpoint			q,				// the query point
	ANNdist				sqRad,			// squared radius search bound
	double				eps)			// the error bound
{
	ANNGlobal global;
	global.ANNkdFRDim = dim;					// copy arguments to static equivs
	global.ANNkdFRQ = q;
	global.ANNkdFRSqRad = sqRad;
	global.ANNkdFRPts = pts;
	global.ANNkdFRPtsVisited = 0;				// initialize count of points visited
	global.ANNkdFRPtsInRange = 0;				// ...and points in the range

	global.ANNkdFRMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

	//global.ANNkdFRPointMK = new ANNmin_k(k);	// create set for closest k points
										// search starting at the root
	
	//----------------------------------------------------------------------
	//	Limit on number of points visited
	//		We have an option for terminating the search early if the
	//		number of points visited exceeds some threshold.  If the
	//		threshold is 0 (its default)  this means there is no limit
	//		and the algorithm applies its normal termination condition.
	//		This is for applications where there are real time constraints
	//		on the running time of the algorithm.
	//----------------------------------------------------------------------
	
	ANNVisit visit;
	visit.ANNmaxPtsVisited = 0;	// maximum number of pts visited
	root->ann_FR_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim),visit,global);
/*
	for (int i = 0; i < k; i++) {		// extract the k-th closest points
		if (dd != NULL)
			dd[i] = global.ANNkdFRPointMK->ith_smallest_key(i);
		if (nn_idx != NULL)
			nn_idx[i] = global.ANNkdFRPointMK->ith_smallest_info(i);
	}
*/
	//delete global.ANNkdFRPointMK;				// deallocate closest point set
	//return global.ANNkdFRPtsInRange;			// return final point count
	return global.indices;
}

//----------------------------------------------------------------------
//	kd_split::ann_FR_search - search a splitting node
//		Note: This routine is similar in structure to the standard kNN
//		search.  It visits the subtree that is closer to the query point
//		first.  For fixed-radius search, there is no benefit in visiting
//		one subtree before the other, but we maintain the same basic
//		code structure for the sake of uniformity.
//----------------------------------------------------------------------

void ANNkd_split::ann_FR_search(ANNdist box_dist, ANNVisit &visit, ANNGlobal& global )
{
										// check dist calc term condition
	if (visit.ANNmaxPtsVisited != 0 && visit.ANNmaxPtsVisited > visit.ANNmaxPtsVisited) return;

										// distance to cutting plane
	ANNcoord cut_diff = global.ANNkdFRQ[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		child[ANN_LO]->ann_FR_search(box_dist,visit,global);// visit closer child first

		ANNcoord box_diff = cd_bnds[ANN_LO] - global.ANNkdFRQ[cut_dim];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
										// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

										// visit further child if in range
		if (box_dist * global.ANNkdFRMaxErr <= global.ANNkdFRSqRad)
			child[ANN_HI]->ann_FR_search(box_dist,visit,global);

	}
	else {								// right of cutting plane
		child[ANN_HI]->ann_FR_search(box_dist,visit,global);// visit closer child first

		ANNcoord box_diff = global.ANNkdFRQ[cut_dim] - cd_bnds[ANN_HI];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
										// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

										// visit further child if close enough
		if (box_dist * global.ANNkdFRMaxErr <= global.ANNkdFRSqRad)
			child[ANN_LO]->ann_FR_search(box_dist,visit,global);

	}
	ANN_FLOP(13)						// increment floating ops
	ANN_SPL(1)							// one more splitting node visited
}

//----------------------------------------------------------------------
//	kd_leaf::ann_FR_search - search points in a leaf node
//		Note: The unreadability of this code is the result of
//		some fine tuning to replace indexing by pointer operations.
//----------------------------------------------------------------------

void ANNkd_leaf::ann_FR_search(ANNdist box_dist, ANNVisit& visit, ANNGlobal& global )
{
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNcoord t;
	register int d;

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		pp = global.ANNkdFRPts[bkt[i]];		// first coord of next data point
		qq = global.ANNkdFRQ;					// first coord of query point
		dist = 0;

		for(d = 0; d < global.ANNkdFRDim; d++) {
			ANN_COORD(1)				// one more coordinate hit
			ANN_FLOP(5)					// increment floating ops

			t = *(qq++) - *(pp++);		// compute length and adv coordinate
										// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > global.ANNkdFRSqRad) {
				break;
			}
		}

		if (d >= global.ANNkdFRDim &&					// among the k best?
		   (ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
												// add it to the list
			//global.ANNkdFRPointMK->insert(dist, bkt[i]);
			//global.ANNkdFRPtsInRange++;				// increment point count
			global.indices.push_back(bkt[i]);
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
	ANN_PTS(n_pts)						// increment points visited
	global.ANNkdFRPtsVisited += n_pts;			// increment number of points visited
}
