// Copyright (c) 2011, Regents of the University of Utah
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "kdTree.H"
#include <float.h>
#include <algorithm>
bool IntDoublePairCompare (const IntDoublePair &p1, const IntDoublePair &p2) {
	return (p1.d() < p2.d());
};

KDTree::KDTree(const std::vector<SlVector3> &pts) {
	std::vector<IntDoublePair> xSortedList;
	std::vector<IntDoublePair> ySortedList;
	std::vector<IntDoublePair> zSortedList;
	std::vector<IntDoublePair> scratch;

	xSortedList.resize(pts.size());
	ySortedList.resize(pts.size());
	zSortedList.resize(pts.size());
	scratch.resize(pts.size());

	std::vector<IntDoublePair>::iterator xi = xSortedList.begin();
	std::vector<IntDoublePair>::iterator yi = ySortedList.begin();
	std::vector<IntDoublePair>::iterator zi = zSortedList.begin();
	std::vector<SlVector3>::const_iterator pi = pts.begin();
	unsigned int i;

	for (i=0; i<pts.size(); i++, xi++, yi++, zi++, pi++) {
		xi->i() = i;
		yi->i() = i;
		zi->i() = i;
		xi->d() = (*pi)[0];
		yi->d() = (*pi)[1];
		zi->d() = (*pi)[2];
	}

	sort(xSortedList.begin(), xSortedList.end(), IntDoublePairCompare);
	sort(ySortedList.begin(), ySortedList.end(), IntDoublePairCompare);
	sort(zSortedList.begin(), zSortedList.end(), IntDoublePairCompare);

	xsplit(pts, xSortedList, ySortedList, zSortedList, scratch, 0, pts.size());

	tree = new int[pts.size()];
	int *ip = tree;

	for (xi = xSortedList.begin(), i=0; i<pts.size(); i++, ip++, xi++) {
		(*ip) = xi->i();
	}
	npts = pts.size();

	distCache = new double[npts];
}

// The split routines pick the median point in the given sorted array
// and then partition all three arrays so that the set of points to each side
// of the median are identical.  This is slightly complicated to deal
// with duplicated coordinates.

void KDTree::xsplit(const std::vector<SlVector3> &pts,
							 std::vector<IntDoublePair> &xSortedList, 
							 std::vector<IntDoublePair> &ySortedList, 
							 std::vector<IntDoublePair> &zSortedList, 
							 std::vector<IntDoublePair> &scratch, 
							 int start, int end) {
	if (start >= end-1) return;
	int median = start+(end-start)/2;
	int medid = xSortedList[median].i();
	double medval = pts[medid][0];
	int medlistindex = -1;
	int i,j;
	
	// Y
	for (i=start, j=0; i<end; i++) {
		const double &d = pts[ySortedList[i].i()][0];
		if (d == medval) {
			if (ySortedList[i].i() == medid) {
				medlistindex = i;
				continue;
			} else {
				for (int k=median-1; k >= start; k--) {
					if (ySortedList[i].i() == xSortedList[k].i()) {
						scratch[j++] = ySortedList[i];
						break;
					}
				}
			}
		}	else if (d < medval) {
			scratch[j++] = ySortedList[i];
		}
	}

	scratch[j++] = ySortedList[medlistindex];

	for (i=start; i<end; i++) {
		const double &d = pts[ySortedList[i].i()][0];
		if (d == medval) {
			if (ySortedList[i].i() == medid) continue;
			else {
				for (int k=median+1; k < end; k++) {
					if (ySortedList[i].i() == xSortedList[k].i()) {
						scratch[j++] = ySortedList[i];
						break;
					}
				}
			}
		} else if (d > medval) {
			scratch[j++] = ySortedList[i];
		}
	}

	
	for (i=start, j=0; i<end; i++, j++) {
		ySortedList[i] = scratch[j];
	}

	// Z
	for (i=start, j=0; i<end; i++) {
		const double &d = pts[zSortedList[i].i()][0];
		if (d == medval) {
			if (zSortedList[i].i() == medid) {
				medlistindex = i;
				continue;
			} else {
				for (int k=median-1; k >= start; k--) {
					if (zSortedList[i].i() == xSortedList[k].i()) {
						scratch[j++] = zSortedList[i];
						break;
					}
				}
			}
		}	else if (d < medval) {
			scratch[j++] = zSortedList[i];
		}
	}

	scratch[j++] = zSortedList[medlistindex];

	for (i=start; i<end; i++) {
		const double &d = pts[zSortedList[i].i()][0];
		if (d == medval) {
			if (zSortedList[i].i() == medid) continue;
			else {
				for (int k=median+1; k < end; k++) {
					if (zSortedList[i].i() == xSortedList[k].i()) {
						scratch[j++] = zSortedList[i];
						break;
					}
				}
			}
		} else if (d > medval) {
			scratch[j++] = zSortedList[i];
		}
	}
	
	for (i=start, j=0; i<end; i++, j++) {
		zSortedList[i] = scratch[j];
	}


	ysplit(pts, xSortedList, ySortedList, zSortedList, scratch, start, median);
	ysplit(pts, xSortedList, ySortedList, zSortedList, scratch, median+1, end);
}

void KDTree::ysplit(const std::vector<SlVector3> &pts,
							 std::vector<IntDoublePair> &xSortedList, 
							 std::vector<IntDoublePair> &ySortedList, 
							 std::vector<IntDoublePair> &zSortedList, 
							 std::vector<IntDoublePair> &scratch, 
							 int start, int end) {
	if (start >= end-1) return;
	int median = start+(end-start)/2;
	int medid = ySortedList[median].i();
	double medval = pts[medid][1];
	int medlistindex = -1;
	int i,j;

	// X
	for (i=start, j=0; i<end; i++) {
		const double &d = pts[xSortedList[i].i()][1];
		if (d == medval) {
			if (xSortedList[i].i() == medid) {
				medlistindex = i;
				continue;
			} else {
				for (int k=median-1; k >= start; k--) {
					if (xSortedList[i].i() == ySortedList[k].i()) {
						scratch[j++] = xSortedList[i];
						break;
					}
				}
			}
		}	else if (d < medval) {
			scratch[j++] = xSortedList[i];
		}
	}

	scratch[j++] = xSortedList[medlistindex];

	for (i=start; i<end; i++) {
		const double &d = pts[xSortedList[i].i()][1];
		if (d == medval) {
			if (xSortedList[i].i() == medid) continue;
			else {
				for (int k=median+1; k < end; k++) {
					if (xSortedList[i].i() == ySortedList[k].i()) {
						scratch[j++] = xSortedList[i];
						break;
					}
				}
			}
		} else if (d > medval) {
			scratch[j++] = xSortedList[i];
		}
	}

	
	for (i=start, j=0; i<end; i++, j++) {
		xSortedList[i] = scratch[j];
	}

	// Z
	for (i=start, j=0; i<end; i++) {
		const double &d = pts[zSortedList[i].i()][1];
		if (d == medval) {
			if (zSortedList[i].i() == medid) {
				medlistindex = i;
				continue;
			} else {
				for (int k=median-1; k >= start; k--) {
					if (zSortedList[i].i() == ySortedList[k].i()) {
						scratch[j++] = zSortedList[i];
						break;
					}
				}
			}
		}	else if (d < medval) {
			scratch[j++] = zSortedList[i];
		}
	}

	scratch[j++] = zSortedList[medlistindex];

	for (i=start; i<end; i++) {
		const double &d = pts[zSortedList[i].i()][1];
		if (d == medval) {
			if (zSortedList[i].i() == medid) continue;
			else {
				for (int k=median+1; k < end; k++) {
					if (zSortedList[i].i() == ySortedList[k].i()) {
						scratch[j++] = zSortedList[i];
						break;
					}
				}
			}
		} else if (d > medval) {
			scratch[j++] = zSortedList[i];
		}
	}
	
	for (i=start, j=0; i<end; i++, j++) {
		zSortedList[i] = scratch[j];
	}

	zsplit(pts, xSortedList, ySortedList, zSortedList, scratch, start, median);
	zsplit(pts, xSortedList, ySortedList, zSortedList, scratch, median+1, end);
}

void KDTree::zsplit(const std::vector<SlVector3> &pts,
							 std::vector<IntDoublePair> &xSortedList, 
							 std::vector<IntDoublePair> &ySortedList, 
							 std::vector<IntDoublePair> &zSortedList, 
							 std::vector<IntDoublePair> &scratch, 
							 int start, int end) {
	if (start >= end-1) return;
	int median = start+(end-start)/2;
	int medid = zSortedList[median].i();
	double medval = pts[medid][2];
	int medlistindex = -1;
	int i,j;
	
	// X
	for (i=start, j=0; i<end; i++) {
		const double &d = pts[xSortedList[i].i()][2];
		if (d == medval) {
			if (xSortedList[i].i() == medid) {
				medlistindex = i;
				continue;
			} else {
				for (int k=median-1; k >= start; k--) {
					if (xSortedList[i].i() == zSortedList[k].i()) {
						scratch[j++] = xSortedList[i];
						break;
					}
				}
			}
		}	else if (d < medval) {
			scratch[j++] = xSortedList[i];
		}
	}

	scratch[j++] = xSortedList[medlistindex];

	for (i=start; i<end; i++) {
		const double &d = pts[xSortedList[i].i()][2];
		if (d == medval) {
			if (xSortedList[i].i() == medid) continue;
			else {
				for (int k=median+1; k < end; k++) {
					if (xSortedList[i].i() == zSortedList[k].i()) {
						scratch[j++] = xSortedList[i];
						break;
					}
				}
			}
		} else if (d > medval) {
			scratch[j++] = xSortedList[i];
		}
	}

	
	for (i=start, j=0; i<end; i++, j++) {
		xSortedList[i] = scratch[j];
	}

	// Y
	for (i=start, j=0; i<end; i++) {
		const double &d = pts[ySortedList[i].i()][2];
		if (d == medval) {
			if (ySortedList[i].i() == medid) {
				medlistindex = i;
				continue;
			} else {
				for (int k=median-1; k >= start; k--) {
					if (ySortedList[i].i() == zSortedList[k].i()) {
						scratch[j++] = ySortedList[i];
						break;
					}
				}
			}
		}	else if (d < medval) {
			scratch[j++] = ySortedList[i];
		}
	}

	scratch[j++] = ySortedList[medlistindex];

	for (i=start; i<end; i++) {
		const double &d = pts[ySortedList[i].i()][2];
		if (d == medval) {
			if (ySortedList[i].i() == medid) continue;
			else {
				for (int k=median+1; k < end; k++) {
					if (ySortedList[i].i() == zSortedList[k].i()) {
						scratch[j++] = ySortedList[i];
						break;
					}
				}
			}
		} else if (d > medval) {
			scratch[j++] = ySortedList[i];
		}
	}

	
	for (i=start, j=0; i<end; i++, j++) {
		ySortedList[i] = scratch[j];
	}


	xsplit(pts, xSortedList, ySortedList, zSortedList, scratch, start, median);
	xsplit(pts, xSortedList, ySortedList, zSortedList, scratch, median+1, end);
}

void checkHeap(std::vector<int> &heap, const double *distCache) {
	for (unsigned int i=0; i<heap.size(); i++) {
		if (2*i+1 < heap.size())
			if (distCache[heap[2*i+1]] > distCache[heap[i]]) std::cout<<"broken heap "<<heap[i]<<" "<<heap[2*i]<<" "<<std::endl;
		if (2*i+2 < heap.size())
			if (distCache[heap[2*i+2]] > distCache[heap[i]]) std::cout<<"broken heap "<<heap[i]<<" "<<heap[2*i+1]<<" "<<std::endl;
	}
}

void heapRemove (std::vector<int> &heap, const double *distCache) {
	heap.front() = heap.back();
	unsigned int i=0;
	unsigned int swapIndex = 0;
	heap.resize(heap.size()-1);
	while (true) {
		if (2*i+1 >= heap.size()) return;
		if (distCache[heap[2*i+1]] > distCache[heap[i]]) swapIndex = 2*i+1;
		if (2*i+2 >= heap.size()) return;
		if (distCache[heap[2*i+2]] > distCache[heap[swapIndex]]) swapIndex = 2*i+2;
		if (swapIndex == i) return;
		unsigned int tmp = heap[i];
		heap[i] = heap[swapIndex];
		heap[swapIndex] = tmp;
		i = swapIndex;
	}
	//checkHeap(heap, distCache);
}

void heapAdd (std::vector<int> &heap, const double *distCache, int x) {
	int i = heap.size();
	heap.push_back(x);
	double di = distCache[x];
	
	while (true) {
		int j = (i-1)/2;
		if (j < 0) return;
		double dj = distCache[heap[j]];
		if (dj >= di) return;
		int tmp = heap[i];
		heap[i] = heap[j];
		heap[j] = tmp;
		if (j == 0) return;
		i = j;
	}
	//checkHeap(heap, distCache);
}

void KDTree::neighbors(const std::vector<SlVector3> &pts, const SlVector3 &x, int num, double r, std::vector<int> &neighbors) {
	neighbors.clear();
	if (num > 0 && r > 0) {
		double r2 = r*r;
		neighbors.reserve(num);
		neighborsRecurse(pts, x, num, r, r2, neighbors, 0, npts, 0);
	} else if (num > 0) {
		neighbors.reserve(num);
		neighborsRecurse(pts, x, num, neighbors, 0, npts, 0);
	} else if (r > 0) {
		double r2 = r*r;
		neighborsRecurse(pts, x, r, r2, neighbors, 0, npts, 0);
	} else {
		std::cout<<"either num or r must be > zero"<<std::endl;
	}
}


int KDTree::neighbor(const std::vector<SlVector3> &pts, const SlVector3 &x, double r) {
	int neighbor = -1;
	double ndist = DBL_MAX;
	if (r > 0) {
		double r2 = r*r;
		neighborRecurse(pts, x, r, r2, neighbor, ndist, 0, npts, 0);
	} else {
		neighborRecurse(pts, x, neighbor, ndist, 0, npts, 0);
	}
	return neighbor;
}

void KDTree::neighborsRecurse(const std::vector<SlVector3> &pts, const SlVector3 &x, unsigned int num, double r, double r2, std::vector<int> &neighbors, int start, int end, int plane) {
	if (start >= end) return;
	int median = start+(end-start)/2;
	int medid = tree[median];
	double dist = x[plane]-pts[medid][plane];
	double dist2 = dist*dist;

	if (dist2 < r2) {
		bool searchBoth = false;
		if (neighbors.size() < num) {
			searchBoth = true;
			distCache[medid] = sqrMag(pts[medid] - x);
			if (distCache[medid] < r2) {
				heapAdd(neighbors, distCache, medid);
			}
		} else if (dist2 < distCache[neighbors.front()]) {
			searchBoth = true;
			distCache[medid] = sqrMag(pts[medid] - x);
			if (distCache[medid] < distCache[neighbors.front()]) {
				heapRemove(neighbors, distCache);
				heapAdd(neighbors, distCache, medid);
			}
		}
		
		if (dist < 0) {
			neighborsRecurse(pts, x, num, r, r2, neighbors, start, median, (plane+1)%3);
			if (searchBoth) neighborsRecurse(pts, x, num, r, r2, neighbors, median+1, end, (plane+1)%3);
		} else {
			neighborsRecurse(pts, x, num, r, r2, neighbors, median+1, end, (plane+1)%3);
			if (searchBoth) neighborsRecurse(pts, x, num, r, r2, neighbors, start, median, (plane+1)%3);
		}
	} else {
		if (dist < 0) {
			neighborsRecurse(pts, x, num, r, r2, neighbors, start, median, (plane+1)%3);
		} else {
			neighborsRecurse(pts, x, num, r, r2, neighbors, median+1, end, (plane+1)%3);
		}
	}
}

void KDTree::neighborsRecurse(const std::vector<SlVector3> &pts, const SlVector3 &x, unsigned int num, std::vector<int> &neighbors, int start, int end, int plane) {
	if (start >= end) return;
	int median = start+(end-start)/2;
	int medid = tree[median];
	double dist = x[plane]-pts[medid][plane];
	double dist2 = dist*dist;

	bool searchBoth = false;
	if (neighbors.size() < num) {
		searchBoth = true;
		distCache[medid] = sqrMag(pts[medid] - x);
		heapAdd(neighbors, distCache, medid);
	} else if (dist2 < distCache[neighbors.front()]) {
		searchBoth = true;
		distCache[medid] = sqrMag(pts[medid] - x);
		if (distCache[medid] < distCache[neighbors.front()]) {
			heapRemove(neighbors, distCache);
			heapAdd(neighbors, distCache, medid);
		}
	}
		
	if (dist < 0) {
		neighborsRecurse(pts, x, num, neighbors, start, median, (plane+1)%3);
		if (searchBoth) neighborsRecurse(pts, x, num, neighbors, median+1, end, (plane+1)%3);
	} else {
		neighborsRecurse(pts, x, num, neighbors, median+1, end, (plane+1)%3);
		if (searchBoth) neighborsRecurse(pts, x, num, neighbors, start, median, (plane+1)%3);
	}
}

void KDTree::neighborsRecurse(const std::vector<SlVector3> &pts, const SlVector3 &x, double r, double r2, std::vector<int> &neighbors, int start, int end, int plane) {
	if (start >= end) return;
	int median = start+(end-start)/2;
	int medid = tree[median];
	double dist = x[plane]-pts[medid][plane];
	double dist2 = dist*dist;

	if (dist2 < r2) {
		if (sqrMag(pts[medid] - x) < r2) neighbors.push_back(medid);
		
		if (dist < 0) {
			neighborsRecurse(pts, x, r, r2, neighbors, start, median, (plane+1)%3);
			neighborsRecurse(pts, x, r, r2, neighbors, median+1, end, (plane+1)%3);
		} else {
			neighborsRecurse(pts, x, r, r2, neighbors, median+1, end, (plane+1)%3);
			neighborsRecurse(pts, x, r, r2, neighbors, start, median, (plane+1)%3);
		}
	} else {
		if (dist < 0) {
			neighborsRecurse(pts, x, r, r2, neighbors, start, median, (plane+1)%3);
		} else {
			neighborsRecurse(pts, x, r, r2, neighbors, median+1, end, (plane+1)%3);
		}
	}
}

void KDTree::neighborRecurse(const std::vector<SlVector3> &pts, const SlVector3 &x, double r, double r2, 
														 int &neighbor, double &ndist, int start, int end, int plane) {
	if (start >= end) return;
	int median = start+(end-start)/2;
	int medid = tree[median];
	
	double dist = x[plane]-pts[medid][plane];
	double dist2 = dist*dist;

	if (dist2 < r2) {
		bool searchBoth = false;
		if (neighbor == -1) {
			searchBoth = true;
			double d = sqrMag(pts[medid] - x);
			if (d < r2) {
				neighbor = medid;
				ndist = d;
			}
		} else if (dist2 < distCache[neighbor]) {
			searchBoth = true;
			double d = sqrMag(pts[medid] - x);
			if (d < ndist) {
				neighbor = medid;
				ndist = d;
			}
		}
		
		if (dist < 0) {
			neighborRecurse(pts, x, r, r2, neighbor, ndist, start, median, (plane+1)%3);
			if (searchBoth) neighborRecurse(pts, x, r, r2, neighbor, ndist, median+1, end, (plane+1)%3);
		} else {
			neighborRecurse(pts, x, r, r2, neighbor, ndist, median+1, end, (plane+1)%3);
			if (searchBoth) neighborRecurse(pts, x, r, r2, neighbor, ndist, start, median, (plane+1)%3);
		}
	} else {
		if (dist < 0) {
			neighborRecurse(pts, x, r, r2, neighbor, ndist, start, median, (plane+1)%3);
		} else {
			neighborRecurse(pts, x, r, r2, neighbor, ndist, median+1, end, (plane+1)%3);
		}
	}
}

void KDTree::neighborRecurse(const std::vector<SlVector3> &pts, const SlVector3 &x,  
														 int &neighbor, double &ndist, int start, int end, int plane) {
	
	if (start >= end) return;
	int median = start+(end-start)/2;
	int medid = tree[median];
	
	double dist = x[plane]-pts[medid][plane];
	double dist2 = dist*dist;

	bool searchBoth = false;
	if (neighbor == -1) {
		searchBoth = true;
		double d = sqrMag(pts[medid] - x);
		neighbor = medid;
		ndist = d;
	} else if (dist2 < ndist) {
		searchBoth = true;
		double d = sqrMag(pts[medid] - x);
		if (d < ndist) {
			neighbor = medid;
			ndist = d;
		}
	}
	
	if (dist < 0) {
		neighborRecurse(pts, x, neighbor, ndist, start, median, (plane+1)%3);
		if (searchBoth) neighborRecurse(pts, x, neighbor, ndist, median+1, end, (plane+1)%3);
	} else {
		neighborRecurse(pts, x, neighbor, ndist, median+1, end, (plane+1)%3);
		if (searchBoth) neighborRecurse(pts, x, neighbor, ndist, start, median, (plane+1)%3);
	}
}
