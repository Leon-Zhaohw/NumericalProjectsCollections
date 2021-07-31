
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

#ifndef SMOOTHINGFILTER_is_defined
#define SMOOTHINGFILTER_is_defined

#include <openvdb/openvdb.h>
#include <openvdb/math/Stencils.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/util/NullInterrupter.h>

class SmoothingFilter {
public:
  	typedef typename openvdb::util::NullInterrupter InterruptT;
    typedef boost::shared_ptr<SmoothingFilter>      Ptr;
	  typedef openvdb::DoubleGrid                    GridType;
    typedef typename GridType::TreeType            TreeType;
    typedef typename TreeType::LeafNodeType        LeafType;
    typedef typename LeafType::ValueType           ValueType;
	  typedef typename openvdb::tree::LeafManager<TreeType>   LeafManagerType;
    typedef typename LeafManagerType::RangeType    RangeType;
    typedef typename LeafManagerType::BufferType   BufferType;
    typedef typename LeafType::ValueOnIter         VoxelIterTOn;
    typedef typename LeafType::ValueOffIter        VoxelIterTOff;
    typedef typename LeafType::ValueAllIter        VoxelIterTAll;


    SmoothingFilter(GridType::Ptr grid, int numBuffers, InterruptT* interrupt = NULL):
        mGrid(grid),
				mTask(0),
				mDx(grid->voxelSize()[0]),
				mDelta(0.0),
				mGrainSize(1),
				mInterrupter(interrupt),
				mLeafs(new LeafManagerType(grid->tree(), numBuffers)),
				accessor(mGrid->getAccessor()),
				stencil(*mGrid,mDx)
	  { }
    

    SmoothingFilter(const SmoothingFilter& other):
        mGrid(other.mGrid),
        mTask(other.mTask),
        mDx(other.mDx),
        mDelta(other.mDelta),
        mGrainSize(other.mGrainSize),
        mInterrupter(other.mInterrupter),
        mLeafs(other.mLeafs),
        accessor(other.accessor),
        stencil(other.stencil)
    { }

	  SmoothingFilter(SmoothingFilter& other, tbb::split):
        mGrid(other.mGrid),
        mTask(other.mTask),
        mDx(other.mDx),
        mDelta(0.0),
        mGrainSize(other.mGrainSize),
        mInterrupter(other.mInterrupter),
        mLeafs(other.mLeafs),
        accessor(other.accessor),
        stencil(other.stencil)
    { }

    
    virtual ~SmoothingFilter() {};

    
    /// Used internally by tbb::parallel_for ()
    void operator()(const RangeType& r) const
    {
        if (mTask) mTask(const_cast<SmoothingFilter*>(this), r);
        else OPENVDB_THROW(openvdb::ValueError, "task is undefined - call offset(), etc");
    }
    
    void operator()(const RangeType& r)
    {
        if (mTask) mTask(this, r);
        else OPENVDB_THROW(openvdb::ValueError, "task is undefined - call offset(), etc");
    }

    void join(const SmoothingFilter& other) { mDelta += other.mDelta;}


    void setGrainSize(int flag);

	  void initializeLevelSet(double cMin, double cMax);
    void initializeLevelSetVariableRadius(GridType &minGrid, GridType &maxGrid);

    void laplacianSmoothing(double dt, int ntimesteps);
	  void biharmonicSmoothing(double dt, int ntimesteps, int redistanceFrequency);
	  void biharmonicSmoothing(double dt, int ntimesteps, int redistanceFrequency, double cgThreshold, int cgMaxIters, bool verbose);
    void redistance(int buf);
    void extractSurface(char *outfname);

private:
    typedef typename boost::function<void(SmoothingFilter*, const RangeType&)> FuncType;
    GridType::Ptr mGrid;
    FuncType mTask;
    ValueType mDx;
    ValueType mDelta;
    int mGrainSize;
    InterruptT* mInterrupter;
    LeafManagerType* mLeafs;
	  openvdb::DoubleGrid::Accessor accessor;
	  openvdb::math::GradStencil<openvdb::DoubleGrid> stencil;

    SmoothingFilter& operator=(const SmoothingFilter& other) { return *this; }

	  void cook(int swapBuffer);
    bool wasInterrupted();

	  // Initialization Methods
    void initialize(double cMin, double cMax);
	  void initializeVariableRadius(GridType &minGrid, GridType &maxGrid);
    void doInitialize(const RangeType &range, double cMin, double cMax);
    void doInitializeVariableRadius(const RangeType &range, GridType &minGrid, GridType &maxGrid);

  	// Grid Activation Methods
	  void setOnOffRead(double updateBand);
    void setOnOffAux1(double updateBand);
    void setOnOffAux1WithConstraints(double updateBand);
    void doSetOnOffRead(const RangeType&, double updateBand);
    void doSetOnOffAux1(const RangeType&, double updateBand);
    void doSetOnOffAux1WithConstraints(const RangeType&, double updateBand);

	  // Explicit Integration Methods
	  void laplacianStep(double dt);
    void biharmonicStep(double dt);
    void computeLaplacianAndGradient();
    void doLaplacianStep(const RangeType &range, double dt);
    void doBiharmonicStep(const RangeType&, double dt);
    void doComputeLaplacianAndGradient(const RangeType&);

	  // Implicit Integration Methods
  	bool CG(double dt, double &tol, int &cgIters, double threshold, int cgMaxIters, bool verbose);
	  void initializeImplicitBiharmonic();
    void zeroBuffers();
  	void computeLaplacian();
    void applyBiharmonic(double dt, int sourceBuffer);
    void computeInitialResidual();
    void CGStep(double alpha);
    void CGDirection(double sigma);
    void applyConstraints(bool &constraintsViolated);
    void innerProduct(int xid, int yid);
    void copySolutionToReadInBand();

	  void doInitializeImplicitBiharmonic(const RangeType &range);
    void doZeroBuffers(const RangeType &range);
    void doComputeLaplacian(const RangeType& range);
    void doApplyBiharmonic(const RangeType& range, double dt, int sourceBuffer);
    void doComputeInitialResidual(const RangeType &range);
    void doCGStep(const RangeType &range, double alpha);
    void doCGDirection(const RangeType &range, double sigma);
  	void doApplyConstraints(const RangeType &range, bool &constraintsViolated);
    void doInnerProduct(const RangeType &range, int xid, int yid);
    void doCopySolutionToReadInBand(const RangeType &range);

  	// Redistancing Methods
    void redistanceInterface();
    void redistanceSweep();
    void copyAllValuesFromReadToAux();
    void copyLockedValuesFromAuxToRead();
    void copyUnlockedValuesFromAuxToRead();
	  void fastSweepingInterfaceHelper(size_t &n, VoxelIterTOn &vIterOn, openvdb::Index &cI, double &x, double &y);
    void doRedistanceInterface(const RangeType&);
    void doRedistanceSweep(const RangeType&);
	  void doCopyAllValuesFromReadToAux(const RangeType &range);
    void doCopyLockedValuesFromAuxToRead(const RangeType &range);
    void doCopyUnlockedValuesFromAuxToRead(const RangeType &range);
};


#endif
