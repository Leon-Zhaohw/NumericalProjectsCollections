
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



#include "smoothingFilter.H"
#include <openvdb/tools/VolumeToMesh.h>

using namespace openvdb;

enum Buffers { Read = -1, MinPhi = 1, MaxPhi = 2, Laplacian = 3, GradMag = 4, Aux1 = 4, Aux2 = 5,
							 Q = 5, X = 6, R = 7, D = 8, Activeset = 9};

inline int sign(float x) {
    if (x > 0.0) return 1;
    else return -1;
}

inline double absmin (double x, double y) {
	if (std::abs(x) < std::abs(y)) return x;
    return y;
}

void
SmoothingFilter::setGrainSize(int flag)
{
    this->mGrainSize = flag;
    return;
}


void
SmoothingFilter::cook(int swapBuffer)
{
    if (mInterrupter) mInterrupter->start("Cook");
#if 0
		(*this)(mLeafs->getRange());
#else
    switch(mGrainSize) {
    case 0:
        (*this)(mLeafs->getRange());
        break;
    case 1:
        tbb::parallel_for(mLeafs->getRange(), *this);
        break;
    case 2:
        tbb::parallel_reduce(mLeafs->getRange(), *this);
        break;
    }
#endif
    if(swapBuffer > 0) mLeafs->swapLeafBuffer(swapBuffer, mGrainSize==0);
    if (mInterrupter) mInterrupter->end();
    return;
}


bool
SmoothingFilter::wasInterrupted()
{
    if (util::wasInterrupted(mInterrupter)) {
        tbb::task::self().cancel_group_execution();
        return true;
    }
    return false;
}


// Initialization methods

// Initialize min and max constraint buffer and the initial surface in the read buffer
void
SmoothingFilter::initializeLevelSet(double rMin, double rMax)
{
    initialize(rMin, rMax);
    redistance(Buffers::MinPhi);
    redistance(Buffers::MaxPhi);
    redistance(Buffers::Read);
}


// Same as above with variable radii
void
SmoothingFilter::initializeLevelSetVariableRadius(GridType &minGrid, GridType &maxGrid)
{
    initializeVariableRadius(minGrid, maxGrid);
    redistance(Buffers::MinPhi);
    redistance(Buffers::MaxPhi);
    redistance(Buffers::Read);
}

// tbb wrapper to call the doInitialize methods
void
SmoothingFilter::initialize(double rMin, double rMax)
{
    if (mInterrupter) mInterrupter->start("Initialize");
    mTask = boost::bind(&SmoothingFilter::doInitialize, _1, _2, rMin, rMax);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::initializeVariableRadius(GridType &minGrid, GridType &maxGrid)
{
    if (mInterrupter) mInterrupter->start("InitializeVariableRadius");
    mTask = boost::bind(&SmoothingFilter::doInitializeVariableRadius, _1, _2, boost::ref(minGrid), boost::ref(maxGrid));
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}

// set the values in the buffers.
void
SmoothingFilter::doInitialize(const RangeType& range, double rMin, double rMax)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& minBuffer = mLeafs->getBuffer(n, Buffers::MinPhi);
        BufferType& maxBuffer = mLeafs->getBuffer(n, Buffers::MaxPhi);
        for (VoxelIterTAll vIter = mLeafs->leaf(n).beginValueAll(); vIter; ++vIter) {
            Index c = vIter.pos();
            const double phi = vIter.getValue();
            minBuffer.setValue(c, phi-rMin);
            maxBuffer.setValue(c, phi-rMax);
        }
    }
}


void
SmoothingFilter::doInitializeVariableRadius(const RangeType& range, GridType &minGrid, GridType &maxGrid)
{
    this->wasInterrupted();
    DoubleGrid::Accessor minAccessor = minGrid.getAccessor();
    DoubleGrid::Accessor maxAccessor = maxGrid.getAccessor();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& minBuffer = mLeafs->getBuffer(n, Buffers::MinPhi);
        BufferType& maxBuffer = mLeafs->getBuffer(n, Buffers::MaxPhi);
        for (VoxelIterTAll vIter = mLeafs->leaf(n).beginValueAll(); vIter; ++vIter) {
            Index c = vIter.pos();
            minBuffer.setValue(c, minAccessor.getValue(vIter.getCoord()));
            maxBuffer.setValue(c, maxAccessor.getValue(vIter.getCoord()));
        }
    }
}


// Grid Activation Methods
// Wrappers to call the doSetOnOff* routines
void
SmoothingFilter::setOnOffRead(double updateBand)
{
    if (mInterrupter) mInterrupter->start("SetonoffRead");
    mTask = boost::bind(&SmoothingFilter::doSetOnOffRead, _1, _2, updateBand);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::setOnOffAux1(double updateBand)
{
    if (mInterrupter) mInterrupter->start("SetOnOffAux1");
    mTask = boost::bind(&SmoothingFilter::doSetOnOffAux1, _1, _2, updateBand);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::setOnOffAux1WithConstraints(double updateBand)
{
    if (mInterrupter) mInterrupter->start("SetOnOffAux1WithConstraints");
    mTask = boost::bind(&SmoothingFilter::doSetOnOffAux1WithConstraints, _1, _2, updateBand);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


// Set voxels with read buffer values less than update band "On" all other "off"
void
SmoothingFilter::doSetOnOffRead(const RangeType& range, double updateBand)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        for (VoxelIterTAll vIterAll = mLeafs->leaf(n).beginValueAll(); vIterAll; ++vIterAll) {
            if(fabs(vIterAll.getValue()) <= updateBand) vIterAll.setValueOn();
            else vIterAll.setValueOff();
        }
    }
}


// Set voxels with Aux1 buffer values less than update band "On" all other "off"
void
SmoothingFilter::doSetOnOffAux1(const RangeType& range, double updateBand)
{
    this->wasInterrupted();
		int on = 0, off = 0;
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& bufferid = mLeafs->getBuffer(n, Buffers::Aux1);
        for (VoxelIterTAll vIterAll = mLeafs->leaf(n).beginValueAll(); vIterAll; ++vIterAll) {
					if(fabs(bufferid.getValue(vIterAll.pos())) <= updateBand) {vIterAll.setValueOn(); on++;}
					else {vIterAll.setValueOff(); off++;}
        }
    }
		//std::cout<<"on = "<<on<<" off = "<<off<<std::endl;
}


// Set voxels with (Aux1 buffer values less than update band) && (not constrained) "On" all other "off"
void
SmoothingFilter::doSetOnOffAux1WithConstraints(const RangeType& range, double updateBand)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& bufferid = mLeafs->getBuffer(n, Buffers::Aux1);
        BufferType& activeSet = mLeafs->getBuffer(n, Buffers::Activeset);
        for (VoxelIterTAll vIterAll = mLeafs->leaf(n).beginValueAll(); vIterAll; ++vIterAll) {
            Index c = vIterAll.pos();
            if(fabs(bufferid.getValue(c)) <= updateBand && activeSet.getValue(c) == 0.0) vIterAll.setValueOn();
            else vIterAll.setValueOff();
        }
    }
    return;
}




// Explicit Integration Methods

// take ntimesteps of laplacian smoothing
void
SmoothingFilter::laplacianSmoothing(double dt, int ntimesteps) 
{
    double updateBand = 4*mDx;
    for (int i = 0; i < ntimesteps; i++) {
        setOnOffRead(updateBand);
        laplacianStep(dt); // this automatically swaps the result into the read buffer
    }
		redistance(Buffers::Read);
}


// take n timesteps of biharmonic smoothing, redistancing ever redistanceFrequency timesteps
void SmoothingFilter::biharmonicSmoothing(double dt, int ntimesteps, int redistanceFrequency) 
{
    double updateBandBiharmonic = 3*mDx, updateBandLaplacian = 4*mDx;

    for (int i = 0; i < ntimesteps; i++) {
        setOnOffRead(updateBandLaplacian);
        computeLaplacianAndGradient();

				setOnOffRead(updateBandBiharmonic);
				mLeafs->swapLeafBuffer(Buffers::Laplacian, mGrainSize==0);
        biharmonicStep(dt); // this automatically swaps the result into the read buffer

        if(i%redistanceFrequency == 0 && i > 0) {
					  redistance(Buffers::Read);
        }
    }
}


// tbb wrappers
// Note that the first two methods place the result in the read buffer
void
SmoothingFilter::laplacianStep(double dt)
{
    if (mInterrupter) mInterrupter->start("LaplacianStep");
    mTask = boost::bind(&SmoothingFilter::doLaplacianStep, _1, _2, dt);
    this->cook(Buffers::Laplacian);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::biharmonicStep(double dt)
{
    if (mInterrupter) mInterrupter->start("BiharmonicStep");
    mTask = boost::bind(&SmoothingFilter::doBiharmonicStep, _1, _2, dt);
    //this->cook(Buffers::Biharmonic);
    this->cook(Buffers::Aux2);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::computeLaplacianAndGradient()
{
    if (mInterrupter) mInterrupter->start("ComputeLaplacianAndGradient");
    mTask = boost::bind(&SmoothingFilter::doComputeLaplacian, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


// take a single step of laplacian smoothing.  The result is placed int he "laplacian" buffer.  (It is swapped to the read buffer above)
void
SmoothingFilter::doLaplacianStep(const RangeType& range, double dt)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& lapBuffer = mLeafs->getBuffer(n, Buffers::Laplacian);
        BufferType& minBuffer = mLeafs->getBuffer(n, Buffers::MinPhi);
        BufferType& maxBuffer = mLeafs->getBuffer(n, Buffers::MaxPhi);

        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            Index c = vIter.pos();
            stencil.moveTo(vIter);
            const double P = stencil.getValue();
            const double L = stencil.laplacian();
            const double minPhi = minBuffer.getValue(c);
            const double maxPhi = maxBuffer.getValue(c);
            math::Vec3<ValueType>  grad = stencil.gradient();
            double updatedPhi = P + L*dt*grad.length();
            updatedPhi = std::min(updatedPhi, minPhi);
            updatedPhi = std::max(updatedPhi, maxPhi);
            lapBuffer.setValue(c, updatedPhi);
        }
    }
}


// take a single step of laplacian smoothing.  The result is placed int the "Biharmonic" buffer.  (It is swapped to the read buffer above)
void
SmoothingFilter::doBiharmonicStep(const RangeType& range, double dt)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
  			// Laplacian buffer currently has the previous phi values, read buffer has laplacian values
  			BufferType& lapBuffer = mLeafs->getBuffer(n, Buffers::Laplacian); 
        BufferType& gradBuffer = mLeafs->getBuffer(n, Buffers::GradMag);
				BufferType& auxBuffer = mLeafs->getBuffer(n, Buffers::Aux2);
				//BufferType& biBuffer = mLeafs->getBuffer(n, Buffers::Biharmonic);
        BufferType& minBuffer = mLeafs->getBuffer(n, Buffers::MinPhi);
        BufferType& maxBuffer = mLeafs->getBuffer(n, Buffers::MaxPhi);

				// Iterate over "On" voxels
        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            stencil.moveTo(vIter); //laplacian
            Index c = vIter.pos();
            const double P = lapBuffer.getValue(c); //phi value
            const double gradMag = gradBuffer.getValue(c);
            const double minPhi = minBuffer.getValue(c);
            const double maxPhi = maxBuffer.getValue(c);
            const double B = stencil.laplacian(); //biharmonic
            double updatedPhi = P - B*dt*gradMag;
            updatedPhi = std::min(updatedPhi, minPhi);
            updatedPhi = std::max(updatedPhi, maxPhi);
            //biBuffer.setValue(c, updatedPhi);
            auxBuffer.setValue(c, updatedPhi);
        }
    }

		// Copy old phi values to "Off" voxels
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& lapBuffer = mLeafs->getBuffer(n, Buffers::Laplacian);
        //BufferType& biBuffer = mLeafs->getBuffer(n, Buffers::Biharmonic);
        BufferType& auxBuffer = mLeafs->getBuffer(n, Buffers::Aux2);
        for (VoxelIterTOff vOffIter = mLeafs->leaf(n).beginValueOff(); vOffIter; ++vOffIter) {
            auxBuffer.setValue(vOffIter.pos(), lapBuffer.getValue(vOffIter.pos()));
        }
    }
}


// Compute the laplacian and magnitude of the gradient for later use in doBiharmonicStep()
void
SmoothingFilter::doComputeLaplacianAndGradient(const RangeType& range)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& lapBuffer = mLeafs->getBuffer(n, Buffers::Laplacian);
        BufferType& gradBuffer = mLeafs->getBuffer(n, Buffers::GradMag);
        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            Index c = vIter.pos();
            stencil.moveTo(vIter);
            const double L = stencil.laplacian();
            lapBuffer.setValue(c, L);
            const double G = stencil.gradient().length();
            gradBuffer.setValue(c, G);
        }
    }
}


// implicit integration methods

// main driver routine, performs ntimesteps steps of size dt.  Each step calls CG until it converges or uses too many iterations.
// CG automatically exits if constraints are violated, and this routine automatically restarts it with the new constraints.
void
SmoothingFilter::biharmonicSmoothing(double dt, int ntimesteps, int redistanceFrequency, double cgThreshold, int cgMaxIters, bool verbose)
{
	  for (int i = 0; i < ntimesteps; i++) {
        int cgIters = 0;
        double tol = -1.0;
        initializeImplicitBiharmonic();

				while (!CG(dt, tol, cgIters, cgThreshold, cgMaxIters, verbose)) {
					// iterate until true
				}
				
        copySolutionToReadInBand();

        if(i % redistanceFrequency == 0) redistance(Buffers::Read);

        if(verbose) std::cout<<"At "<<i<<" th timestep"<<" "<<cgIters<<" iterations."<<std::endl;
    }
}


// Main CG routine
bool
SmoothingFilter::CG(double dt, double &tol, int &cgIters, double cgThreshold, int cgMaxIters, bool verbose)
{
	  // initialization stuff.  Apply the operator to the current field and get the initial residual
	  double absoluteThreshold = 1E-05;
    double updatebandL = 4*mDx, updatebandB = 3*mDx;
    zeroBuffers();

    mLeafs->swapLeafBuffer(Buffers::X, mGrainSize==0);
    setOnOffAux1(updatebandL);
    computeLaplacian();
    mLeafs->swapLeafBuffer(Buffers::X, mGrainSize==0);

    mLeafs->swapLeafBuffer(Buffers::Laplacian, mGrainSize==0);
    setOnOffAux1WithConstraints(updatebandB);
    applyBiharmonic(dt, Buffers::X);
    mLeafs->swapLeafBuffer(Buffers::Laplacian, mGrainSize==0);

    computeInitialResidual();

    mDelta = 0.0;
    innerProduct(Buffers::R, Buffers::R);
    double deltaNew = mDelta;

		// we get the tolerance from the first call to CG for each timestep
    if(tol < 0.0) tol = cgThreshold*cgThreshold*deltaNew;

    while(cgIters++ < cgMaxIters && deltaNew > tol && deltaNew > absoluteThreshold) {
  			// apply operator
			  mLeafs->swapLeafBuffer(Buffers::D, mGrainSize==0);
        setOnOffAux1(updatebandL);
        computeLaplacian();
        mLeafs->swapLeafBuffer(Buffers::D, mGrainSize==0);

        mLeafs->swapLeafBuffer(Buffers::Laplacian, mGrainSize==0);
        setOnOffAux1WithConstraints(updatebandB);
        applyBiharmonic(dt, Buffers::D);
        mLeafs->swapLeafBuffer(Buffers::Laplacian, mGrainSize==0);

				// compute stepsize (alpha) and take a step
        mDelta = 0.0;
        innerProduct(Buffers::D, Buffers::Q);
        double alpha = deltaNew/mDelta;
        CGStep(alpha);

				// check if the solution violates any constraints
				bool constraintsViolated = false;
        setOnOffAux1WithConstraints(updatebandB);
				applyConstraints(constraintsViolated);
				if (constraintsViolated) return false;

				// update residual
        double deltaOld=deltaNew;
        setOnOffAux1WithConstraints(updatebandB);
        mDelta = 0.0;
        innerProduct(Buffers::R, Buffers::R);
        deltaNew = mDelta;

				// compute next conjugate direction
        double beta = deltaNew/deltaOld;
        setOnOffAux1WithConstraints(updatebandB);
        CGDirection(beta);

				if(verbose) std::cout<<"At iteration "<<cgIters<<": deltaNew = "<<deltaNew<<"."<<std::endl;

    }
    return true;
}


// wrapper functions.  None of these do any automatic swaps.
void
SmoothingFilter::initializeImplicitBiharmonic()
{
    if (mInterrupter) mInterrupter->start("InitializeImplicitBiharmonic");
    mTask = boost::bind(&SmoothingFilter::doInitializeImplicitBiharmonic, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::zeroBuffers()
{
    if (mInterrupter) mInterrupter->start("ZeroBuffers");
    mTask = boost::bind(&SmoothingFilter::doZeroBuffers, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::computeLaplacian()
{
    if (mInterrupter) mInterrupter->start("ComputeLaplacian");
    mTask = boost::bind(&SmoothingFilter::doComputeLaplacian, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::applyBiharmonic(double dt, int sourceBuffer)
{
    if (mInterrupter) mInterrupter->start("ApplyBiharmonic");
    mTask = boost::bind(&SmoothingFilter::doApplyBiharmonic, _1, _2, dt, sourceBuffer);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::computeInitialResidual()
{
    if (mInterrupter) mInterrupter->start("ComputeInitialResidual");
    mTask = boost::bind(&SmoothingFilter::doComputeInitialResidual, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::CGStep(double alpha) 
{
    if (mInterrupter) mInterrupter->start("CGStep");
    mTask = boost::bind(&SmoothingFilter::doCGStep, _1, _2, alpha);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::CGDirection(double sigma) 
{
    if (mInterrupter) mInterrupter->start("CGDirection");
    mTask = boost::bind(&SmoothingFilter::doCGDirection, _1, _2, sigma);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::applyConstraints(bool &constraintsViolated)
{
    if (mInterrupter) mInterrupter->start("ApplyConstraints");
    mTask = boost::bind(&SmoothingFilter::doApplyConstraints, _1, _2, boost::ref(constraintsViolated));
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


// parallel reduce
void
SmoothingFilter::innerProduct(int xid, int yid)
{
    if (mInterrupter) mInterrupter->start("InnerProduct");
    mTask = boost::bind(&SmoothingFilter::doInnerProduct, _1, _2, xid, yid);
    int tmp = mGrainSize;
    mGrainSize = 2;
    this->cook(-1);
    mGrainSize = tmp;
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::copySolutionToReadInBand()
{
    if (mInterrupter) mInterrupter->start("CopySolutionToReadInBand");
    mTask = boost::bind(&SmoothingFilter::doCopySolutionToReadInBand, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


// This initializes the solution (X) to the current phi field and creates a backup in Aux1.
// It also zeroes out the laplacian and activeset.
// There may be some redundant/unnecessary zeroing here.
void
SmoothingFilter::doInitializeImplicitBiharmonic(const RangeType &range) 
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& lapBuffer = mLeafs->getBuffer(n, Buffers::Laplacian);
        BufferType& xBuffer = mLeafs->getBuffer(n, Buffers::X);
        BufferType& auxBuffer = mLeafs->getBuffer(n, Buffers::Aux1);
        BufferType& activeSet = mLeafs->getBuffer(n, Buffers::Activeset);

        for (VoxelIterTAll iter = mLeafs->leaf(n).beginValueAll(); iter; ++iter) {
            Index c = iter.pos();
            double val=iter.getValue();
            auxBuffer.setValue(c, val);
            xBuffer.setValue(c, val);
            lapBuffer.setValue(c, 0.0);
            activeSet.setValue(c, 0.0);
        }
    }
}


// This is called at the beginning of CG.  It zeros out several buffers: Laplacian, D, R, and Q.
// There may be some redundant/unnecessary zeroing here.
void
SmoothingFilter::doZeroBuffers(const RangeType &range)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& lapBuffer = mLeafs->getBuffer(n, Buffers::Laplacian);
        BufferType& dBuffer = mLeafs->getBuffer(n, Buffers::D);
        BufferType& rBuffer = mLeafs->getBuffer(n, Buffers::R);
        BufferType& qBuffer = mLeafs->getBuffer(n, Buffers::Q);

        for (VoxelIterTAll iter = mLeafs->leaf(n).beginValueAll(); iter; ++iter) {
            Index c = iter.pos();
            lapBuffer.setValue(c, 0.0);
            dBuffer.setValue(c, 0.0);
            rBuffer.setValue(c, 0.0);
            qBuffer.setValue(c, 0.0);
        }
    }
}


// This computes the intial residual after applying the operator to the
// initial solution.
// There may be some redundant/unnecessary zeroing here.
void
SmoothingFilter::doComputeInitialResidual(const RangeType& range)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& auxBuffer = mLeafs->getBuffer(n, Buffers::Aux1);
        BufferType& dBuffer = mLeafs->getBuffer(n, Buffers::D);
        BufferType& rBuffer = mLeafs->getBuffer(n, Buffers::R);
        BufferType& qBuffer = mLeafs->getBuffer(n, Buffers::Q);
        BufferType& lapBuffer = mLeafs->getBuffer(n, Buffers::Laplacian);

        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            Index c = vIter.pos();
            double val = auxBuffer.getValue(c)-qBuffer.getValue(c); // phi_0 - q
            dBuffer.setValue(c, val);
            rBuffer.setValue(c, val);
            qBuffer.setValue(c, 0.0);
            lapBuffer.setValue(c, 0.0);
        }
        for (VoxelIterTOff vIter = mLeafs->leaf(n).beginValueOff(); vIter; ++vIter) {
            Index c = vIter.pos();
            lapBuffer.setValue(c, 0.0);
            qBuffer.setValue(c, 0.0);
        }
    }
}


// compute the laplacian and store it in the laplacian buffer
void
SmoothingFilter::doComputeLaplacian(const RangeType& range)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
			  BufferType& lapBuffer = mLeafs->getBuffer(n, Buffers::Laplacian);
        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            stencil.moveTo(vIter);
            const double L = stencil.laplacian();
            lapBuffer.setValue(vIter.pos(), L);
        }
    }

    return;
}


// compute the full operator (I + dt*biharmonic) and store it in Q.
// We may apply this operator to either X or D, so we take a parameter to tell us which sourceBuffer
// for the identity part of the operator.  The rest is computed from the laplacian, which is in the read buffer.
void
SmoothingFilter::doApplyBiharmonic(const RangeType& range, double dt, int sourceBuffer)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& srcBuffer = mLeafs->getBuffer(n, sourceBuffer);
        BufferType& qBuffer = mLeafs->getBuffer(n, Buffers::Q);
        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            Index c = vIter.pos();
            stencil.moveTo(vIter); //laplacian
            qBuffer.setValue(c, srcBuffer.getValue(c)+dt*stencil.laplacian());
        }
    }
    return;
}


// Take a step of length alpha in direction D and add it to X.
// Also, update the residual by a step of length alpha in direction Q.
void 
SmoothingFilter::doCGStep(const RangeType &range, double alpha) 
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& dBuffer = mLeafs->getBuffer(n, Buffers::D);
        BufferType& rBuffer = mLeafs->getBuffer(n, Buffers::R);
        BufferType& xBuffer = mLeafs->getBuffer(n, Buffers::X);
        BufferType& qBuffer = mLeafs->getBuffer(n, Buffers::Q);

        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            Index c = vIter.pos();
            double xval = xBuffer.getValue(c);
						double rval = rBuffer.getValue(c);
            xBuffer.setValue(c, xval + alpha*dBuffer.getValue(c));
            rBuffer.setValue(c, rval - alpha*qBuffer.getValue(c));
        }
    }

    return;
}


// compute a new direction for cg
void 
SmoothingFilter::doCGDirection(const RangeType &range, double beta)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& dBuffer = mLeafs->getBuffer(n, Buffers::D);
        BufferType& rBuffer = mLeafs->getBuffer(n, Buffers::R);
        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            Index c = vIter.pos();
            double pval = dBuffer.getValue(c);
            dBuffer.setValue(c, rBuffer.getValue(c) + beta*pval);
        }
    }
    return;
}


// Apply the min/max constraints.  If any are violated update constraintsViolated.
void
SmoothingFilter::doApplyConstraints(const RangeType &range, bool &constraintsViolated) 
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& xBuffer = mLeafs->getBuffer(n, Buffers::X);
        BufferType& minBuffer = mLeafs->getBuffer(n, Buffers::MinPhi);
        BufferType& maxBuffer = mLeafs->getBuffer(n, Buffers::MaxPhi);
        BufferType& activeSet = mLeafs->getBuffer(n, Buffers::Activeset);

        for (VoxelIterTOn vIter = mLeafs->leaf(n).beginValueOn(); vIter; ++vIter) {
            Index c = vIter.pos();
            double xBufferVal = xBuffer.getValue(c);
            double minPhi = minBuffer.getValue(c);
            double maxPhi = maxBuffer.getValue(c);

            if(xBufferVal > minPhi) {
							  constraintsViolated = true;
                activeSet.setValue(c, 1.0);
                xBuffer.setValue(c, minPhi);
            } else if(xBufferVal < maxPhi) {
							  constraintsViolated = true;
							  activeSet.setValue(c, 1.0);
                xBuffer.setValue(c, maxPhi);
            }
        }
    }
    return;
}


// compute an inner product between buffers xid and yid.
void
SmoothingFilter::doInnerProduct(const RangeType &range, int xid, int yid)
{
    this->wasInterrupted();
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& xBuffer = mLeafs->getBuffer(n, xid);
        BufferType& yBuffer = mLeafs->getBuffer(n, yid);
        for (VoxelIterTOn iter = mLeafs->leaf(n).beginValueOn(); iter; ++iter) {
					  Index c = iter.pos();
            double xval = xBuffer.getValue(c);
						double yval = yBuffer.getValue(c);
            mDelta += xval * yval;
        }
    }
    return;
}


// Copy the solution (X) into the read buffer inside of the updateBand.  This
// effectively updates the solution to the next time level.
void 
SmoothingFilter::doCopySolutionToReadInBand(const RangeType &range)
{
    this->wasInterrupted();
    double updateband = 3*mDx;
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& bufferx = mLeafs->getBuffer(n, Buffers::X);
        BufferType& bufferdup = mLeafs->getBuffer(n, Buffers::Aux1);
        for (VoxelIterTAll iter = mLeafs->leaf(n).beginValueAll(); iter; ++iter) {
            if(fabs(bufferdup.getValue(iter.pos())) <= updateband) iter.setValue(bufferx.getValue(iter.pos()));
        }
    }
    return;
}


// code to extract surface and write a triangulation to file.
typedef boost::scoped_array<openvdb::Vec3s> PointList;
typedef boost::scoped_array<openvdb::tools::PolygonPool> PolygonPoolList;

void
SmoothingFilter::extractSurface(char *outfname) 
{
    openvdb::tools::VolumeToMesh mesher(0.0);
    mesher.operator()<openvdb::DoubleGrid>(*mGrid);
    std::ofstream out(outfname, std::ios::out);

    PointList& p = mesher.pointList();
    for (unsigned int i = 0; i < mesher.pointListSize(); i++) {
        out<<"v "<<p[i][0]<<" "<<p[i][1]<<" "<<p[i][2]<<std::endl;
    }

    openvdb::tools::PolygonPoolList& polygonPoolList = mesher.polygonPoolList();
    for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
        const openvdb::tools::PolygonPool& polygons = polygonPoolList[n];
        for (size_t i = 0, I = polygons.numQuads(); i < I; ++i) {
            const openvdb::Vec4I& quad = polygons.quad(i);
            out<<"f "<<(quad[0]+1)<<" "<<(quad[3]+1)<<" "<<(quad[2]+1)<<std::endl;
            out<<"f "<<(quad[0]+1)<<" "<<(quad[2]+1)<<" "<<(quad[1]+1)<<std::endl;
        }
    }
}


// Begin Redistancing code
// Buffer Aux2 will hold the field we are redistancing.  We lock values next to the interface by setting them to OFF.

void SmoothingFilter::redistance(int buf)
{
    if (mInterrupter) mInterrupter->start("Redistance");
    if(buf > 0) mLeafs->swapLeafBuffer(buf, mGrainSize==0);

    copyAllValuesFromReadToAux();
    redistanceInterface(); 
    copyLockedValuesFromAuxToRead(); // these are the only values that changed in aux
		
    for (int i = 0; i < 8; i++) {
			  redistanceSweep(); // sweeps over read, placing new values in aux
        copyUnlockedValuesFromAuxToRead(); // copies new values back to read
				//mLeafs->swapLeafBuffer(Buffers::Aux2, mGrainsSize==0);
    }

    if(buf > 0)
        mLeafs->swapLeafBuffer(buf, mGrainSize==0);
    if (mInterrupter) mInterrupter->end();
}


// wrapper routines
void
SmoothingFilter::redistanceInterface()
{
    if (mInterrupter) mInterrupter->start("RedistanceInterface");
    mTask = boost::bind(&SmoothingFilter::doRedistanceInterface, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();

}


void
SmoothingFilter::redistanceSweep()
{
    if (mInterrupter) mInterrupter->start("RedistanceSweep");
    mTask = boost::bind(&SmoothingFilter::doRedistanceSweep, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();

}

void
SmoothingFilter::copyAllValuesFromReadToAux() 
{
    if (mInterrupter) mInterrupter->start("copyAllValuesFromReadToAux");
    mTask = boost::bind(&SmoothingFilter::doCopyAllValuesFromReadToAux, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::copyLockedValuesFromAuxToRead() 
{
    if (mInterrupter) mInterrupter->start("copyLockedValuesFromAuxToRead");
    mTask = boost::bind(&SmoothingFilter::doCopyLockedValuesFromAuxToRead, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


void
SmoothingFilter::copyUnlockedValuesFromAuxToRead()
{
    if (mInterrupter) mInterrupter->start("copyUnlockedValuesFromAuxToRead");
    mTask = boost::bind(&SmoothingFilter::doCopyUnlockedValuesFromAuxToRead, _1, _2);
    this->cook(-1);
    if (mInterrupter) mInterrupter->end();
}


// fastSweepingInterfaceHelper helps initialize values next to the interface
void
SmoothingFilter::fastSweepingInterfaceHelper(size_t &n, VoxelIterTOn &vIterOn, Index &cI, 
                          double &x, double &y)
{
    this->wasInterrupted();
    double pointnineh = 0.9*mDx;
    if (sign(x) != sign(y)) {
			  // If the gradient is close to one, make it one.  Otherwise copy the old value.
			  if (fabs(x)+fabs(y)>pointnineh) {
            mLeafs->getBuffer(n, Buffers::Aux2).setValue(cI, absmin(mLeafs->getBuffer(n, Buffers::Aux2).getValue(cI), mDx*x/fabs(x-y)));
        } else {
					  mLeafs->getBuffer(n, Buffers::Aux2).setValue(cI, absmin(mLeafs->getBuffer(n, Buffers::Aux2).getValue(cI), x));
        }
				vIterOn.setValueOff();
    }
}


// just redistances values next to the interface, looks in all six directions for a sign change
void
SmoothingFilter::doRedistanceInterface(const RangeType &range)
{
    this->wasInterrupted();
    VoxelIterTOn vIterOn;
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        for (vIterOn = mLeafs->leaf(n).beginValueOn(); vIterOn; ++vIterOn) {
            Coord c = vIterOn.getCoord();
            Index cI = vIterOn.pos();
            double x = vIterOn.getValue();
            double y = accessor.getValue(c.offsetBy(1, 0, 0));
            fastSweepingInterfaceHelper(n, vIterOn, cI, x, y);
            y = accessor.getValue(c.offsetBy(0, 1, 0));
            fastSweepingInterfaceHelper(n, vIterOn, cI, x, y);
            y = accessor.getValue(c.offsetBy(0, 0, 1));
            fastSweepingInterfaceHelper(n, vIterOn, cI, x, y);
            y = accessor.getValue(c.offsetBy(-1, 0, 0));
            fastSweepingInterfaceHelper(n, vIterOn, cI, x, y);
            y = accessor.getValue(c.offsetBy(0, -1, 0));
            fastSweepingInterfaceHelper(n, vIterOn, cI, x, y);
            y = accessor.getValue(c.offsetBy(0, 0, -1));
            fastSweepingInterfaceHelper(n, vIterOn, cI, x, y);
        }
    }
}


// sweeps over the entire grid updating unlocked values.
void
SmoothingFilter::doRedistanceSweep(const RangeType &range)
{
    this->wasInterrupted();
    size_t n, e=range.end();
    VoxelIterTOn vIterOn;
    double background =  mGrid->background();
    for (n=range.begin();n != e; ++n) {
        for (vIterOn = mLeafs->leaf(n).beginValueOn(); vIterOn; ++vIterOn) {
            Coord cC = vIterOn.getCoord();
            double a = std::min<double>(fabs(accessor.getValue(cC.offsetBy(-1, 0, 0))), fabs(accessor.getValue(cC.offsetBy(1, 0, 0))));
            double b = std::min<double>(fabs(accessor.getValue(cC.offsetBy(0, -1, 0))), fabs(accessor.getValue(cC.offsetBy(0, 1, 0))));
            double c = std::min<double>(fabs(accessor.getValue(cC.offsetBy(0, 0, -1))), fabs(accessor.getValue(cC.offsetBy(0, 0, 1))));
            double temp;
            if(a > b) { temp=a; a=b; b=temp; }
            if(a > c) { temp=a; a=c; c=temp; }
            if(b > c) { temp=b; b=c; c=temp; }
            if(math::isExactlyEqual(a, background)) continue;
            double x = a + mDx;
            if (x > b) {
                x = 0.5*(a+b+math::Sqrt(2*math::Pow2(mDx)-math::Pow2(a-b)));
                if (x > c) {
                    x = (a+b+c+math::Sqrt(3*math::Pow2(mDx)-math::Pow2(a-b)-math::Pow2(b-c)-math::Pow2(c-a)))/3.0;
                }
            }
            mLeafs->getBuffer(n, Buffers::Aux2).setValue(vIterOn.pos(), sign(vIterOn.getValue())*x);

        }
    }
}


// copy all values from the read buffer to Aux2
void
SmoothingFilter::doCopyAllValuesFromReadToAux(const RangeType &range)
{
    this->wasInterrupted();
    VoxelIterTAll vIterAll;
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
			  BufferType& auxBuffer = mLeafs->getBuffer(n, Buffers::Aux2);		
        for (vIterAll = mLeafs->leaf(n).beginValueAll(); vIterAll; ++vIterAll) {
            vIterAll.setValueOn();
            auxBuffer.setValue(vIterAll.pos(), vIterAll.getValue());
        }
    }
}


// copy the locked values from Aux2 to the read buffer
void
SmoothingFilter::doCopyLockedValuesFromAuxToRead(const RangeType &range)
{
    this->wasInterrupted();
    VoxelIterTOff vIterOff;
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& auxBuffer = mLeafs->getBuffer(n, Buffers::Aux2);
        for (vIterOff = mLeafs->leaf(n).beginValueOff(); vIterOff; ++vIterOff) {
            vIterOff.setValue(auxBuffer.getValue(vIterOff.pos()));
        }
    }
}


// copy (updated) unlocked values from Aux2 to the read buffer
void
SmoothingFilter::doCopyUnlockedValuesFromAuxToRead(const RangeType &range) 
{
    this->wasInterrupted();
    VoxelIterTOn vIterOn;
    for (size_t n=range.begin(), e=range.end(); n != e; ++n) {
        BufferType& auxBuffer = mLeafs->getBuffer(n, Buffers::Aux2);
        for (vIterOn = mLeafs->leaf(n).beginValueOn(); vIterOn; ++vIterOn) {
            vIterOn.setValue(auxBuffer.getValue(vIterOn.pos()));
        }
    }
}


