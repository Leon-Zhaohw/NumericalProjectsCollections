///////////////////////////////////////////////////////////////////////////////////////////////////
// I. LICENSE CONDITIONS
//
// Copyright (c) 2018 by Disney-Pixar
//
// Permission is hereby granted to use this software solely for non-commercial applications
// and purposes including academic or industrial research, evaluation and not-for-profit media
// production. All other rights are retained by Pixar. For use for or in connection with
// commercial applications and purposes, including without limitation in or in connection with
// software products offered for sale or for-profit media production, please contact Pixar at
// tech-licensing@pixar.com.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
// NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TET_NEWTON_SOLVER_H
#define TET_NEWTON_SOLVER_H

#include "EigenCG.h"
#include "LineSearches.h"
#include "Settings.h"

namespace CubeSim
{

class Material;

struct TetNewtonResult final
{
    Scalar forceResidual;
    int newtonIterationCount;
    int cgIterationCount;
    bool succeeded;
};

template<class T>
class TetNewtonSolver final
{

public:

    // Resizes the internal storage. Must be called before Solve.
    void Resize(const int ndofs)
    {
        _forces.conservativeResize(ndofs);
        _minimizedStep.conservativeResize(ndofs);
        _delta.conservativeResize(ndofs);
        // _H.conservativeResize(ndofs, ndofs);
    }

    SparseMatrixX& GetH()
    {
        return _H;
    }

    // Initializes the linear solve. Must be called before Solve.
    void InitLinearSolver(const int maxIters, const Scalar& tol)
    {
        _cgSolver.setMaxIterations(maxIters);
        _cgSolver.setTolerance(tol);
    }

    void ComputeDelta(const SparseMatrixX& A, const VectorX& b, VectorX& x, int& cgIters)
    {
        _cgSolver.compute(A);
        x = _cgSolver.solve(b);
        cgIters = int(_cgSolver.iterations());
    }

    TetNewtonResult Solve(const bool printStatus, const Material& material, const Scalar& forceTol, const int minIters, const int maxIters, T& mesh)
    {
        if (printStatus)
        {
            std::cout << "Newton Solve:" << std::endl;
        }

        mesh.ComputeForce(material, _forces);

        TetNewtonResult results;
        results.newtonIterationCount = 0;
        results.cgIterationCount = 0;
        results.forceResidual = _forces.dot(_forces); // forces.lpNorm<Eigen::Infinity>();
        results.succeeded = false;

        if (printStatus)
        {
            std::cout << "    Initial residual: " << results.forceResidual << std::endl;
        }

        while ((results.newtonIterationCount < minIters) || (results.newtonIterationCount < maxIters && results.forceResidual > forceTol))
        {
            if (printStatus)
            {
                std::cout << "    Iteration " << results.newtonIterationCount << std::endl;
            }
            mesh.ComputeHessian(material, _H);

            int numCgItersInStep = 0;
            ComputeDelta(_H, _forces, _delta, numCgItersInStep);
            results.cgIterationCount += numCgItersInStep;
            if (printStatus)
            {
                std::cout << "        CG Iterations: " << results.cgIterationCount << std::endl;
            }

            Scalar alphaMin;
            Scalar Umin;
            LineSearches::MinimizeInSearchDirection(mesh, material, _delta, _minimizedStep, alphaMin, Umin);

            if (printStatus)
            {
                std::cout << "        Alpha: " << alphaMin << std::endl;
                std::cout << "        Energy: " << Umin << std::endl;
            }

            mesh.PerturbSimulatedVertices(_minimizedStep);
            mesh.SwapStateCaches();
            mesh.UpdateCurrentCachedState();

            mesh.ComputeForce(material, _forces);
            results.forceResidual = _forces.dot(_forces); // forces.lpNorm<Eigen::Infinity>();
            if (printStatus)
            {
                std::cout << "        Residual: " << results.forceResidual << std::endl;
            }
            results.newtonIterationCount++;
        }

        results.succeeded = results.forceResidual <= forceTol;
        return results;
    }

private:

    VectorX _forces;
    VectorX _minimizedStep;
    VectorX _delta;
    SparseMatrixX _H;

    Eigen::ConjugateGradient<SparseMatrixX, Eigen::Lower|Eigen::Upper> _cgSolver;

};

}

#endif
