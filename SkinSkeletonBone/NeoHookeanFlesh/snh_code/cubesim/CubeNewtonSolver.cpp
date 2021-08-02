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

#include "CubeNewtonSolver.h"

#include "CubeMesh.h"
#include "Material.h"
#include "CGInterface.h"

#include <iostream>

namespace CubeSim
{

CubeNewtonSolver::CubeNewtonSolver(const int ndofs)
: _u(ndofs)
, _uOld(ndofs)
, _uDelta(ndofs)
, _f(ndofs)
, _K(ndofs, ndofs)
{}

SparseMatrixX& CubeNewtonSolver::GetK()
{
    return _K;
}

QuasistaticSolveResult CubeNewtonSolver::SolveQuasistatic(const bool printStatus, const int maxIterations, const Scalar& tol, const Material& material, const CubeMesh& cubeMesh, CGInterface& cg, DynamicState& dynamicState)
{
    cubeMesh.ComputeDisplacements(dynamicState, _u);

    int newtonIterations = 0;
    int totalCGIterations = 0;
    Scalar residual = std::numeric_limits<Scalar>::max();

    while (true)
    {
        if (printStatus)
        {
            std::cout << "  Iteration: " << newtonIterations << std::endl;
        }

        // Build (minus) the right hand side
        cubeMesh.ComputeForces(dynamicState, material, _f);

        // Check the force residual for termination
        residual = _f.dot(_f);
        if (printStatus)
        {
            std::cout << "    Residual: " << residual << std::endl;
        }
        if (residual < tol)
        {
            break;
        }

        // If we've maxed out the number of acceptable iterations
        if (newtonIterations >= maxIterations)
        {
            break;
        }

        // Compute the stiffness matrix
        if (printStatus)
        {
            std::cout << "    Computing K ... " << std::flush;
        }
        cubeMesh.ComputeStiffnessMatrixSparse(material, dynamicState, _K);
        if (printStatus)
        {
            std::cout << "done. " << std::endl;
        }

        // Solve the system
        if (printStatus)
        {
            std::cout << "    Solving CG ... " << std::flush;
        }
        _uDelta.setZero();
        cg.Solve(_K, -_f, _uDelta);
        if (printStatus)
        {
            std::cout << "done." << std::endl;
            std::cout << "      Iterations: " << cg.Iterations() << " \t error: " << cg.Error() << std::endl;
        }
        totalCGIterations += cg.Iterations();

        if (printStatus && !cg.LastSolveSucceeded())
        {
            std::cout << RED_C << "    CG failed to converge with residual: " << (_K * _uDelta + _f).lpNorm<Eigen::Infinity>() << NO_C << std::endl;
        }

        // Line search
        {
            if (printStatus)
            {
                std::cout << "    Minimizing in search direction ... " << std::flush;
            }
            Scalar alphaMin = _MinimizeInSearchDirection(_uDelta, material, cubeMesh, dynamicState);
            if (printStatus)
            {
                std::cout << "done." << std::endl;
                std::cout << "      Alpha: " << alphaMin << std::endl;
            }
            if (alphaMin == 0.0)
            {
                if (printStatus)
                {
                    std::cout << RED_C << "    Newton failed to make progress, halting solve." << NO_C << std::endl;
                }
                return QuasistaticSolveResult(false, residual, newtonIterations, totalCGIterations);
            }
            _u += alphaMin * _uDelta;
            cubeMesh.ScatterDisplacements(_u, dynamicState);
        }

        newtonIterations++;
    }

    if (newtonIterations == maxIterations && residual >= tol)
    {
        if (printStatus)
        {
            std::cout << RED_C << "  Newton did not converge." << NO_C << std::endl;
        }
        return QuasistaticSolveResult(false, residual, newtonIterations, totalCGIterations);
    }

    if (printStatus)
    {
        std::cout << GREEN_C << "  Newton iterations: " << newtonIterations << NO_C << std::endl;
    }
    return QuasistaticSolveResult(true, residual, newtonIterations, totalCGIterations);
}

Scalar CubeNewtonSolver::_MinimizeInSearchDirection(const VectorX& delta, const Material& material, const CubeMesh& cubeMesh, DynamicState& dynamicState)
{
    _uOld = _u;

    // Middle point in bracket
    Scalar alphaMiddle = 0.0;
    Scalar Umiddle = cubeMesh.ComputeEnergy(dynamicState, material);

    // Right bracket on minimum
    _u = _uOld + delta;
    cubeMesh.ScatterDisplacements(_u, dynamicState);
    Scalar alphaRight = 1.0;
    Scalar Uright = cubeMesh.ComputeEnergy(dynamicState, material);

    // Left bracket on minimum
    Scalar alphaLeft;
    Scalar Uleft;

    if (Umiddle > Uright)
    {
        Uleft = Umiddle;
        alphaLeft = alphaMiddle;
        Umiddle = Uright;
        alphaMiddle = alphaRight;

        // Search for an alpha with energy greater than the middle
        while (Uright <= Umiddle)
        {
            alphaRight *= 2.0;
            _u = _uOld + alphaRight * delta;
            cubeMesh.ScatterDisplacements(_u, dynamicState);
            Uright = cubeMesh.ComputeEnergy(dynamicState, material);
        }
    }
    else if (Umiddle < Uright)
    {
        // Uright = Umiddle;
        // alphaRight = alphaMiddle;
        // Umiddle = Uleft;
        // alphaMiddle = alphaLeft;

        alphaLeft = -1.0;
        _u = _uOld + alphaLeft * delta;
        cubeMesh.ScatterDisplacements(_u, dynamicState);
        Uleft = cubeMesh.ComputeEnergy(dynamicState, material);

        // Search for an alpha with energy greater than the middle
        while (Uleft <= Umiddle)
        {
            alphaLeft *= 2.0;
            _u = _uOld + alphaLeft * delta;
            cubeMesh.ScatterDisplacements(_u, dynamicState);
            Uleft = cubeMesh.ComputeEnergy(dynamicState, material);
        }
    }
    //if (Umiddle == Uright)
    else
    {
        std::cerr << "Step yielded same energy... this corner case isn't coded up yet." << std::endl;
        return 0.0;
        // std::exit(EXIT_FAILURE);
    }

    // Search for a minimum in the bracket
    constexpr Scalar tol = 1.0e-4;

    while ((alphaRight - alphaLeft) > tol)
    {
        assert(alphaLeft <= alphaMiddle);
        assert(alphaRight >= alphaMiddle);
        assert(Uleft >= Umiddle);
        assert(Uright >= Umiddle);

        // Recurse in the right subinterval
        if ((alphaRight - alphaMiddle) > (alphaMiddle - alphaLeft))
        {
            const Scalar alpha_new = 0.5 * (alphaMiddle + alphaRight);
            // Scalar alpha_new = alphaLeft + (alphaRight - alphaMiddle);
            // if (alpha_new == alphaMiddle) { alpha_new = 0.5 * (alphaMiddle + alphaRight); }
            assert(alpha_new > alphaMiddle); assert(alpha_new < alphaRight);
            _u = _uOld + alpha_new * delta;
            cubeMesh.ScatterDisplacements(_u, dynamicState);
            const Scalar U_new = cubeMesh.ComputeEnergy(dynamicState, material);
            // If this is the new energy minimum, tighten the left bound
            if (U_new < Umiddle)
            {
                alphaLeft = alphaMiddle;
                Uleft = Umiddle;
                alphaMiddle = alpha_new;
                Umiddle = U_new;
                assert(Uleft >= Umiddle);
                assert(Uright >= Umiddle);
            }
            // Otherwise, tighten the right bound
            else
            {
                alphaRight = alpha_new;
                Uright = U_new;
                assert(Uleft >= Umiddle);
                assert(Uright >= Umiddle);
            }
        }
        // Recurse in the left subinterval
        else
        {
            const Scalar alpha_new = 0.5 * (alphaLeft + alphaMiddle);
            // Scalar alpha_new = alphaRight - (alphaMiddle - alphaLeft);
            // if (alpha_new == alphaMiddle) { alpha_new = 0.5 * (alphaLeft + alphaMiddle); }
            assert(alpha_new > alphaLeft); assert(alpha_new < alphaMiddle);
            _u = _uOld + alpha_new * delta;
            cubeMesh.ScatterDisplacements(_u, dynamicState);
            const Scalar U_new = cubeMesh.ComputeEnergy(dynamicState, material);
            // If this is a new energy minimum, tighten the right bound
            if (U_new < Umiddle)
            {
                alphaRight = alphaMiddle;
                Uright = Umiddle;
                alphaMiddle = alpha_new;
                Umiddle = U_new;
                assert(Uleft >= Umiddle);
                assert(Uright >= Umiddle);
            }
            // Otherwise, tighten the left bound
            else
            {
                alphaLeft = alpha_new;
                Uleft = U_new;
                assert(Uleft >= Umiddle);
                assert(Uright >= Umiddle);
            }
        }
    }

    _u = _uOld;
    cubeMesh.ScatterDisplacements(_u, dynamicState);

    return alphaMiddle;
}

}
