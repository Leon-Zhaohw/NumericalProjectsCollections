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

#ifndef LINE_SEARCHES_H
#define LINE_SEARCHES_H

#include "Material.h"

namespace CubeSim
{

namespace LineSearches
{
    // Finds a local minimum in the search direction
    template<class T>
    void MinimizeInSearchDirection(const T& mesh, const Material& material, const VectorX& direction, VectorX& minimizedStep, Scalar& alphaMin, Scalar& Umin)
    {
        // Width of the bracket
        constexpr Scalar tol = 1.0e-4;

        minimizedStep.resize(direction.size());

        // Middle point in bracket
        Scalar alpha_middle = 0.0;
        Scalar U_middle = mesh.ComputePotentialEnergy(material);

        // Right bracket on minimum
        Scalar alpha_right = 1.0;
        mesh.UpdatePerturbedCachedState(direction);
        Scalar U_right = mesh.ComputePerturbedPotentialEnergy(material);

        // Left bracket on minimum
        Scalar alpha_left;
        Scalar U_left;

        if (U_middle > U_right)
        {
            U_left = U_middle;
            alpha_left = alpha_middle;
            U_middle = U_right;
            alpha_middle = alpha_right;

            // Search for an alpha with energy greater than the middle
            while (U_right <= U_middle)
            {
                alpha_right *= 2.0;
                minimizedStep = alpha_right * direction;
                mesh.UpdatePerturbedCachedState(minimizedStep);
                U_right = mesh.ComputePerturbedPotentialEnergy(material);
            }
        }
        else if (U_middle < U_right)
        {
            alpha_left = -1.0;
            minimizedStep = alpha_left * direction;
            mesh.UpdatePerturbedCachedState(minimizedStep);
            U_left = mesh.ComputePerturbedPotentialEnergy(material);

            // Search for an alpha with energy greater than the middle
            while (U_left <= U_middle)
            {
                alpha_left *= 2.0;
                minimizedStep = alpha_left * direction;
                mesh.UpdatePerturbedCachedState(minimizedStep);
                U_left = mesh.ComputePerturbedPotentialEnergy(material);
            }
        }
        else // (U_middle == U_right)
        {
            assert(U_middle == U_right);
            // Try the other direction
            alpha_left = -1.0;
            minimizedStep = alpha_left * direction;
            mesh.UpdatePerturbedCachedState(minimizedStep);
            U_left = mesh.ComputePerturbedPotentialEnergy(material);
            if (U_left != U_middle)
            {
                if (U_left > U_middle)
                {
                    std::cerr << "Right gives same energy...";
                    std::cerr << "and left step has greater energy... this corner case isn't coded up yet." << std::endl;
                    std::cerr << "U_left, U_center, U_right: " << U_left << "    " << U_middle << "    " << U_right << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                else // U_left < U_middle
                {
                    // Shift the values to the right
                    U_right = U_middle;
                    alpha_right = alpha_middle;
                    U_middle = U_left;
                    alpha_middle = alpha_left;

                    // Search for an alpha with energy greater than the middle
                    alpha_left = -2.0;
                    minimizedStep = alpha_left * direction;
                    mesh.UpdatePerturbedCachedState(minimizedStep);
                    U_left = mesh.ComputePerturbedPotentialEnergy(material);
                    while (U_left <= U_middle)
                    {
                        alpha_left *= 2.0;
                        minimizedStep = alpha_left * direction;
                        mesh.UpdatePerturbedCachedState(minimizedStep);
                        U_left = mesh.ComputePerturbedPotentialEnergy(material);
                    }
                }
            }
            // U_left == U_middle == U_right
            // else
            // {
            //     std::cout << "and left step gives same energy. Probably a flat energy landscape. Searching anyway." << std::endl;
            // }
        }

        while ((alpha_right - alpha_left) > tol)
        {
            assert(alpha_left <= alpha_middle);
            assert(alpha_right >= alpha_middle);
            assert(U_left >= U_middle);
            assert(U_right >= U_middle);

            // Recurse in the right subinterval
            if ((alpha_right - alpha_middle) > (alpha_middle - alpha_left))
            {
                const Scalar alpha_new = 0.5 * (alpha_middle + alpha_right);
                assert(alpha_new > alpha_middle); assert(alpha_new < alpha_right);
                minimizedStep = alpha_new * direction;
                mesh.UpdatePerturbedCachedState(minimizedStep);
                const Scalar U_new = mesh.ComputePerturbedPotentialEnergy(material);
                // If this is the new energy minimum, tighten the left bound
                if (U_new < U_middle)
                {
                    alpha_left = alpha_middle;
                    U_left = U_middle;
                    alpha_middle = alpha_new;
                    U_middle = U_new;
                    assert(U_left >= U_middle);
                    assert(U_right >= U_middle);
                }
                // Otherwise, tighten the right bound
                else
                {
                    alpha_right = alpha_new;
                    U_right = U_new;
                    assert(U_left >= U_middle);
                    assert(U_right >= U_middle);
                }
            }
            // Recurse in the left subinterval
            else
            {
                const Scalar alpha_new = 0.5 * (alpha_left + alpha_middle);
                assert(alpha_new > alpha_left); assert(alpha_new < alpha_middle);
                minimizedStep = alpha_new * direction;
                mesh.UpdatePerturbedCachedState(minimizedStep);
                const Scalar U_new = mesh.ComputePerturbedPotentialEnergy(material);
                // If this is a new energy minimum, tighten the right bound
                if (U_new < U_middle)
                {
                    alpha_right = alpha_middle;
                    U_right = U_middle;
                    alpha_middle = alpha_new;
                    U_middle = U_new;
                    assert(U_left >= U_middle);
                    assert(U_right >= U_middle);
                }
                // Otherwise, tighten the left bound
                else
                {
                    alpha_left = alpha_new;
                    U_left = U_new;
                    assert(U_left >= U_middle);
                    assert(U_right >= U_middle);
                }
            }
        }

        alphaMin = alpha_middle;
        Umin = U_middle;
        minimizedStep = alpha_middle * direction;
    }
}

}

#endif
