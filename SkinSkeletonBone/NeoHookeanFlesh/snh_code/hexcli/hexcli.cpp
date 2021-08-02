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

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "cubesim/CGInterface.h"
#include "cubesim/CubeMesh.h"
#include "cubesim/CubeNewtonSolver.h"
#include "cubesim/InitialMesh.h"
#include "cubesim/StableNeoHookean.h"

static std::vector<std::array<int,3>> g_surfaceTriangles;

static int g_nvx;
static int g_nvy;
static int g_nvz;

static CubeSim::CubeMesh g_cubeMesh;
static CubeSim::DynamicState g_dynamicState;
static std::unique_ptr<CubeSim::Material> g_material;
static CubeSim::CubeNewtonSolver g_newtonSolver;
static std::unique_ptr<CubeSim::CGInterface> g_cg;

static void QuasistaticSolve()
{
    constexpr bool printStatus = false;

    constexpr int maxIters = 50;
    constexpr CubeSim::Scalar tol = 1.0e-6;

    const CubeSim::QuasistaticSolveResult result = g_newtonSolver.SolveQuasistatic(printStatus, maxIters, tol, *g_material, g_cubeMesh, *g_cg, g_dynamicState);

    if (!result.success)
    {
        std::cerr << CubeSim::RED_C << "Warning, quasistatic solve failed." << CubeSim::NO_C << std::endl;
    }
}

static void UpdateBoundaryConditions(const int stepNum, const CubeSim::Scalar& stepDelta)
{
    const CubeSim::Scalar newNegativeBoundary = -1.0 - stepNum * stepDelta;
    for (int k = 0; k < g_nvz; k++)
    {
        constexpr int j = 0;
        for (int i = 0; i < g_nvx; i++)
        {
            const int vidx = i + j * g_nvy + k * g_nvx * g_nvy;
            g_dynamicState.vertices[static_cast<unsigned long>(vidx)].y() = newNegativeBoundary;
        }
    }

    const CubeSim::Scalar newPositiveBoundary = 1.0 + stepNum * stepDelta;
    for (int k = 0; k < g_nvz; k++)
    {
        const int j = g_nvy - 1;
        for (int i = 0; i < g_nvx; i++)
        {
            const int vidx = i + j * g_nvy + k * g_nvx * g_nvy;
            g_dynamicState.vertices[static_cast<unsigned long>(vidx)].y() = newPositiveBoundary;
        }
    }
}

static std::string GenerateOutputName(const std::string& dirName, const int stepNum)
{
    char buffer[1024];
    sprintf(buffer, "%03i.obj", stepNum);
    return dirName + "/" + std::string(buffer);
}

static void WriteMesh(std::ofstream& os)
{
    // Output the vertices
    for (const CubeSim::Vector3& p : g_dynamicState.vertices)
    {
        os << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    // Output the faces
    for (const std::array<int,3>& f : g_surfaceTriangles)
    {
        os << "f " << f[0] + 1 << " " << f[1] + 1 << " " << f[2] + 1 << std::endl;
    }
}

static int SaveMesh(const std::string& fileName)
{
    std::ofstream outputStream(fileName);
    if (!outputStream.is_open())
    {
        std::cerr << "Error, failed to open: " << fileName << std::endl;
        return EXIT_FAILURE;
    }
    WriteMesh(outputStream);
    return EXIT_SUCCESS;
}

static int RunSimulation(const std::string& outputDirectory)
{
    constexpr int numSteps = 25;
    constexpr CubeSim::Scalar stepSize = 0.1;

    if (!outputDirectory.empty())
    {
        const std::string fileName = GenerateOutputName(outputDirectory, 0);
        std::cout << "Saving mesh to: " << fileName << std::endl;
        if (SaveMesh(fileName) != EXIT_SUCCESS)
        {
            return EXIT_FAILURE;
        }
    }

    for (int stepNum = 1; stepNum <= numSteps; stepNum++)
    {
        std::cout << "Stepping to " << stepNum << " of " << numSteps << std::endl;
        UpdateBoundaryConditions(stepNum, stepSize);
        QuasistaticSolve();
        if (!outputDirectory.empty())
        {
            const std::string fileName = GenerateOutputName(outputDirectory, stepNum);
            std::cout << "Saving mesh to: " << fileName << std::endl;
            if (SaveMesh(fileName) != EXIT_SUCCESS)
            {
                return EXIT_FAILURE;
            }
        }
    }

    return EXIT_SUCCESS;
}

template<typename T>
bool ParseString(const std::string& str, T& val)
{
    std::stringstream sstrm(str);
    sstrm >> val;
    return !sstrm.fail();
}

static CubeSim::InitialMesh GenerateTestMesh(const int resolution, const int nv, const int nc)
{
    CubeSim::InitialMesh initialMesh;

    // Create the vertices
    const CubeSim::Scalar delta = 2.0 / CubeSim::Scalar(resolution);
    initialMesh.vertsInit.resize(3, nv);
    const Eigen::Vector3d origin(-1.0, -1.0, -1.0);
    for (int k = 0; k < g_nvz; k++)
    {
        for (int j = 0; j < g_nvy; j++)
        {
            for (int i = 0; i < g_nvx; i++)
            {
                const int vrtIdx = i + j * g_nvy + k * g_nvx * g_nvy;
                initialMesh.vertsInit.col(vrtIdx) = origin + delta * CubeSim::Vector3(i, j, k);
            }
        }
    }
    initialMesh.vertsRest = initialMesh.vertsInit;

    // Kinematically script two ends of the cube
    initialMesh.kinematic.resize(static_cast<unsigned long>(nv), false);
    for (int k = 0; k < g_nvz; k++)
    {
        constexpr int j = 0;
        for (int i = 0; i < g_nvx; i++)
        {
            const int vidx = i + j * g_nvy + k * g_nvx * g_nvy;
            initialMesh.kinematic[static_cast<unsigned long>(vidx)] = true;
        }
    }
    for (int k = 0; k < g_nvz; k++)
    {
        const int j = g_nvy - 1;
        for (int i = 0; i < g_nvx; i++)
        {
            const int vidx = i + j * g_nvy + k * g_nvx * g_nvy;
            initialMesh.kinematic[static_cast<unsigned long>(vidx)] = true;
        }
    }

    // Create the cubes
    initialMesh.cubes.resize(8, nc);
    for (int k = 0; k < g_nvz - 1; k++)
    {
        for (int j = 0; j < g_nvy - 1; j++)
        {
            for (int i = 0; i < g_nvx - 1; i++)
            {
                const int c0 = i + j * g_nvy + k * g_nvx * g_nvy;
                const int c1 = c0 + 1;
                const int c2 = c0 + g_nvx;
                const int c3 = c1 + g_nvx;
                const int c4 = c0 + g_nvx * g_nvy;
                const int c5 = c1 + g_nvx * g_nvy;
                const int c6 = c2 + g_nvx * g_nvy;
                const int c7 = c3 + g_nvx * g_nvy;
                const int cube_idx = i + j * (g_nvy - 1) + k * (g_nvx - 1) * (g_nvy - 1);
                initialMesh.cubes.col(cube_idx) << c0, c1, c2, c3, c4, c5, c6, c7;
            }
        }
    }

    return initialMesh;
}

static std::vector<std::array<int,3>> ExtractSurfaceTriangles(const int nvx, const int nvy, const int nvz)
{
    std::vector<std::array<int,3>> surfaceTriangles;

    // -y side
    for (int k = 0; k < nvz - 1; k++)
    {
        constexpr int j = 0;
        for (int i = 0; i < nvx - 1; i++)
        {
            const int c0 = i + j * nvy + k * nvx * nvy;
            const int c1 = c0 + 1;
            const int c4 = c0 + nvx * nvy;
            const int c5 = c1 + nvx * nvy;
            surfaceTriangles.emplace_back(std::array<int,3>{c0, c1, c4});
            surfaceTriangles.emplace_back(std::array<int,3>{c4, c1, c5});
        }
    }
    // +y side
    for (int k = 0; k < nvz - 1; k++)
    {
        const int j = nvy - 2;
        for (int i = 0; i < nvx - 1; i++)
        {
            const int c0 = i + j * nvy + k * nvx * nvy;
            const int c1 = c0 + 1;
            const int c2 = c0 + nvx;
            const int c3 = c1 + nvx;
            const int c6 = c2 + nvx * nvy;
            const int c7 = c3 + nvx * nvy;
            surfaceTriangles.emplace_back(std::array<int,3>{c2, c6, c3});
            surfaceTriangles.emplace_back(std::array<int,3>{c3, c6, c7});
        }
    }
    // -z side
    {
        constexpr int k = 0;
        for (int j = 0; j < nvy - 1; j++)
        {
            for (int i = 0; i < nvx - 1; i++)
            {
                const int c0 = i + j * nvy + k * nvx * nvy;
                const int c1 = c0 + 1;
                const int c2 = c0 + nvx;
                const int c3 = c1 + nvx;
                surfaceTriangles.emplace_back(std::array<int,3>{c0, c2, c1});
                surfaceTriangles.emplace_back(std::array<int,3>{c1, c2, c3});
            }
        }
    }
    // +z side
    {
        const int k = nvz - 2;
        for (int j = 0; j < nvy - 1; j++)
        {
            for (int i = 0; i < nvx - 1; i++)
            {
                const int c0 = i + j * nvy + k * nvx * nvy;
                const int c1 = c0 + 1;
                const int c2 = c0 + nvx;
                const int c3 = c1 + nvx;
                const int c4 = c0 + nvx * nvy;
                const int c5 = c1 + nvx * nvy;
                const int c6 = c2 + nvx * nvy;
                const int c7 = c3 + nvx * nvy;
                surfaceTriangles.emplace_back(std::array<int,3>{c4, c5, c6});
                surfaceTriangles.emplace_back(std::array<int,3>{c5, c7, c6});
            }
        }
    }
    // -x side
    for (int k = 0; k < nvz - 1; k++)
    {
        for (int j = 0; j < nvy - 1; j++)
        {
            constexpr int i = 0;
            const int c0 = i + j * nvy + k * nvx * nvy;
            // const int c1 = c0 + 1;
            const int c2 = c0 + nvx;
            // const int c3 = c1 + nvx;
            const int c4 = c0 + nvx * nvy;
            // const int c5 = c1 + nvx * nvy;
            const int c6 = c2 + nvx * nvy;
            // const int c7 = c3 + nvx * nvy;
            surfaceTriangles.emplace_back(std::array<int,3>{c0, c6, c2});
            surfaceTriangles.emplace_back(std::array<int,3>{c0, c4, c6});
        }
    }
    // +x side
    for (int k = 0; k < nvz - 1; k++)
    {
        for (int j = 0; j < nvy - 1; j++)
        {
            const int i = nvx - 2;
            {
                const int c0 = i + j * nvy + k * nvx * nvy;
                const int c1 = c0 + 1;
                const int c3 = c1 + nvx;
                const int c5 = c1 + nvx * nvy;
                const int c7 = c3 + nvx * nvy;
                surfaceTriangles.emplace_back(std::array<int,3>{c1, c3, c7});
                surfaceTriangles.emplace_back(std::array<int,3>{c1, c7, c5});
            }
        }
    }

    return surfaceTriangles;
}

int main(int argc, char** argv)
{
    if (argc != 5 && argc != 6)
    {
        std::cerr << "Usage: " << argv[0] << " resolution material_model mu lambda [output_directory]" << std::endl;
        return EXIT_FAILURE;
    }

    // Read the resolution
    int resolution;
    if (!ParseString(argv[1], resolution))
    {
        std::cerr << "Error, failed to read resolution parameter. Please provide an integer." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Resolution: " << resolution << std::endl;
    if (resolution <= 0)
    {
        std::cerr << "Error, resolution parameter must be positive." << std::endl;
        return EXIT_FAILURE;
    }

    // Read the material model
    const std::string materialModel = argv[2];
    std::cout << "Material model: " << materialModel << std::endl;
    if (materialModel != "stable_neo_hookean")
    {
        std::cerr << "Error, unsupported material model \"" << materialModel << "\" requested." << std::endl;
        std::cerr << "Valid options are: stable_neo_hookean" << std::endl;
        return EXIT_FAILURE;
    }

    // Read mu
    CubeSim::Scalar mu;
    if (!ParseString(argv[3], mu))
    {
        std::cerr << "Error, failed to read mu. Please provide a non-negative scalar." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "mu: " << mu << std::endl;
    if (mu < 0)
    {
        std::cerr << "Error, mu parameter must be non-negative." << std::endl;
        return EXIT_FAILURE;
    }

    // Read lambda
    CubeSim::Scalar lambda;
    if (!ParseString(argv[4], lambda))
    {
        std::cerr << "Error, failed to read lambda. Please provide a non-negative scalar." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "lambda: " << lambda << std::endl;
    if (lambda < 0)
    {
        std::cerr << "Error, lambda parameter must be non-negative." << std::endl;
        return EXIT_FAILURE;
    }

    // If present, read the output directory
    std::string outputDirectory;
    if (argc == 6)
    {
        outputDirectory = argv[5];
    }

    g_nvx = 1 + resolution;
    g_nvy = 1 + resolution;
    g_nvz = 1 + resolution;
    const int nv = g_nvx * g_nvy * g_nvz;
    std::cout << "Vertex count: " << nv << std::endl;
    const int nc = (g_nvx - 1) * (g_nvy - 1) * (g_nvz - 1);
    std::cout << "Cube count: " << nc << std::endl;

    if (materialModel == "stable_neo_hookean")
    {
        g_material.reset(new CubeSim::StableNeoHookean(mu, lambda));
    }
    else
    {
        std::cerr << "Invalid materila model. This is a bug." << std::endl;
        return EXIT_FAILURE;
    }

    const CubeSim::InitialMesh initialMesh = GenerateTestMesh(resolution, nv, nc);
    g_cubeMesh = CubeSim::CubeMesh(initialMesh);
    g_cubeMesh.CreateInitialDynamicState(initialMesh, g_dynamicState);

    // Generate a mesh to visualize the surface
    g_surfaceTriangles = ExtractSurfaceTriangles(g_nvx, g_nvy, g_nvz);

    g_newtonSolver = CubeSim::CubeNewtonSolver(g_cubeMesh.GetNumSimulatedDOFs());
    g_cubeMesh.CreateInitialSparseMatrix(g_newtonSolver.GetK());

    {
        constexpr CubeSim::Scalar cgTol = 1.0e-6;
        const int maxCGIters = 2 * g_cubeMesh.GetNumSimulatedDOFs();
        g_cg = CubeSim::GenerateCGSolver(CubeSim::CGSolverType::Eigen, cgTol, maxCGIters);
    }

    return RunSimulation(outputDirectory);
}
