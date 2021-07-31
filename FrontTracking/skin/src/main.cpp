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


#include <getopt.h>
#include "particleIO.H"

#ifdef OPENVDB
#include <openvdb/openvdb.h>
#include <tbb/task_scheduler_init.h>
#include "smoothingFilter.H"

int main(int argc, char **argv) {
    openvdb::initialize();
    static int verboseFlag = 0, explicitFlag = 0, vrFlag = 0, helpFlag = 0;
    int lnTimesteps = 15, bnTimesteps = -1, redistanceFrequency = -1, cgMaxIters=1500, threads=-1;
    double rmin = -DBL_MAX, rmax = -DBL_MAX, rinit = -DBL_MAX,
			  ldt = -DBL_MAX, bdt = -DBL_MAX, cgThreshold = 0.0001;
    double rratio = 2;
    static struct option long_options[] = {
        {"help", no_argument, &helpFlag, 1},
        {"verbose", no_argument, &verboseFlag, 1},
        {"explicit", no_argument, &explicitFlag, 1},
        {"threads", required_argument, 0, 't'},
        {"variable_radius", no_argument, &vrFlag, 1},
        {"vr", no_argument, &vrFlag, 1},
        {"rratio", required_argument, 0, 'r'},
        {"rmin", required_argument, 0, 'm'},
        {"rmax", required_argument, 0, 'M'},
        {"rinit", required_argument, 0, 'i'},
        {"redist_freq", required_argument, 0, 1},
        {"ldt", required_argument, 0, 2},
				{"ltimesteps", required_argument, 0, 3},
				{"bdt", required_argument, 0, 4},
				{"btimesteps", required_argument, 0, 5},
        {"cg_threshold", required_argument, 0, 6},
        {"cg_max_iter", required_argument, 0, 7},
        {0, 0, 0, 0}
    };
    while(1)
    {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hvet:r:m:M:i:",
                             long_options, &option_index);
        if (c == -1) break;
        switch (c) {
        case 0:
            if (long_options[option_index].flag != 0)
                break;
            break;

        case 'h':
            helpFlag = 1;
            break;

        case 'v':
            verboseFlag = 1;
            break;

        case 'e':
            explicitFlag = 1;
            break;

        case 't':
            threads = atoi(optarg);
            break;

        case 'r':
            rratio = atof(optarg);
            break;

        case 'm':
            rmin = atof(optarg);
            break;

        case 'M':
            rmax = atof(optarg);
            break;

        case 'i':
            rinit = atof(optarg);
            break;

        case 1:
            redistanceFrequency = atoi(optarg);
            break;

        case 2:
            ldt = atof(optarg);
            break;

        case 3:
            lnTimesteps = atoi(optarg);
            break;

        case 4:
            bdt = atof(optarg);
            break;

        case 5:
            bnTimesteps = atoi(optarg);
            break;

        case 6:
            cgThreshold = atof(optarg);
            break;

        case 7:
            cgMaxIters = atoi(optarg);
            break;

        case '?':
            break;

        default:
            abort ();
        }
    }

 
    if (helpFlag || argc == 1)
    {
        std::cout<<"Welcome! "<<std::endl;
        std::cout<<"Usage: "<<argv[0]<<" grid_spacing input_file output_file"<<std::endl;
        std::cout<<"Available switches are......."<<std::endl;
        std::cout<<"-h/--help -> display this message "<<std::endl;
        std::cout<<"-v/--verbose -> output information to terminal "<<std::endl;
				std::cout<<"-t/--threads number -> Number of threads (default is automatic)"<<std::endl;
 		    std::cout<<"-e/--explicit -> Use explicit integration "<<std::endl;
        std::cout<<"--vr/--variable_radius -> The particles have variable radiuses,\n\tin this case r_min, r_max, and r_init are multipliers for the\n\tindividual particle radii, you will want to use -m and -M with this option"<<std::endl;
        std::cout<<"-r/--rratio number -> r_max =  number * r_min, Default value is 2\n\t(over-riden by -M)"<<std::endl;
        std::cout<<"-m/--rmin number -> Minimum radius,\n\tDefault value is (sqrt(3) * grid_spacing)"<<std::endl;
        std::cout<<"-M/--rmax number -> Maximum radius, Default value is (2 * rmin)"<<std::endl;
        std::cout<<"-i/--rinit number -> r_init = number * r_min, otherwise\n\tr_init = 0.5 * (r_min + r_max) "<<std::endl;
        std::cout<<"--redist_freq number -> Frequency of redistancing, Default value is 1 (implicit) or 50 (explicit)"<<std::endl;
        std::cout<<"--ldt number -> Timestep for Laplacian smoothing\n\t Deafult value is 0.1*h^2"<<std::endl;
				std::cout<<"--ltimesteps number -> Number of Laplacian smoothing timesteps,\n\tDefault value is 15"<<std::endl;
        std::cout<<"--bdt number -> Timestep for biharmonic smoothing,\n\tDefault value is 20*h^4 (implicit) or 0.01*h^4 (explicit)"<<std::endl;
        std::cout<<"--btimestesps number -> Number of biharmonic smoothing timesteps,\n\tDefault value is 5 (implicit) or 500 (explicit)"<<std::endl;
		    std::cout<<"--cg_threshold number -> Threshold for CG convergence, Default value is .0001"<<std::endl;
				std::cout<<"--cg_max_iter number -> Upper limit for CG iterations, Default value is 1500"<<std::endl;
        exit(1);
    }

    tbb::task_scheduler_init *init;
    if (threads < 0) init = new tbb::task_scheduler_init(tbb::task_scheduler_init::automatic);
    else init = new tbb::task_scheduler_init(threads);

    double h = atof(argv[optind++]);
    char *infname =  argv[optind++];
    char *outfname = argv[optind++];

    if (rmin == -DBL_MAX) rmin = 2*0.86603 * h;
    if (rmax == -DBL_MAX) rmax = rratio*rmin;
    if (rinit == -DBL_MAX) rinit = 0.5*(rmin+rmax);
    if (ldt == -DBL_MAX) ldt = 0.1*h*h;
    if (explicitFlag) {
			  if (redistanceFrequency == -1) redistanceFrequency = 50;
			  if (bnTimesteps == -1) bnTimesteps = 500;
			  if (bdt == -DBL_MAX) bdt = 0.01*h*h*h*h;
		} else {
			  if (redistanceFrequency == -1) redistanceFrequency = 1;
			  if (bnTimesteps == -1) bnTimesteps = 5;
			  if (bdt == -DBL_MAX) bdt = 20*h*h*h*h;
		}

    const double voxelSize = h;
    const int halfWidth = 5 + ceil(std::max((rmax-rinit)/h, (rinit-rmin)/h));
    openvdb::DoubleGrid::Ptr grid, maxgrid, mingrid;
    grid = openvdb::createLevelSet<openvdb::DoubleGrid>(voxelSize, halfWidth);

    if (vrFlag) {
        maxgrid = openvdb::createLevelSet<openvdb::DoubleGrid>(voxelSize, halfWidth);
        mingrid = openvdb::createLevelSet<openvdb::DoubleGrid>(voxelSize, halfWidth);
				rasterizeVariableRadius(grid, maxgrid, mingrid, rinit, rmin, rmax, infname);
    } else {
			  rasterize(grid, rinit, infname);
		}

    SmoothingFilter *filter;
    if (!explicitFlag) filter = new SmoothingFilter(grid, 9);
    else filter = new SmoothingFilter(grid, 5);

    if (vrFlag) filter->initializeLevelSetVariableRadius(*mingrid, *maxgrid);
    else filter->initializeLevelSet((rmin-rinit), (rmax-rinit));

    filter->laplacianSmoothing(ldt, lnTimesteps);

    if (!explicitFlag) {
        filter->biharmonicSmoothing(bdt, bnTimesteps, redistanceFrequency, cgThreshold, cgMaxIters, verboseFlag);
    } else {
        filter->biharmonicSmoothing(bdt, bnTimesteps, redistanceFrequency);
		}

    filter->extractSurface(outfname);

    return 0;
}

#else //!OPENVDB

#include "smoothingGrid.H"
#include "marchingCube.H"
#include <float.h>


void writeObj (char *outfname, const SmoothingGrid &grid) {
    std::vector<SlTri> triangles; // triangles in polygon mesh
    std::vector<SlVector3> meshPts; // points in polygon mesh
    std::vector<SlVector3> normals; // normals (at points) in polygon mesh

    MarchingCube mc(grid.nx, grid.ny, grid.nz, grid.h, grid.bbMin);
    mc.buildTriangleMesh(grid.phi, triangles, meshPts);

    std::ofstream out(outfname, std::ios::out);
    std::vector<SlVector3>::const_iterator p;
    std::vector<SlTri>::const_iterator t;
    std::vector<SlVector3>::const_iterator n;
    for (p = meshPts.begin(); p != meshPts.end(); p++)
        out<<"v "<<(*p)[0]<<" "<<(*p)[1]<<" "<<(*p)[2]<<std::endl;
    for (t = triangles.begin(); t != triangles.end(); t++)
        out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;
}


int main(int argc, char **argv) {
    unsigned int flags = 0;
    static int verboseFlag = 0, explicitFlag = 0, helpFlag = 0, vrFlag = 0, naFlag = 0, vaFlag = 0;
    int lnTimesteps = -1, bnTimesteps = -1, redistanceFrequency = -1, cgMaxIters = 1500;
    double rmin = -DBL_MAX, rmax = -DBL_MAX, rinit = -DBL_MAX, velGain = 1.0,
  			ldt = -DBL_MAX, bdt = -DBL_MAX, cgThreshold = 0.0001;
    double rratio = 2, maxStretch = 4;

    static struct option long_options[] = {
        {"help", no_argument, &helpFlag, 1},
        {"verbose", no_argument, &verboseFlag, 1},
        {"explicit", no_argument, &explicitFlag, 1},
        {"variable_radius", no_argument, &vrFlag, 1},
        {"vr", no_argument, &vrFlag, 1},
        {"rratio", required_argument, 0, 'r'},
        {"rmin", required_argument, 0, 'm'},
        {"rmax", required_argument, 0, 'M'},
        {"rinit", required_argument, 0, 'i'},
        {"redist_freq", required_argument, 0, 1},
        {"ldt", required_argument, 0, 2},
				{"ltimesteps", required_argument, 0, 3},
				{"bdt", required_argument, 0, 4},
				{"btimesteps", required_argument, 0, 5},
        {"cg_threshold", required_argument, 0, 6},
        {"cg_max_iter", required_argument, 0, 7},
				{"neighbor_anisotropy", no_argument, &naFlag, 1},
				{"na", no_argument, &naFlag, 1},
				{"velocity_anisotropy", required_argument, 0, 8},
				{"va", required_argument, 0, 8},
				{"max_stretch", required_argument, 0, 9},
        {0, 0, 0, 0}
    };
    while(1)
    {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hvet:r:m:M:i:",
                             long_options, &option_index);
        if (c == -1) break;
        switch (c) {

        case 0:
            if (long_options[option_index].flag != 0)
                break;
            break;

        case 'h':
            helpFlag = 1;
            break;

        case 'v':
            verboseFlag = 1;
            break;

        case 'e':
            explicitFlag = 1;
            break;

        case 'r':
            rratio = atof(optarg);
            break;

        case 'm':
            rmin = atof(optarg);
            break;

        case 'M':
            rmax = atof(optarg);
            break;

        case 'i':
            rinit = atof(optarg);
            break;

        case 1:
            redistanceFrequency = atoi(optarg);
            break;

        case 2:
            ldt = atof(optarg);
            break;

        case 3:
            lnTimesteps = atoi(optarg);
            break;

        case 4:
            bdt = atof(optarg);
            break;

        case 5:
            bnTimesteps = atoi(optarg);
            break;

        case 6:
            cgThreshold = atof(optarg);
            break;

        case 7:
            cgMaxIters = atoi(optarg);
            break;

        case 8:
					  vaFlag = 1;
						velGain = atof(optarg);
            break;

				case 9:
					  maxStretch = atof(optarg);
						break;

        default:
            abort ();
        }
    }


    if (helpFlag || argc == 1) {
        std::cout<<"Welcome! "<<std::endl;
        std::cout<<"Usage: "<<argv[0]<<" grid_spacing input_file output_file"<<std::endl;
        std::cout<<"Available switches are......."<<std::endl;
        std::cout<<"-h/--help -> display this message "<<std::endl;
        std::cout<<"-v/--verbose -> output information to terminal "<<std::endl;
 		    std::cout<<"-e/--explicit -> Use explicit integration "<<std::endl;
        std::cout<<"--vr/--variable_radius -> The particles have variable radiuses,\n\tin this case r_min, r_max, and r_init are multipliers for the\n\tindividual particle radii, you will want to use -m and -M with this option"<<std::endl;
        std::cout<<"-r/--rratio number -> r_max =  number * r_min, Default value is 2\n\t(over-riden by -M)"<<std::endl;
        std::cout<<"-m/--rmin number -> Minimum radius,\n\tDefault value is (sqrt(3) * grid_spacing)"<<std::endl;
        std::cout<<"-M/--rmax number -> Maximum radius, Default value is (2 * rmin)"<<std::endl;
        std::cout<<"-i/--rinit number -> r_init = number * r_min, otherwise\n\tr_init = 0.5 * (r_min + r_max) "<<std::endl;
        std::cout<<"--redist_freq number -> Frequency of redistancing, Default value is 1 (implicit) or 50 (explicit)"<<std::endl;
        std::cout<<"--ldt number -> Timestep for Laplacian smoothing\n\t Deafult value is 0.1*h^2"<<std::endl;
				std::cout<<"--ltimesteps number -> Number of Laplacian smoothing timesteps,\n\tDefault value is 15"<<std::endl;
        std::cout<<"--bdt number -> Timestep for biharmonic smoothing,\n\tDefault value is 20*h^4 (implicit) or 0.01*h^4 (explicit)"<<std::endl;
        std::cout<<"--btimestesps number -> Number of biharmonic smoothing timesteps,\n\tDefault value is 5 (implicit) or 500 (explicit)"<<std::endl;
		    std::cout<<"--cg_threshold number -> Threshold for CG convergence, Default value is .0001"<<std::endl;
				std::cout<<"--cg_max_iter number -> Upper limit for CG iterations, Default value is 1500"<<std::endl;
        std::cout<<"--na/--neighbor_anisotropy -> Turn on neighborhood based anisotropy "<<std::endl;
        std::cout<<"--va/--velocity_anisotropy number -> Turn on velocity-based anisotropy,\n\tthe number is a gain on the amount of anisotropy,\n\tlarger values lead to more anisotropy "<<std::endl;
        std::cout<<"--max_stretch number -> maximum amount of anisotropy (condition number of G),\n\tDefault value is 4 "<<std::endl;
        exit(1);
    }

    double h = atof(argv[optind++]);
    char *infname =  argv[optind++];
    char *outfname = argv[optind++];

    if (rmin == -DBL_MAX) {
			rmin = 1.73295 * h;  // sqrt(3)*h
			// there could be a stretch up to this amount, we want to make sure that particles
			// all touch at least one grid point, so make the radius bigger...
			if (naFlag) {
        rmin *= sqrt(maxStretch);
			} else if (vaFlag) {
        rmin *= cbrt(maxStretch);
			}
		}

    if (rmax == -DBL_MAX) rmax = rratio*rmin;
    if (rinit == -DBL_MAX) rinit = 0.5*(rmin+rmax);
    if (ldt == -DBL_MAX) ldt = 0.1*h*h;
    if (explicitFlag) {
			  if (redistanceFrequency == -1) redistanceFrequency = 50;
			  if (bnTimesteps == -1) bnTimesteps = 500;
			  if (bdt == -DBL_MAX) bdt = 0.01*h*h*h*h;
		} else {
			  if (redistanceFrequency == -1) redistanceFrequency = 1;
			  if (bnTimesteps == -1) bnTimesteps = 5;
			  if (bdt == -DBL_MAX) bdt = 20*h*h*h*h;
		}
		
		if (verboseFlag) flags |= SmoothingGrid::VERBOSE;
		if (naFlag) flags |= SmoothingGrid::NEIGHBOR_ANISOTROPY;
		if (vaFlag) flags |= SmoothingGrid::VELOCITY_ANISOTROPY;

    std::vector<SlVector3> particles, velocities;
    std::vector<double> radii;
    readfile(infname, particles, radii, velocities);

    SmoothingGrid grid(h, rmin, rmax, rinit, velGain, maxStretch, flags, particles, radii, velocities);

		grid.doLaplacianSmoothing(ldt, lnTimesteps, redistanceFrequency);

    if (!explicitFlag) {
        grid.doBiharmonicSmoothing(bdt, bnTimesteps, redistanceFrequency, cgThreshold, cgMaxIters);
    } else {
        grid.doBiharmonicSmoothing(bdt, bnTimesteps, redistanceFrequency);
    }

    writeObj(outfname, grid);

    return 0;
}

#endif // OPENVDB
