//
// Source code for the paper
// 
// D. Heck and T. Schloemer and O. Deussen, "Blue Noise Sampling with
// Controlled Aliasing", ACM Trans. Graph., 2013, in press
//
// Copyright (C) 2012,2013 Daniel Heck and Thomas Schloemer
//
#include "param.hh"
#include "common.hh"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cfloat>
#include <signal.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int nbins = -1;                 // default: number of points
bool randomize_force = false;
int maxattempts = 5;
float Tmax = 1, Tmin = 1e-3, Tfactor = 0.9;

// A flag set on SIGKILL to stop optimization
bool abort_flag = false;

float CalcGradients(const std::vector<Point> &pts, 
        const Curve &rdf, const Curve &target,
        std::vector<Point> &gradient) {
    Curve force(nbins, 0, 0.5f);

    for (int j=0; j<force.size(); j++)
        force[j] = force.dx * (rdf[j] - target[j]);
    for (int j=1; j<force.size(); j++)
        force[j] += force[j-1];
    for (int j=1; j<force.size(); j++) {
        float x = force.ToX(j);
        force[j] /= (pts.size() * x * x);
        // if (force[j] < 0) force[j] = 0;
    }
    force.Write("force.dat");

    float maxforce = 0.0f;
    gradient.resize(pts.size());
    
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned int i=0; i<pts.size(); i++) {
        Point grad;
        for (unsigned int j = 0; j < pts.size(); j++) {
            if (i == j)
                continue;

            float dx = pts[j].x - pts[i].x;
            float dy = pts[j].y - pts[i].y;
            dx += (dx < -0.5f) - (dx > 0.5f);
            dy += (dy < -0.5f) - (dy > 0.5f);
            float dist2 = dx * dx + dy * dy;

            // RDF is only reliable up to 0.5, if we go higher, the periodic
            // boundary conditions cause anisotropies in the final spectrum
            if (dist2 > 0.49f * 0.49f || dist2 == 0) continue;

            float dist = sqrtf(dist2);
            float f = force.At(dist);// / dist;
            grad.x -= dx * f;
            grad.y -= dy * f;
        }
        float ff = grad.x*grad.x + grad.y*grad.y;
        if (ff > maxforce)
            maxforce = ff;
        gradient[i] = grad;
    }
}
    return sqrtf(maxforce);
}

std::vector<Point> MovePoints(const std::vector<Point> &pts,
        const std::vector<Point> &gradient, float stepsize)
{
    std::vector<Point> newpts(pts.size());
    for (unsigned i=0; i<pts.size(); i++) {
        Point force = gradient[i];
        if (randomize_force) {
            force.x *= drand48();
            force.y *= drand48();
        }
        float x = pts[i].x + stepsize * force.x;
        float y = pts[i].y + stepsize * force.y;

        while (x >= 1) x -= 1;
        while (y >= 1) y -= 1;
        while (x < 0) x += 1;
        while (y < 0) y += 1;

        newpts[i].x = x;
        newpts[i].y = y;
    }
    return newpts;
}


float CalcEnergy(const Curve &a, const Curve &b) {
    Curve c(a);
    for (int i=0; i<c.size(); i++) 
        c[i] = powf(c[i] - b[i], 2.0f);
    return sqrtf(Integrate(c));
//    fprintf(stderr, "e = %g\n", x);
//    return x;
}

float MainOptimization(const std::vector<Point> &pts,
        Curve &target,
        float smoothing,
        std::vector<Point> &output)
{

    std::vector<Point> current = pts, best = pts;
    std::vector<Point> gradient;
    Curve rdf = CalcRDF(nbins, current.size(), &current[0].x, smoothing);
    float bestenergy = CalcEnergy(rdf, target);
    float T = Tmax;             // temperature
    int attempts = 0;
    while (!abort_flag && T >= Tmin) {
        // Calculate gradients and move points
        float maxgrad = CalcGradients(current, rdf, target, gradient);
        float stepsize = T / (sqrt(current.size()) * maxgrad);
        current = MovePoints(current, gradient, stepsize);

        // Perform iteration control
        rdf = CalcRDF(nbins, current.size(), &current[0].x, smoothing);
        float energy = CalcEnergy(rdf, target);
        attempts++;
        if (energy < bestenergy) {
            attempts = 0;
            best = current;
            bestenergy = energy;
        } else if (energy > bestenergy * 1.2) {
            attempts = maxattempts;
            current = best;
        }
        printf("%1.5g %1.5g\t%5.5g\t%g\n", 
                T, stepsize, energy, bestenergy);
        if (attempts >= maxattempts) {
            attempts = 0;
            T *= Tfactor;
        }
    }
    output = best;
    return bestenergy;
}


void FunctionStep(float critx, Curve &c, int npts) {
    float x0 = critx * sqrt(2/(sqrt(3)*npts));
    for (int i=0; i<c.size(); i++) {
        c[i] = (c.ToX(i) >= x0) ? 1 : 0;
    }
}

void FunctionJinc(float critx, Curve &c, int npts) {
    float x0 = critx / sqrt(2/(sqrt(3)*npts));
    for (int i=0; i<c.size(); i++) {
        float xx0 = c.ToX(i);
    
        float a = 0.0;
        const int N=20;
        for (int j=0; j<N; j++) {
            float x = xx0 + j*c.dx/N;
            if (x < 1e-5) a += 1-M_PI*x0*x0/npts;
            else {
                float xx = 2*M_PI*x0*x;
                a += 1 - 2*M_PI*x0*x0/npts * jn(1, xx) / xx;
            }
        }
        c[i] = a / N;
    }
}

void FunctionPeak(float critx, float peaky, Curve &c, int npts) {
    for (int i=0; i<c.size(); i++) {
        c[i] = (c.ToX(i) > critx) ? 1 : 0;
        c[i] = ((int) c.ToX(i) == (int) critx) ? peaky : c[i];
    }
}


void Usage() {
    std::cerr << "Usage: targetrdf [options]\n"
              << "  --help                          show this message\n"
              << "  --out file                      output file\n"
              << "  --dda file                      output RDF in DDA form, Zhou et al.-style"
              << "\nTarget RDF\n"
              << "  --reference pointset            determine RDF from point set\n"
              << "  --steprp                        generate step power\n"
              << "  --steprdf                       generate step RDF\n"
              << "  --peakrp                        generate peak power\n"
              << "  --crit freq\n"
              << "  --peak power                    height of peak\n"
              << "  --peaksmooth sigma (default 0)  Gaussian applied to peak power\n"
              << "  --rdf file                      load desired RDF from file\n"
              << "  --spectrum file                 load desired spectrum from file\n"
              << "  --rp-as-rdf file                load spectrum from file, use as RDF\n"
              << "\nOptimization Parameters\n"
              << "  --nbins n (default #points)     histogram size for RDF\n"
              << "  --smoothing sigma (default 8)   Gaussian applied to all RDFs\n"
              << "  --randomize                     randomly scale force\n"
              << "\nInput RDF\n"
              << "  --in file                       input point set\n"
              << "  --npts n (default 4096)         number of random points\n"
              << "  --seed num                      random seed\n"
        ;
}

static void sighandler(int sig) {
    fprintf(stderr, "Aborting...\n");
    abort_flag = true;
}

enum AnalyticCurve { 
    RDF_Step,
    RDF_Jinc,
    RP_Peak
};

int main(int argc, char **argv) {
    ParamList params;
    params.Define("in", "");
    params.Define("out", "targetrdf.txt");
    params.Define("dda", "targetdda.pfm");
    params.Define("reference", "");
    params.Define("rdf", "");
    params.Define("rp-as-rdf", "");
    params.Define("spectrum", "");
    params.Define("seed", "1");
    params.Define("nbins", "-1");
    params.Define("smoothing", "8");
    params.Define("randomize", "false");
    params.Define("npts", "4096");
    params.Define("steprp", "false");
    params.Define("steprdf", "false");
    params.Define("peakrp", "false");
    params.Define("crit", "");
    params.Define("peak", "");
    params.Define("peaksmooth", "0");
    params.Define("help", "false");

    std::vector<std::string> args;
    params.Parse(argc, argv, args);

    std::string infile = params.GetString("in");
    std::string outfile = params.GetString("out");
    std::string ddafile = params.GetString("dda");
    std::string reference = params.GetString("reference");
    std::string rdffile = params.GetString("rdf");
    std::string spectrum = params.GetString("spectrum");
    std::string rp_as_rdf = params.GetString("rp-as-rdf");
    nbins = params.GetInt("nbins");
    randomize_force = params.GetBool("randomize");
    int npts = params.GetInt("npts");
    long int seed = params.GetInt("seed");
    float smoothing = params.GetFloat("smoothing");

    AnalyticCurve analyticCurve = RDF_Jinc;
    float critfreq = 0.0f, peakpower = 1.0f, peaksmooth = 0.0f;

    if (params.GetBool("steprp")) {
        analyticCurve = RDF_Jinc;
        critfreq = params.GetFloat("crit", 0.606f);
    }
    if (params.GetBool("steprdf")) {
        analyticCurve = RDF_Step;
        critfreq = params.GetFloat("crit", 1.0f);
    }
    if (params.GetBool("peakrp")) {
        analyticCurve = RP_Peak;
        critfreq = params.GetFloat("crit", 36.5f);
        peakpower = params.GetFloat("peak", 1.0f);
        peaksmooth = params.GetFloat("peaksmooth", 0.0f);
    }
    bool show_usage = params.GetBool("help");

    if (!args.empty()) {
        std::cerr << "Unexpected argument '" << args[0] << "'\n";
        show_usage = true;
    }
    if (const Param *p = params.UnusedOption()) {
        std::cerr << "Unknown option '" << p->name << "'\n";
        show_usage = true;
    }

    if (!reference.empty() && !rdffile.empty()) {
        std::cerr << "--rdf and --reference are mutually exclusive.\n";
        show_usage = true;
    }
  
    if (show_usage) {
        Usage();
        exit(1);
    }

    std::vector<Point> pts, result;
    srand48(seed);
    if (infile.empty()) {
        fprintf(stderr, "No input file specified, using %d random points\n", npts);
        pts.resize(npts);
        for (int i=0; i<npts; i++) {
            pts[i].x = drand48();
            pts[i].y = drand48();
        }
    } else
        LoadPoints(infile, pts);
    result.resize(pts.size());

    if (nbins == -1)
        nbins = pts.size();

    printf("Using %d bins, smoothing=%g\n", nbins, smoothing);

    // Prepare target RDF
    Curve target(nbins, 0, 0.5f);
    if (!reference.empty()) {
        std::vector<Point> referencepts;
        LoadPoints(reference, referencepts);
        target = CalcRDF(nbins, referencepts.size(), &referencepts[0].x, 
                smoothing);
    } else if (!rdffile.empty()) {
        Curve tmp = Curve::Read(rdffile.c_str());
        target = FilterGauss(nbins, tmp, smoothing);
    } else if (!spectrum.empty()) {
        Curve tmp  = Curve::Read(spectrum.c_str());
        target = Power2RDF(npts, tmp, nbins, 0, 0.5f, smoothing);
    } else if (!rp_as_rdf.empty()) {
        Curve tmp = Curve::Read(rp_as_rdf.c_str());
        for (int i=0; i<nbins; i++) {
            int j = tmp.ToIndex(target.ToX(i) * npts);
            if (j >= 0 && j < tmp.size())
                target[i] = tmp[j];
            else
                target[i] = 1;
        }
    } else if (analyticCurve == RDF_Jinc) {
        FunctionJinc(critfreq, target, pts.size());
        Curve tmp(nbins, 0, npts / 2);
        tmp = RDF2Power(npts, target, nbins, 0, npts / 2);
        tmp.Write("target_rp.dat");
    } else if (analyticCurve == RDF_Step) {
        FunctionStep(critfreq, target, pts.size());
    } else if (analyticCurve == RP_Peak) {
        Curve tmp(nbins, 0, npts / 2);
        FunctionPeak(critfreq, peakpower, tmp, pts.size());
        if (peaksmooth > 0.f)
            tmp = FilterGauss(nbins, tmp, peaksmooth);
        tmp.Write("target_rp.dat");
        target = Power2RDF(npts, tmp, nbins, 0, 0.5f, smoothing);
    }
    target.Write("target_rdf.dat");
    printf("Target RDF written to 'target_rdf.dat'\n");
    
    // std::string fname("cccvt_dda.pfm");
    // int width, height;
    // float* data = NULL;
    // ReadPFM(width, height, data, fname.c_str());
    // printf("Read %s of size %dx%d\n", fname.c_str(), width, height);
    // for (int x = 0; x < width; ++x)
    //     for (int y = 0; y < height; ++y)
    //         printf("%d, %d: %f\n", x, y, data[x + y*width]);
    // if (data) delete[] data;
    
    if (!ddafile.empty()) {
        // 128, 0.125
        // 256, 0.175
        // 512, 0.225
        const int size = 256;
        const float range = 0.175f;  // Guessed empirically
        std::vector<float> dda(size * size);
        for (int x = 0; x < size; ++x) {
            for (int y = 0; y < size; ++y) {
                float dx = (x - size/2) * range / (size/2);
                float dy = (y - size/2) * range / (size/2);
                float dist = sqrtf(dx*dx + dy*dy);
                dda[x + y*size] = target.At(dist) / (size/2); // Scaling by size/2 guessed
            }
        }
        fprintf(stderr, "Writing '%s'\n", ddafile.c_str());
        WritePFM(size, size, &dda[0], ddafile.c_str());
    }

    fprintf(stderr, "Starting optimization\n");
    signal(SIGINT, sighandler);

    result = pts;
    timeval t0, t;
    
    gettimeofday(&t0, NULL);
    float bestd = MainOptimization(pts, target, smoothing, result);
    gettimeofday(&t, NULL);
    
    long long delta = (long long)(t.tv_sec  - t0.tv_sec ) * 1000ll +
                      (long long)(t.tv_usec - t0.tv_usec) / 1000ll;
    fprintf(stderr, "- finished: %5.5g, time: %3.fs -\n", bestd,
            delta / 1000.0);
    fprintf(stderr, "Writing '%s'\n", outfile.c_str());
    WritePoints(outfile, result);
}
