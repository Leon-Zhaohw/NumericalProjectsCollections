/*
 * Post optimize the point set over different ranks
 * 2017-01-23: First implemented by Abdalla Ahmed
 *
 */

#include "geometry.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <complex>
#include <fstream>
#include <string>

typedef std::string String;

struct Tile {
    std::vector<Points> p;                                                  // Point coordinates at different ranks
    std::vector<Lists> list;                                                // Lists of belonging sites of different points in different ranks
    TList children;                                                         // ids of sub-tiles
    TList order;
    void read(std::fstream &file, int nChildren) {
        int nRanks = nChildren - 1;
        p.resize(nRanks);
        list.resize(nRanks);
        for (int pointNo = 0; pointNo < nRanks; pointNo++) {
            p[pointNo].resize(nRanks - pointNo + 1);                        // The last point is reserved for the point of the corresponding child point
            list[pointNo].resize(nRanks - pointNo);                                  // For each point and every rank we maintain a list.
        }
        double x, y;
        file >> x >> y;
        p[0][0] = Point(x, y);                                              // The optimal position of first point with single point per tile.
        children.resize(nChildren);
        for (int i = 0; i < nChildren; i++) {
            file >> children[i];
        }
        order.resize(nChildren);
        for (int i = 0; i < nChildren; i++) {
            file >> order[i];
        }
    };
};

struct Site {                                                                   // This is to handle the Delaunay triangulation.
    Point p;
    VH9 vh;                                                                     // Handles all vertices in triangulation
    unsigned X, Y;                                                              // Stratum coordinates
    unsigned tileID;
    unsigned pointNo;                                                           // Order of this point in the tile
    unsigned rank;                                                              // Rank at which this site appears
    bool isStable, isUnstable;                                                  // Used during optimization to skip processing a site if neither it nor neighbors move.
    TList *list;                                                                // List containing this site; that is, list of same tileID, pointNo, and rank.
    inline Vector getDisplacement() {
        return p - Point(X, Y);
    };
};

class TileSet {
private:
    TList idList;
    int one;                                                                    // Width of domain
    int idCount;                                                                // Distinct ID's in 1D
    int tileCount;                                                              // number of IDs in 2D; = idCount * idCount
    unsigned scale;                                                             // Self-similarity scale in 1D
    int nChildren;                                                              // Number of children tiles; = scale * scale
    int nRanks;                                                                 // Number of ranks (densities) within the octave; = nChildren - 1
    std::vector<Matrix> matrix;                                                 // Transformation matrix to children tiles coordinates
    std::vector<double> dhex, avg_area, target_rf, target_rc;                   // Spatial measures over different ranks
    double max_area_deviation;                                                  // Target maximum deviation of cell areas (capacity)
    bool isAdaptive = false;                                                    // Whether the target use is adaptive sampling; then we maintain coherence across ranks and octaves,
    std::vector<double> weights;                                                //   by using a filter with these weights to smooth the point locations over ranks and octaves
    bool allStable;                                                             // To implement termination criteria of optimization
    int stableCount;                                                            // To track optimization
    double maxShift;
    double maxStep;
    std::vector<Tile> tiles;                                                    // The list of tiles in the set
    std::vector<TList> g;                                                       // Lists of sites carrying same ID and rank.
    std::vector<DT> dts;                                                        // A list of Delaunay triangulations, one for each rank.
    std::vector<Site> s;                                                        // A list of sites (points) in all triangulations
    int total;                                                                  // Total number of sites
    std::vector<unsigned> N, beginID, endID;                                    // Number of sites, begin, and end site index in each triangulation

    inline Point childPoint(unsigned tileID, unsigned pointNo);                 // Location, in parent's coordinates, of the point of the pointNo'th child in the list of children
    inline double toroidallinearDist(double x1, double x2) const {              // 1D Nearest distance between replicas of two points
        double dx = x1 - x2;                                                    // Find distance in primary period
        if (dx > one/2) dx -= one;                                              // If larger than half the period length another replica is certainly nearer
        if (dx < -one/2) dx += one;                                             // Same, but opposite ordering of points
        return dx;
    };
    inline double toroidalSqDist(const Point &p1, const Point &p2) {            // 2D Nearest distance; square to avoid costly square root operation
        double dx = toroidallinearDist(p1.x(), p2.x());
        double dy = toroidallinearDist(p1.y(), p2.y());
        return dx * dx + dy * dy;
    };
    inline Point replica(Point p, int i) {                                      // Find one of the 9 replicas of a point
        i = (i+4) % 9;                                                          // We make the middle replica at index 0
        double x = p.x() + (i%3 - 1) * one;                                     // Add -one, 0, or one to x
        double y = p.y() + (i/3 - 1) * one;                                     // Same for y
        return Point(x, y);
    };
    Point marginBL, marginTR;                                                   // Points within these margins affect points in the pattern
    void setSite(unsigned siteID, int tileID, int pointNo, int rank, int X, int Y);
    void moveSite(int siteID, Point p);
    inline void moveSite(int siteID, Vector& shift);
    inline void markUnstable(Site &site);
    void moveAll(unsigned tileID, unsigned pointNo, unsigned rank, Vector shift, bool clamp = true);
    inline void moveAll(unsigned tileID, unsigned pointNo, unsigned rank, Point target, bool clamp = true);
    Vector getDisplacement(unsigned tileID, unsigned pointNo, unsigned rank);
    void updateStability();
    void capacity(unsigned siteID);
    void conflict(unsigned siteID);
    void coverage(unsigned siteID);
    void farthestPoint(unsigned siteID);
    void lloyd(unsigned siteID);
    void lloyd();
    void conflict();
    void coverage();
    String outputPath = "";                                                     // Path for output files; could be set to a folder in /tmp
    void setupPointSets();
    inline void smooth(Points &p);
public:
    TileSet(std::fstream &file);
    void setOutputPath(String v) { outputPath = v; }
    void setdmin(double v);
    void setrc(double v);
    void setMaxDrift(double v);
    void set_max_area_deviation(double v);
    void set_maxStep(double v);
    void setSmoothingKernel(double standardDeviation);
    void printText(String baseName);
    void serial(String seq);
    void parallel(String seq);
    void latinize();
    void setAllUnstable();
    int getStableCount() { return stableCount; };
    bool isAllStable() { return allStable; };
    void printTable(String fileName);
};

TileSet::TileSet(std::fstream &file) {
    file >> scale >> idCount;
    file >> one;
    idList.resize(one);
    for (int i = 0; i < one; i++) file >> idList[i];
    nChildren = scale * scale;
    nRanks = nChildren - 1;
    tileCount = idCount * idCount;
    tiles.resize(tileCount);
    for (int id = 0; id < tiles.size(); id++) {
        tiles[id].read(file, nChildren);
        if (file.eof()) {
            fprintf(stderr, "Failed to load tiles\n");
            exit(1);
        }
    }

    /************* Populate the children transformation matrix: ***************/
    double step = 1.0 / scale;
    matrix.resize(nChildren);
    for (int Y = 0; Y < scale; Y++) {
        for (int X = 0; X < scale; X++) {
            matrix[Y * scale + X] = {step, step * X, step * Y};
        }
    }
    fprintf(stderr, "Loaded %4d tiles\n", tiles.size());

    /********************** Initialize point locations ************************/
    for (int tileID = 0; tileID < tileCount; tileID++) {
        for (int pointNo = 0; pointNo < nRanks; pointNo++) {
            for (int rank = pointNo; rank < nRanks; rank++) {
                if (pointNo + rank == 0) continue;                              // Skip first point in first rank, it is already set.
                Points &p = tiles[tileID].p[pointNo];                           // Select the list of points
                p[rank - pointNo] = childPoint(tileID, pointNo);
            }
        }
    }
    setupPointSets();
}

Point TileSet::childPoint(unsigned tileID, unsigned pointNo) {                  // Coordinates of rank-0 point of the pointNo'th child tile, relative to parent's origin
    unsigned childNo = tiles[tileID].order[pointNo];
    unsigned childID = tiles[tileID].children[childNo];
    return transform(tiles[childID].p[0][0], matrix[childNo]);
}

void TileSet::setupPointSets() {
    maxStep = one;                                                              // unbound maximum step.
    double margin = 4;                                                          // Margin for toroidal domain.
    marginBL = {-margin, -margin};                                              // Bottom-left of primary period + margin
    marginTR = {(one * scale) + margin, (one * scale) + margin};                // Top-right. In our convention BL is included, TR is excluded
    total = one * one * (nRanks * (1 + nRanks) / 2);                            // Using sum of arithmetic progression
    fprintf(stderr, "Total %d sites\n", total);
    s.resize(total);
    N.resize(nRanks);
    beginID.resize(nRanks);
    endID.resize(nRanks);
    dts.resize(nRanks);
    dhex.resize(nRanks);
    target_rf.resize(nRanks);
    target_rc.resize(nRanks);
    avg_area.resize(nRanks);
    int siteID = 0;
    for (int rank = 0; rank < nRanks; rank++) {
        int pointCount = rank + 1;
        N[rank] = one * one * pointCount;
        beginID[rank] = siteID;
        avg_area[rank] = 1.0 / pointCount;                                      // Tile area is 1, hence average cell area is inversely proportional to number of points in the tile
        dhex[rank] = sqrt(avg_area[rank]) * sqrt(2/sqrt(3));                    // dhex scales in accordance with the square root of average area
        for (int Y = 0; Y < one; Y++) {
            for (int X = 0; X < one; X++) {
                unsigned tileID = idList[Y] * idCount + idList[X];
                for (int pointNo = 0; pointNo < pointCount; pointNo++) {
                    unsigned childNo = tiles[tileID].order[pointNo];
                    unsigned childID = tiles[tileID].children[childNo];
                    setSite(siteID++, tileID, pointNo, rank, X, Y);
                }
            }
        }
        endID[rank] = siteID;
    }
    setdmin(0.87);                                                              // Default settings
    setrc(0.67);
    set_max_area_deviation(0.03);
}

void TileSet::setSite(unsigned siteID, int tileID, int pointNo, int rank, int X, int Y) {    // Set location of the indexed point in pattern and insert it in triangulation
    Site &site = s[siteID];
    site.tileID = tileID;
    site.pointNo = pointNo;
    site.rank = rank;
    site.X = X;
    site.Y = Y;
    site.list = &(tiles[tileID].list[pointNo][rank - pointNo]);
    site.list->push_back(siteID);
    Point p = tiles[tileID].p[pointNo][rank-pointNo] + Vector(X, Y);
    site.p = p;                                                             // Save a handy copy of point coordinates
    DT &dt = dts[rank];
    for (int i = 0; i < 9; i++) {                                           // We loop through the 9 replicas,
        if (isInRect(replica(p, i), marginBL, marginTR)) {                  // if the original (not the shifted) location of a replica is within margin
            site.vh[i] = dt.insert(replica(p, i));                      // insert replica in triangulation and keep handle to it
            site.vh[i]->info() = siteID;
        }
        else site.vh[i] = NULL;
    }
    int childNo = tiles[tileID].order[pointNo];
};

void TileSet::moveSite(int siteID, Point p) {
    Site &site = s[siteID];
    markUnstable(site);                                                         // The site should be marked unstable, whether it's moved or not.
    if (
        (int(floor(p.x())) != site.X) ||
        (int(floor(p.y())) != site.Y)
    ) return;                                                                   // Only move point if target is inside the stratum.
    double l = (p - site.p).squared_length();
    maxShift = std::max(maxShift, l);
    site.p = p;                                                                 // Save a handy copy of updated point coordinates
    for (int i = 0; i < 9; i++) {
        if (site.vh[i] != NULL) {
            site.vh[i] = dts[site.rank].move(site.vh[i], replica(p, i));
        }
    }
    markUnstable(site);                                                         // To mark the new neighbors unstable.
}

void TileSet::moveSite(int siteID, Vector& shift) {
    moveSite(siteID, s[siteID].p + shift);
}

void TileSet::markUnstable(Site &site) {
    site.isUnstable = true;                                                     // Mark the site unstable
    VC vc = dts[site.rank].incident_vertices(site.vh[0]), done(vc);
    do {                                                                        // Mark neighbors unstable
        s[ vc->info() ].isUnstable = true;
    } while (++vc != done);
    allStable = false;                                                          // Mark the whole point set unstable
}

void TileSet::updateStability() {
    stableCount = 0;
    for (int siteID = 0; siteID < total; siteID++) {
        s[siteID].isStable = !s[siteID].isUnstable;
        if (s[siteID].isStable) stableCount++;
    }
}

void TileSet::moveAll(unsigned tileID, unsigned pointNo, unsigned rank, Vector shift, bool doClamp) {
    Vector currentDisplacement = getDisplacement(tileID, pointNo, rank);
    Vector targetDisplacement = currentDisplacement + shift;
    TList &list = tiles[tileID].list[pointNo][rank - pointNo];
    for (int i = 0; i < list.size(); i++) {                                     // Loop through all sites in the group
        int siteID = list[i];                                                   // Retrieve ID of site.
        moveSite(siteID, shift);                                                // Apply shift to each site
    }
    tiles[tileID].p[pointNo][rank-pointNo] = CGAL::ORIGIN + getDisplacement(tileID, pointNo, rank);
}

void TileSet::moveAll(unsigned tileID, unsigned pointNo, unsigned rank, Point target, bool clamp) {
    Point current = CGAL::ORIGIN + getDisplacement(tileID, pointNo, rank);
    moveAll(tileID, pointNo, rank, target - current, clamp);
}

Vector TileSet::getDisplacement(unsigned tileID, unsigned pointNo, unsigned rank) {
    int siteID = tiles[tileID].list[pointNo][rank - pointNo][0];                // Representative site
    return s[siteID].getDisplacement();
}

void TileSet::farthestPoint(unsigned siteID) {
    Site &site = s[siteID];
    if (site.list->size() > 1) conflict(siteID);                                // FPO does not accept averaging; use conflict instead for multi-context ID's
    Points neighbors;                                                           // We use an auxiliary Delaunay triangulation for this purpose
    VC vc = dts[site.rank].incident_vertices(site.vh[0]), done(vc);
    do {                                                                        // Populate the list of neighbors
        neighbors.push_back(vc->point());
    } while (++vc != done);
    DT localDT(neighbors.begin(), neighbors.end());
    Point fp = {-1.0, -1.0};
    double max = 0;
    DT::Finite_faces_iterator fit = localDT.finite_faces_begin();
    for ( ; fit != localDT.finite_faces_end(); fit++) {                         // Iterated through faces of the locally set up Delaunay triangulation
        Point c = localDT.circumcenter(fit);                                    // Circumcenters of these faces are the candidate farthest points
        if (localDT.is_infinite(localDT.locate(c))) continue;                   // A point in infinite face of local DT means it is outside its convex hull; skip
        double rr = CGAL::squared_distance(c, fit->vertex(0)->point());         // Squared circumradius
        if (rr > max) {                                                         // If larger than previous candidate
            fp = c;                                                             // Make this the candidate
            max = rr;                                                           // Update largest circumradius
        }
    }
    double wt = 1.0 / tiles[site.tileID].list[site.pointNo][site.rank - site.pointNo].size();
    moveAll(site.tileID, site.pointNo, site.rank, wt * (fp - site.p), false);
}

void TileSet::conflict(unsigned siteID) {
    Site &site = s[siteID];
    DT &dt = dts[site.rank];
    double &dmin = target_rf[site.rank];
    Point &p = site.p;
    Vector shift[30];
    int nID[30];
    int m = 0;
    VC vc = dt.incident_vertices(site.vh[0]), done(vc);
    do {
        Vector edge = vc->point() - p;
        double l = VL(edge);
        if (l < dmin) {
            shift[m] = (1.001 * dmin/l - 1) * edge;
            nID[m] = vc->info();
        } else nID[m] = -1;
        m++;
    } while (++vc != done);
    for (int i = 0; i < m; i++) {
        if (nID[i] >= 0) {
            moveAll(s[ nID[i] ].tileID, s[ nID[i] ].pointNo, s[ nID[i] ].rank, shift[i]);
        }
    }
}

void TileSet::coverage(unsigned siteID) {
    Site &site = s[siteID];
    DT &dt = dts[site.rank];
    double &rc = target_rc[site.rank];
    Point &p = site.p;
    int nID[30];
    double scale[30];
    Vector edge[30];
    int m = 0;
    FC fc = dt.incident_faces(site.vh[0]), done(fc);
    VC vc = dt.incident_vertices(site.vh[0], fc);
    do {
        vc++;
        edge[m] = site.p - vc->point();
        nID[m] = vc->info();
        if (triangleType (fc) <= 0) {
            Point c = dt.circumcenter(fc);
            double l = VL(c - site.p);
            scale[m] = rc / l;
        } else scale[m] = 2;                                                    // > 1
        m++;
    } while (++fc != done);
    for (int i = 0; i < m; i++) {
        double scl = std::min(scale[i], scale[(i+1)%m]);
        if (scl < 1) {
            Vector shift = (1.001 - scl) * edge[i];
            moveAll(s[ nID[i] ].tileID, s[ nID[i] ].pointNo, s[ nID[i] ].rank, shift);
        }
    }
}

void TileSet::capacity(unsigned siteID) {
    Site &site = s[siteID];
    DT &dt = dts[site.rank];
    Point &p = site.p;
    double &A = avg_area[site.rank];
    double d[30], el[30], pressure;
    double sum_w = 0;
    double a = 0;                                                               // Area of Voronoi cell
    double XProduct;
    FC fc2 = dt.incident_faces(site.vh[0]), fc1(fc2++);                         // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
    VC vc = dt.incident_vertices(site.vh[0], fc2), done(vc);                    // The vertex sharing fc1 anf fc2 with s[siteID].vh[0]
    int m = 0;                                                                  // Number of neighbors
    Vector dir[30];                                                             // Direction vectors to neighbors
    int nID[30];                                                                // Id's of neighbors. We can't use the circulator for updating
    do {
        Point c1 = dt.circumcenter(fc1), c2 = dt.circumcenter(fc2);             // Circumcenters of faces are endpoints of Voronoi cell edge
        XProduct = c1.x() * c2.y() - c1.y() * c2.x();
        a += XProduct;                                                          // Accumulate areas
        el[m] = sqrt((c2 - c1).squared_length());                               // Length of Voronoi edge
        dir[m] = (vc->point() - site.p);
        d[m] = sqrt(dir[m].squared_length());                                   // Distance to neighbor (= 2 x distance to Voronoi edge)
        dir[m] = dir[m] / d[m];                                                 // Normalize direction vector
        nID[m] = vc->info();
        ++fc1;
        ++fc2;
        ++m;
    } while (++vc != done);
    a /= 2;
    double dA = a - A;                                                          // Required expansion or contraction
    if (fabs(dA / A) > max_area_deviation) {
        for (int j = 0; j < m; j++) {
            sum_w += el[j] * el[j];
        }
        pressure = -2 * dA / sum_w;
        for (int j = 0; j < m; j++) {                                           // Loop again through neighbors to give each an appropriately sized nudge
            Vector force = pressure * el[j] * dir[j];
            moveAll(s[ nID[j] ].tileID, s[ nID[j] ].pointNo, s[ nID[j] ].rank, force);
        }
    }
}

void TileSet::lloyd(unsigned siteID) {
    Site &site = s[siteID];
    Vector shift = centroid(dts[site.rank], site.vh[0]);
    moveAll(site.tileID, site.pointNo, site.rank, shift);
}

void TileSet::serial(String seq) {
    maxShift = 0;
    for (int siteID = 0; siteID < total; siteID++)
        s[siteID].isUnstable = false;                                           // Assume all points will be found stable.
    allStable = true;                                                           // Assume the point set will be found stable.
    std::vector<unsigned> order;
    for (unsigned rank = 0; rank < nRanks; rank++) {
        order.resize(N[rank]);
        shuffle(order);
        for (int k = 0; k < seq.length(); k++) {
            for (int i = 0; i < N[rank]; i++) {
                unsigned siteID = beginID[rank] + order[i];
                if (s[siteID].isStable) continue;                               // Skip stable sites
                switch (seq[k]) {
                    case 'a': coverage(siteID); break;
                    case 'b': conflict(siteID); break;
                    case 'c': capacity(siteID); break;
                    case 'f': farthestPoint(siteID); break;
                    case 'l': lloyd(siteID); break;
                }
            }
        }
    }
    if (isAdaptive) {                                                           // Smooth the point locations between ranks and between octaves
        for (int tileID = 0; tileID < tileCount; tileID++) {                    // Update child points across the tile-set to ensure inter-octave coherence
            for (int pointNo = 0; pointNo < nRanks; pointNo++) {
                Points &p = tiles[tileID].p[pointNo];                           // Retrieve the list of point locations
                p.back() = childPoint(tileID, pointNo);                         // The last location is the location of the respective child point
            }
        }
        for (int tileID = 0; tileID < tileCount; tileID++) {                    // Iterate again over all tiles,
            for (int pointNo = 0; pointNo < nRanks; pointNo++) {                //  and all points in the octave
                Points &p = tiles[tileID].p[pointNo];                           // Retrieve the list of point, including child's point at the back
                smooth(p);                                                      // Smooth the point locations to enforce intra-octave coherence
                for (int rank = pointNo; rank < nRanks; rank++) {               // Iterate through the ranks,
                    moveAll(tileID, pointNo, rank, p[rank - pointNo], false); //   and update the point locations in the triangulation
                }
            }
        }
    }
    updateStability();
    maxShift = sqrt(maxShift);
}

void TileSet::lloyd() {
    Vectors v(nRanks * tileCount);
    std::vector<unsigned> count(nRanks * tileCount);
    for (int rank = 1; rank < nRanks; rank++) {
        int nPoints = rank + 1;
        int nDistinct = nPoints * tileCount;
        for (int i = 0; i < nDistinct; i++) {
            v[i] = Vector(0, 0);
            count[i] = 0;
        }
        for (int siteID = beginID[rank]; siteID < endID[rank]; siteID++) {
            unsigned index = s[siteID].pointNo * tileCount + s[siteID].tileID;
            v[index] = v[index] + centroid(dts[rank], s[siteID].vh[0]);
            count[index]++;
        }
        for (int pointNo = 0; pointNo < nPoints; pointNo++) {
            for (int tileID = 0; tileID < tileCount; tileID++) {
                unsigned index = pointNo * tileCount + tileID;
                moveAll(tileID, pointNo, rank, v[index] / count[index]);
            }
        }
    }
}

void TileSet::conflict() {
    Vectors v(nRanks * tileCount);
    std::vector<unsigned> count(nRanks * tileCount);
    for (int rank = 0; rank < nRanks; rank++) {
        int nPoints = rank + 1;
        int nDistinct = tileCount * nPoints;
        for (int i = 0; i < nDistinct; i++) {
            v[i] = Vector(0, 0);
            count[i] = 0;
        }
        double dmin = target_rf[rank];
        double dmin2 = dmin * dmin;
        for (int siteID = beginID[rank]; siteID < endID[rank]; siteID++) {
            Site &site = s[siteID];
            VC vc = dts[rank].incident_vertices(site.vh[0]), done(vc);
            do {
                Vector edge = site.p - vc->point();
                double ll = LL(edge);
                if (ll < dmin2) {
                    double l = sqrt(ll);
                    unsigned index = site.pointNo * tileCount + site.tileID;
                    v[index] = v[index] + 0.5 * (1.001 * dmin/l - 1) * edge;    // 0.5 Because we push both ends of the edge
                    count[index]++;
                }
            } while (++vc != done);
        }
        for (int pointNo = 0; pointNo < nPoints; pointNo++) {
            for (int tileID = 0; tileID < tileCount; tileID++) {
                unsigned index = pointNo * tileCount + tileID;
                moveAll(tileID, pointNo, rank, v[index] / count[index]);
            }
        }
    }
}

void TileSet::coverage() {
    Vectors v(nRanks * tileCount);
    std::vector<unsigned> count(nRanks * tileCount);
    for (int rank = 0; rank < nRanks; rank++) {
        int nPoints = rank + 1;
        int nDistinct = tileCount * nPoints;
        TList count(nDistinct, 0);
        Vectors v(nDistinct);
        double &rc = target_rc[rank];
        double rc2 = rc * rc;
        DT &dt = dts[rank];
        for (int i = 0; i < tileCount; i++) v[i] = Vector(0, 0);
        Point BL(0, 0), TR(one, one);                                           // We will consider only faces in one period
        DT::Finite_faces_iterator fit = dt.finite_faces_begin();
        for ( ; fit != dt.finite_faces_end(); fit++) {                          // Iterate through all finite faces in triangulation
            Point c = dt.circumcenter(fit);
            if (!isInRect(c, BL, TR)) continue;
            double rr = CGAL::squared_distance(c, fit->vertex(0)->point());     // Farthest distance in current face
            if (rr > rc2) {
                double r = sqrt(rr);
                for (int j = 0; j < 3; j++) {
                    Point p = fit->vertex(j)->point();
                    unsigned &siteID = fit->vertex(j)->info();
                    unsigned index = s[siteID].pointNo * tileCount + s[siteID].tileID;
                    v[index] = v[index] + (1.001 - rc/r) * (c - p);
                    count[index]++;
                }
            }
        }
        for (int pointNo = 0; pointNo < nPoints; pointNo++) {
            for (int tileID = 0; tileID < tileCount; tileID++) {
                unsigned index = pointNo * tileCount + tileID;
                moveAll(tileID, pointNo, rank, v[index] / count[index]);
            }
        }
    }
}

void TileSet::parallel(String seq) {
    maxShift = 0;
    for (int siteID = 0; siteID < total; siteID++) {
        s[siteID].isUnstable = false;                                           // Assume all points will be found stable.
    }
    allStable = true;                                                           // Assume the point set will be found stable.
    for (int k = 0; k < seq.length(); k++) {
        switch (seq[k]) {
            case 'a': coverage(); break;
            case 'b': conflict(); break;
            case 'l': lloyd(); break;
        }
    }
    if (isAdaptive) {                                                           // Smooth the point locations between ranks and between octaves
        for (int tileID = 0; tileID < tileCount; tileID++) {                    // Update child points across the tile-set to ensure inter-octave coherence
            for (int pointNo = 0; pointNo < nRanks; pointNo++) {
                Points &p = tiles[tileID].p[pointNo];                           // Retrieve the list of point locations
                p.back() = childPoint(tileID, pointNo);                         // The last location is the location of the respective child point
            }
        }
        for (int tileID = 0; tileID < tileCount; tileID++) {                    // Iterate again over all tiles,
            for (int pointNo = 0; pointNo < nRanks; pointNo++) {                //  and all points in the octave
                Points &p = tiles[tileID].p[pointNo];                           // Retrieve the list of point, including child's point at the back
                smooth(p);                                                      // Smooth the point locations to enforce intra-octave coherence
                for (int rank = pointNo; rank < nRanks; rank++) {               // Iterate through the ranks,
                    moveAll(tileID, pointNo, rank, p[rank - pointNo], false); //   and update the point locations in the triangulation
                }
            }
        }
    }
    updateStability();
    maxShift = sqrt(maxShift);
}

void TileSet::latinize() {
    std::vector<double> x, y;
    std::vector<unsigned> order;
    for (int rank = 0; rank < nRanks; rank++) {
        x.resize(N[rank]);
        y.resize(N[rank]);
        order.resize(N[rank]);
        for (int i = 0; i < N[rank]; i++) {
            int siteID = beginID[rank] + i;
            x[i] = s[siteID].p.x();
            y[i] = s[siteID].p.y();
        }
        int m = (rank + 1) * one;                                               // Number of points per row or column
        double step = 1.0 / m;
        sortIndexes(x, order);
        for (int i = 0; i < N[rank]; i++) {
            x[ order[i] ] = step * (i % m) + rnd(step);
        }
        sortIndexes(y, order);
        for (int i = 0; i < N[rank]; i++) {
            y[ order[i] ] = step * (i % m) + rnd(step);
        }
        for (int tileID = 0; tileID < tileCount; tileID++) {
            for (int pointNo = 0; pointNo <= rank; pointNo++) {
                double avgx(0), avgy(0);
                TList &list = tiles[tileID].list[pointNo][rank - pointNo];
                for (int i = 0; i < list.size(); i++) {
                    int siteID = list[i];
                    avgx += x[siteID - beginID[rank]];
                    avgy += y[siteID - beginID[rank]];
                }
                avgx /= list.size();
                avgy /= list.size();
                moveAll(tileID, pointNo, rank, Point(avgx, avgy), false);
            }
        }
    }
    setAllUnstable();
}

void TileSet::smooth(Points &p) {
    int n = p.size();
    std::vector<double> x(n, 0.0), y(n, 0.0);
    for (int i = 0; i < n; i++) {
        double sum_w(0);
        for (int j = 0; j < n; j++) {
            int offset = std::abs(i - j);
            sum_w += weights[offset];
            x[i] += weights[offset] * p[j].x();
            y[i] += weights[offset] * p[j].y();
        }
        x[i] /= sum_w;
        y[i] /= sum_w;
    }
    for (int i = 0; i < n; i++) {
        p[i] = Point(x[i], y[i]);
    }
}

void TileSet::printText(String baseName) {                                      // Generate a text printout
    char buffer[20];
    int nRanks = nChildren - 1;
    for (int rank = 0; rank < nRanks; rank++) {
        sprintf(buffer, "%02d.txt", rank);
        String fileName = outputPath + baseName + buffer;
        FILE *file = fopen(fileName.c_str(), "w");
        fprintf(file, "%d\n", N[rank]);
        for (unsigned siteID = beginID[rank]; siteID < endID[rank]; siteID++) {
            Point &p = s[siteID].p;
            fprintf(file, "%0.16f %0.16f\n", (1./one) * p.x(), (1./one) * p.y());
        }
        fclose(file);
    }
}

void TileSet::printTable(String fileName) {
    const char *fullFileName = (outputPath + fileName).c_str();                 // Include path and/or prefixes
    FILE *file = fopen(fullFileName, "w");
    if (!file) {
        fprintf(stderr, "Failed to open %s\n", fullFileName);
        exit(1);
    }
    fprintf(file, "%d %d\n", scale, idCount);
    fprintf(file, "%d", one);
    for (int i = 0; i < one; i++) {
        fprintf(file, " %d", idList[i]);
    }
    fprintf(file, "\n");
    int nRanks = nChildren - 1;
    for (int tileID = 0; tileID < tileCount; tileID++) {
        Tile &tile = tiles[tileID];
        for (int i = 0; i < scale * scale; i++) {                               // ID's of children
            fprintf(file, " %4d", tile.children[i]);
        }
        for (int i = 0; i < scale * scale; i++) {                               // Order of children
            fprintf(file, " %2d", tile.order[i]);
        }
        for (int rank = 0; rank < nRanks; rank++) {
            for (int pointNo = 0; pointNo <= rank; pointNo++) {
                Vector d = getDisplacement(tileID, pointNo, rank);
                fprintf(file, "%9.6lf %9.6lf ", d.x(), d.y());
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void TileSet::setAllUnstable() {
    for (int siteID = 0; siteID < total; siteID++) s[siteID].isStable = false;
    allStable = false;
}

void TileSet::setdmin(double v) {
    for (int rank = 0; rank < nRanks; rank++) {
        target_rf[rank] = v * dhex[rank];
    }
}

void TileSet::setrc(double v) {
    for (int rank = 0; rank < nRanks; rank++) {
        target_rc[rank] = v * dhex[rank];
    }
}

void TileSet::set_max_area_deviation(double v) { max_area_deviation = v; }

void TileSet::set_maxStep(double v) { maxStep = v; }

void TileSet::setSmoothingKernel(double SD) {
    isAdaptive = true;
    weights.resize(nChildren);
    for (int i = 0; i < nChildren; i++) {
        weights[i] = exp(-0.5 * i * i / (SD * SD));
        fprintf(stderr, "weights[%2d] = %5.3f\n", i, weights[i]);
    }
}

const char *USAGE_MESSAGE = "Usage: %s [options] <dataFileName>\n"
    "Options:\n"
    "-n <number of iterations>\n"
    "-p <path/prefix of output files>\n"
    "-d <target dmin for spring>\n"
    "-v <normalized maximum deviation of cell areas>\n"
    "-r <target coverage radius>\n"
    "-S <serial optimization sequence>\n"
    "-P <parallel optimization sequence>\n"
    "-T <max step>\n"                                                           // Used to be "relative to dhex", but for regular grid I made it relative to unit
    "-s <standard deviation of smoothing kernel>\n"
    "-c <latinization cycle>"
;

int main(int argc,char **argv) {
    srand(time(NULL)); srand48(time(NULL));                                     // Random seeds to random number generators (int and real)
    int opt;                                                                    // For use by getopt, the command line options utility
    int iterations = 1;                                                         // Number of iterations to apply
    char *fileNamePrefix = (char *)"";                                          // Current working directory is default output directory
    double target_rf = 0.87;                                                         // Target conflict radius
    double rc = 0.67;                                                           // Target coverage radius
    double max_area_deviation = -1;                                                         // Target maximum area deviation in Voronoi cells
    double maxStep = -1;                                                        // A constraint to steps in optimization
    String serialSeq = "";
    String parallelSeq = "";
    int latinizeCycle = -1;
    double smoothingKernelSD = -1;
    while ((opt = getopt(argc, argv, "n:p:S:P:d:r:v:m:T:s:c:")) != -1) {     // Modify default settings with command line options
        switch (opt) {
            case 'n': iterations = atoi(optarg); break;
            case 'p': fileNamePrefix = optarg; break;
            case 'd': target_rf = atof(optarg); break;
            case 'r': rc = atof(optarg); break;
            case 'v': max_area_deviation = atof(optarg); break;
            case 'S': serialSeq = optarg; break;
            case 'P': parallelSeq = optarg; break;
            case 'T': maxStep = atof(optarg); break;
            case 's': smoothingKernelSD = atof(optarg); break;
            case 'c': latinizeCycle = atoi(optarg); break;
            default: fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
        }
    }
    if (optind > argc - 1) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
    }
    String dataFileName = argv[optind];
    std::fstream file(dataFileName);

    TileSet tileSet(file);
    tileSet.setOutputPath(fileNamePrefix);
    tileSet.setdmin(target_rf);
    tileSet.setrc(rc);
    tileSet.set_max_area_deviation(max_area_deviation);
    if (maxStep > 0.) tileSet.set_maxStep(maxStep);
    if (smoothingKernelSD > 0.) tileSet.setSmoothingKernel(smoothingKernelSD);

    int iteration = 0;
    tileSet.setAllUnstable();
    clock_t tOpt0 = clock();
    for ( ; iteration < iterations; iteration++) {                                      // The main optimization loop
        tileSet.serial(serialSeq);
        tileSet.parallel(parallelSeq);
        if ((latinizeCycle > 0) && (iteration % latinizeCycle == 0))
            tileSet.latinize();
        fprintf(stderr, "Iteration %4d, stable: %6d\n", iteration, tileSet.getStableCount());
    }
    clock_t tOpt1 = clock();
    double totalTime = (double)(tOpt1 - tOpt0) / CLOCKS_PER_SEC;
    fprintf(stderr, "Total optimization time: %10.6fs\n", totalTime);
    tileSet.printText("rank-");
    tileSet.printTable("post-optimized.dat");
}
