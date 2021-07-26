/*
 * Optimize ART tilesets
 * 2015-12-29: Created by Abdalla Ahmed
 * 2017-04-26: This version, revised for inclusion with the paper
 */

#include "geometry.h"
#include "words.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <complex>

struct Tile {                                                                   // This record describes a distinct template tile
    Point p;                                                                    // Location of sample point in the tile
    std::vector<unsigned> children;                                             // ids of sub-tiles
    std::vector<unsigned> order;
    std::vector<unsigned> list;                                                 // Lists of instances in the tiling
};

struct Site {                                                                   // This record describes an instance tile in the tiling, and handles triangulation
    unsigned tileID;                                                            // ID of template tile
    Point p;                                                                    // Global coordinates of sample point
    VH vh[9];                                                                   // Handles to up to 9 replicas to maintain toroidal domain
    bool isStable, isUnstable;                                                  // Used to speed up optimization by skipping stable sites
};

class TileSet {                                                                 // Actually this is a representative periodic tiling that we use for optimization
private:
    String S;                                                                   // The constituent string of bits for the optimization environment
    int idWidth;                                                                // Width of identification strings.
    Strings idStrings;                                                          // Distinct length-idWidth sub-strings (factors) of S corresponding to the IDs
    int idCount;                                                                // Number of distinct IDs in 1D
    std::vector<unsigned> idList;                                               // 1D list of consecutive IDs
    std::vector<Tile> tiles;                                                    // List of distinct (template) tiles
    int nTiles;                                                                 // number of distinct tiles, equals idCount x idCount
    int ONE;                                                                    // Length of S, width domain. Instead of using unit domain we prefer 1 unit area per sample.
    int leftMargin, rightMargin;                                                // Number of bits on each side of the tile in the ID string.
    unsigned scale = 0;                                                         // Self-similarity scale in 1D
    unsigned nChildren;                                                         // Number of children tiles; equals scale * scale
    inline Point childPoint(unsigned tileID, unsigned childNo);                 // Location, in parent's coordinates, of the point of the indicated child (in natural order)
    std::vector<Matrix> matrix;                                                 // Transformation matrices to child tiles coordinates
    unsigned morphCount;                                                        // How many times to morph S for next level; equals log_2(scale);
    double target_rf;                                                           // Target conflict radius
    double target_rc;                                                           // Target coverage radius
    double max_area_deviation;                                                  // Target maximum deviation of Voronoi cell areas. Mean area is 1.
    bool allStable;                                                             // To implement termination criteria
    int stableCount;                                                            // Number of stable list in last iteration
    double largestShift;                                                        // Largest shift in last iteration
    double maxStep;                                                             // Maximum allowed shift in ONE iteration
    DT dt;                                                                      // We maintain a Delaunay triangulation for optimization
    int nSites;                                                                 // Total number of sites; that is ONE x ONE.
    std::vector<Site> s;                                                        // List of instance tiles and associated points
    inline double toroidallinearDist(const double &x1, const double &x2);       // 1D Nearest distance between replicas of two points
    inline double toroidalSqDist(const Point &p1, const Point &p2);             // 2D Nearest distance; square to avoid costly square root operation
    inline Point replica(const Point &p, int i);                                // Find ONE of the 9 replicas of a point
    Point marginBL, marginTR;                                                   // Margins to make a toroidal domain
    void moveSite(unsigned siteID, const Point &p);                             // Move the sample point of an instance tile to a given location
    inline void moveSite(unsigned siteID, const Vector& shift);                 // Displace the sample point by a given shift
    void moveAll(unsigned tileID, Vector shift, bool clamp = true);             // Displace the sample points in all tiles with a given ID
    void moveAll(unsigned tileID, const Point &target, bool clamp = true);      // Reposition the sample points in all tiles with a given ID
    Point getDisplacement(int tileID);                                          // Retrieve the relative sample location in an instance tile
    void markUnstable(unsigned siteID);                                         // Mark a site unstable so that it is considered in the following iteration
    void setupPointSet();                                                       // Setup the point set for optimization
    void coverage(unsigned siteID);                                             // Coverage optimization for ONE site, for serial optimization
    void conflict(unsigned siteID);                                             // Conflict optimization for ONE site, for serial optimization
    void capacity(unsigned siteID);                                             // Capacity optimization for ONE site, for serial optimization
    void farthestPoint(unsigned siteID);                                        // Farthest point optimization for ONE site, for serial optimization
    void lloyed(unsigned siteID);                                               // Lloyd optimization for ONE site, for serial optimization
    String outputPath = "";                                                     // Path for output files; e.g. a sub-folder under /tmp
    double energy(int siteID, int seqNo, int start, int end);                   // Energy function for void-and-cluster ranking method.
    inline Tile &T(int siteID) { return tiles[ s[siteID].tileID ]; };           // A handy shortcut for retrieving the template tile of a given instance
public:
    TileSet(String S, int idWidth, String path = "");                           // Constructor. We output many files for testing purposes
    void set_target_rf(double value) { target_rf = value * dhex; };
    void set_target_rc(double value) { target_rc = value * dhex; };
    void set_max_area_deviation(double value) { max_area_deviation = value; }
    void set_maxStep(double value) { maxStep = value; }
    void setScale(unsigned value);
    double getMaxDisplacement();
    void printVectors(String fileName);                                         // Print a list of sample point locations in distinct tiles
    void loadVectors(const char *fileName);                                     // Load a list of sample point locations in distinct tiles
    void printText(String fileName);                                            // Print a list of global sample locations in [0,1) x [0,1)
    void jitter(double weight);                                                 // Apply a random jitter to the sample locations in each distinct tile; for use after loading a file
    double serial(String seq);                                                  // Optimize the point set serially, ONE by ONE.
    void latinize();                                                            // Latinize the point set
    void setAllUnstable();                                                      // Mark all the list unstable so that they are considered in the following iteration
    bool isAllStable() { return allStable; };                                   // Whether all list are stable in last iteration; use to terminate optimization
    void snap(String fileName);                                                 // Make the tile set self-coincident
    void assignRanks(String fileName);                                          // Assign sequencing order of child tiles
    void printData(String fileName);                                            // Print the complete data of the tile set
    void printStatistics(FILE *file, int serialNo);                             // Display running statistics during optimization
    void printStatisticsHeader(FILE *file);                                     // Print the header line for the statistics
};

TileSet::TileSet(String S, int idWidth, String path) {                          // Constructor
    this->S = S;
    this->idWidth = idWidth;
    outputPath = path;
    ONE = S.length();
    idList.resize(ONE);                                                         // 1D Lists of factor IDs at each position
    leftMargin = idWidth - 1;                                                   // For convenience we make the whole idString on the left; though it is better centered
    rightMargin = idWidth - leftMargin - 1;
    String X = S.substr(ONE - leftMargin) + S + S.substr(0, rightMargin);       // Make the string periodic
    for (int i = 0; i < ONE; i++) {                                             // Loop through all 1D coordinates
        idList[i] = findOrInsert(X.substr(i, idWidth), idStrings);              // Find the idString if known, or define a new one.
    }
    idCount = idStrings.size();
    nTiles = idCount * idCount;
    fprintf(stderr, "idCount = %d, nTiles = %d\n", idCount, nTiles);
    scale = 0;                                                                  // No recursion by default
    tiles.resize(nTiles);
    for (unsigned tileID = 0; tileID < nTiles; tileID++) {
        tiles[tileID].p = Point(drand48(), drand48());                          // Initialize sample locations randomly
    }
    double margin = 2.0;                                                        // Margin for toroidal domain, two tiles wide should be sufficient
    marginBL = Point(-margin, -margin);                                         // Bottom-left of primary period + margin
    marginTR = Point(ONE + margin, ONE + margin);                               // Top-right. In our convention BL is included, TR is excluded
    nSites = ONE * ONE;                                                         // Total number of sites
    setupPointSet();
    // Initialize optimization parameters
    maxStep = ONE;                                                              // unbound maximum step.
    set_target_rf(0.87);
    set_target_rc(0.67);
    set_max_area_deviation(0.03);
}

void TileSet::setupPointSet() {
    dt.clear();                                                                 // Clear the Delaunay triangulation to rebuild from scratch, used for loading
    for (int tileID = 0; tileID < nTiles; tileID++) {                           // Clear list of instances; even though it doesn't change
        tiles[tileID].list.clear();
    }
    s.resize(nSites);
    fprintf(stderr, "Constructing %d x %d = %d sites\n", ONE, ONE, nSites);
    for (unsigned Y = 0; Y < ONE; Y++) {
        for (unsigned X = 0; X < ONE; X++) {
            unsigned siteID = Y * ONE + X;
            unsigned tileID = idList[Y] * idCount + idList[X];
            s[siteID].tileID = tileID;
            tiles[tileID].list.push_back(siteID);                               // Add site to the list of the template tile; needed to coordinate the displacements
            Point p = tiles[tileID].p + Vector(X, Y);
            s[siteID].p = p;                                                    // Save a handy copy of point coordinates
            for (int i = 0; i < 9; i++) {                                       // We loop through the 9 replicas,
                if (isInRect(replica(p, i), marginBL, marginTR)) {              // if the location of a replica is within margin
                    s[siteID].vh[i] = dt.insert(replica(p, i));                 // insert replica in triangulation and keep handle to it
                    s[siteID].vh[i]->info() = siteID;                           // Point the vertex back to the site
                }
                else s[siteID].vh[i] = NULL;
            }
        }
    }
}

double TileSet::toroidallinearDist(const double &x1, const double &x2) {        // 1D Nearest distance between replicas of two points
    double dx = x1 - x2;                                                        // Find distance in primary period
    if (dx > ONE/2) dx -= ONE;                                                  // If larger than half the period length another replica is certainly nearer
    if (dx < -ONE/2) dx += ONE;                                                 // Same, but opposite ordering of points
    return dx;
}

double TileSet::toroidalSqDist(const Point &p1, const Point &p2) {              // 2D Nearest distance; square to avoid costly square root operation
    double dx = toroidallinearDist(p1.x(), p2.x());
    double dy = toroidallinearDist(p1.y(), p2.y());
    return dx * dx + dy * dy;
}

Point TileSet::replica(const Point &p, int i) {                                 // Find ONE of the 9 replicas of a point
    i = (i+4) % 9;                                                              // We make the middle replica at index 0
    double x = p.x() + (i%3 - 1) * ONE;                                         // Add -ONE, 0, or ONE to x
    double y = p.y() + (i/3 - 1) * ONE;                                         // Same for y
    return Point(x, y);
}

void TileSet::setScale(unsigned value) {
    if (value <= 1) return;
    if (value & (value-1)) {
        fprintf(stderr, "Subdivision scale must be a power of 2\n");
        exit(1);
    }
    scale = value;
    nChildren = scale * scale;
    typedef std::vector<unsigned> List;
    // Morph IDs in 1D
    std::vector<List> childID(idCount);
    int totalLength = idWidth * scale;                                          // Total length of morphed ID string
    int requiredLength = idWidth + scale - 1;                                   // Required length to contain children idStrings
    morphCount = round(log2(scale));
    for (unsigned parentID = 0; parentID < idCount; parentID++) {
        String morphed = THMorph(idStrings[parentID], morphCount);
        morphed = morphed.substr(totalLength - requiredLength);
        for (int i = 0; i < scale; i++) {
            unsigned id = findOrAbort(morphed.substr(i, idWidth), idStrings);
            childID[parentID].push_back(id);
        }
    }
    // Now compute 2D morphing
    for (unsigned idy = 0; idy < idCount; idy++) {
        for (unsigned idx = 0; idx < idCount; idx++) {
            unsigned tileID = idy * idCount + idx;
            for (unsigned i = 0; i < scale; i++) {
                for (unsigned j = 0; j < scale; j++) {
                    unsigned id = childID[idy][i] * idCount + childID[idx][j];
                    tiles[tileID].children.push_back(id);
                }
            }
        }
    }
    /************* Populate the children transformation matrix: ***************/
    double step = 1.0 / scale;
    matrix.resize(nChildren);
    for (int i = 0; i < scale; i++) {
        for (int j = 0; j < scale; j++) {
            matrix[i * scale + j] = {step, step * j, step * i};
        }
    }
    // Resize ordering lists
    for (unsigned tileID = 0; tileID < nTiles; tileID++) {
        tiles[tileID].order.resize(nChildren);
    }
}

Point TileSet::childPoint(unsigned tileID, unsigned childNo) {                  // Coordinates of the sample point of indicated child relative to parent's origin
    unsigned childID = tiles[tileID].children[childNo];
    //Point p = getDisplacement(childID);
    //return transform(p, matrix[childNo]);
    return transform(tiles[childID].p, matrix[childNo]);
}

void TileSet::assignRanks(String fileName) {
    snap("");                                                                   // make the set self-similar.
    double inv = 1.0 / scale;
    for (unsigned tileID = 0; tileID < nTiles; tileID++) {
        std::vector<unsigned> &order = tiles[tileID].order;
        shuffle(order);                                                         // Initialize to a random order
        Point &p = tiles[tileID].p;
        unsigned i = unsigned(scale * p.y());                                   // Row of children containing parent's sample point
        unsigned j = unsigned(scale * p.x());                                   // Column ~
        unsigned heir = i * scale + j;                                          // sequence number of heir child
        for (int i = 0; i < nChildren; i++) {
            if (order[i] == heir) {
                std::swap(order[0], order[i]);                                  // Assign order[0] to heir
                break;
            }
        }
    }
    std::vector<unsigned> order(nTiles);
    for (int seqNo = 1; seqNo < nChildren/2; seqNo++) {
        fprintf(stderr, "SeqNo: %3d\n", seqNo);
        for (int iteration = 0; iteration < 200; iteration++) {
            bool allStable = true;
            int stableCount(0);
            shuffle(order);
            for (int i = 0; i < nTiles; i++) {
                int tileID = order[i];
                double minEnergy(10e12);
                int minSeqNo;
                for (int candidate = seqNo; candidate < nChildren; candidate++) {
                    double a(0);
                    std::vector<unsigned> &list = tiles[tileID].list;
                    for (int j = 0; j < list.size(); j++) {
                        int siteID = list[j];
                        double e = energy(siteID, candidate, 0, seqNo);
                        a += e;
                    }
                    a /= list.size();
                    if (a < minEnergy) {
                        minEnergy = a;
                        minSeqNo = candidate;
                    }
                }
                if (minSeqNo == seqNo) stableCount++;
                if (minSeqNo != seqNo) {
                    allStable = false;
                    std::swap(tiles[tileID].order[seqNo], tiles[tileID].order[minSeqNo]);
                }
            }
            fprintf(stderr, "Iteration %3d, #stable: %5d\n", iteration, stableCount);
            if (allStable) break;
        }
    }
    for (int seqNo = nChildren/2; seqNo < nChildren - 1; seqNo++) {
        fprintf(stderr, "SeqNo: %3d\n", seqNo);
        for (int iteration = 0; iteration < 200; iteration++) {
            bool allStable = true;
            int stableCount(0);
            shuffle(order);
            for (int i = 0; i < nTiles; i++) {
                int tileID = order[i];
                double selectedEnergy(10e12);
                int selectedSeqNo;
                for (int candidate = seqNo; candidate < nChildren; candidate++) {
                    double a(0);
                    std::vector<unsigned> &list = tiles[tileID].list;
                    for (int j = 0; j < list.size(); j++) {
                        int siteID = tiles[tileID].list[j];
                        double e = energy(siteID, candidate, seqNo + 1, nChildren - 1);
                        a -= e;
                    }
                    a /= list.size();
                    if (a < selectedEnergy) {
                        selectedEnergy = a;
                        selectedSeqNo = candidate;
                    }
                }
                if (selectedSeqNo == seqNo) stableCount++;
                if (selectedSeqNo != seqNo) {
                    allStable = false;
                    std::swap(tiles[tileID].order[seqNo], tiles[tileID].order[selectedSeqNo]);
                }
            }
            fprintf(stderr, "Iteration %3d, #stable: %5d\n", iteration, stableCount);
            if (allStable) break;
        }
    }
}

const double R = 4.;                                                            // Size of neighborhood considered in ranking

double filter(double rr) {
    if (rr > R * R) return 0;
    /***************************************************************************
     * Original setting in Ulickney's void-and-cluster is:                     *
     *      base = exp(-1/(2*sigma^2));                                        *
     * where sigma = 1.5 / scale;                                              *
     * But we experimentally found this setting to work well.                  *
     **************************************************************************/
    const double base = 0.0001;
    return pow(base, rr);
}

double TileSet::energy(int refSID, int refSeqNo, int startSeqNo, int endSeqNo) {
    int X0 = refSID % ONE;
    int Y0 = refSID / ONE;
    int refChildNo = T(refSID).order[refSeqNo];
    Point ref = childPoint(s[refSID].tileID, refChildNo) + Vector(X0, Y0);
    double result(0.0);
    for (int i = -R; i <= R; i++) {
        int Y = (Y0 + ONE + i) % ONE;
        for (int j = -R; j <= R; j++) {
            int X = (X0 + j + ONE) % ONE;
            int siteID = Y * ONE + X;
            for (int seqNo = startSeqNo; seqNo <= endSeqNo; seqNo++) {
                int childNo = T(siteID).order[seqNo];
                if ((i == 0) && (j == 0) && (childNo == refChildNo)) continue;      // Skip the point itself
                Point p = childPoint(s[siteID].tileID, childNo) + Vector(X, Y);
                double rr = toroidalSqDist(ref, p);
                result += filter(rr);
            }
        }
    }
    return result;
}

void TileSet::snap(String fileName) {
    if (scale <= 1) return;
    std::vector<unsigned> i(nTiles), j(nTiles), targetID(nTiles);               // row, column, and ID of nominated heir child
    for (unsigned tileID = 0; tileID < nTiles; tileID++) {
        Point parentPoint = getDisplacement(tileID);
        double min = 1.0;
        for (unsigned childNo = 0; childNo < nChildren; childNo++) {
            Point p = childPoint(tileID, childNo);
            double sqDist = LL(p - parentPoint);
            if (sqDist < min) {
                min = sqDist;
                targetID[tileID] = tiles[tileID].children[childNo];
                i[tileID] = childNo / scale;
                j[tileID] = childNo % scale;
            }
        }
    }
    unsigned maxDepth = 31 / morphCount;                                        // We need morphCount bits for each level
    double inv = 1.0 / (1lu << (maxDepth * morphCount));
    for (int tileID = 0; tileID < nTiles; tileID++) {
        unsigned x(0), y(0);
        int id = tileID;
        for (int depth = 0; depth < maxDepth; depth++) {
            x = (x << morphCount) + j[id];
            y = (y << morphCount) + i[id];
            id = targetID[id];
        }
        Point target = Point(inv * x, inv * y);
        moveAll(tileID, target, false);
    }
    if (!fileName.empty()) printVectors(fileName);
}

void TileSet::latinize() {
    std::vector<double> x(nSites), y(nSites);
    for (int siteID = 0; siteID < nSites; siteID++) {
        x[siteID] = s[siteID].p.x();
        y[siteID] = s[siteID].p.y();
    }
    double step = 1.0 / ONE;
    std::vector<unsigned> order(nSites);
    sortIndexes(x, order);
    for (int i = 0; i < nSites; i++) {
        int siteID = order[i];
        x[siteID] = step * (i % ONE) + drand48() * step;
    }
    sortIndexes(y, order);
    for (int i = 0; i < nSites; i++) {
        int siteID = order[i];
        y[siteID] = step * (i % ONE) + drand48() * step;
    }
    for (int tileID = 0; tileID < nTiles; tileID++) {
        double avgx(0), avgy(0);
        std::vector<unsigned> &list = tiles[tileID].list;
        for (int i = 0; i < list.size(); i++) {
            int siteID = tiles[tileID].list[i];
            avgx += x[siteID];
            avgy += y[siteID];
        }
        avgx /= list.size();
        avgy /= list.size();
        Point target = Point(avgx, avgy);
        moveAll(tileID, target, false);
    }
}

Point TileSet::getDisplacement(int tileID) {
    int siteID = tiles[tileID].list[0];                                                  // Index of the first point in the group block; to represent the group.
    int X = siteID % ONE;
    int Y = siteID / ONE;
    return s[siteID].p - Vector(X, Y);
}

void TileSet::printVectors(String fileName) {                                   // Print shift vectors of groups; for use by rollout program
    const char *fullFileName = (outputPath + fileName).c_str();                 // Include path and/or prefixes
    FILE *file = fopen(fullFileName, "w");
    if (!file) {
        fprintf(stderr, "Failed to open %s\n", fullFileName);
        exit(1);
    }
    double max = 0;                                                             // Maximum distance traveled by points after relaxation
    for (int tileID = 0; tileID < nTiles; tileID++) {                                        // Iterate through groups
        Point &p = tiles[tileID].p;
        fprintf(file, "%lf %lf\n", p.x(), p.y());                               // Print it out
    }
    fclose(file);
}

void TileSet::printData(String fileName) {                                      // For use by rollout programs
    const char *fullFileName = (outputPath + fileName).c_str();                 // Include path and/or prefixes
    FILE *file = fopen(fullFileName, "w");
    if (!file) {
        fprintf(stderr, "Failed to open %s\n", fullFileName);
        exit(1);
    }
    assignRanks("");
    fprintf(file, "%d %d\n", scale, idCount);
    fprintf(file, "%d", ONE);
    for (int i = 0; i < ONE; i++) {
        fprintf(file, " %d", idList[i]);                                        // A periodic list of IDs that the roll out program may use
    }
    fprintf(file, "\n");
    for (int tileID = 0; tileID < nTiles; tileID++) {                           // Iterate through tiles
        Point p = getDisplacement(tileID);
        fprintf(file, "%12.9lf %12.9lf", p.x(), p.y());                         // Location of sample point
        for (int i = 0; i < nChildren; i++) {
            fprintf(file, " %4d", tiles[tileID].children[i]);                   // IDs of children
        }
        for (int i = 0; i < nChildren; i++) {
            fprintf(file, " %2d", tiles[tileID].order[i]);                      // Order of children
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void TileSet::loadVectors(const char *fileName) {                               // Load shift vectors saved by printVectors; useful for resuming optimization
    FILE *vectorsFile = fopen(fileName, "r");
    if (!vectorsFile) {
        fprintf(stderr, "Failed to open %s\n", fileName);
        exit(1);
    }
    for (int tileID = 0; tileID < nTiles; tileID++) {                                    // Loop through all groups (in the same order as for saving)
        double dx, dy;
        fscanf(vectorsFile, " %lf %lf", &dx, &dy);                              // Read a shift vector entry
        if (!feof(vectorsFile)) {                                               // If successfully read
            tiles[tileID].p = Point(dx, dy);
        } else {                                                            // End of file reached?
            fprintf(stderr, "Failed to read all vectors from %s\n", fileName);
            exit(1);
        }
    }
    setupPointSet();
}

void TileSet::jitter(double weight) {                                           // Apply a random jitter to each template tile
    for (int tileID = 0; tileID < nTiles; tileID++) {
        Point current = getDisplacement(tileID);
        Point random = {drand48(), drand48()};
        moveAll(tileID, weight * (random - current), false);
    }
}

void TileSet::printText(String fileName) {                                      // Generate a text printout
    const char *fullFileName = (outputPath + fileName).c_str();
    FILE *file = fopen(fullFileName, "w");
    if (!file) {
        fprintf(stderr, "Failed to open %s\n", fullFileName);
        exit(1);
    }
    fprintf(file, "%d\n", nSites);                                              // Print number of point is first line; a common convention
    double inv = 1. / ONE;
    for (int i = 0; i < nSites; i++) {                                          // Loop through all points
        Point &p = s[i].p;                                                      // Normalize and wrap back to unit torus
        fprintf(file, "%0.16f %0.16f\n", inv * p.x(), inv * p.y());
    }
    fclose(file);
}

void TileSet::markUnstable(unsigned siteID) {
    s[siteID].isUnstable = true;                                                   // Mark the site unstable
    VC vc = dt.incident_vertices(s[siteID].vh[0]), done(vc);
    do {                                                                        // Mark neighbors unstable
        s[ vc->info() ].isUnstable = true;
    } while (++vc != done);
    allStable = false;                                                          // Mark the whole point set unstable
}

void TileSet::moveSite(unsigned siteID, const Point &p) {                       // Adjust location of indexed point ; we left it double to cater for potential extensions(in t-domain) and update in triangulation
    markUnstable(siteID);                                                       // The site should be marked unstable, whether it's moved or not.
    if (
        (int(floor(p.x())) != siteID % ONE) ||
        (int(floor(p.y())) != siteID / ONE)
    ) return;                                                                   // Only move point if target is inside the stratum.
    double l = (p - s[siteID].p).squared_length();
    largestShift = std::max(largestShift, l);
    s[siteID].p = p;                                                            // Save a handy copy of updated point coordinates
    for (int i = 0; i < 9; i++) {
        if (s[siteID].vh[i] != NULL)
            s[siteID].vh[i] = dt.move(s[siteID].vh[i], replica(p, i));
    }
    markUnstable(siteID);                                                       // To mark the new neighbors unstable.
}

void TileSet::moveSite(unsigned siteID, const Vector &shift) {
    moveSite(siteID, s[siteID].p + shift);
}

void TileSet::moveAll(unsigned tileID, Vector shift, bool doClamp) {
    if (doClamp) clamp(shift, maxStep);
    std::vector<unsigned> &list = tiles[tileID].list;
    for (int i = 0; i < list.size(); i++) {                                     // Loop through all instances
        int siteID = list[i];                                                   // Retrieve site ID
        moveSite(siteID, shift);                                                // Apply shift to each site
    }
    tiles[tileID].p = getDisplacement(tileID);                                  // Update the sample location in the template tile
}

void TileSet::moveAll(unsigned tileID, const Point &target, bool doClamp) {
    Point current = getDisplacement(tileID);
    moveAll(tileID, target - current, doClamp);
}

void TileSet::capacity(unsigned siteID) {                                                      // Immediately update a site and neighbors
    double d[20], el[20], pressure;
    double sum_w = 0;
    double a = 0;                                                               // Area of Voronoi cell
    double XProduct;
    FC fc2 = dt.incident_faces(s[siteID].vh[0]), fc1(fc2++);                         // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
    VC vc = dt.incident_vertices(s[siteID].vh[0], fc2), done(vc);                    // The vertex sharing fc1 anf fc2 with s[siteID].vh[0]
    int m = 0;                                                                  // Number of neighbors
    Vector dir[20];                                                             // Direction vectors to neighbors
    int nID[20];                                                                // Id's of neighbors. We can't use the circulator for updating
    do {
        Point c1 = dt.circumcenter(fc1), c2 = dt.circumcenter(fc2);             // Circumcenters of faces are endpoints of Voronoi cell edge
        XProduct = c1.x() * c2.y() - c1.y() * c2.x();
        a += XProduct;                                                          // Accumulate areas
        el[m] = sqrt((c2 - c1).squared_length());                               // Length of Voronoi edge
        dir[m] = (vc->point() - s[siteID].p);
        d[m] = sqrt(dir[m].squared_length());                                   // Distance to neighbor (= 2 x distance to Voronoi edge)
        dir[m] = dir[m] / d[m];                                                 // Normalize direction vector
        nID[m] = vc->info();
        ++fc1;
        ++fc2;
        ++m;
    } while (++vc != done);
    a /= 2;
    double dA = a - 1;                                                          // Required expansion or contraction
    if (fabs(dA) > max_area_deviation) {
        for (int j = 0; j < m; j++) {
            sum_w += el[j] * el[j];
        }
        pressure = -2 * dA / sum_w;
        for (int j = 0; j < m; j++) {                                               // Loop again through neighbors to give each an appropriately sized nudge
            Vector force = pressure * el[j] * dir[j];
            moveAll(s[ nID[j] ].tileID, force);
        }
    }
}

void TileSet::conflict(unsigned siteID) {
    Point &p = s[siteID].p;
    Vector shift[30];
    int nID[30];
    int m = 0;
    VC vc = dt.incident_vertices(s[siteID].vh[0]), done(vc);
    do {
        Vector edge = vc->point() - p;
        double l = VL(edge);
        if (l < target_rf) {
            shift[m] = (1.001 * target_rf/l - 1) * edge;
            nID[m] = vc->info();
        } else nID[m] = -1;
        m++;
    } while (++vc != done);
    /************************************ Test ********************************/
    //if (m == 0) farthestPoint(siteID);
    /**************************************************************************/
    for (int i = 0; i < m; i++) {
        if (nID[i] >= 0) {
            moveAll(s[ nID[i] ].tileID, shift[i]);
        }
    }
}

void TileSet::coverage (unsigned siteID) {                                             // Apply push and pull shifts to optimze conflict and coverage, respectivly
    int nID[30];
    double scale[30];
    Vector edge[30];
    int m = 0;
    FC fc = dt.incident_faces(s[siteID].vh[0]), done(fc);
    VC vc = dt.incident_vertices(s[siteID].vh[0], fc);
    do {
        vc++;
        edge[m] = s[siteID].p - vc->point();
        nID[m] = vc->info();
        if (triangleType (fc) <= 0) {
            Point c = dt.circumcenter(fc);
            double l = VL(c - s[siteID].p);
            scale[m] = target_rc / l;
        } else scale[m] = 2;                                                    // > 1
        m++;
    } while (++fc != done);
    for (int i = 0; i < m; i++) {
        double scl = std::min(scale[i], scale[(i+1)%m]);
        if (scl < 1) {
            Vector shift = (1.001 - scl) * edge[i];
            moveAll(s[ nID[i] ].tileID, shift);
        }
    }
}

void TileSet::lloyed(unsigned siteID) {
    moveAll(s[siteID].tileID, centroid(dt, s[siteID].vh[0]));
}

void TileSet::farthestPoint(unsigned siteID) {
    if (T(siteID).list.size() > 1) { conflict(siteID); return; }                // FPO is not average-able; use conflict instead if there is more than one instance
    Points neighbors;                                                           // We use an auxiliary Delaunay triangulation for this purpose
    VC vc = dt.incident_vertices(s[siteID].vh[0]), done(vc);
    do {                                                                        // Populate the list of neighbors
        neighbors.push_back(vc->point());
    } while (++vc != done);
    DT localDT(neighbors.begin(), neighbors.end());
    Point fp = {-1.0, -1.0};
    bool fpFound = false;
    double max = 0;
    DT::Finite_faces_iterator fit = localDT.finite_faces_begin();
    for ( ; fit != localDT.finite_faces_end(); fit++) {                         // Iterated through faces of the locally set up Delaunay triangulation
        Point c = localDT.circumcenter(fit);                                    // Circumcenters of these faces are the candidate farthest points
        if (localDT.is_infinite(localDT.locate(c))) continue;                   // A point in infinite face of local DT means it is outside its convex hull; skip
        fpFound = true;
        double rr = CGAL::squared_distance(c, fit->vertex(0)->point());         // Squared circumradius
        if (rr > max) {                                                         // If larger than previous candidate
            fp = c;                                                             // Make this the candidate
            max = rr;                                                           // Update largest circumradius
        }
    }
    moveAll(s[siteID].tileID, fp - s[siteID].p, false);
}

double TileSet::serial(String seq) {
    largestShift = 0;
    for (int siteID = 0; siteID < nSites; siteID++) {
        s[siteID].isUnstable = false;                                           // Assume all points will be found stable.
    }
    allStable = true;                                                           // Assume the point set will be found stable.
    std::vector<unsigned> order(nSites);
    shuffle(order);
    for (int i = 0; i < nSites; i++) {
        int siteID = order[i];                                                     // Translate index according to the random order
        if (s[siteID].isStable) continue;                                        // Skip stable points
        for (int k = 0; k < seq.length(); k++) {
            switch (seq[k]) {
                case 'a': coverage(siteID); break;
                case 'b': conflict(siteID); break;
                case 'c': capacity(siteID); break;
                case 'f': farthestPoint(siteID); break;
                case 'l': lloyed(siteID); break;
            }
        }
    }
    stableCount = 0;
    for (int siteID = 0; siteID < nSites; siteID++) {
        s[siteID].isStable = !s[siteID].isUnstable;
        if (s[siteID].isStable) stableCount++;
    }
    largestShift = sqrt(largestShift);
    return largestShift;                                                            // Return maximum shift which can be used to monitor convergence
}

void TileSet::setAllUnstable() {
    for (int siteID = 0; siteID < nSites; siteID++) {
        s[siteID].isStable = false;
    }
    allStable = false;
}

void TileSet::printStatistics(FILE *output, int serialNo)  {                    // Printout spatial statistics
    // Pre-compute circumcenters that will be used in more than a measure:
    DT::Finite_faces_iterator fit = dt.finite_faces_begin();
    for ( ; fit != dt.finite_faces_end(); fit++) {                              // Iterate through all (finite) faces in triangulation
        fit->info().c = dt.circumcenter(fit);                                   // Circumcenter of face is one end of a Voronoi edge
    }
    // --------- Conflict radius and BOO
    double rf(ONE);                                                             // rf can't go beyond ONE :)
    double rf_avg(0);                                                           // average conflict radius
    double orientorder;                                                         // All BOO computations are taken from publicly available PSA code.
    std::complex<double> acc = 0;
    unsigned long nacc = 0;
    for (int i = 0; i < nSites; i++) {                                          // Loop through all sites
        double rf_local = ONE;
        std::complex<double> localacc = 0;
        VC vc = dt.incident_vertices(s[i].vh[0]), done(vc), next;
        do {
            next = vc; ++next;
            const Point &v1 = vc->point();
            const Point &v2 = next->point();
            // Local rf
            double dist = (v1 - s[i].p).squared_length();
            rf_local = std::min(rf_local, dist);
            // Orientational order
            std::complex<double> c1(v1.x(), v1.y());
            std::complex<double> c2(v2.x(), v2.y());
            localacc += std::polar(1.0, 6.0 * arg(c1 - c2));
            ++nacc;
        } while (++vc != done);
        rf_local = sqrt(rf_local);
        rf = std::min(rf, rf_local);
        rf_avg += rf_local;
        acc += abs(localacc);
    }
    rf /= dhex;
    rf_avg /= nSites * dhex;
    orientorder = abs(acc) / nacc;                                              // I don't quite get why we should average this once rather than average the averages
    // -------- Coverage Radius
    double rc(0), rc_avg(0);
    unsigned count(0);                                                          // Count the number of faces for cross-checking
    Point BL(0, 0), TR(ONE, ONE);                                               // We consider only faces in ONE period
    for (fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); fit++) {  // Iterate through all finite faces in triangulation
        Point &c = fit->info().c;                                               // Retrieve the precomputed circumcenter
        if (isInRect(c, BL, TR)) {                                              // If it is in the bounding box we consider the face
            double r = VL(fit->vertex(0)->point() - c);                         // Farthest distance (coverage radius) in current face
            rc = std::max(rc, r);                                               // If larger than current farthest distance update the latter
            rc_avg += r;
            count++;
        }
    }
    assert(count == 2 * nSites);                                                     // This is a property of the toroidal domain
    rc /= dhex;
    rc_avg /= count * dhex;
    // -------- Valence Frequency
    std::vector<unsigned> histogram(30, 0);                                     // A simple histogram of 30 buckets; 0, 1, and 2 would stay empty
    for (int i = 0; i < nSites; i++) {                                               // Iterate through all sites
        VC vc = dt.incident_vertices(s[i].vh[0]), done(vc);
        int m(0); do { m++; } while (++vc != done);                             // Count number of neighbors
        ++histogram[m];                                                         // Increment the respective histogram bucket
    }
    double N4 = 100.0 * (double)histogram[4] / nSites;                               // Compute percentages
    double N5 = 100.0 * (double)histogram[5] / nSites;                               // 3 and 9+ are very unlikely after optimization
    double N6 = 100.0 * (double)histogram[6] / nSites;
    double N7 = 100.0 * (double)histogram[7] / nSites;
    double N8 = 100.0 * (double)histogram[8] / nSites;

    // -------- capacity variance and CVT energy
    double meanSquaredArea = 0;
    double cvte = 0;
    double max_area_deviation(0);
    for (int i = 0; i < nSites; i++) {
        Point &p = s[i].p;
        double area = 0;                                                        // Area of Voronoi cell
        FC fc2 = dt.incident_faces(s[i].vh[0]), fc1(fc2++), done(fc1);          // fc1 and fc2 are two consecutive (ccw) faces incident to current vertex
        VC vc = dt.incident_vertices(s[i].vh[0], fc2);                          // The vertex sharing fc1 anf fc2 with v[i].vh
        do {
            Point &c1 = fc1->info().c, &c2 = fc2->info().c;                     // Circumcenters of faces are the endpoints of Voronoi cell edge
            double XProduct = c1.x() * c2.y() - c1.y() * c2.x();
            area += XProduct;                                                   // Accumulate areas
            Point mid = p + 0.5 * (vc->point() - p);                            // Mid point to neighbor site
            cvte += cvtEnergy(p, c1, c2, mid);
            ++fc2;
            ++vc;
        } while (++fc1 != done);
        area *= 0.5;
        max_area_deviation = std::max(max_area_deviation, fabs(area - 1));                              // 1 is the average capacity
        meanSquaredArea += area * area;
    }
    cvte /= nSites;
    double areaVariance = (meanSquaredArea / nSites) - 1;                            // E(a^2) - (E(a)) ^ 2; where E(a) = 1;
    // -------- Print the record:
    fprintf(
        output,
        "%5d  "                                                                 // Iteration number, supplied
        "%5d "                                                                  // Stable count
        "%9.6f "                                                                // Largest shift
        "%09.7f  %09.7f "                                                       // Standard and maximum deviation of capacity
        "{%4.1f, %4.1f, %4.1f, %4.1f, %4.1f}  "                                 // Valence histogram
        "%0.3f  %0.3f  "                                                        // rc rc_avg
        "%0.3f  %0.3f  "                                                        // rf rf_avg
        "%0.3f  "                                                               // BOO
        "%0.4f\n",                                                              // CVT energy
        serialNo,
        stableCount,
        largestShift,
        sqrt(areaVariance), max_area_deviation,
            N4, N5, N6, N7, N8,
            rc, rc_avg,
            rf, rf_avg,
            orientorder,
            cvte
    );
}

void TileSet::printStatisticsHeader(FILE *file) {
    fprintf(file,
            "    i stable  maxShift    sdA        mdA      valence 4,5,6,7,8 per cent     rcmax  rcavg  rfmin  rfavg   BOO   Ecvt\n");
    /********    0      0  0.295730 0.0721181  0.3593222 { 0.5, 21.6, 56.2, 21.0,  0.8}  0.806  0.616  0.645  0.829  0.398  0.0851 */
}

const char *USAGE_MESSAGE = "Usage: %s [options] <bitString> <idWidth>\n"
    "Options:\n"
    "-n <number of iterations>\n"
    "-p <path/prefix of output files>\n"
    "-S <scale of self-similarity>\n"
    "-l <path to vectors table to load>\n"
    "-j <size of pre-optimization jitter (in grid scale)>\n"
    "-d <target conflict radius>\n"
    "-r <target coverage radius>\n"
    "-v <target maximum deviation of cell areas>\n"
    "-T <max step in one iteration>\n"
    "-q <optimization sequence>\n"
    "-c <snapping cycle>\n"
    "-C <latinize cycle>\n"
    "Optimization sequences:\n"
    "  a: coverage, b: conflict, c: capacity, f: FPO, l: Lloyd\n"
;

int main(int argc,char **argv) {
    srand(time(NULL)); srand48(time(NULL));                                     // Random seeds to random number generators (int and real)
    int opt;                                                                    // For use by getopt, the command line options utility
    int iterations = 0;                                                         // Number of iterations to apply
    char *fileNamePrefix = (char *)"";                                          // Current working directory is default output directory
    double jitter = 0.;                                                         // Size of jitter applied prior to optimization; a value around 1 is good
    char *vectorsFileName = NULL;                                               // Path to file for loading a shift-ve    for (unsigned tileID = 0; tileID < nTiles; tileID++) {
    unsigned selfSimilarityScale = 0;
    double target_rf = 0.87;                                                         // Target dmin for use by our spring relaxation
    double target_rc = 0.67;                                                           // Target coverage radius
    double max_area_deviation = -1;                                                         // Target maximum area deviation in Voronoi cells
    double maxStep = -1;                                                        // A constraint to steps in optimization
    String optimizationSequence = "lf";
    int snappingCycle = 0;
    int latinize_cycle = 0;
    while ((opt = getopt(argc, argv, "n:p:q:l:j:d:r:v:T:S:c:C:")) != -1) {     // Modify default settings with command line options
        switch (opt) {
            case 'n': iterations = atoi(optarg); break;
            case 'S': selfSimilarityScale = atoi(optarg); break;
            case 'p': fileNamePrefix = optarg; break;
            case 'l': vectorsFileName = optarg; break;
            case 'j': jitter = atof(optarg); break;
            case 'd': target_rf = atof(optarg); break;
            case 'r': target_rc = atof(optarg); break;
            case 'v': max_area_deviation = atof(optarg); break;
            case 'q': optimizationSequence = optarg; break;
            case 'T': maxStep = atof(optarg); break;
            case 'c': snappingCycle = atoi(optarg); break;
            case 'C': latinize_cycle = atoi(optarg); break;
            default: fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
        }
    }
    if (optind > argc - 2) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
    }
    String bitString = argv[optind];
    int idWidth = atoi(argv[optind+1]);
    TileSet tileSet(bitString, idWidth, fileNamePrefix);             // Create the pattern object
    tileSet.printStatisticsHeader(stderr);
    if (vectorsFileName != NULL) {                                              // If a vectors table is supplied load and apply it
        tileSet.loadVectors(vectorsFileName);
        tileSet.printStatistics(stderr, -5);
        tileSet.printStatistics(stderr, -4);
        tileSet.printText("loaded.txt");
        if (jitter != 0) {
            tileSet.jitter(jitter);
            tileSet.printStatistics(stderr, -3);
            tileSet.printText("jittered.txt");
        }
    }
    tileSet.set_target_rf(target_rf);                                                           // Set target dmin for spring-based relaxations
    tileSet.set_target_rc(target_rc);
    if (selfSimilarityScale > 0) tileSet.setScale(selfSimilarityScale);
    if (maxStep > 0) tileSet.set_maxStep(maxStep);
    if (max_area_deviation >= 0) tileSet.set_max_area_deviation(max_area_deviation);
    char fName[200];                                                            // A buffer to compile output file names
    int iteration = 0;
    tileSet.setAllUnstable();
    clock_t tOpt0 = clock();
    for ( ; iteration < iterations; iteration++) {                                              // The main optimization loop
        tileSet.serial(optimizationSequence);
        if ((latinize_cycle > 0) && (iteration % latinize_cycle == 0))
            tileSet.latinize();
        if ((snappingCycle > 0) && (iteration % snappingCycle == 0))
            tileSet.snap("");
        tileSet.printStatistics(stderr, iteration);
        if (tileSet.isAllStable()) { break; }
    }
    clock_t tOpt1 = clock();
    double totalTime = (double)(tOpt1 - tOpt0) / CLOCKS_PER_SEC;
    fprintf(stderr, "Total optimization time: %10.6fs\n", totalTime);
    tileSet.printText("before-snapping.txt");
    tileSet.printVectors("before-snapping.dat");
    if (selfSimilarityScale > 0) {
        tileSet.snap("self-similar.dat");
        tileSet.printText("self-similar.txt");
        tileSet.printData("table.dat");
    }
}
