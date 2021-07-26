/* Test recursive tiles set for halftoning
 * 2015-10-16: Created Abdalla Ahmed
 * 2017-01-13: Adapted to use new tile structure ranked by order, not dmin.
 */

#include <algorithm>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <locale.h>

struct Point {
    double x, y;
};

typedef std::vector<Point> Points;
typedef std::vector<int> TList;                                                 // List of indexes

class CImage {
    int w, h;
    int size;                                                                   // largest of w, h
    std::vector<int> bmp;
    int getpixel(int x, int y) {
        if (x >= 0 && x < w && y >= 0 && y < h)
            return bmp[y * w + x];
        return 0;
    }
public:
    CImage(const char *fname, double colorScale);                               // 256 grades, no comments, no errors.
    int width() { return w; };
    int height() { return h; };
    int getDensity(Point p) { return getpixel(size * p.x, size * p.y); };
    int getMaxDensity(Point bottomLeft, Point topRight);
};

CImage::CImage(const char *fileName, double colorScale) {
    FILE *pgmfile = fopen(fileName, "r");
    if (!pgmfile) {
        std::cerr << "Failed to open file " << fileName << "\n";
        exit(1);
    }
    int grades, x, y;
    fscanf(pgmfile, "P2 %d %d %d", &w, &h, &grades);
    if (grades != 255) {
        fprintf(stderr, "Error: incorrect image format\n");
        exit(1);
    }
    size =  std::max(w, h);
    bmp.resize(h * w);
    for (int i = 0; i < h * w; i++) {
        int c;
        fscanf(pgmfile, " %d", &c);
        if (!feof(pgmfile)) {
            int x = i%w, y = (h - 1) - (i / w);
            bmp[y * w + x] = round(colorScale * (grades - c));
        } else {
            fprintf(stderr,
                    "Sorry, failed to read all image information\n");
            exit(1);
        }
    }
    fclose(pgmfile);
}

int CImage::getMaxDensity(Point bottomLeft, Point topRight) {
    unsigned left   = std::max(int(size * bottomLeft.x  ), 0);
    unsigned right  = std::min(int(size * topRight.x + 1), w);
    unsigned bottom = std::max(int(size * bottomLeft.y  ), 0);
    unsigned top    = std::min(int(size * topRight.y + 1), h);
    int max(0);
    for (int y = bottom; y < top; y++) {
        for (int x = left; x < right; x++) {
            max = std::max(max, bmp[y * w + x]);
        }
    }
    return max;
}

struct Matrix {                                                                // Postscript-style matrix, but only with scale and translate
    double scale, tx, ty;
    Matrix concat(const Matrix &transformation) const {
        Matrix result;
        result.scale = transformation.scale * scale;
        result.tx = transformation.tx * scale + tx;
        result.ty = transformation.ty * scale + ty;
        return result;
    };
};

inline Point transform(const Point &p, const Matrix &m) {
    double X = m.scale * p.x + m.tx;
    double Y = m.scale * p.y + m.ty;
    return {X, Y};
};

/********************** Global Variables *******/
int idCount, tileCount;
unsigned ONE;
TList idList;
int scale;
int nChildren, nRanks;
const Matrix IDENTITY = {1, 0, 0};
Matrix *M;
/***********************************************/
struct Tile {
    std::vector<Points> p;                                                      // Point coordinates at different ranks
    TList children;                                                             // ids of sub-tiles
    TList order;
    void read(std::fstream &file) {
        children.resize(nChildren);
        for (int i = 0; i < nChildren; i++) {
            file >> children[i];
        }
        order.resize(nChildren);
        for (int i = 0; i < nChildren; i++) {
            file >> order[i];
        }
        p.resize(nRanks);
        for (int rank = 0; rank < nRanks; rank++) {
            p[rank].resize(rank + 1);
            for (int pointNo = 0; pointNo <= rank; pointNo++) {
                Point &point = p[rank][pointNo];
                file >> point.x >> point.y;
            }
        }
    };
};

std::vector<Tile> tiles;

void loadTileSet(std::fstream &file) {
    file >> scale >> idCount;
    file >> ONE;
    idList.resize(ONE);
    for (int i = 0; i < ONE; i++) file >> idList[i];
    nChildren = scale * scale;
    nRanks = nChildren - 1;
    tileCount = idCount * idCount;
    tiles.resize(tileCount);
    for (int id = 0; id < tiles.size(); id++) {
        tiles[id].read(file);
        if (file.eof()) {
            fprintf(stderr, "Failed to load tiles\n");
            exit(1);
        }
    }
    // Populate the matrix:
    M = new Matrix[nChildren];
    double step = 1.0 / scale;
    for (int Y = 0; Y < scale; Y++) {
        for (int X = 0; X < scale; X++) {
            M[Y * scale + X] = {step, step * X, step * Y};
        }
    }
    fprintf(stderr, "Loaded %4d tiles\n", tiles.size());
}


/********************** Global Variables *******/
bool showtiles = false;
int maxDepth = 10;
int minDepth = 1;

/***********************************************/

Points sample(CImage img, int rootID) {
    Points samples;
    samples.reserve(1000000);
    fprintf(stderr, "Root tile = %d\n", rootID);
    const int MAX_TILES = 100000;
    double ratio = 1.0 / scale;
    const double MAX = 1.0 - 1e-6;
    std::vector<Point> probeSet = {{0.,0.}, {0.,MAX}, {MAX,MAX}, {MAX,0.}};
    struct TRecord {
        int tileID;
        int depth;
        Matrix matrix;
        int tileRank;
    };
    std::vector<TRecord> stack;
    stack.reserve(MAX_TILES);
    stack.push_back({rootID, 0, IDENTITY, 0});
    int counter(0);
    int reachedDepth(0);
    std::cout <<
        "%!PS-Adobe-3.0 EPSF-3.0\n"
        "%%BoundingBox: 0 0 1000 1000\n"
        "1000 dup scale\n"
        "/r " << 0.001 << " def\n"
        "/p {r 0 360 arc fill} def\n"
        "r setlinewidth\n"
        "/tile {gsave 1 1 0 setrgbcolor pop dup rectstroke grestore} def\n"
    ;
    do {
        TRecord record = stack.back();
        stack.pop_back();
        Tile &tile = tiles[record.tileID];
        Matrix &matrix = record.matrix;
        int &depth = record.depth;
        int digit = std::pow(nChildren, depth);                                 // Value of the current (depth'th) digit in base nChildren.
        int &tileRank = record.tileRank;
        reachedDepth = std::max(reachedDepth, depth);
        Point bottomLeft = transform({0., 0.}, matrix);
        Point topRight = transform({1., 1.}, matrix);
        int maxDensity = img.getMaxDensity(bottomLeft, topRight);
        int maxRank = nRanks * digit + tileRank;                                // Highest rank in the current tile
        if (showtiles) {
            printf("%f %f %f %d tile\n", matrix.tx, matrix.ty, matrix.scale, depth);
        }
        if ((maxDensity > maxRank) && (record.depth < maxDepth)) {
            for (int seqNo = 0; seqNo < nChildren; seqNo++) {
                unsigned childNo = tile.order[seqNo];
                unsigned childID = tile.children[childNo];
                Matrix m = matrix.concat( M[childNo] );
                int rank = seqNo * digit + tileRank;
                stack.push_back({childID, depth + 1, m, rank});
            }
        }
        else {
            int nPoints = (maxDensity - tileRank - 1) / digit;
            for (int pointNo = 0; pointNo <= nPoints; pointNo++) {
                counter++;
                Point p = transform(tile.p[nPoints][pointNo], matrix);
                int rank = pointNo * digit + tileRank;
                int density = img.getDensity(p);
                if (density > rank + 1) {
                    samples.push_back(p);
                }
            }
        }
    } while(!stack.empty());
    fprintf(stderr, "Made %5d hits out of %5d generated samples\n", int(samples.size()), counter);
    fprintf(stderr, "Reached depth = %3d\n", reachedDepth);
    return samples;
}

void plotEPS(std::ostream &file, Points &samples) {
    for (int i = 0; i < samples.size(); i++) {
        file << samples[i].x << " " << samples[i].y << " p\n";
    }
    file << "showpage\n";
}

const char *USAGE_MESSAGE = "Usage: %s <pgmfile> <tileSet>\n"
    "Options:\n"
    "-s <value>                     Scales color to add more samples\n"
    "-i <number>                    id of root tile\n"
    "-M <max subdivision depth>\n"
    "-m <min subdivision depth>\n"
    "-t                             plot tiles\n"
;

int main(int argc,char **argv) {
    srand(time(NULL));
    int opt;                                                                    // For use by getopt, the command line options utility
    int rootID = -1;
    double colorScale = 1.0;
    while ((opt = getopt(argc, argv, "s:i:m:M:t")) != -1) {                       // Modify default settings with command line options
        switch (opt) {
            case 's': colorScale = atof(optarg); break;
            case 'i': rootID = atoi(optarg); break;
            case 'M': maxDepth = atoi(optarg); break;
            case 'm': minDepth = atoi(optarg); break;
            case 't': showtiles = true; break;
            default: fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
        }
    }
    if (optind > argc - 2) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]);
        exit(1);
    }
    char *pgmFileName = argv[optind];
    char *tilesFileName = argv[optind + 1];
    CImage img(pgmFileName, colorScale);
    std::fstream file(tilesFileName);
    loadTileSet(file);
    file.close();
    if (rootID < 0) rootID = rand();
    rootID %= tileCount;
    clock_t t0 = clock();
    Points samples = sample(img, rootID);
    clock_t t1 = clock();
    double totalTime = (double)(t1 - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "generated %'d samples in %.6fs (%'d points per second)..\n",
            int(samples.size()), totalTime, (int)(samples.size()/totalTime));
    plotEPS(std::cout, samples);
}
