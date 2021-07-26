/*
 * Generates a sequence of blue noise samples similar to LD sequences
 * Abdalla Ahmed
 * 2016-12-01
 */

#include "self-similar.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <algorithm>

inline Point getSample(unsigned tileID, unsigned sampleNo, Matrix m = IDENTITY) {
    while (sampleNo > 0) {
        unsigned childNo = tiles[tileID].order[ sampleNo & mask ];
        tileID = tiles[tileID].children[childNo];
        m = m.concat( M[childNo] );
        sampleNo = (sampleNo >> shift);
    }
    return transform(tiles[tileID].p, m);
}

void generate(Points &samples) {
    int n = samples.size();
    int tileID = rand() % tileCount;
    for (int i = 0; i < n; i++) {
        samples[i] = getSample(tileID, i);
    }
}

const char *USAGE_MESSAGE = "Usage: %s <data> <n>\n";

int main(int argc,char **argv) {
    srand(time(NULL));
    if (argc != 3) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]);
        exit(1);
    }
    char *fileName = argv[1];
    std::fstream file(fileName);
    loadTileSet(file);
    int n = atoi(argv[2]);
    Points s(n);
    clock_t t0 = clock();
    generate(s);
    clock_t t1 = clock();
    double totalTime = (double)(t1 - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "generated %'d samples in %.6fs (%'d points per second)..\n",
            n, totalTime, (int)(n/totalTime));
    printf("%d\n", n);
    for (int i = 0; i < n; i++) printf("%0.16f %0.16f\n", s[i].x, s[i].y);
}
