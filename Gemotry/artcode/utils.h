/*
 * 2016-02-26: Created by Abdalla Ahmed
 */

#ifndef UTILS_H
#define UTILS_H
#include <algorithm>
#include <vector>
#include <math.h>

typedef std::vector<int>                                        TList;
typedef std::vector<TList>                                      Lists;

std::vector<unsigned> shuffle(const unsigned N) {                               // Return a randomly ordered list; adapted from FPO code
    std::vector<unsigned> list(N);
    for (unsigned i = 0; i < N; i++) list[i] = i;
    for (unsigned i = 0; i < N-1; i++) {
        unsigned r = i + rand() % (N - i);
        std::swap(list[i], list[r]);
    }
    return list;
}

void shuffle(std::vector<unsigned> &list) {                                     // make a randomly ordered list.
    unsigned N = list.size();
    for (unsigned i = 0; i < N; i++) list[i] = i;
    for (unsigned i = 0; i < N-1; i++) {
        unsigned r = i + rand() % (N - i);
        std::swap(list[i], list[r]);
    }
}

int find(double value, const std::vector<double> &th) {                         // Locate a value in a list of intervals. The list is assumed to be sorted.
    int result = -1;
    for (unsigned i = 1; i < th.size(); i++) {
        if (th[i] > value) {
            result = i-1;
            break;
        }
    }
    return result;
}

struct TCompare {
    double *p;
    bool operator() (int i1, int i2) { return (p[i1] < p[i2]); }
};

std::vector<unsigned> sortIndexes(std::vector<double> &x) {
    int count = x.size();
    TCompare cmp;
    cmp.p = x.data();
    std::vector<unsigned> l(count);
    for (int i = 0; i < count; i++) l[i] = i;
    std::sort (l.begin(), l.end(), cmp);
    return l;
}
void sortIndexes(std::vector<double> &x, std::vector<unsigned> &order) {
    int count = x.size();
    TCompare cmp;
    cmp.p = x.data();
    for (int i = 0; i < count; i++) order[i] = i;
    std::sort(order.begin(), order.end(), cmp);
}



inline double rnd(const double &size) {
    return drand48() * size;
}

inline double jtr(const double &size) {
    return (drand48() - 0.5) * size;
}

inline double frac(const double &x) { return x - floor(x); }

#endif

