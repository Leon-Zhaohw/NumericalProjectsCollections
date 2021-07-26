/*
 * A common place for CGAL geometry-related utils
 * Abdalla Ahmed
 * 2015-10-05
 */

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#define LL(x) ((x).squared_length())
#define VL(x) (sqrt((x).squared_length()))

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Point_2                                                  Point;
typedef K::Vector_2                                                 Vector;
struct FInfo {                                                                  // Information stored in faces in DT.
    Point c;                                                                    // Circumcenter, to avoid repeated calculation
};
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Tvb;
typedef CGAL::Triangulation_face_base_with_info_2<FInfo, K>         Tfb;
typedef CGAL::Triangulation_data_structure_2<Tvb,Tfb>               Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      DT;
typedef DT::Vertex_handle                                           VH;
typedef DT::Face_handle                                             FH;
typedef DT::Vertex_circulator                                       VC;
typedef DT::Face_circulator                                         FC;

typedef VH                                                          VH9[9];
typedef std::vector<DT>                                             DTS;
typedef std::vector<Point>                                          Points;
typedef std::vector<Vector>                                         Vectors;

const double dhex = 1.07456993182354;                                           // Point spacing in triangular grid when average capacity is 1; = sqrt(2/(sqrt(3)))


inline double crossProduct(const Vector &v1, const Vector &v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}

inline double crossProduct(const Point &ref, const Point &p1, const Point &p2) {
    return (p1.x() - ref.x()) * (p2.y() - ref.y())
    - (p1.y() - ref.y()) * (p2.x() - ref.x());
}

inline bool isConvex(Point &p1, Point &p2, Point &p3, Point &p4) {
    bool result = true;
    if (
        (crossProduct(p1, p2, p4) < 0) ||
        (crossProduct(p2, p3, p1) < 0) ||
        (crossProduct(p3, p4, p2) < 0) ||
        (crossProduct(p4, p1, p3) < 0)
    ) {
        return false;
    }
    return true;
}

inline Point intersect(
    const Point &p1,
    const Point &p2,
    const Point &p3,
    const Point &p4
) {                                                                             // See https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    const double &x1 = p1.x(), &x2 = p2.x(), &x3 = p3.x(), &x4 = p4.x();
    const double &y1 = p1.y(), &y2 = p2.y(), &y3 = p3.y(), &y4 = p4.y();
    double x12 = x1 - x2, y12 = y1 - y2, x34 = x3 - x4, y34 = y3 - y4;
    double xy12 = x1 * y2 - y1 * x2, xy34 = x3 * y4 - y3 * x4;
    double det = 1.0 / (x12 * y34 - y12 * x34);
    return {det * (xy12 * x34 - xy34 * x12), det * (xy12 * y34 - xy34 * y12)};
}

double area(const Points &p, bool lastEqFirst = false) {
    int n = p.size();
    if (n == 0) return 0;
    double a = 0;
    for (int i = 0; i < n-1; i++) {
        a += p[i].x() * p[i+1].y() - p[i].y() * p[i+1].x();
    }
    if (!lastEqFirst) {
        a += p[n-1].x() * p[0].y() - p[n-1].y() * p[0].x();
    }
    return a / 2.0;
}

Points VoronoiCell(const DT &dt, const VH &vh) {
    Points result, empty;
    FC fc = dt.incident_faces(vh), done(fc);
    do {
        if (dt.is_infinite(fc)) { return empty; }
        result.push_back(dt.circumcenter(fc));
    } while (++fc != done);
    return result;
}

Vector centroid(const DT &dt, const VH &vh) {                                   // Centroid relative to vertex See http://en.wikipedia.org/wiki/Centroid
    double a(0), cx(0), cy(0);                                                  // Cell area (actually twice the area) and centroid coordinates
    double XProduct;                                                            // Cross product of vectors to adjacent vertices
    FC fc = dt.incident_faces(vh), done(fc);
    do {
        Point p1 = dt.circumcenter(fc);
        Point p2 = dt.circumcenter(++fc);
        XProduct = p1.x() * p2.y() - p1.y() * p2.x();
        a += XProduct;                                                          // Accumulate areas
        cx += (p1.x() + p2.x()) * XProduct;
        cy += (p1.y() + p2.y()) * XProduct;
    } while (fc != done);
    cx /= 3.0 * a;
    cy /= 3.0 * a;
    return Point(cx, cy) - vh->point();                                         // Return shift from current position to centroid
};

Vector qcentroid(const DT &dt, const VH &vh) {                                  // Centroid of the neighbors qhull relative to vertex
    double a(0), cx(0), cy(0);                                                  // Cell area (actually twice the area) and centroid coordinates
    double XProduct;                                                            // Cross product of vectors to adjacent vertices
    VC vc = dt.incident_vertices(vh), done(vc);
    do {
        Point p1 = (vc++)->point();
        Point p2 = vc->point();
        XProduct = p1.x() * p2.y() - p1.y() * p2.x();
        a += XProduct;                                                          // Accumulate areas
        cx += (p1.x() + p2.x()) * XProduct;
        cy += (p1.y() + p2.y()) * XProduct;
    } while (vc != done);
    cx /= 3.0 * a;
    cy /= 3.0 * a;
    return Point(cx, cy) - vh->point();                                         // Return shift from current position to centroid
};


inline bool isInRect(const Point &p, const Point &BL, const Point &TR) {        // Check if a point is inside a given rectangle
    return BL.x() <= p.x() && BL.y() <= p.y() &&
    p.x() < TR.x() && p.y() < TR.y();
}

inline double area(Point &p1, Point &p2, Point &p3) {
    return 0.5 * fabs(crossProduct(p1, p2, p3));
}

inline Point midPoint(Point &p1, Point &p2, double ratio = 0.5) {               // A point (halfway) between two.
    return p1 + ratio * (p2 - p1);
}

inline Vector ortho(Point &p1, Point &p2) {                                     // Orthogonal vector, if p1 is origin and p2 is direction x then give dir y.
    return {p1.y() - p2.y(), p2.x() - p1.x()};
}

inline Vector unit(Vector &v) {
    return v / sqrt(v.squared_length());
}

inline Point circumcenter(DT &dt, VH &vh, FC &fc) {
    if (dt.is_infinite(fc)) {
        VC vc = dt.incident_vertices(vh, fc);
        double sgn = 1.0;
        if (dt.is_infinite(vc)) {
            vc++;
            sgn = -1.0;
        }
        Point m = midPoint(vh->point(), vc->point());
        Vector dir = ortho(vh->point(), vc->point());
        return m + 1e6 * sgn * unit(dir);
    }
    else return dt.circumcenter(fc);
}

inline double triangleType(Point &p1, Point &p2, Point &p3) {
    double sq1 = (p3 - p2).squared_length();
    double sq2 = (p1 - p3).squared_length();
    double sq3 = (p2 - p1).squared_length();
    if (sq1 < sq2) std::swap(sq1, sq2);
    if (sq1 < sq3) std::swap(sq1, sq3);
    return sq1 - (sq2 + sq3);                                                   // 0 for right, < 0 for acute, > 0 for obtuse
}

inline double triangleType(FC &fc) {
    return triangleType(
        fc->vertex(0)->point(),
        fc->vertex(1)->point(),
        fc->vertex(2)->point()
    );
}

void clamp(Vector &v, double max) {
    double l = VL(v);
    if (l > max) {
        v = max * (v / l);
    }
}

void clamp2sq(Vector &v, double halfWidth) {
    double largest = std::max(fabs(v.x()), fabs(v.y()));
    if (largest > halfWidth) {
        v = (halfWidth/largest) * v;
    }
}

// See https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
inline bool isInsideEdge(const Point &p, const Point &edgeStart, const Point &edgeEnd) {
    return ( crossProduct(p, edgeStart, edgeEnd) > 0 );
}

Points clipByEdge(Points &p, Point &edgeStart, Point &edgeEnd) {
    Points result;
    int startIndex = p.size() - 1;
    for (int endIndex = 0; endIndex < p.size(); endIndex++) {
        Point start = p[startIndex];
        Point end = p[endIndex];
        if (isInsideEdge(end, edgeStart, edgeEnd)) {
            if (!isInsideEdge(start, edgeStart, edgeEnd)) {
                result.push_back(intersect(start, end, edgeStart, edgeEnd));
            }
            result.push_back(end);
        }
        else if (isInsideEdge(start, edgeStart, edgeEnd)) {
            result.push_back(intersect(start, end, edgeStart, edgeEnd));
        }
        startIndex = endIndex;
    }
    return result;
}

Points crop(const Points &subject, const Points &boundry) {                     // See https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
    Points result = subject;
    Point edgeStart = boundry.back();
    for (int i = 0; i < boundry.size(); i++) {
        Point edgeEnd = boundry[i];
        result = clipByEdge(result, edgeStart, edgeEnd);
        edgeStart = edgeEnd;
    }
    return result;
}


Vector clipDisplacementByEdge(
    const Point &origin,
    const Vector &displacement,
    const Point &edgeStart,
    const Point &edgeEnd
) {
    if (!isInsideEdge(origin, edgeStart, edgeEnd)) {
        return {0, 0};
    }
    Point endPoint = origin + displacement;
    if (!isInsideEdge(origin, edgeStart, edgeEnd)) {
        return displacement;
    }
    endPoint = intersect(origin, endPoint, edgeStart, edgeEnd);
    return endPoint - origin;
}

Vector clampDisplacementByPolygon(
    const Point &origin,
    const Vector &displacement,
    Points polygon
) {
    Vector result = displacement;
    int startIndex = polygon.size() - 1;
    for (int endIndex = 0; endIndex < polygon.size(); endIndex++) {
        Point start = polygon[startIndex];
        Point end = polygon[endIndex];
        result = clipDisplacementByEdge(origin, result, start, end);
        startIndex = endIndex;
    }
    return result;
}

double cvtEnergy(
    const Point &p,                                                             // A point in the set
    const Point &v1,                                                            // A Voronoi vertex
    const Point &v2,                                                            // The following Voronoi vertex in CCW
    const Point &m                                                              // The projection point of height, = mid point to neighbor
) {
    Vector pv1 = v1 - p, pv2 = v2 - p, pm = m - p;
    double A1 = crossProduct(pv1, pm);                                          // Could be negative, but it's appropriate
    double A2 = crossProduct(pm, pv2);                                          // Also becomes negative when appropriate
    A1 *= pv1.squared_length() + 2 * pm.squared_length();
    A2 *= pv2.squared_length() + 2 * pm.squared_length();
    return (A1 + A2) / 24;
}

struct Matrix {                                                                 // Postscript-style matrix, but only with scale and translate
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
    double x = m.scale * p.x() + m.tx;
    double y = m.scale * p.y() + m.ty;
    return Point(x, y);
}



#endif
