#ifndef _R3_
#define _R3_ 1

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592308
#endif

#ifndef FEQ_EPS
#define FEQ_EPS 1e-12
#endif

namespace TML
{
  R3 RotParam(const R3& from, const R3& to, double& angle);
  R3 Rotate(const R3& p, const R3& axis, double angle);
  
  class R3
    {
    protected:
      inline void copy(const R3& v) { x = v.x; y = v.y; z = v.z; }
      
    public:
      // Attributes
      double x, y, z;

      //Constructors
      R3()                             { set(); }
      R3(const R3& v)                  { copy(v); }
      R3(double a, double b, double c) { set(a,b,c); }
      ~R3() { }
      
      //Set
      inline void set()                             { x = y = z = 0.; }
      inline void set(double a, double b, double c) { x = a; y = b; z = c; } 
      inline void set(const R3& v)                  { copy(v); }

      // Length
      inline double length2() { return (x*x + y*y + z*z); }
      inline double length()  { return std::sqrt(length2()); }
      void  normalize()
	{ 
	  double l = length(); 
	  if (l<FEQ_EPS) throw Error("TML::R3::normalize: length is zero");
	  x/=l; y/=l; z/=l;
	}
      
      // Arithmetic Operators
      bool operator==(const R3& v) const
	{
	  double dx=x-v.x; double dy=y-v.y; double dz=z-v.z;
	  return ( (dx*dx + dy*dy + dz*dz) < FEQ_EPS);
	}
      bool operator!=(const R3& v) const
	{
	  double dx=x-v.x; double dy=y-v.y; double dz=z-v.z;
	  return ( (dx*dx + dy*dy + dz*dz) >= FEQ_EPS);
	}
      
      R3& operator=(const R3& v)  { copy(v); return *this; }
      R3& operator+=(const R3& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
      R3& operator-=(const R3& v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
      R3& operator*=(double k)    { x*=k; y*=k; z*=k; return *this; }
      R3& operator/=(double k)    { x/=k; y/=k; z/=k; return *this; }
      
      R3 operator+(const R3& v) const { return R3(x+v.x, y+v.y, z+v.z); }
      R3 operator-(const R3& v) const { return R3(x-v.x, y-v.y, z-v.z); }
      R3 operator*(double k) const    { return R3(k*x, k*y, k*z); }
      R3 operator/(double k) const    { return R3(x/k, y/k, z/k); }
      R3 operator-() const            { return R3(-x,-y,-z); }
    };
  
  inline R3 operator*(double k, const R3 &vi)
    { return R3( k*vi.x, k*vi.y, k*vi.z ); }

  // Vetorial operators
  inline double dot(const R3& a, const R3& b)
    { return ( a.x*b.x + a.y*b.y + a.z*b.z ); }
  inline R3 cross(const R3& a, const R3& b)
    { 
      return R3( a.y*b.z - a.z*b.y,
		 a.z*b.x - a.x*b.z,
		 a.x*b.y - a.y*b.x ); 
    }
  inline Matrix extprod(const R3& a, const R3& b)
    { 
      Matrix e(3,3);
      e[0][0] = a.x*b.x;  e[0][1] = a.x*b.y; e[0][2] = a.x*b.z;
      e[1][0] = a.y*b.x;  e[1][1] = a.y*b.y; e[1][2] = a.y*b.z;
      e[2][0] = a.z*b.x;  e[2][1] = a.z*b.y; e[2][2] = a.z*b.z;
      return Matrix(e);
    }

  // IO
  inline ostream& operator<<(ostream& out, const R3& v)
    { return out << v.x << " " << v.y << " " << v.z; }
  inline istream& operator>>(istream& in, R3& v)
    { return in >> v.x >> v.y >> v.z; }

  // |a-b| ~ 0
  inline bool zero(double a, double b, double eps=FEQ_EPS)
    { return (std::fabs(a-b) < eps); }
}

#endif
