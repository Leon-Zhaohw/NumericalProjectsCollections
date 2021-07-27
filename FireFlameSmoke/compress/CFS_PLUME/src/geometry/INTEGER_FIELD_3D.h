#ifndef INTEGER_FIELD_3D_H
#define INTEGER_FIELD_3D_H

#include "EIGEN.h"
#include "VECTOR.h"
#include "FIELD_3D.h"

using std::vector;

class INTEGER_FIELD_3D {
  public:
    INTEGER_FIELD_3D(); 
    INTEGER_FIELD_3D(const int& xRes, const int& yRes, const int& zRes);
    INTEGER_FIELD_3D(const int* data, const int& xRes, const int& yRes, const int& zRes);
    INTEGER_FIELD_3D(const INTEGER_FIELD_3D& m);
    ~INTEGER_FIELD_3D();
    
    // Accessors
    inline int& operator()(int x, int y, int z) { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _xRes*_yRes  + y * _xRes + x]; 
    } 
    const int operator()(int x, int y, int z) const { 
      assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
      return _data[z * _xRes*_yRes  + y * _xRes + x]; 
    }
    const int* dataConst() const { return _data; }

    inline int& operator[](int x) { return _data[x]; }
    const int operator[](int x) const { return _data[x]; }
    const int xRes() const { return _xRes; }
    const int yRes() const { return _yRes; }
    const int zRes() const { return _zRes; }
    const int totalCells() const { return _totalCells; }
    int* data() { return _data; };
    
    // Return a flattened VECTOR of the array contents
    VECTOR flattened() const;
    VECTOR flattenedRow() const;
    
    // set all entries to zero
    void clear();
    void clear(const vector<int>& nonZeros);
    void clear(const vector<int>& nonZeros, const int size);
    
    // sum of the squares of the entries
    int energySum();

    // Resize to the desired dimensions and set entries to zero
    void resizeAndWipe(int xRes, int yRes, int zRes);

    // Overloaded operators
    
    INTEGER_FIELD_3D& operator*=(const double& alpha);
    INTEGER_FIELD_3D& operator=(const INTEGER_FIELD_3D& A);

    // Get a subfield
    INTEGER_FIELD_3D subfield(const int xBegin, const int xEnd, 
                    const int yBegin, const int yEnd, 
                    const int zBegin, const int zEnd) const;
    // Swap contents with another 3d integer field
    void swapPointers (INTEGER_FIELD_3D& field);

  private:
    int _xRes;
    int _yRes;
    int _zRes;
    int _totalCells;
    int* _data;


};


#endif
