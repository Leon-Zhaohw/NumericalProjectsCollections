#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <iostream>
#include <fftw3.h>

#include "EIGEN.h"
#include "VECTOR3_FIELD_3D.h"
#include "INTEGER_FIELD_3D.h"
#include "MATRIX_COMPRESSION_DATA.h"

////////////////////////////////////////////////////////
// A custom class that tracks non-zero matrix entries
// but avoids resizing
////////////////////////////////////////////////////////
class NONZERO_ENTRIES {
  public:
    NONZERO_ENTRIES() :
      _data(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE),
      _currentEntry(-1),
      _maxEntries(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE) {};
    ~NONZERO_ENTRIES() {};
    void clear() {
      //TIMER functionTimer("NONZERO_ENTRIES clear");
      _currentEntry = -1;
    };
    void clearAll() {
      _currentEntry = -1;
      for (unsigned int x = 0; x < _data.size(); x++)
        _data[x] = 0;
    };
    const int size() const { return _currentEntry + 1; };
    const vector<int>& data() const { return _data; };
    inline int& operator[](int index) {
      assert(index >= 0);
      assert(index < _currentEntry + 1);
      return _data[index];
    };
    inline const int operator[](int index) const {
      assert(index >= 0);
      assert(index < _currentEntry + 1);
      return _data[index];
    }
    void push_back(const int entry) {
      _currentEntry++;
      _data[_currentEntry] = entry;
      assert(_currentEntry < _maxEntries);
    };

  private:
    vector<int> _data;
    int _currentEntry;
    int _maxEntries;
};

////////////////////////////////////////////////////////
// Function signatures
////////////////////////////////////////////////////////

// round a FIELD_3D to an INTEGER_FIELD_3D
void RoundFieldToInt(const FIELD_3D& F, INTEGER_FIELD_3D* castedField);

// cast an INTEGER_FIELD_3D to a FIELD_3D
void CastIntFieldToDouble(const INTEGER_FIELD_3D& F, FIELD_3D* castedField);

// TK: Adding a version that does it on the raw pointer
void CastIntFieldToDouble(const INTEGER_FIELD_3D& F, const int totalCells, double* castedField);
void CastIntFieldToDouble(const INTEGER_FIELD_3D& F, const int totalCells, double* castedField, const vector<int>& nonZeros);
void CastIntFieldToDouble(const INTEGER_FIELD_3D& F, const int totalCells, double* castedField, const NONZERO_ENTRIES& nonZeros);

// form the cumulative sum starting at zero of a passed in integer vector
void ModifiedCumSum(const VectorXi& V, VectorXi* sum);

// build the block indices matrix from the block lengths matrix using
// modified cum sum
void BuildBlockIndicesMatrix(COMPRESSION_DATA* data);

// builds the block indices matrix from the block lengths matrix.
// uses explicitly passed in matrices for debugging purposes!
void BuildBlockIndicesMatrixDebug(const MatrixXi blockLengths, MatrixXi* blockIndices);

// extract the three scalar fields from a vector field
void GetScalarFields(const VECTOR3_FIELD_3D& V, FIELD_3D* X, FIELD_3D* Y, FIELD_3D* Z);

// make an fftw 3d dct plan. direction 1 is forward, -1 is inverse

void Create_DCT_Plan(double* in, int direction, fftw_plan* plan);

// perform a unitary normalization for the forward 3d dct
void DCT_Unitary_Normalize(double* buffer);

// perform a unitary normalization in preparation for the inverse 3d dct
void UndoNormalize(FIELD_3D* F);

// perform a unitary normalization on a flattened-out FIELD_3D in preparation
// for the inverse 3d dct
void UndoNormalizeEigen(VectorXd* F);

// do the corresponding dct based on the plan and direction
void DCT_Smart_Unitary(const fftw_plan& plan, int direction, double* in, FIELD_3D* F);


// do the corresponding dct on a flattened-out FIELD_3D based on the plan and direction
void DCT_Smart_Unitary_Eigen(const fftw_plan& plan, int direction, double* in, VectorXd* F);

// given passed in dimensions, compute the amount of padding required to get to
// dimensions that are evenly divisible by BLOCK_SIZE
void GetPaddings(const VEC3I& v, VEC3I* paddings);

// given passed in field, parse it into a vector of padded 3d blocks of size BLOCK_SIZE
// in row-major order
void GetBlocks(const FIELD_3D& F, vector<FIELD_3D>* blocks);

// given a passed in FIELD_3D, pad it with zeros and  parse it
// into a vector of flattened 8 x 8 x 8 blocks (listed in row-major order
void GetBlocksEigen(const FIELD_3D& F, vector<VectorXd>* blocks);

// reassemble a big field from a row-major list of smaller blocks
void AssimilateBlocks(const VEC3I& dims, const vector<FIELD_3D>& V, FIELD_3D* assimilatedField);

// reconstruct a FIELD_3D with the passed in dims
// from a list of 8 x 8 x 8 flattened blocks
void AssimilateBlocksEigen(const VEC3I& dims, vector<VectorXd>* blocks,
    FIELD_3D* assimilatedField);

// perform a unitary dct on each block of a passed in list. direction 1 is dct,
// -1 is idct
void UnitaryBlockDCT(int direction, vector<FIELD_3D>* blocks);

// perform a unitary dct on each block of a passed in list of flattened-out fields.
// direction 1 is dct, -1 is idct.
void UnitaryBlockDCTEigen(int direction, vector<VectorXd>* blocks);

// build a block diagonal matrix with A's as the kernel for
// a 'count' number of times. inefficient usage of memory
// since it fails to use sparse matrix
void BlockDiagonal(const MatrixXd& A, int count, MatrixXd* B);

// build a sparse block diagonal matrix with A's as the kernel
// for a 'count' number of times
void SparseBlockDiagonal(const MatrixXd& A, int count, SparseMatrix<double>* B);

// given a passed in vec3 field, build a matrix with each column
// corresponding to the x-, y-, and z-components
void BuildXYZMatrix(const VECTOR3_FIELD_3D& V, MatrixXd* A);

// build a transformed vector field using a coordinate transform computed from the
// svd decomposition of the original x-, y-, and z- coordinates
// uses v^T, not v!
void TransformVectorFieldSVD(VectorXd* s, MatrixXd* v, VECTOR3_FIELD_3D* transformedV);

// build a transformed vector field using a coordinate transform computed from
// svd decomposition of the original x-, y-, and z- coordinates. update the
// compression data to store the v matrix and the singular values.
// can be called repeatedly in chain-like fashion.
void TransformVectorFieldSVDCompression(VECTOR3_FIELD_3D* V, COMPRESSION_DATA* data);

// perform the coordinate transform without having
// to recompute the SVD by reading in from a cached 3 x 3
// v matrix. used for the projection trick later on.
void TransformVectorFieldSVDCached(Matrix3d* v, VECTOR3_FIELD_3D* V);

// undo the effects of a previous svd coordinate transform on a vector field
void UntransformVectorFieldSVD(const MatrixXd& v, VECTOR3_FIELD_3D* transformedV);

// normalize the block to a resolution of nBits based on the DC component.
// update the sList.
void PreprocessBlock(FIELD_3D* F, int blockNumber, int col, COMPRESSION_DATA* data);

// Binary search to find the appropriate gamma given
// desired percent threshold within maxIterations. Prints
// out information as it goes.
///////////////////////////////////////////////////////
void TuneGammaVerbose(const FIELD_3D& F, int blockNumber, int col,
    COMPRESSION_DATA* data, FIELD_3D* damp);

// do a binary search to find the appropriate gamma given the desired percent
// energy accuracy and max iterations. the variable damp will be rewritten to the
// desired damping array. updates gamaList.
void TuneGamma(const FIELD_3D& F, int blockNumber, int col, COMPRESSION_DATA* data,
    FIELD_3D* damp);

// version that uses fastPow instead of the slower pow
void TuneGammaFastPow(const FIELD_3D& F, int blockNumber, int col, COMPRESSION_DATA* data,
    FIELD_3D* damp);

// quantizes the gamma values to quarter-integers
void TuneGammaQuantized(const FIELD_3D& F, int blockNumber, int col, COMPRESSION_DATA* data,
    FIELD_3D* damp);

// simply sets gamma equal to zero for no damping. for
// debug purposes only.
void TuneGammaDebug(const FIELD_3D& F, int blockNumber, int col,
    COMPRESSION_DATA* data, FIELD_3D* damp);

// takes a passed in FIELD_3D (which is intended to be
// the result of a DCT post-preprocess). calculates the best gamma value
// for a damping array. then damps by that array and quantizes the result to an integer.
// stores the value of gamma for the damping.
void EncodeBlock(const FIELD_3D& F, int blockNumber, int col, COMPRESSION_DATA* data,
    INTEGER_FIELD_3D* quantized);

// takes a passed in FIELD_3D (which is intended to be
// the result of a DCT post-preprocess). for debug purposes
// only, does not perform any damping---gamma is set to 0.
void EncodeBlockDebug(const FIELD_3D& F, int blockNumber, int col, COMPRESSION_DATA* data,
    INTEGER_FIELD_3D* quantized);

// takes a passed in INTEGER_FIELD_3D (which is intended to be run-length
// decoded and unzigzagged) at a particular blockNumber and column of the matrix.
// undoes the effects of damping/quantization as best as it can.
/*
void DecodeBlock(const INTEGER_FIELD_3D& intBlock, int blockNumber, int col,
    const DECOMPRESSION_DATA& decompression_data, FIELD_3D* decoded);
*/

// performs the same operations as DecodeBlock, but with passed in compression data
// rather than passed in decompression data. due to const poisoning, compression
// data cannot be marked const, but is treated as such.
//void DecodeBlockWithCompressionData(const INTEGER_FIELD_3D& intBlock,
//  int blockNumber, int col, COMPRESSION_DATA* data, FIELD_3D* decoded);

// TK: This actually only needs the raw pointer for "decoded", and then we avoid an
// expensive copy at the end.
void DecodeBlockWithCompressionData(const INTEGER_FIELD_3D& intBlock,
  int blockNumber, int col, COMPRESSION_DATA* data, Real* decoded);
void DecodeBlockWithCompressionDataSparse(const INTEGER_FIELD_3D& intBlock,
  int blockNumber, int col, COMPRESSION_DATA* data, Real* decoded, const NONZERO_ENTRIES& nonZeros);
void DecodeBlockWithCompressionDataSparseStackless();
void DecodeBlockWithCompressionDataSparseQuantized(const INTEGER_FIELD_3D& intBlock,
  int blockNumber, int col, COMPRESSION_DATA* data, Real* decoded, const NONZERO_ENTRIES& nonZeros);

// flattens an INTEGER_FIELD_3D through a zig-zag scan
// into a VectorXi. Since the scan always follows the same order,
// we precompute the zigzag scan array, pass it
// as a parameter, and then just do an index lookup
void ZigzagFlatten(const INTEGER_FIELD_3D& F, const INTEGER_FIELD_3D& zigzagArray,
    VectorXi* zigzagged);

// unflattens a VectorXi into an INTEGER_FIELD_3D
void ZigzagUnflatten(const VectorXi& V, const INTEGER_FIELD_3D& zigzagArray,
    INTEGER_FIELD_3D* unflattened);

// helper function to run-length encoding to find runs
int FindRun(const VectorXi& V, int i);

// given a zigzagged integer vector, write it to a binary
// file via run-length encoding. updates the blockLengths matrix
void RunLengthEncodeBinary(const char* filename, int blockNumber, int col,
    const VectorXi& zigzaggedArray, COMPRESSION_DATA* compression_data);


// decode a run-length encoded binary file and fill
// a VectorXi with the contents.
void RunLengthDecodeBinary(int* allData, int blockNumber, int col,
    COMPRESSION_DATA* compression_data, VectorXi* parsedData);

// TK: do a version that does the zigzag unflatten at the same time
void RunLengthDecodeBinaryInPlace(int* allData, int blockNumber, int col,
    const INTEGER_FIELD_3D& reverseZigzag,
    COMPRESSION_DATA* compression_data,
    INTEGER_FIELD_3D& parsedDataField);
void RunLengthDecodeBinaryInPlaceSparse(int* allData, int blockNumber, int col,
    const INTEGER_FIELD_3D& reverseZigzag,
    COMPRESSION_DATA* compression_data,
    INTEGER_FIELD_3D& parsedDataField, NONZERO_ENTRIES& nonZeros);
void RunLengthDecodeBinaryInPlaceSparseStackless();

// takes an input FIELD_3D which is the result of
// an SVD coordinate transform, compresses it according
// to the general scheme, and writes it to a binary file
void CompressAndWriteField(const char* filename, const FIELD_3D& F, int col,
    COMPRESSION_DATA* compression_data);

// takes an input FIELD_3D at a particular matrix column
// which is the result of an SVD coordinate transform, compresses
// it according to the general scheme, and writes it to a binary file.
// meant to be called in a chain so that the binary file
// continues to grow. for debugging, gamma is set to zero everywhere!
void CompressAndWriteFieldDebug(const char* filename, const FIELD_3D& F, int col,
    COMPRESSION_DATA* compression_data);

// print four different percents for how far along each column we are
void PrintProgress(int col, int numCols);

// generate the header information for the encoded binary file
void WriteMetaData(const char* filename, const COMPRESSION_DATA& compression_data);

// write the singular values and V matrices to a binary file
void WriteSVDData(const char* filename, COMPRESSION_DATA* data);

// delete a binary file if it already exists
void DeleteIfExists(const char* filename);

// compress all of the scalar field components
// of a matrix (which represents a vector field) and write them to
// a binary file. applies svd coordinate transform first
void CompressAndWriteMatrixComponents(const char* filename, const MatrixXd& U,
      COMPRESSION_DATA* data0, COMPRESSION_DATA* data1, COMPRESSION_DATA* data2);

// compress all of the scalar field components
// of a matrix (which represents a vector field) and write them to
// a binary file. applies svd coordinate transform first.
// uses gamma as zero everywhere for debugging!
void CompressAndWriteMatrixComponentsDebug(const char* filename, const MatrixXd& U,
      COMPRESSION_DATA* data0, COMPRESSION_DATA* data1, COMPRESSION_DATA* data2);

// reads from a binary file of the SVD data and sets the initializations
// inside compression data
void ReadSVDData(const char* filename, COMPRESSION_DATA* data);

// reads from a binary file into a buffer and sets initializations
// inside compression data
int* ReadBinaryFileToMemory(const char* filename,
    COMPRESSION_DATA* data);

// decode a particular scalar component of the corresponding column of the
// matrix
void DecodeScalarField(COMPRESSION_DATA* compression_data, int* allData,
    int col, FIELD_3D* decoded);

// decode an entire scalar field of a particular column from matrix compression data
// *without* going back to the spatial domain (or the SVD transform). leave them
// in a list of blocks as well
void DecodeScalarFieldEigen(COMPRESSION_DATA* compression_data, int* allData,
    int col, vector<VectorXd>* decoded);

// try doing it all sparsely, i.e. skipping all the zeros
void DecodeScalarFieldEigenSparse(COMPRESSION_DATA* compression_data, int* allData,
    int col, vector<VectorXd>* decoded);
void DecodeScalarFieldEigenSparse(COMPRESSION_DATA* compression_data, int* allData,
    int col, vector<VectorXd>* decoded, vector<NONZERO_ENTRIES>& nonZeros);
void DecodeScalarFieldEigenSparseStackless(COMPRESSION_DATA* compression_data, int* allData,
    int col, vector<VectorXd>* decoded, vector<NONZERO_ENTRIES>& nonZeros);

// uses DecodeScalarField three times to form a vector field, and then undoes
// the SVD transform to recover the original vector field corresponding to col
void DecodeVectorField(MATRIX_COMPRESSION_DATA* data, int col,
    VECTOR3_FIELD_3D* decoded);

// uses DecodeVectorField on each column to recover a lossy entire matrix
void DecodeMatrix(MATRIX_COMPRESSION_DATA* data, MatrixXd* decoded);

// compute the block-wise dot product between two lists and sum them into one
// large dot product
double GetDotProductSum(const vector<VectorXd>& Vlist, const vector<VectorXd>& Wlist);

// helper function for frequency domain projection. transforms
// V first by the SVD and then DCT. fills up the three vectors
// with the three components.
void TransformSVDAndDCT(int col, const VECTOR3_FIELD_3D& V,
    MATRIX_COMPRESSION_DATA* U_data,
    vector<VectorXd>* Xpart, vector<VectorXd>* Ypart, vector<VectorXd>* Zpart);

// helper function for frequency domain projection. transforms
// V by the DCT. fills up the three vectors
// with the three components.
void TransformDCT(const VECTOR3_FIELD_3D& V,
    MATRIX_COMPRESSION_DATA* U_data, vector<VectorXd>* Xpart, vector<VectorXd>* Ypart,
    vector<VectorXd>* Zpart);

// do peeled compressed projection naively in the regular spatial domain
void PeeledCompressedProject(const VECTOR3_FIELD_3D& V, MATRIX_COMPRESSION_DATA* U_data,
    VectorXd* q);

// projection, implemented in the frequency domain
void PeeledCompressedProjectTransform(const VECTOR3_FIELD_3D& V,
    MATRIX_COMPRESSION_DATA* U_data, VectorXd* q);

// projection, implemented in the frequency domain. assumes no SVD!
void PeeledCompressedProjectTransformNoSVD(const VECTOR3_FIELD_3D& V,
    MATRIX_COMPRESSION_DATA* U_data, VectorXd* q);

// set zeros at the places where we have artificially padded
void SetZeroPadding(vector<VectorXd>* blocks, COMPRESSION_DATA* data);

// unproject the reduced coordinate into the peeled cells in this field
// using compression data
void PeeledCompressedUnproject(MATRIX_COMPRESSION_DATA* U_data, const VectorXd& q,
    VECTOR3_FIELD_3D* V);

// unproject the reduced coordinate into the peeled cells in this field
// using compression data. stays in the frequency domain until the end
void PeeledCompressedUnprojectTransform(MATRIX_COMPRESSION_DATA* U_data, const VectorXd& q,
    VECTOR3_FIELD_3D* V);

// scale a vector<VectorXd> each by the same scalar
void ScaleVectorEigen(double alpha, vector<VectorXd>* V);

// increment each element of a vector<VectorXd> correspondingly
// by another vector<VectorXd>
void AddVectorEigen(const vector<VectorXd> V, vector<VectorXd>* W);

// given a row number and the dimensions, computes
// which block number we need for the decoder. populates
// blockIndex with the corresponding value as well.
void ComputeBlockNumber(int row, const VEC3I& dims, int* blockNumber, int* blockIndex);


// given a (row, col), decode the vector
// at that cell from the lossy matrix. assumes row
// is divisible by 3 since it is the start of the vector.
void DecodeFromRowCol(int row, int col, MATRIX_COMPRESSION_DATA* data, Vector3d* cell);


// given a start row, computes the 3 x numCols submatrix
// given the compression data. assumes start row is
// divisible by 3!
void GetSubmatrix(int startRow, MATRIX_COMPRESSION_DATA* data, MatrixXd* submatrix);

// given a start row, computes the 3 x numCols submatrix
// given the compression data. assumes start row is
// divisible by 3! assumes no SVD!
void GetSubmatrixNoSVD(int startRow, MATRIX_COMPRESSION_DATA* data, MatrixXd* submatrix);
void GetSubmatrixNoSVDSparse(int startRow, MATRIX_COMPRESSION_DATA* data, MatrixXd* submatrix);

// clear out all the non-zeros from the passed in blocks
void clearNonZeros(vector<VectorXd>& blocks, const vector<NONZERO_ENTRIES>& allNonZeros);

#endif


