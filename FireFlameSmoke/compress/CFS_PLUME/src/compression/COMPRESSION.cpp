#include <iostream>
#include <sys/stat.h>
#include "COMPRESSION.h"

using std::vector;
using std::cout;
using std::endl;

const double DCT_NORMALIZE = 1.0 / sqrt( 8 * BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE );
const double SQRT_ONEHALF = 1.0/sqrt(2.0);
const double SQRT_TWO = sqrt(2.0);

////////////////////////////////////////////////////////
// global vars, just to see if pushing args onto
// the stack is getting expensive
////////////////////////////////////////////////////////
COMPRESSION_DATA* _compression_data = NULL;
int* _allData;
int _blockNumber;
int _col;
//INTEGER_FIELD_3D* _reverseZigzag = NULL;
INTEGER_FIELD_3D* _unflattened = NULL;
vector<NONZERO_ENTRIES>* _allNonZeros = NULL;
vector<VectorXd>* _decoded = NULL;

////////////////////////////////////////////////////////
// Function Implementations
////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
// cast a FIELD_3D to an INTEGER_FIELD_3D by rounding
// *NOTE* an in-place method for doing the same operation
// is now supported in FIELD_3D
////////////////////////////////////////////////////////
void RoundFieldToInt(const FIELD_3D& F, INTEGER_FIELD_3D* castedField) {
  TIMER functionTimer(__FUNCTION__);

  assert( F.totalCells() == castedField->totalCells() );

  for (int i = 0; i < F.totalCells(); i++) {
    (*castedField)[i] = rint(F[i]);
  }
}

////////////////////////////////////////////////////////
// cast an INTEGER_FIELD_3D to a FIELD_3D
////////////////////////////////////////////////////////
void CastIntFieldToDouble(const INTEGER_FIELD_3D& F, const int totalCells, double* castedField, const NONZERO_ENTRIES& nonZeros)
{
  TIMER functionTimer(__FUNCTION__);

  for (int x = 0; x < nonZeros.size(); x++)
  {
    const int i = nonZeros[x];
    //assert(i < F.totalCells());
    //assert(i >= 0);
    castedField[i] = F[i];
  }
}

////////////////////////////////////////////////////////
// cast an INTEGER_FIELD_3D to a FIELD_3D
////////////////////////////////////////////////////////
void CastIntFieldToDouble(const INTEGER_FIELD_3D& F, const int totalCells, double* castedField, const vector<int>& nonZeros)
{
  TIMER functionTimer(__FUNCTION__);

  for (int x = 0; x < nonZeros.size(); x++)
  {
    const int i = nonZeros[x];
    assert(i < F.totalCells());
    assert(i >= 0);
    castedField[i] = F[i];
  }
}


////////////////////////////////////////////////////////
// cast an INTEGER_FIELD_3D to a FIELD_3D
////////////////////////////////////////////////////////
void CastIntFieldToDouble(const INTEGER_FIELD_3D& F, const int totalCells, double* castedField) {
  TIMER functionTimer(__FUNCTION__);
  // DEBUG
  //static int zeros = 0;
  //static int total = 0;
  //total += totalCells;

  for (int i = 0; i < totalCells; i++) {
    //castedField[i] = (double)F[i];

    // TK: faster? only a little bit ...
    const int entry = F[i];
    castedField[i] = entry ? static_cast<double>(entry) : 0.0;

    //if (entry == 0) zeros++;
  }

  // DEBUG
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << zeros << " zeros of " << total << " entries (" << 100.0 * (float)zeros / total << "%)" << endl;
}

////////////////////////////////////////////////////////
// cast an INTEGER_FIELD_3D to a FIELD_3D
////////////////////////////////////////////////////////
void CastIntFieldToDouble(const INTEGER_FIELD_3D& F, FIELD_3D* castedField) {
  TIMER functionTimer(__FUNCTION__);

  assert( F.totalCells() == castedField->totalCells() );

  for (int i = 0; i < F.totalCells(); i++) {
    (*castedField)[i] = (double)F[i];
  }
}

////////////////////////////////////////////////////////
// operates on a VectorXi and fills another VectorXi
// with its cumulative sum starting at zero and omitting the last
// entry. e.g. if the input vector was (1, 2, 3, 4),
// the result would be (0, 1, 3, 6)
////////////////////////////////////////////////////////
void ModifiedCumSum(const VectorXi& V, VectorXi* sum)
{
  TIMER functionTimer(__FUNCTION__);

  // wipe the output and set it to the appropriate size
  sum->setZero(V.size());
  int accumulation = 0;

  // note the loop starts offset at 1
  for (int i = 1; i < V.size(); i++) {

    // accumulate the previous value
    accumulation += V[i - 1];
    (*sum)[i] = accumulation;
  }

}


////////////////////////////////////////////////////////
// returns the 3 component scalar fields
// from a passed in vector field
////////////////////////////////////////////////////////
void GetScalarFields(const VECTOR3_FIELD_3D& V, FIELD_3D* X, FIELD_3D* Y, FIELD_3D* Z) {
  //TIMER functionTimer(__FUNCTION__);
  *X = V.scalarField(0);
  *Y = V.scalarField(1);
  *Z = V.scalarField(2);
}

// ZigzagFlattned/Unflatten


////////////////////////////////////////////////////////
// Given a passed in buffer and a 'direction'
// (1 for forward, -1 for inverse),
// we return an fftw plan for doing an in-place 3d dct
// which is linked to the in buffer
////////////////////////////////////////////////////////
void Create_DCT_Plan(double* in, int direction, fftw_plan* plan) {
  TIMER functionTimer(__FUNCTION__);

  // direction is 1 for a forward transform, -1 for a backward transform
  assert( direction == 1 || direction == -1 );

  int xRes = BLOCK_SIZE;
  int yRes = BLOCK_SIZE;
  int zRes = BLOCK_SIZE;

  fftw_r2r_kind kind;
  if (direction == 1) {
    kind = FFTW_REDFT10;
  }
  else {
    kind = FFTW_REDFT01;
  }
  // 'in' appears twice since it is in-place
  *plan = fftw_plan_r2r_3d(zRes, yRes, xRes, in, in, kind, kind, kind, FFTW_MEASURE);
}

////////////////////////////////////////////////////////
// perform a unitary normalization on the passed in
// buffer of a field
////////////////////////////////////////////////////////
void DCT_Unitary_Normalize(double* buffer)
{
  TIMER functionTimer(__FUNCTION__);

  int totalCells = BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE;
  int slabSize = BLOCK_SIZE * BLOCK_SIZE;

  for (int i = 0; i < totalCells; i++) {
    buffer[i] *= DCT_NORMALIZE;
  }

  for (int z = 0; z < BLOCK_SIZE; z++) {
    for (int y = 0; y < BLOCK_SIZE; y++) {
      buffer[z * slabSize + y * BLOCK_SIZE] *= SQRT_ONEHALF;
    }
  }

  for (int y = 0; y < BLOCK_SIZE; y++) {
    for (int x = 0; x < BLOCK_SIZE; x++) {
      buffer[y * BLOCK_SIZE + x] *= SQRT_ONEHALF;
    }
  }

  for (int z = 0; z < BLOCK_SIZE; z++) {
    for (int x = 0; x < BLOCK_SIZE; x++) {
      buffer[z * slabSize + x] *= SQRT_ONEHALF;
    }
  }
}

////////////////////////////////////////////////////////
// undo the unitary normalization prior to doing an
// fftw-style idct
////////////////////////////////////////////////////////
void UndoNormalize(FIELD_3D* F)
{
  TIMER functionTimer(__FUNCTION__);
  assert( F->xRes() == BLOCK_SIZE && F->yRes() == BLOCK_SIZE && F->zRes() == BLOCK_SIZE );

  int totalCells = BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE;
  int slabSize = BLOCK_SIZE * BLOCK_SIZE;
  double* buffer = F->data();

  for (int i = 0; i < totalCells; i++) {
    buffer[i] *= DCT_NORMALIZE;
  }

  for (int z = 0; z < BLOCK_SIZE; z++) {
    for (int y = 0; y < BLOCK_SIZE; y++) {
      buffer[z * slabSize + y * BLOCK_SIZE] *= SQRT_TWO;
    }
  }

  for (int y = 0; y < BLOCK_SIZE; y++) {
    for (int x = 0; x < BLOCK_SIZE; x++) {
      buffer[y * BLOCK_SIZE + x] *= SQRT_TWO;
    }
  }

  for (int z = 0; z < BLOCK_SIZE; z++) {
    for (int x = 0; x < BLOCK_SIZE; x++) {
      buffer[z * slabSize + x] *= SQRT_TWO;
    }
  }
}
////////////////////////////////////////////////////////
// undo the unitary normalization on a flattened
// FIELD_3D  prior to doing an fftw-style idct
////////////////////////////////////////////////////////
void UndoNormalizeEigen(VectorXd* F)
{
  TIMER functionTimer(__FUNCTION__);

  int totalCells = BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE;
  assert( F->size() == totalCells );

  int slabSize = BLOCK_SIZE * BLOCK_SIZE;

  for (int i = 0; i < totalCells; i++) {
    (*F)[i] *= DCT_NORMALIZE;
  }

  for (int z = 0; z < BLOCK_SIZE; z++) {
    for (int y = 0; y < BLOCK_SIZE; y++) {
      (*F)[z * slabSize + y * BLOCK_SIZE] *= SQRT_TWO;
    }
  }

  for (int y = 0; y < BLOCK_SIZE; y++) {
    for (int x = 0; x < BLOCK_SIZE; x++) {
      (*F)[y * BLOCK_SIZE + x] *= SQRT_TWO;
    }
  }

  for (int z = 0; z < BLOCK_SIZE; z++) {
    for (int x = 0; x < BLOCK_SIZE; x++) {
      (*F)[z * slabSize + x] *= SQRT_TWO;
    }
  }
}
////////////////////////////////////////////////////////
// given a passed in FIELD_3D, fftw plan, and
// corresponding 'in' buffer, performs the corresponding
// transform on the field. this one is *unitary normalization*
////////////////////////////////////////////////////////
void DCT_Smart_Unitary(const fftw_plan& plan, int direction, double* in, FIELD_3D* F)
{
  //TIMER functionTimer(__FUNCTION__);

  int xRes = F->xRes();
  int yRes = F->yRes();
  int zRes = F->zRes();
  int totalCells = xRes * yRes * zRes;

  assert ( xRes == BLOCK_SIZE && yRes == BLOCK_SIZE && zRes == BLOCK_SIZE );

  if (direction == -1) { // inverse transform; need to pre-normalize!
    UndoNormalize(F);
  }

  // fill the 'in' buffer
  memcpy(in, F->data(), totalCells * sizeof(double));


  //TIMER fftTimer("fftw execute");
  fftw_execute(plan);
  //fftTimer.stop();

  // 'in' is now overwritten to the result of the transform

  if (direction == 1) { // forward transform; need to post-normalize!
    DCT_Unitary_Normalize(in);
  }

  // rewrite F's data with the new contents of in
  memcpy(F->data(), in, totalCells * sizeof(double));
}

////////////////////////////////////////////////////////
// given a passed in flattened-out FIELD_3D, fftw plan, and
// corresponding 'in' buffer, performs the corresponding
// transform on the field. this one is *unitary normalization*
////////////////////////////////////////////////////////
void DCT_Smart_Unitary_Eigen(const fftw_plan& plan, int direction, double* in, VectorXd* F)
{
  TIMER functionTimer(__FUNCTION__);

  int totalCells = F->size();

  assert ( totalCells == BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE );

  if (direction == -1) { // inverse transform; need to pre-normalize!
    UndoNormalizeEigen(F);
  }

  // fill the 'in' buffer
  memcpy(in, F->data(), totalCells * sizeof(double));

  //TIMER fftTimer("fftw execute");
  fftw_execute(plan);
  //fftTimer.stop();
  // 'in' is now overwritten to the result of the transform

  if (direction == 1) { // forward transform; need to post-normalize!
    DCT_Unitary_Normalize(in);
  }

  // rewrite F's data with the new contents of in
  memcpy(F->data(), in, totalCells * sizeof(double));
}



////////////////////////////////////////////////////////
// given passed in dimensions, computes how much we
// have to pad by in each dimension to reach the next
// multiple of BLOCK_SIZE for even block subdivision
////////////////////////////////////////////////////////

void GetPaddings(const VEC3I& v, VEC3I* paddings)
{
  //TIMER functionTimer(__FUNCTION__);
  int xRes = v[0];
  int yRes = v[1];
  int zRes = v[2];
  int xPadding = (BLOCK_SIZE - (xRes % BLOCK_SIZE)) % BLOCK_SIZE;     // how far are you from the next multiple of 8?
  int yPadding = (BLOCK_SIZE - (yRes % BLOCK_SIZE)) % BLOCK_SIZE;
  int zPadding = (BLOCK_SIZE - (zRes % BLOCK_SIZE)) % BLOCK_SIZE;
  (*paddings)[0] = xPadding;
  (*paddings)[1] = yPadding;
  (*paddings)[2] = zPadding;
}

////////////////////////////////////////////////////////
// given a passed in FIELD_3D, pad it and  parse it
// into a vector of 8 x 8 x 8 blocks (listed in row-major order)
////////////////////////////////////////////////////////

void GetBlocks(const FIELD_3D& F, vector<FIELD_3D>* blocks)
{
  TIMER functionTimer(__FUNCTION__);

  int xRes = F.xRes();
  int yRes = F.yRes();
  int zRes = F.zRes();
  VEC3I v(xRes, yRes, zRes);

  VEC3I paddings(0, 0, 0);
  // fill these in with the appropriate paddings
  GetPaddings(v, &paddings);

  FIELD_3D F_padded = F.pad_xyz(paddings);

  // update the resolutions to the padded ones
  xRes = F_padded.xRes();
  yRes = F_padded.yRes();
  zRes = F_padded.zRes();

  // sanity check that our padder had the desired effect
  assert(xRes % BLOCK_SIZE == 0);
  assert(yRes % BLOCK_SIZE == 0);
  assert(zRes % BLOCK_SIZE == 0);

  for (int z = 0; z < zRes/BLOCK_SIZE; z++) {
    for (int y = 0; y < yRes/BLOCK_SIZE; y++) {
      for (int x = 0; x < xRes/BLOCK_SIZE; x++) {
        blocks->push_back(F_padded.subfield(BLOCK_SIZE*x, BLOCK_SIZE*(x+1),
            BLOCK_SIZE*y, BLOCK_SIZE*(y+1), BLOCK_SIZE*z, BLOCK_SIZE*(z+1)));
      }
    }
  }
}
////////////////////////////////////////////////////////
// given a passed in FIELD_3D, pad it with zeros and  parse it
// into a vector of flattened 8 x 8 x 8 blocks (listed in row-major order)
////////////////////////////////////////////////////////

void GetBlocksEigen(const FIELD_3D& F, vector<VectorXd>* blocks)
{
  TIMER functionTimer(__FUNCTION__);

  int xRes = F.xRes();
  int yRes = F.yRes();
  int zRes = F.zRes();
  VEC3I v(xRes, yRes, zRes);

  VEC3I paddings(0, 0, 0);
  // fill these in with the appropriate paddings
  GetPaddings(v, &paddings);

  // call zero pad rather than continuous value pad
  FIELD_3D F_padded = F.zeroPad_xyz(paddings);

  // update the resolutions to the padded ones
  xRes = F_padded.xRes();
  yRes = F_padded.yRes();
  zRes = F_padded.zRes();

  // sanity check that our padder had the desired effect
  assert(xRes % BLOCK_SIZE == 0);
  assert(yRes % BLOCK_SIZE == 0);
  assert(zRes % BLOCK_SIZE == 0);

  // resize blocks appropriately
  blocks->resize(xRes/BLOCK_SIZE * yRes/BLOCK_SIZE * zRes/BLOCK_SIZE);

  int index = 0;
  for (int z = 0; z < zRes/BLOCK_SIZE; z++) {
    for (int y = 0; y < yRes/BLOCK_SIZE; y++) {
      for (int x = 0; x < xRes/BLOCK_SIZE; x++, index++) {
        // add the flattened out block to the list
        (*blocks)[index] = (F_padded.subfield(BLOCK_SIZE*x, BLOCK_SIZE*(x+1),
            BLOCK_SIZE*y, BLOCK_SIZE*(y+1), BLOCK_SIZE*z, BLOCK_SIZE*(z+1)).flattenedEigen());
      }
    }
  }
}
////////////////////////////////////////////////////////
// reconstruct a FIELD_3D with the passed in dims
// from a list of 8 x 8 x 8 blocks
////////////////////////////////////////////////////////

void AssimilateBlocks(const VEC3I& dims, const vector<FIELD_3D>& V, FIELD_3D* assimilatedField)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = dims[0];
  const int yRes = dims[1];
  const int zRes = dims[2];

  assert( xRes % BLOCK_SIZE == 0 && yRes % BLOCK_SIZE == 0 && zRes % BLOCK_SIZE == 0 );
  assert( xRes == assimilatedField->xRes() && yRes == assimilatedField->yRes() && zRes == assimilatedField->zRes() );

  for (int z = 0; z < zRes; z++) {
    for (int y = 0; y < yRes; y++) {
      for (int x = 0; x < xRes; x++) {
        int index = (x/BLOCK_SIZE) + (y/BLOCK_SIZE) * (xRes/BLOCK_SIZE) + (z/BLOCK_SIZE) * (xRes/BLOCK_SIZE) * (yRes/BLOCK_SIZE);     // warning, evil integer division happening!
        (*assimilatedField)(x, y, z) = V[index](x % BLOCK_SIZE, y % BLOCK_SIZE, z % BLOCK_SIZE);
      }
    }
  }

}

////////////////////////////////////////////////////////
// reconstruct a FIELD_3D with the passed in dims
// from a list of 8 x 8 x 8 flattened blocks
////////////////////////////////////////////////////////
void AssimilateBlocksEigen(const VEC3I& dims, vector<VectorXd>* V,
    FIELD_3D* assimilatedField)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = dims[0];
  const int yRes = dims[1];
  const int zRes = dims[2];

  assert( xRes % BLOCK_SIZE == 0 && yRes % BLOCK_SIZE == 0 && zRes % BLOCK_SIZE == 0 );
  assert( xRes == assimilatedField->xRes() && yRes == assimilatedField->yRes() && zRes == assimilatedField->zRes() );
  assert( V->size() == (xRes / BLOCK_SIZE) * (yRes / BLOCK_SIZE) * (zRes / BLOCK_SIZE) );

  int slabSize = BLOCK_SIZE * BLOCK_SIZE;
  for (int z = 0; z < zRes; z++) {
    for (int y = 0; y < yRes; y++) {
      for (int x = 0; x < xRes; x++) {
        int index = (x/BLOCK_SIZE) + (y/BLOCK_SIZE) * (xRes/BLOCK_SIZE) + (z/BLOCK_SIZE) * (xRes/BLOCK_SIZE) * (yRes/BLOCK_SIZE);     // warning, evil integer division happening!
        (*assimilatedField)(x, y, z) = (*V)[index]( (z % BLOCK_SIZE) * slabSize + (y % BLOCK_SIZE) * BLOCK_SIZE + (x % BLOCK_SIZE) );
      }
    }
  }
}
////////////////////////////////////////////////////////
// performs a UNITARY dct/idct on each individual block of a passed in
// vector of blocks. direction 1 is dct, -1 is idct.
////////////////////////////////////////////////////////

void UnitaryBlockDCT(int direction, vector<FIELD_3D>* blocks)
{
  TIMER functionTimer(__FUNCTION__);

  // allocate a buffer for the size of an block-size block block
  double* in = (double*) fftw_malloc(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * sizeof(double));

  // make the appropriate plan
  fftw_plan plan;
  Create_DCT_Plan(in, direction, &plan);

  for (auto itr = blocks->begin(); itr != blocks->end(); ++itr) {
    // take the transform at *itr (which is a FIELD_3D)
    // and overwrite its contents
    DCT_Smart_Unitary(plan, direction, in, &(*itr));
  }

  fftw_free(in);
  fftw_destroy_plan(plan);
  fftw_cleanup();
}
////////////////////////////////////////////////////////
// performs a UNITARY dct/idct on each individual flattened
// block of a passed in vector of blocks. direction 1 is dct, -1 is idct.
////////////////////////////////////////////////////////

void UnitaryBlockDCTEigen(int direction, vector<VectorXd>* blocks)
{
  TIMER functionTimer(__FUNCTION__);

  // allocate a buffer for the size of an block-size block
  double* in = (double*) fftw_malloc(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * sizeof(double));

  // make the appropriate plan
  fftw_plan plan;
  Create_DCT_Plan(in, direction, &plan);

  for (auto itr = blocks->begin(); itr != blocks->end(); ++itr) {
    // take the transform at *itr (which is a FIELD_3D)
    // and overwrite its contents
    DCT_Smart_Unitary_Eigen(plan, direction, in, &(*itr));
  }

  fftw_free(in);
  fftw_destroy_plan(plan);
  fftw_cleanup();
}

////////////////////////////////////////////////////////
// build a block diagonal matrix with repeating copies of the passed
// in matrix A along the diagonal. very poor memory usage!
////////////////////////////////////////////////////////
void BlockDiagonal(const MatrixXd& A, int count, MatrixXd* B)
{
  TIMER functionTimer(__FUNCTION__);

  *B = MatrixXd::Zero(A.rows() * count, A.cols() * count);

  // **************************************************
  // check the memory consumption
  double GB = A.rows() * count * A.cols() * count * sizeof(double) / pow(2.0, 30);
  cout << "B is consuming " << GB << " GB of memory!" << endl;
  // **************************************************

  for (int i = 0; i < count; i++) {
    B->block(i * A.rows(), i * A.cols(), A.rows(), A.cols()) = A;
  }
}

////////////////////////////////////////////////////////
// build a block diagonal matrix with repeating copies of the passed
// in matrix A along the diagonal. assumes for simplicity
// that A is 3 x 3 for the SVD case
////////////////////////////////////////////////////////
void SparseBlockDiagonal(const MatrixXd& A, int count, SparseMatrix<double>* B)
{
  TIMER functionTimer(__FUNCTION__);

  assert( A.rows() == 3 && A.cols() == 3 );

  typedef Eigen::Triplet<double> T;
  vector<T> tripletList;

  int nonzeros = A.rows() * A.cols() * count;
  tripletList.reserve(nonzeros);
  for (int i = 0; i < A.rows() * count; i++) {
    if (i % 3 == 0) {
      tripletList.push_back(T(i, i,     A(0, 0)));
      tripletList.push_back(T(i, i + 1, A(0, 1)));
      tripletList.push_back(T(i, i + 2, A(0, 2)));
      }
    else if (i % 3 == 1) {
      tripletList.push_back(T(i, i - 1, A(1, 0)));
      tripletList.push_back(T(i, i,     A(1, 1)));
      tripletList.push_back(T(i, i + 1, A(1, 2)));
    }
    else { // i % 3 == 2
      tripletList.push_back(T(i, i - 2, A(2, 0)));
      tripletList.push_back(T(i, i - 1, A(2, 1)));
      tripletList.push_back(T(i, i,     A(2, 2)));
    }
  }

  *B = SparseMatrix<double>(A.rows() * count, A.cols() * count);
  B->setFromTriplets(tripletList.begin(), tripletList.end());

}

////////////////////////////////////////////////////////
// given a passed in vec3 field, build a matrix with 3 columns,
// one for each of the x, y, and z components
////////////////////////////////////////////////////////
void BuildXYZMatrix(const VECTOR3_FIELD_3D& V, MatrixXd* A)
{
  TIMER functionTimer(__FUNCTION__);

  int N = V.xRes() * V.yRes() * V.zRes();
  *A = MatrixXd::Zero(N, 3);

  FIELD_3D Vx, Vy, Vz;
  GetScalarFields(V, &Vx, &Vy, &Vz);

  A->col(0) = Vx.flattenedEigen();
  A->col(1) = Vy.flattenedEigen();
  A->col(2) = Vz.flattenedEigen();
}

////////////////////////////////////////////////////////
// find a new 3d coordinate system for a vector field using svd and transform into it
// fill s with the singular values and v^T with the coordinate transform matrix
////////////////////////////////////////////////////////
void TransformVectorFieldSVD(VectorXd* s, MatrixXd* v, VECTOR3_FIELD_3D* V)
{
  TIMER functionTimer(__FUNCTION__);

  MatrixXd xyzMatrix;
  BuildXYZMatrix(*V, &xyzMatrix);
  JacobiSVD<MatrixXd> svd(xyzMatrix, ComputeThinU | ComputeThinV);
  *s = svd.singularValues();
  *v = svd.matrixV();

  int count = V->xRes() * V->yRes() * V->zRes();

  SparseMatrix<double> B;
  SparseBlockDiagonal(v->transpose(), count, &B);

  VectorXd transformProduct = B * V->flattenedEigen();
  memcpy(V->data(), transformProduct.data(), 3 * count * sizeof(double));
}

////////////////////////////////////////////////////////
// perform the coordinate transform without having
// to recompute the SVD by reading in from a cached 3 x 3
// v matrix
////////////////////////////////////////////////////////
void TransformVectorFieldSVDCached(Matrix3d* v, VECTOR3_FIELD_3D* V)
{
  TIMER functionTimer(__FUNCTION__);

  // collapse V into 3 columns by coordinate
  MatrixXd xyzMatrix;
  BuildXYZMatrix(*V, &xyzMatrix);

  // build the sparse block diagonal matrix with V^T on the diagonal
  int count = V->xRes() * V->yRes() * V->zRes();
  SparseMatrix<double> B;
  SparseBlockDiagonal(v->transpose(), count, &B);

  // compute the product and copy it into the result
  VectorXd transformProduct = B * V->flattenedEigen();
  memcpy(V->data(), transformProduct.data(), 3 * count * sizeof(double));
}
////////////////////////////////////////////////////////
// find a new 3d coordinate system for a vector field using svd and transform into it
// fill s with the singular values and v^T with the coordinate transform matrix.
// update the compression data to account for the transform matrix and its
// corresponding singular values.
////////////////////////////////////////////////////////
void TransformVectorFieldSVDCompression(VECTOR3_FIELD_3D* V, COMPRESSION_DATA* data)
{
  TIMER functionTimer(__FUNCTION__);

  // build the N x 3 matrix from V
  MatrixXd xyzMatrix;
  BuildXYZMatrix(*V, &xyzMatrix);

  // compute the thin svd
  JacobiSVD<MatrixXd> svd(xyzMatrix, ComputeThinU | ComputeThinV);

  // fetch the data to be updated
  int numCols = data->get_numCols();
  vector<Vector3d>* singularList = data->get_singularList();
  vector<Matrix3d>* vList = data->get_vList();

  // if it's the first time calling TransformVectorFieldSVD in a chain,
  // preallocate
  if (singularList->size() <= 0 && vList->size() <= 0) {
    singularList->reserve(numCols);
    vList->reserve(numCols);
  }

  // update the compression data
  singularList->push_back(svd.singularValues());
  vList->push_back(svd.matrixV());

  // build the sparse block diagonal transform matrix
  int count = V->xRes() * V->yRes() * V->zRes();
  SparseMatrix<double> B;
  SparseBlockDiagonal(svd.matrixV().transpose(), count, &B);

  // compute the transformation using a matrix-vector multiply
  VectorXd transformProduct = B * V->flattenedEigen();

  // copy the result into V (in-place)
  memcpy(V->data(), transformProduct.data(), 3 * count * sizeof(double));
}

////////////////////////////////////////////////////////
// undo the effects of a previous svd coordinate transformation using the passed
// in v matrix
////////////////////////////////////////////////////////
void UntransformVectorFieldSVD(const MatrixXd& v, VECTOR3_FIELD_3D* transformedV)
{
  TIMER functionTimer(__FUNCTION__);

  MatrixXd xyzMatrix;
  BuildXYZMatrix(*transformedV, &xyzMatrix);

  int count = transformedV->xRes() * transformedV->yRes() * transformedV->zRes();

  SparseMatrix<double> B;
  SparseBlockDiagonal(v, count, &B);

  VectorXd transformProduct = B * transformedV->flattenedEigen();
  memcpy(transformedV->data(), transformProduct.data(), 3 * count * sizeof(double));
}

////////////////////////////////////////////////////////
// Normalize the block contents to a resolution of
// nBits based on the DC component. Update the sList.
////////////////////////////////////////////////////////
void PreprocessBlock(FIELD_3D* F, int blockNumber, int col, COMPRESSION_DATA* data)
{
  TIMER functionTimer(__FUNCTION__);
  int nBits = data->get_nBits();

  // *****************************************************
  // debugging test--let's use absMax instead of DC component
  // *****************************************************
  double Fmax = F->absMax();
  double s = (pow(2, nBits - 1) - 1) / Fmax;
  // normalize so that the DC component is at 2 ^ {nBits - 1} - 1
  // double s = (pow(2, nBits - 1) - 1) / (*F)[0];
  (*F) *= s;

  // fetch data for updating sList
  MatrixXd* sListMatrix = data->get_sListMatrix();
  int numBlocks = data->get_numBlocks();
  int numCols = data->get_numCols();

  // if it's the first time PreprocessBlock is called in a chain, resize
  if (sListMatrix->cols() <= 0) {
    sListMatrix->setZero(numBlocks, numCols);
  }

  // update sList
  (*sListMatrix)(blockNumber, col) = s;
}

////////////////////////////////////////////////////////
// Binary search to find the appropriate gamma given
// desired percent threshold within maxIterations. Prints
// out information as it goes.
////////////////////////////////////////////////////////
void TuneGammaVerbose(const FIELD_3D& F, int blockNumber, int col,
    COMPRESSION_DATA* data, FIELD_3D* damp)
{
  TIMER functionTimer(__FUNCTION__);

  // fetch parameters from data
  int nBits = data->get_nBits();
  int maxIterations = data->get_maxIterations();
  double percent = data->get_percent();

  double lower = 0.0;
  // QUESTION: how should we define upper?
  double upper = nBits;
  cout << "Upper: " << upper << endl;
  // arbitrarily set epsilon to be 0.5%
  double epsilon = 0.005;
  double gamma = 0.5 * (upper + lower);
  cout << "Initial gamma: " << gamma << endl;
  damp->toPower(gamma);
  cout << "Initial damping array: " << endl;
  cout << damp->flattened() << endl;

  cout << "F: " << endl;
  cout << F.flattened() << endl;
  // the total amount of energy in the Fourier space
  double totalEnergy = F.sumSq();
  FIELD_3D damped = ( F / (*damp) );
  damped.roundInt();
  cout << "Damped block: " << endl;
  cout << damped.flattened() << endl;

  cout << "Undamped block: " << endl;
  cout << ((*damp) * damped).flattened() << endl;

  double energyDiff = abs(totalEnergy - ( (*damp) * damped ).sumSq());
  cout << "Absolute energy difference: " << energyDiff << endl;
  double percentEnergy = 1.0 - (energyDiff / totalEnergy);
  cout << "Initial percent: " << percentEnergy << endl;
  int iterations = 0;

  while ( abs( percent - percentEnergy ) > epsilon && iterations < maxIterations) {

    if (percentEnergy < percent) { // too much damping; need to lower gamma
      upper = gamma;
      gamma = 0.5 * (upper + lower);

      // to the power of 1 / upper brings it back to the vanilla state,
      // from which we raise it to the new gamma
      damp->toPower(gamma / upper);
    }

    else { // not enough damping; need to increase gamma
      lower = gamma;
      gamma = 0.5 * (upper + lower);

      // to the power of 1 / lower brings it back to the vanilla state,
      // from which we raise it to the new gamma
      damp->toPower(gamma / lower);
    }

    cout << "New damping array: " << endl;
    cout << damp->flattened() << endl;
    // update percentEnergy
    damped = ( F / (*damp) );
    damped.roundInt();
    cout << "New damped block: " << endl;
    cout << damped.flattened() << endl;
    cout << "New undamped block: " << endl;
    cout << ((*damp) * damped).flattened() << endl;
    energyDiff = abs(totalEnergy - ( (*damp) * damped ).sumSq());
    percentEnergy =  1.0 - (energyDiff / totalEnergy);
    cout << "New percent energy: " << percentEnergy << endl;
    cout << "New gamma: " << gamma << endl;
    iterations++;
    cout << "Next iteration: " << iterations << endl;
  }

  cout << "Took " << iterations << " iterations to compute gamma!\n";
  cout << "Percent Energy ended up at : " << percentEnergy << endl;
  cout << "Gamma ended up at: " << gamma << endl;

  // fetch data to update gammaList
  MatrixXd* gammaListMatrix = data->get_gammaListMatrix();
  int numBlocks = data->get_numBlocks();
  int numCols = data->get_numCols();

  // if it's the first time TuneGamma is called in a chain, resize
  if (gammaListMatrix->cols() <= 0) {
    gammaListMatrix->setZero(numBlocks, numCols);
  }

  (*gammaListMatrix)(blockNumber, col) = gamma;

}

////////////////////////////////////////////////////////
// Binary search to find the appropriate gamma given
// desired percent threshold within maxIterations.
///////////////////////////////////////////////////////
void TuneGamma(const FIELD_3D& F, int blockNumber, int col,
    COMPRESSION_DATA* data, FIELD_3D* damp)
{
  TIMER functionTimer(__FUNCTION__);

  // fetch parameters from data
  int nBits = data->get_nBits();
  int maxIterations = data->get_maxIterations();
  double percent = data->get_percent();

  double lower = 0.0;
  // QUESTION: how should we define upper?
  double upper = nBits;
  // arbitrarily set epsilon to be 0.5%
  double epsilon = 0.005;
  double gamma = 0.5 * (upper + lower);
  damp->toPower(gamma);

  // the total amount of energy in the Fourier space
  double totalEnergy = F.sumSq();
  FIELD_3D damped = ( F / (*damp) );
  damped.roundInt();

  double energyDiff = abs(totalEnergy - ( (*damp) * damped ).sumSq());
  double percentEnergy = 1.0 - (energyDiff / totalEnergy);
  int iterations = 0;

  while ( abs( percent - percentEnergy ) > epsilon && iterations < maxIterations) {

    if (percentEnergy < percent) { // too much damping; need to lower gamma
      upper = gamma;
      gamma = 0.5 * (upper + lower);

      // to the power of 1 / upper brings it back to the vanilla state,
      // from which we raise it to the new gamma
      damp->toPower(gamma / upper);
    }

    else { // not enough damping; need to increase gamma
      lower = gamma;
      gamma = 0.5 * (upper + lower);

      // to the power of 1 / lower brings it back to the vanilla state,
      // from which we raise it to the new gamma
      damp->toPower(gamma / lower);
    }

    // update percentEnergy
    damped = ( F / (*damp) );
    damped.roundInt();
    energyDiff = abs(totalEnergy - ( (*damp) * damped ).sumSq());
    percentEnergy =  1.0 - (energyDiff / totalEnergy);
    iterations++;
  }


  // fetch data to update gammaList
  MatrixXd* gammaListMatrix = data->get_gammaListMatrix();
  int numBlocks = data->get_numBlocks();
  int numCols = data->get_numCols();

  // if it's the first time TuneGamma is called in a chain, resize
  if (gammaListMatrix->cols() <= 0) {
    gammaListMatrix->setZero(numBlocks, numCols);
  }

  (*gammaListMatrix)(blockNumber, col) = gamma;

}

////////////////////////////////////////////////////////
// Binary search to find the appropriate gamma given
// desired percent threshold within maxIterations.
///////////////////////////////////////////////////////
void TuneGammaFastPow(const FIELD_3D& F, int blockNumber, int col,
    COMPRESSION_DATA* data, FIELD_3D* damp)
{
  TIMER functionTimer(__FUNCTION__);

  // fetch parameters from data
  int nBits = data->get_nBits();
  int maxIterations = data->get_maxIterations();
  double percent = data->get_percent();

  // cache the original damping matrix
  FIELD_3D vanilla(*damp);

  double lower = 0.0;
  // QUESTION: how should we define upper?
  double upper = nBits;
  // arbitrarily set epsilon to be very small
  double epsilon = 0.0001;
  // double gamma = 0.5 * (upper + lower);
  double gamma = 0.0;
  damp->toFastPower(gamma);

  // the total amount of energy in the Fourier space
  double totalEnergy = F.sumSq();
  FIELD_3D damped = ( F / (*damp) );
  damped.roundInt();

  double energyDiff = abs(totalEnergy - ( (*damp) * damped ).sumSq());
  double percentEnergy = 1.0 - (energyDiff / totalEnergy);
  int iterations = 0;

  while ( abs( percent - percentEnergy ) > epsilon && iterations < maxIterations) {

    if (percentEnergy < percent) { // too much damping; need to lower gamma
      upper = gamma;
      gamma = 0.5 * (upper + lower);

      (*damp) = vanilla;
      damp->toFastPower(gamma);
    }

    else { // not enough damping; need to increase gamma
      lower = gamma;
      gamma = 0.5 * (upper + lower);

      (*damp) = vanilla;
      damp->toFastPower(gamma);
    }

    // update percentEnergy
    damped = ( F / (*damp) );
    damped.roundInt();
    energyDiff = abs(totalEnergy - ( (*damp) * damped ).sumSq());
    percentEnergy =  1.0 - (energyDiff / totalEnergy);
    iterations++;
  }


  // fetch data to update gammaList
  MatrixXd* gammaListMatrix = data->get_gammaListMatrix();
  int numBlocks = data->get_numBlocks();
  int numCols = data->get_numCols();

  // if it's the first time TuneGamma is called in a chain, resize
  if (gammaListMatrix->cols() <= 0) {
    gammaListMatrix->setZero(numBlocks, numCols);
  }

  (*gammaListMatrix)(blockNumber, col) = gamma;

}

////////////////////////////////////////////////////////
// Binary search to find the appropriate gamma given
// desired percent threshold within maxIterations.
// Quantizes gamma values to quarter-integers in the range
// 0 to nBits.
///////////////////////////////////////////////////////
void TuneGammaQuantized(const FIELD_3D& F, int blockNumber, int col,
    COMPRESSION_DATA* data, FIELD_3D* damp)
{
  TIMER functionTimer(__FUNCTION__);

  // fetch parameters from data
  int nBits = data->get_nBits();
  int maxIterations = data->get_maxIterations();
  double percent = data->get_percent();

  // cache the original damping matrix
  FIELD_3D vanilla(*damp);

  double lower = 0.0;
  // QUESTION: how should we define upper?
  double upper = nBits;
  // arbitrarily set epsilon to be very small
  double epsilon = 0.000001;
  // double gamma = 0.5 * (upper + lower);
  double gamma = 0.0;
  damp->toFastPower(gamma);

  // the total amount of energy in the Fourier space
  double totalEnergy = F.sumSq();
  FIELD_3D damped = ( F / (*damp) );
  damped.roundInt();

  double energyDiff = abs(totalEnergy - ( (*damp) * damped ).sumSq());
  double percentEnergy = 1.0 - (energyDiff / totalEnergy);
  int iterations = 0;

  while ( abs( percent - percentEnergy ) > epsilon && iterations < maxIterations) {

    if (percentEnergy < percent) { // too much damping; need to lower gamma
      upper = gamma;
      gamma = 0.5 * (upper + lower);

      (*damp) = vanilla;
      damp->toFastPower(gamma);
    }

    else { // not enough damping; need to increase gamma
      lower = gamma;
      gamma = 0.5 * (upper + lower);

      (*damp) = vanilla;
      damp->toFastPower(gamma);
    }

    // update percentEnergy
    damped = ( F / (*damp) );
    damped.roundInt();
    energyDiff = abs(totalEnergy - ( (*damp) * damped ).sumSq());
    percentEnergy =  1.0 - (energyDiff / totalEnergy);
    iterations++;
  }

  // fetch data to update gammaList
  MatrixXd* gammaListMatrix = data->get_gammaListMatrix();
  int numBlocks = data->get_numBlocks();
  int numCols = data->get_numCols();

  // if it's the first time TuneGamma is called in a chain, resize
  if (gammaListMatrix->cols() <= 0) {
    gammaListMatrix->setZero(numBlocks, numCols);
  }

  // quantize gamma to the nearest quarter-integer
  gamma = rint(4 * gamma) / 4.0;
  (*gammaListMatrix)(blockNumber, col) = gamma;

  // re-compute damp in case gamma changed
  (*damp) = vanilla;
  damp->toFastPower(gamma);

}
////////////////////////////////////////////////////////
// simply sets gamma equal to zero for no damping. for
// debug purposes only.
///////////////////////////////////////////////////////
void TuneGammaDebug(const FIELD_3D& F, int blockNumber, int col,
    COMPRESSION_DATA* data, FIELD_3D* damp)
{
  TIMER functionTimer(__FUNCTION__);

  // set gamma to zero always, for debugging
  double gamma = 0.0;
  damp->toPower(gamma);

  // fetch data to update gammaList
  MatrixXd* gammaListMatrix = data->get_gammaListMatrix();
  int numBlocks = data->get_numBlocks();
  int numCols = data->get_numCols();

  // if it's the first time TuneGammaDebug is called in a chain, resize
  if (gammaListMatrix->cols() <= 0) {
    gammaListMatrix->setZero(numBlocks, numCols);
  }

  (*gammaListMatrix)(blockNumber, col) = gamma;

}

////////////////////////////////////////////////////////
// takes a passed in FIELD_3D (which is intended to be
// the result of a DCT post-preprocess). calculates the best gamma value
// for a damping array. then damps by that array and
// quantizes the result to an integer. stores the
// value of gamma for the damping.
////////////////////////////////////////////////////////
void EncodeBlock(const FIELD_3D& F, int blockNumber, int col, COMPRESSION_DATA* data,
    INTEGER_FIELD_3D* quantized)
{

  TIMER functionTimer(__FUNCTION__);

  // size the return value appropriately
  quantized->resizeAndWipe(F.xRes(), F.yRes(), F.zRes());

  // grab the pre-cached vanila damping array
  const FIELD_3D& dampingArray = data->get_dampingArray();

  // make a copy for modification during TuneGamma
  FIELD_3D damp = dampingArray;

  int numBlocks = data->get_numBlocks();
  assert(blockNumber >= 0 && blockNumber < numBlocks);

  // finds best gamma given the percent. updates gammaList
  // and updates damp
  if (FIELD_3D::usingFastPow()) {
    // ADJ: switching to quantized gamma!
    //TuneGammaFastPow(F, blockNumber, col, data, &damp);
    TuneGammaQuantized(F, blockNumber, col, data, &damp);
  }
  else {
    TuneGamma(F, blockNumber, col, data, &damp);
  }

  // fill the return value with rounded damped entries
  RoundFieldToInt( (F / damp), quantized );

}

////////////////////////////////////////////////////////
// takes a passed in FIELD_3D (which is intended to be
// the result of a DCT post-preprocess). for debug purposes
// only, does not perform any damping---gamma is set to 0.
////////////////////////////////////////////////////////
void EncodeBlockDebug(const FIELD_3D& F, int blockNumber, int col, COMPRESSION_DATA* data,
    INTEGER_FIELD_3D* quantized)
{
  TIMER functionTimer(__FUNCTION__);

  // size the return value appropriately
  quantized->resizeAndWipe(F.xRes(), F.yRes(), F.zRes());

  // grab the pre-cached vanila damping array
  const FIELD_3D& dampingArray = data->get_dampingArray();

  // make a copy for modification during TuneGamma
  FIELD_3D damp = dampingArray;

  int numBlocks = data->get_numBlocks();
  assert(blockNumber >= 0 && blockNumber < numBlocks);

  // sets gamma to zero for debugging purposes.
  // updates damp to be all 1s
  TuneGammaDebug(F, blockNumber, col, data, &damp);

  // fill the return value with rounded damped entries
  RoundFieldToInt( (F / damp), quantized );

}
////////////////////////////////////////////////////////
// takes a passed in INTEGER_FIELD_3D (which is inteneded to
// be run-length decoded and unzigzagged) corresponding to
// a particulary blockNumber and column of the matrix. undoes
// the effects of damping and quantization as best it can.
////////////////////////////////////////////////////////
/*
void DecodeBlock(const INTEGER_FIELD_3D& intBlock, int blockNumber, int col,
    const DECOMPRESSION_DATA& decompression_data, FIELD_3D* decoded)
{

  TIMER functionTimer(__FUNCTION__);

  int numBlocks = decompression_data.get_numBlocks();
  // make sure we are not accessing an invalid block
  assert( (blockNumber >= 0) && (blockNumber < numBlocks) );

  // we use u, v, w rather than x, y , z to indicate the spatial frequency domain
  const int uRes = intBlock.xRes();
  const int vRes = intBlock.yRes();
  const int wRes = intBlock.zRes();
  // size the decoded block appropriately and fill it with the block data
  decoded->resizeAndWipe(uRes, vRes, wRes);
  CastIntFieldToDouble(intBlock, decoded);

  // use the appropriate scale factor to decode
  const MatrixXd& sListMatrix = decompression_data.get_sListMatrix();
  double s = sListMatrix(blockNumber, col);
  double sInv = 1.0 / s;

  const MatrixXd& gammaListMatrix = decompression_data.get_gammaListMatrix();
  double gamma = gammaListMatrix(blockNumber, col);

  // dequantize by inverting the scaling by s and contracting by the
  // appropriate gamma-modulated damping array
  const FIELD_3D& dampingArray = decompression_data.get_dampingArray();
  FIELD_3D damp = dampingArray;
  damp.toPower(gamma);

  // undo the dampings and preprocess
  (*decoded) *= dampingArray;
  (*decoded) *= sInv;

}
*/

////////////////////////////////////////////////////////
// does the same operations as DecodeBlock, but with a passed
// in compression data parameter rather than decompression data
// due to const poisoning, compression data cannot be marked const,
// but nonetheless it is not modified.
////////////////////////////////////////////////////////
void DecodeBlockWithCompressionData(const INTEGER_FIELD_3D& intBlock,
  int blockNumber, int col, COMPRESSION_DATA* data, Real* decoded)
  // TK: This actually only needs the raw pointer for "decoded", and then we avoid an
  // expensive copy at the end.
  //int blockNumber, int col, COMPRESSION_DATA* data, FIELD_3D* decoded)
{
  TIMER functionTimer(__FUNCTION__);

  int numBlocks = data->get_numBlocks();
  // make sure we are not accessing an invalid block
  assert( (blockNumber >= 0) && (blockNumber < numBlocks) );

  TIMER castTimer("Decode Cast");
  // TK: decoded is always the same size, so have the caller size it
  /*
  // we use u, v, w rather than x, y , z to indicate the spatial frequency domain
  const int uRes = intBlock.xRes();
  const int vRes = intBlock.yRes();
  const int wRes = intBlock.zRes();
  // size the decoded block appropriately and fill it with the block data
  decoded->resizeAndWipe(uRes, vRes, wRes);
  */
  //CastIntFieldToDouble(intBlock, decoded);
  CastIntFieldToDouble(intBlock, intBlock.totalCells(), decoded);

  // use the appropriate scale factor to decode
  MatrixXd* sList = data->get_sListMatrix();
  double s = (*sList)(blockNumber, col);
  double sInv = 1.0 / s;
  MatrixXd* gammaList = data->get_gammaListMatrix();
  double gamma = (*gammaList)(blockNumber, col);

  // dequantize by inverting the scaling by s and contracting by the
  // appropriate gamma-modulated damping array
  TIMER dampingTimer("Decode Damping Copy");
  const FIELD_3D& dampingArray = data->get_dampingArray();
  // TK: Skip the allocation every time, just do a memcpy. Using a static
  // will probably cause OpenMP to misbehave later.
  //FIELD_3D damp = dampingArray;
  static FIELD_3D damp(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  memcpy(damp.data(), dampingArray.dataConst(), damp.totalCells() * sizeof(Real));
  damp.toPower(gamma);

  // undo the dampings and preprocess
  TIMER applyPowerTimer("Decode Apply Power");
  //(*decoded) *= damp;
  //(*decoded) *= sInv;

  // TK: loop through the data just once
  /*
  Real* decodedData = decoded->data();
  for (int x = 0; x < decoded->totalCells(); x++)
    decodedData[x] *= damp[x] * sInv;
  */
  // TK: do it on the raw pointer instead
  for (int x = 0; x < intBlock.totalCells(); x++)
    decoded[x] *= damp[x] * sInv;
}

////////////////////////////////////////////////////////
// does the same operations as DecodeBlock, but with a passed
// in compression data parameter rather than decompression data
// due to const poisoning, compression data cannot be marked const,
// but nonetheless it is not modified.
////////////////////////////////////////////////////////
void DecodeBlockWithCompressionDataSparseStackless()
{
  TIMER functionTimer(__FUNCTION__);

  const int numBlocks = _compression_data->get_numBlocks();
  // make sure we are not accessing an invalid block
  assert( (_blockNumber >= 0) && (_blockNumber < numBlocks) );

  double* decoded = (*_decoded)[_blockNumber].data();
  NONZERO_ENTRIES& nonZeros = (*_allNonZeros)[_blockNumber];
  CastIntFieldToDouble((*_unflattened), (*_unflattened).totalCells(), decoded, nonZeros);

  // use the appropriate scale factor to decode
  MatrixXd* sList = _compression_data->get_sListMatrix();
  //double s = (*sList)(_blockNumber, _col);
  //double sInv = 1.0 / s;
  const double sInv = 1.0 / (*sList)(_blockNumber, _col);
  MatrixXd* gammaList = _compression_data->get_gammaListMatrix();
  const double gamma = (*gammaList)(_blockNumber, _col);

  // dequantize by inverting the scaling by s and contracting by the
  // appropriate gamma-modulated damping array
  //
  // TK: Lots of time can be spent here, but it appears to be because of
  // the pow call, not because the sparsity is not fully exploited
  TIMER dampingTimer("Decode Damping Copy");
  const FIELD_3D& dampingArray = _compression_data->get_dampingArray();
  for (unsigned int x = 0; x < nonZeros.size(); x++)
  {
    const int i = nonZeros[x];
    decoded[i] *= FIELD_3D::fastPow(dampingArray[i],gamma) * sInv;
  }
}

////////////////////////////////////////////////////////
// does the same operations as DecodeBlock, but with a passed
// in compression data parameter rather than decompression data
// due to const poisoning, compression data cannot be marked const,
// but nonetheless it is not modified.
////////////////////////////////////////////////////////
void DecodeBlockWithCompressionDataSparse(const INTEGER_FIELD_3D& intBlock,
  int blockNumber, int col, COMPRESSION_DATA* data, Real* decoded, const NONZERO_ENTRIES& nonZeros)
{
  TIMER functionTimer(__FUNCTION__);
  // DEBUG
  puts("Called the non-quantized decode block!");

  int numBlocks = data->get_numBlocks();
  // make sure we are not accessing an invalid block
  assert( (blockNumber >= 0) && (blockNumber < numBlocks) );

  //TIMER castTimer("Decode Cast");
  //CastIntFieldToDouble(intBlock, intBlock.totalCells(), decoded);
  //FIELD_3D dummy(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  //CastIntFieldToDouble(intBlock, intBlock.totalCells(), dummy.data(), nonZeros);
  CastIntFieldToDouble(intBlock, intBlock.totalCells(), decoded, nonZeros);

  // use the appropriate scale factor to decode
  MatrixXd* sList = data->get_sListMatrix();
  const double sInv = 1.0 / (*sList)(blockNumber, col);
  MatrixXd* gammaList = data->get_gammaListMatrix();
  const double gamma = (*gammaList)(blockNumber, col);

  // dequantize by inverting the scaling by s and contracting by the
  // appropriate gamma-modulated damping array
  const FIELD_3D& dampingArray = data->get_dampingArray();
  for (unsigned int x = 0; x < nonZeros.size(); x++)
  {
    const int i = nonZeros[x];
    decoded[i] *= FIELD_3D::fastPow(dampingArray[i],gamma) * sInv;
  }
}

////////////////////////////////////////////////////////
// does the same operations as DecodeBlock, but with a passed
// in compression data parameter rather than decompression data
// due to const poisoning, compression data cannot be marked const,
// but nonetheless it is not modified.
// assumes we have encoded using TuneGammaQuantized!
////////////////////////////////////////////////////////
void DecodeBlockWithCompressionDataSparseQuantized(const INTEGER_FIELD_3D& intBlock,
  int blockNumber, int col, COMPRESSION_DATA* data, Real* decoded, const NONZERO_ENTRIES& nonZeros)
{
  TIMER functionTimer(__FUNCTION__);

  int numBlocks = data->get_numBlocks();
  // make sure we are not accessing an invalid block
  assert( (blockNumber >= 0) && (blockNumber < numBlocks) );

  CastIntFieldToDouble(intBlock, intBlock.totalCells(), decoded, nonZeros);

  // use the appropriate scale factor to decode
  MatrixXd* sList = data->get_sListMatrix();
  const double sInv = 1.0 / (*sList)(blockNumber, col);
  MatrixXd* gammaList = data->get_gammaListMatrix();
  const double gamma = (*gammaList)(blockNumber, col);

  // dequantize by inverting the scaling by s and contracting by the
  // appropriate gamma-modulated damping array
  //const FIELD_3D& dampingArray = data->get_dampingArray();
  bool arrayListBuilt = data->get_arrayListBuilt();
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << "array list built is: " << arrayListBuilt << endl;
  //if (!arrayListBuilt) { 
    //puts("Have to rebuild array list...");
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    data->set_dampingArrayList(); 
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //}
  const vector<FIELD_3D>& dampingArrayList = data->get_dampingArrayList();

  // gamma is a quarter-integer, so rescaling gives the appropriate index 
  // in the list
  int index = 4 * gamma;
  for (unsigned int x = 0; x < nonZeros.size(); x++)
  {
    const int i = nonZeros[x];
    //decoded[i] *= FIELD_3D::fastPow(dampingArray[i],gamma) * sInv;
    decoded[i] *= dampingArrayList[index][i] * sInv;
  }
}

////////////////////////////////////////////////////////
// given a zigzagged integer buffer, write it to a binary
// file via run-length encoding. updates the blockLengthsMatrix.
////////////////////////////////////////////////////////
void RunLengthEncodeBinary(const char* filename, int blockNumber, int col,
    const VectorXi& zigzaggedArray, COMPRESSION_DATA* compression_data)
{
  TIMER functionTimer(__FUNCTION__);

  FILE* pFile;
  // open a file in append mode since we will call this function repeatedly
  pFile = fopen(filename, "ab+");
  if (pFile == NULL) {
    perror ("Error opening file.");
  }
  else {

    // we use a C++ vector container for our data since we don't know
    // a priori how many entries it will have once encoded
    vector<int> dataList;
    // reserve plenty of space just to be on the safe side
    dataList.reserve(2 * BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

    // initialize variables
    int data = 0;
    int runLength = 0;
    int encodedLength = 0;

    // assuming BLOCK_SIZE blocks
    int length = BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE;

    for (int i = 0; i < length; i++) {

      data = zigzaggedArray[i];

      // add the value once no matter what
      dataList.push_back(data);
      encodedLength++;

      // we already wrote one value, so runLength starts from 1
      runLength = 1;

      // if the current value and the next value agree, increment the run length.
      // don't allow the next value to go out of bounds
      while ( (i + 1 < length) && (zigzaggedArray[i] == zigzaggedArray[i + 1]) ) {
        runLength++;
        i++;
      }

      // we don't bother to write run lengths for singletons
      if (runLength > 1) {
        // use a single repeated value as an 'escape' to indicate a run
        dataList.push_back(data);

        // push the runLength to the data vector
        dataList.push_back(runLength);

        encodedLength += 2;
      }


    }

    // the size of the dataList vector is how long the encoded block will be
    // int encodedLength = dataList.size();

    // fetch the blockLengthsMatrix for updating
    MatrixXi* blockLengthsMatrix = compression_data->get_blockLengthsMatrix();
    int numBlocks = compression_data->get_numBlocks();
    int numCols = compression_data->get_numCols();

    // if the matrix isn't yet allocated, prellocate
    if (blockLengthsMatrix->cols() <= 0) {
      blockLengthsMatrix->setZero(numBlocks, numCols);
    }

    // update the appropriate entry
    (*blockLengthsMatrix)(blockNumber, col) = encodedLength;

    fwrite(dataList.data(), sizeof(int), encodedLength, pFile);
    // this write assumes that C++ vectors are stored in contiguous memory!

    fclose(pFile);
  }
}

////////////////////////////////////////////////////////
// decode a run-length encoded binary file and fill
// a VectorXi with the contents.
////////////////////////////////////////////////////////
void RunLengthDecodeBinary(int* allData, int blockNumber, int col,
    COMPRESSION_DATA* compression_data, VectorXi* parsedData)
{

  TIMER functionTimer(__FUNCTION__);

  MatrixXi* blockLengthsMatrix = compression_data->get_blockLengthsMatrix();
  int compressedBlockSize = (*blockLengthsMatrix)(blockNumber, col);
  assert(compressedBlockSize >= 0 &&
      compressedBlockSize <= 2 * BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);


  MatrixXi* blockIndicesMatrix = compression_data->get_blockIndicesMatrix();
  int blockIndex = (*blockIndicesMatrix)(blockNumber, col);


  /*
  // ***************************************************
  // for debugging

    cout << "compressedBlockSize: " << compressedBlockSize << endl;

    VECTOR block(compressedBlockSize);
    for (int i = 0; i < block.size(); i++) {
      block[i] = allData[blockIndex + i];
    }
    cout << "blockNumber " << blockNumber << ", column " << col << ": " << endl;
    cout << block << endl;

  // ***************************************************
  */


  int i = 0;
  int runLength = 1;

  // TK: This is always the same size, so caller should set
  //parsedData->resize(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);
  parsedData->setZero();
  int* data = parsedData->data();
  int j = 0;

  TIMER whileTimer("Run length while timer");
  while (i < compressedBlockSize) {

    // write the value once
    data[j] = allData[blockIndex + i];


    /*
    // ***************************************************
    // purely for debugging---the zeroth entry of each block
    // should decode to the largest signed nBit integer!
    if (j == 0) {
      int nBits = compression_data->get_nBits();
      // assert( abs(data[j]) == pow(2, nBits - 1) - 1 );
      if ( abs(data[j]) != pow(2, nBits - 1) - 1 ) {
        cout << "DC was mapped to: " << data[j] << endl;
      }
    }
    // ***************************************************
    */


    if ( (i + 1 < compressedBlockSize) &&
        allData[blockIndex + i] == allData[blockIndex + i + 1]) {
      i += 2;
      runLength = allData[blockIndex + i];

      assert(runLength > 1 && runLength <= BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

      // TK: skip it if it's zero
      if (allData[blockIndex + i - 2] != 0)
      {
        // TK: cache this? doesn't seem to make a big difference
        //int* start = data + j + 1;
        //std::fill(start, start + runLength - 1, allData[blockIndex + i - 2]);
        std::fill(data + j + 1, data + j + 1 + runLength - 1, allData[blockIndex + i - 2]);
      }
      j += (runLength - 1);
    }

    i++;
    j++;
  }
}

////////////////////////////////////////////////////////
// decode a run-length encoded binary file and fill
// a VectorXi with the contents.
////////////////////////////////////////////////////////
void RunLengthDecodeBinaryInPlace(int* allData, int blockNumber, int col,
    const INTEGER_FIELD_3D& reverseZigzag,
    COMPRESSION_DATA* compression_data,
    INTEGER_FIELD_3D& parsedDataField)
{
  TIMER functionTimer(__FUNCTION__);
  //MatrixXi* blockLengthsMatrix = compression_data->get_blockLengthsMatrix();
  const int compressedBlockSize = (*compression_data->get_blockLengthsMatrix())(blockNumber, col);
  //assert(compressedBlockSize >= 0 && compressedBlockSize <= 2 * BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  //MatrixXi* blockIndicesMatrix = compression_data->get_blockIndicesMatrix();
  const int blockIndex = (*compression_data->get_blockIndicesMatrix())(blockNumber, col);

  int i = 0;
  //int runLength = 1;
  int j = 0;

  //TIMER whileTimer("Run length in-place while timer");
  while (i < compressedBlockSize)
  {
    const int index = blockIndex + i;
    const int value = allData[index];
    const int next = allData[index + 1];
    const int runLength = allData[index + 2];

    if (value == 0 && next == 0)
    {
      i += 3;
      j += runLength;
    }
    else
    {
      parsedDataField[reverseZigzag[j]] = value;
      if ((i + 1 < compressedBlockSize) && value == next)
      {
        for (int x = 0; x < runLength - 1; x++)
          parsedDataField[reverseZigzag[j + 1 + x]] = value;

        i += 2;
        j += (runLength - 1);
      }
      i++;
      j++;
    }
  }
}

////////////////////////////////////////////////////////
// decode a run-length encoded binary file and fill
// a VectorXi with the contents.
////////////////////////////////////////////////////////
void RunLengthDecodeBinaryInPlaceSparse(int* allData, int blockNumber, int col,
    const INTEGER_FIELD_3D& reverseZigzag,
    COMPRESSION_DATA* compression_data,
    INTEGER_FIELD_3D& parsedDataField,
    NONZERO_ENTRIES& nonZeros)
    //vector<int>& nonZeros)
{
  TIMER functionTimer(__FUNCTION__);
  //MatrixXi* blockLengthsMatrix = compression_data->get_blockLengthsMatrix();
  const int compressedBlockSize = (*compression_data->get_blockLengthsMatrix())(blockNumber, col);
  //assert(compressedBlockSize >= 0 && compressedBlockSize <= 2 * BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  //MatrixXi* blockIndicesMatrix = compression_data->get_blockIndicesMatrix();
  const int blockIndex = (*compression_data->get_blockIndicesMatrix())(blockNumber, col);

  int i = 0;
  //int runLength = 1;
  int j = 0;

  //TIMER whileTimer("Run length in-place while timer");
  while (i < compressedBlockSize)
  {
    const int index = blockIndex + i;
    const int value = allData[index];
    const int next = allData[index + 1];
    const int runLength = allData[index + 2];

    if (value == 0 && next == 0)
    {
      i += 3;
      j += runLength;
    }
    else
    {
      parsedDataField[reverseZigzag[j]] = value;
      if (value != 0)
        nonZeros.push_back(reverseZigzag[j]);
      if ((i + 1 < compressedBlockSize) && value == next)
      {
        for (int x = 0; x < runLength - 1; x++)
        {
          parsedDataField[reverseZigzag[j + 1 + x]] = value;

          if (value != 0)
            nonZeros.push_back(reverseZigzag[j + 1 + x]);
        }

        i += 2;
        j += (runLength - 1);
      }
      i++;
      j++;
    }
  }
}

////////////////////////////////////////////////////////
// decode a run-length encoded binary file and fill
// a VectorXi with the contents.
////////////////////////////////////////////////////////
void RunLengthDecodeBinaryInPlaceSparseStackless()
{
  TIMER functionTimer(__FUNCTION__);
  //MatrixXi* blockLengthsMatrix = compression_data->get_blockLengthsMatrix();
  const int compressedBlockSize = (*_compression_data->get_blockLengthsMatrix())(_blockNumber, _col);
  //assert(compressedBlockSize >= 0 && compressedBlockSize <= 2 * BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  //MatrixXi* blockIndicesMatrix = compression_data->get_blockIndicesMatrix();
  const int blockIndex = (*_compression_data->get_blockIndicesMatrix())(_blockNumber, _col);

  const INTEGER_FIELD_3D& reverseZigzag = _compression_data->get_reverseZigzag();
  NONZERO_ENTRIES& nonZeros = (*_allNonZeros)[_blockNumber];

  int i = 0;
  //int runLength = 1;
  int j = 0;

  //TIMER whileTimer("Run length in-place while timer");
  while (i < compressedBlockSize)
  {
    const int index = blockIndex + i;
    const int value = _allData[index];
    const int next = _allData[index + 1];
    const int runLength = _allData[index + 2];

    if (value == 0 && next == 0)
    {
      i += 3;
      j += runLength;
    }
    else
    {
      (*_unflattened)[reverseZigzag[j]] = value;
      if (value != 0)
        nonZeros.push_back(reverseZigzag[j]);
      if ((i + 1 < compressedBlockSize) && value == next)
      {
        for (int x = 0; x < runLength - 1; x++)
        {
          (*_unflattened)[reverseZigzag[j + 1 + x]] = value;

          if (value != 0)
            nonZeros.push_back(reverseZigzag[j + 1 + x]);
        }

        i += 2;
        j += (runLength - 1);
      }
      i++;
      j++;
    }
  }
}

////////////////////////////////////////////////////////
// Flattends an INTEGER_FIELD_3D through a zig-zag scan
// into a VectorXi. Since the scan always follows the same order,
// we precompute the zigzag scan array, pass it
// as a parameter, and then just do an index lookup
////////////////////////////////////////////////////////
void ZigzagFlatten(const INTEGER_FIELD_3D& F, const INTEGER_FIELD_3D& zigzagArray,
    VectorXi* zigzagged)
{

  TIMER functionTimer(__FUNCTION__);

  int xRes = F.xRes();
  int yRes = F.yRes();
  int zRes = F.zRes();
  assert(xRes == BLOCK_SIZE && yRes == BLOCK_SIZE && zRes == BLOCK_SIZE);

  zigzagged->resize(xRes * yRes * zRes);
  int totalCells = BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE;

  for (int i = 0; i < totalCells; i++) {
    int index = zigzagArray[i];
    (*zigzagged)[index] = F[i];
  }

}

////////////////////////////////////////////////////////
// Unflattens a VectorXi into an INTEGER_FIELD_3D.
// uses precomputed zigzagArray and simple lookups.
////////////////////////////////////////////////////////
void ZigzagUnflatten(const VectorXi& V, const INTEGER_FIELD_3D& zigzagArray,
    INTEGER_FIELD_3D* unflattened)
{
  TIMER functionTimer(__FUNCTION__);

  assert(V.size() == BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  // TK: unflattened is always the same size, so have caller size it
  //unflattened->resizeAndWipe(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);

  const int totalCells = BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE;
  for (int i = 0; i < totalCells; i++) {
    //int index = zigzagArray[i];
    //(*unflattened)[i] = V[index];
    // TK: Single-liner gives the compiler a hint?
    (*unflattened)[i] = V[zigzagArray[i]];
  }

}

////////////////////////////////////////////////////////
// reads from a binary file the SVD data and sets the
// values inside compression data
////////////////////////////////////////////////////////
void ReadSVDData(const char* filename, COMPRESSION_DATA* data)
{

  FILE* pFile;
  pFile = fopen(filename, "rb");

  if (pFile == NULL) {
    perror("Error opening SVD data");
    exit(EXIT_FAILURE);
  }

  // fetch number of columns to preallocate the lists
  int numCols = data->get_numCols();
  vector<Vector3d>* singularList = data->get_singularList();
  vector<Matrix3d>* vList        = data->get_vList();

  // preallocate the correct size
  singularList->resize(numCols);
  vList->resize(numCols);

  // for each column, first read the singular values, then read
  // the values in the V matrix. there are 3 singular values,
  // and the V matrices are all 3 x 3.
  for (int col = 0; col < numCols; col++) {

    fread((*singularList)[col].data(), sizeof(double), 3, pFile);

    fread((*vList)[col].data(), sizeof(double), 3 * 3, pFile);
  }

  fclose(pFile);
}

////////////////////////////////////////////////////////
// reads from a binary file into a buffer, and sets
// important initializations inside compression data
////////////////////////////////////////////////////////

int* ReadBinaryFileToMemory(const char* filename,
    COMPRESSION_DATA* data)
{
  TIMER functionTimer(__FUNCTION__);

  // initialize what we will return
  int* allData = NULL;

  FILE* pFile;

  pFile = fopen(filename, "rb");
  if (pFile == NULL) {
    perror("Error opening file.");
    exit(EXIT_FAILURE);
  }

  else {

    // build the damping array and zigzag arrays
    data->set_dampingArray();
    data->set_zigzagArray();

    // read nBits and set it
    int nBits;
    fread(&nBits, 1, sizeof(int), pFile);
    data->set_nBits(nBits);
    // cout << "nBits: " << nBits << endl;

    // set the damping list array 
    // must be done *after* reading in nBits!
    data->set_dampingArrayList();

    // read dims, numCols, and numBlocks
    int xRes, yRes, zRes;
    fread(&xRes, 1, sizeof(int), pFile);
    fread(&yRes, 1, sizeof(int), pFile);
    fread(&zRes, 1, sizeof(int), pFile);
    VEC3I dims(xRes, yRes, zRes);
    data->set_dims(dims);
    // cout << "dims: " << dims << endl;

    int numCols, numBlocks;
    fread(&numCols, 1, sizeof(int), pFile);
    fread(&numBlocks, 1, sizeof(int), pFile);
    // set the decompression data accordingly
    data->set_numCols(numCols);
    data->set_numBlocks(numBlocks);
    // cout << "numCols: " << numCols << endl;
    // cout << "numBlocks: " << numBlocks << endl;

    // read in the sListMatrix and set the data
    int blocksXcols = numBlocks * numCols;
    MatrixXd* sListMatrix = data->get_sListMatrix();
    sListMatrix->resize(numBlocks, numCols);
    fread(sListMatrix->data(), blocksXcols, sizeof(double), pFile);

    // do the same for gammaListMatrix
    MatrixXd* gammaListMatrix = data->get_gammaListMatrix();
    gammaListMatrix->resize(numBlocks, numCols);
    fread(gammaListMatrix->data(), blocksXcols, sizeof(double), pFile);

    // do the same for the blockLengthsMatrix, except the data are ints
    MatrixXi* blockLengthsMatrix = data->get_blockLengthsMatrix();
    blockLengthsMatrix->resize(numBlocks, numCols);
    fread(blockLengthsMatrix->data(), blocksXcols, sizeof(int), pFile);

    // cout << "blockLengthsMatrix, column 0: " << endl;
    VectorXi blockLengths0 = blockLengthsMatrix->col(0);
    // cout << EIGEN::convertInt(blockLengths0) << endl;

    // store the total length of all blocks to be able to
    // read in the full compressed data later
    int totalLength = blockLengthsMatrix->sum();
    // cout << "totalLength: " << totalLength << endl;

    // read in blockIndicesMatrix
    MatrixXi* blockIndicesMatrix = data->get_blockIndicesMatrix();
    blockIndicesMatrix->resize(numBlocks, numCols);
    fread(blockIndicesMatrix->data(), blocksXcols, sizeof(int), pFile);

    // cout << "blockIndicesMatrix, column 0: " << endl;
    // cout << EIGEN::convertInt(blockIndicesMatrix->col(0)) << endl;

    // finally, read in the full compressed data
    allData = (int*) malloc(totalLength * sizeof(int));
    if (allData == NULL) {
      perror("Malloc failed to allocate allData!");
      exit(EXIT_FAILURE);
    }

    fread(allData, totalLength, sizeof(int), pFile);
  }

  return allData;
}

////////////////////////////////////////////////////////
// deletes a file if it already exists
////////////////////////////////////////////////////////
void DeleteIfExists(const char* filename)
{
  TIMER functionTimer(__FUNCTION__);
    struct stat buf;                   // dummy to pass in to stat()
    if (stat(filename, &buf) == 0) {   // if a file named 'filename' already exists
      cout << filename << " exists; deleting it first..." << endl;
      if ( remove(filename) != 0 ) {
        perror("Error deleting file.");
        return;
      }
      else {
        cout << filename << " successfully deleted." << endl;
        return;
      }
    }
    else {                               // file does not yet exist
      cout << filename << " does not yet exist; safely opening in append mode." << endl;
    }
    return;
  }


////////////////////////////////////////////////////////
// takes an input FIELD_3D at a particular matrix column
// which is the result of an SVD coordinate transform, compresses
// it according to the general scheme, and writes it to a binary file.
// meant to be called in a chain so that the binary file
// continues to grow.
////////////////////////////////////////////////////////

void CompressAndWriteField(const char* filename, const FIELD_3D& F, int col,
    COMPRESSION_DATA* compression_data)
{
  TIMER functionTimer(__FUNCTION__);

  // fetch some compression data
  int numBlocks = compression_data->get_numBlocks();
  int numCols = compression_data->get_numCols();

  const INTEGER_FIELD_3D& zigzagArray = compression_data->get_zigzagArray();
  MatrixXi* blockLengthsMatrix = compression_data->get_blockLengthsMatrix();

  // if it's the first time calling this routine in a chain, preallocate
  // the matrix
  if (blockLengthsMatrix->cols() <= 0) {
    blockLengthsMatrix->resize(numBlocks, numCols);
  }

  // subdivide F into blocks
  vector<FIELD_3D> blocks;
  GetBlocks(F, &blocks);

  // do the forward transform
  UnitaryBlockDCT(1, &blocks);

  // initialize the relevant variables before looping through all the blocks
  VectorXi blockLengths(numBlocks);
  INTEGER_FIELD_3D intEncoded_i;
  VectorXi zigzagged_i;

  // loop through the blocks and apply the encoding procedure
  for (int i = 0; i < numBlocks; i++) {

    // rescales data and updates sList
    PreprocessBlock(&(blocks[i]), i, col, compression_data);

    // performs quantization and damping. updates gammaList
    EncodeBlock(blocks[i], i, col, compression_data, &intEncoded_i);

    // do the zigzag scan for run-length encoding
    ZigzagFlatten(intEncoded_i, zigzagArray, &zigzagged_i);

    // performs run-length encoding. updates blockLengthsMatrix. since
    // it opens 'filename' in append mode, it can be called in a chain
    RunLengthEncodeBinary(filename, i, col, zigzagged_i, compression_data);
  }

}

////////////////////////////////////////////////////////
// takes an input FIELD_3D at a particular matrix column
// which is the result of an SVD coordinate transform, compresses
// it according to the general scheme, and writes it to a binary file.
// meant to be called in a chain so that the binary file
// continues to grow. for debugging, gamma is set to zero everywhere!
////////////////////////////////////////////////////////

void CompressAndWriteFieldDebug(const char* filename, const FIELD_3D& F, int col,
    COMPRESSION_DATA* compression_data)
{
  TIMER functionTimer(__FUNCTION__);

  // fetch some compression data
  int numBlocks = compression_data->get_numBlocks();
  int numCols = compression_data->get_numCols();

  const INTEGER_FIELD_3D& zigzagArray = compression_data->get_zigzagArray();
  MatrixXi* blockLengthsMatrix = compression_data->get_blockLengthsMatrix();

  // if it's the first time calling this routine in a chain, preallocate
  // the matrix
  if (blockLengthsMatrix->cols() <= 0) {
    blockLengthsMatrix->resize(numBlocks, numCols);
  }

  // subdivide F into blocks
  vector<FIELD_3D> blocks;
  GetBlocks(F, &blocks);

  // do the forward transform
  UnitaryBlockDCT(1, &blocks);

  // initialize the relevant variables before looping through all the blocks
  VectorXi blockLengths(numBlocks);
  INTEGER_FIELD_3D intEncoded_i;
  VectorXi zigzagged_i;

  // loop through the blocks and apply the encoding procedure
  for (int i = 0; i < numBlocks; i++) {

    // rescales data and updates sList
    PreprocessBlock(&(blocks[i]), i, col, compression_data);

    // performs quantization and damping. updates gammaList
    EncodeBlockDebug(blocks[i], i, col, compression_data, &intEncoded_i);

    // do the zigzag scan for run-length encoding
    ZigzagFlatten(intEncoded_i, zigzagArray, &zigzagged_i);

    // performs run-length encoding. updates blockLengthsMatrix. since
    // it opens 'filename' in append mode, it can be called in a chain
    RunLengthEncodeBinary(filename, i, col, zigzagged_i, compression_data);
  }

}

////////////////////////////////////////////////////////
// build the block indices matrix from the block lengths matrix
////////////////////////////////////////////////////////
void BuildBlockIndicesMatrix(COMPRESSION_DATA* data)
{
  MatrixXi* blockLengths = data->get_blockLengthsMatrix();
  MatrixXi* blockIndices = data->get_blockIndicesMatrix();

  // copy the block lengths to start
  *blockIndices = (*blockLengths);

  // flatten column-wise into a vector
  blockIndices->resize(blockLengths->rows() * blockLengths->cols(), 1);

  // compute the cumulative sum
  VectorXi sum;
  ModifiedCumSum(blockIndices->col(0), &sum);

  // copy back into block indices
  memcpy(blockIndices->data(), sum.data(), sum.size() * sizeof(int));

  // reshape appropriately
  blockIndices->resize(blockLengths->rows(), blockLengths->cols());

}

////////////////////////////////////////////////////////
// build the block indices matrix from the block lengths matrix.
// uses explicitly passed in matrices for debugging!
////////////////////////////////////////////////////////
void BuildBlockIndicesMatrixDebug(const MatrixXi blockLengths, MatrixXi* blockIndices)
{
  TIMER functionTimer(__FUNCTION__);

  // copy the block lengths to start
  *blockIndices = blockLengths;

  // flatten column-wise into a vector
  blockIndices->resize(blockLengths.rows() * blockLengths.cols(), 1);

  // compute the cumulative sum
  VectorXi sum;
  ModifiedCumSum(blockIndices->col(0), &sum);

  // copy back into block indices
  memcpy(blockIndices->data(), sum.data(), sum.size() * sizeof(int));

  // reshape appropriately
  blockIndices->resize(blockLengths.rows(), blockLengths.cols());

}

////////////////////////////////////////////////////////
// given a row number and the dimensions, computes
// which block number we need for the decoder. populates
// blockIndex with the corresponding value as well.
////////////////////////////////////////////////////////

void ComputeBlockNumber(int row, const VEC3I& dims, int* blockNumber, int* blockIndex)
{
  TIMER functionTimer(__FUNCTION__);

  int xRes = dims[0];
  int yRes = dims[1];
  int zRes = dims[2];

  assert( row >= 0 && row < 3 * xRes * yRes * zRes);

  // evil integer division!
  int index = row / 3;
  int z = index / (xRes * yRes);         // index = (xRes * yRes) * z + remainder1
  int rem1 = index - (xRes * yRes * z);
  int y = rem1 / xRes;                   // rem1  = xRes * y          + remainder2
  int rem2 = rem1 - xRes * y;
  int x = rem2;

  int u = x % BLOCK_SIZE;
  int v = y % BLOCK_SIZE;
  int w = z % BLOCK_SIZE;
  // return the blockIndex
  *blockIndex = u + BLOCK_SIZE * v + BLOCK_SIZE * BLOCK_SIZE * w;

  // sanity check!
  assert(index == z * xRes * yRes + y * xRes + x);

  // get the padded resolutions
  VEC3I paddings;
  GetPaddings(dims, &paddings);
  xRes += paddings[0];
  yRes += paddings[1];
  zRes += paddings[2];

  // more evil integer division to compute the block number!
  *blockNumber = x/BLOCK_SIZE + (y/BLOCK_SIZE * (xRes/BLOCK_SIZE)) + (z/BLOCK_SIZE * (xRes/BLOCK_SIZE) * (yRes/BLOCK_SIZE));

}

////////////////////////////////////////////////////////
// given a (row, col), decode the vector
// at that cell from the lossy matrix. assumes row
// is divisible by 3 since it is the start of the vector.
////////////////////////////////////////////////////////

void DecodeFromRowCol(int row, int col, MATRIX_COMPRESSION_DATA* data, Vector3d* cell)
{
  TIMER functionTimer(__FUNCTION__);

  assert ( row % 3 == 0 );

  // using X is arbitrary---it will be the same for all three
  COMPRESSION_DATA* dataX = data->get_compression_dataX();

  // initialize useful data
  const VEC3I& dims = dataX->get_dims();
  int numBlocks = dataX->get_numBlocks();
  const INTEGER_FIELD_3D& zigzagArray = dataX->get_zigzagArray();

  int* allDataX = data->get_dataX();
  int* allDataY = data->get_dataY();
  int* allDataZ = data->get_dataZ();
  COMPRESSION_DATA* dataY = data->get_compression_dataY();
  COMPRESSION_DATA* dataZ = data->get_compression_dataZ();

  // compute the block number and index
  int blockIndex, blockNumber;
  ComputeBlockNumber(row, dims, &blockNumber, &blockIndex);
  cout << "block number: " << blockNumber << endl;
  cout << "block index: " << blockIndex << endl;

  // container for the encoded blocks
  vector<FIELD_3D> blocks(numBlocks);

  // variable to store decoded run-length blocks
  // TK: this is always the same size, so set it once
  VectorXi runLengthDecoded(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  // vector prior to undoing the svd
  Vector3d preSVD;

  /*
   *****************************************************************
   * X coordinate
   *****************************************************************
  */

  // decode the run-length scheme
  RunLengthDecodeBinary(allDataX, blockNumber, col, dataX, &runLengthDecoded);

  // undo the zigzag scan
  INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

  // undo the scaling from the quantizer
  FIELD_3D decodedBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  //DecodeBlockWithCompressionData(unflattened, blockNumber, col, dataX, &decodedBlock);
  // TK: Do it on the raw pointer instead
  DecodeBlockWithCompressionData(unflattened, blockNumber, col, dataX, decodedBlock.data());

  // inverse DCT
  const fftw_plan& plan = data->get_plan();
  double* in = data->get_dct_in();
  int direction = -1;
  DCT_Smart_Unitary(plan, direction, in, &decodedBlock);
  preSVD[0] = decodedBlock[blockIndex];

  /*
   *****************************************************************
   * Y coordinate
   *****************************************************************
  */

  // decode the run-length scheme
  RunLengthDecodeBinary(allDataY, blockNumber, col, dataY, &runLengthDecoded);

  // undo the zigzag scan
  ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

  // undo the scaling from the quantizer
  // TK: Do it on the raw pointer instead
  //DecodeBlockWithCompressionData(unflattened, blockNumber, col, dataY, &decodedBlock);
  DecodeBlockWithCompressionData(unflattened, blockNumber, col, dataY, decodedBlock.data());

  // inverse DCT
  DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

  preSVD[1] = decodedBlock[blockIndex];

  /*
   *****************************************************************
   * Z coordinate
   *****************************************************************
  */

  // decode the run-length scheme
  RunLengthDecodeBinary(allDataZ, blockNumber, col, dataZ, &runLengthDecoded);

  // undo the zigzag scan
  ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

  // undo the scaling from the quantizer
  // TK: Do it on the raw pointer instead
  //DecodeBlockWithCompressionData(unflattened, blockNumber, col, dataZ, &decodedBlock);
  DecodeBlockWithCompressionData(unflattened, blockNumber, col, dataZ, decodedBlock.data());

  // inverse DCT
  DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

  preSVD[2] = decodedBlock[blockIndex];

  // fetch the V matrix to undo the SVD
  vector<Matrix3d>* vList = dataX->get_vList();
  Matrix3d V = (*vList)[col];
  *cell = V * preSVD;

}


////////////////////////////////////////////////////////
// given a start row, computes the 3 x numCols submatrix
// given the compression data. assumes start row is
// divisible by 3!
////////////////////////////////////////////////////////
void GetSubmatrix(int startRow, MATRIX_COMPRESSION_DATA* data, MatrixXd* submatrix)
{
  TIMER functionTimer(__FUNCTION__);

  assert ( startRow % 3 == 0 );

  // start fetching useful parameters
  COMPRESSION_DATA* compression_dataX = data->get_compression_dataX();
  const VEC3I& dims = compression_dataX->get_dims();
  int numCols = compression_dataX->get_numCols();

  // info for undoing the SVD transformation. using X is arbitrary and by convention
  vector<Matrix3d>* vList = compression_dataX->get_vList();

  // if submatrix is not yet allocated, resize it to be 3 x numCols
  if (submatrix->cols() <= 0) {
    submatrix->resize(3, numCols);
  }

  // compute the block number and index
  int blockNumber, blockIndex;
  ComputeBlockNumber(startRow, dims, &blockNumber, &blockIndex);

  // compare the cache to the current
  int cachedBlockNumber = data->get_cachedBlockNumber();
  cout << "cachedBlockNumber: " << cachedBlockNumber << endl;
  if (blockNumber == cachedBlockNumber) { // if we've already decoded this block
    //TIMER cacheTimer("cached block");

    cout << "Used cache!" << endl;

    // load the previously decoded data
    vector<FIELD_3D>* cachedBlocksX = data->get_cachedBlocksX();
    vector<FIELD_3D>* cachedBlocksY = data->get_cachedBlocksY();
    vector<FIELD_3D>* cachedBlocksZ = data->get_cachedBlocksZ();

    // SVD transformation matrix container
    Matrix3d V_i;

    for (int i = 0; i < numCols; i++) {

      // a dummy container for each column of the matrix
      Vector3d col_i;
      col_i[0] = (*cachedBlocksX)[i][blockIndex];
      col_i[1] = (*cachedBlocksY)[i][blockIndex];
      col_i[2] = (*cachedBlocksZ)[i][blockIndex];

      // undo the SVD transformation
      V_i = (*vList)[i];
      Vector3d postSVD_i = V_i * col_i;

      // set the column
      submatrix->col(i) = postSVD_i;
    }
  return;
  }

  else { // no cache; have to compute it from scratch
    cout << "Didn't use cache!" << endl;
    TIMER uncachedTimer("uncached block");

    // initialize useful data
    int numBlocks = compression_dataX->get_numBlocks();
    const INTEGER_FIELD_3D& zigzagArray = compression_dataX->get_zigzagArray();

    int* allDataX = data->get_dataX();
    int* allDataY = data->get_dataY();
    int* allDataZ = data->get_dataZ();
    COMPRESSION_DATA* compression_dataY = data->get_compression_dataY();
    COMPRESSION_DATA* compression_dataZ = data->get_compression_dataZ();

    // container for the encoded blocks
    vector<FIELD_3D> blocks(numBlocks);

    // variable to store decoded run-length blocks
    // TK: this is always the same size, so set it once
    VectorXi runLengthDecoded(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

    // vector prior to undoing the svd
    Vector3d preSVD;

    // data for doing the inverse DCT
    const fftw_plan& plan = data->get_plan();
    double* in = data->get_dct_in();
    int direction = -1;


    // pointers to the cache
    vector<FIELD_3D>* cachedBlocksX = data->get_cachedBlocksX();
    vector<FIELD_3D>* cachedBlocksY = data->get_cachedBlocksY();
    vector<FIELD_3D>* cachedBlocksZ = data->get_cachedBlocksZ();

    for (int i = 0; i < numCols; i++) {

      // a container for column i
      Vector3d preSVD_i;
      /*
       *****************************************************************
       * X coordinate
       *****************************************************************
      */

      // decode the run-length scheme
      RunLengthDecodeBinary(allDataX, blockNumber, i, compression_dataX, &runLengthDecoded);

      // undo the zigzag scan
      INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
      ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

      // undo the scaling from the quantizer
      FIELD_3D decodedBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
      // TK: Do it on the raw pointer instead
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataX, &decodedBlock);
      DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataX, decodedBlock.data());

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[0] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksX)[i] = decodedBlock;

      /*
       *****************************************************************
       * Y coordinate
       *****************************************************************
      */

      // decode the run-length scheme
      RunLengthDecodeBinary(allDataY, blockNumber, i, compression_dataY, &runLengthDecoded);

      // undo the zigzag scan
      ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

      // undo the scaling from the quantizer
      // TK: Do it on the raw pointer instead
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataY, &decodedBlock);
      DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataY, decodedBlock.data());

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[1] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksY)[i] = decodedBlock;

      /*
       *****************************************************************
       * Z coordinate
       *****************************************************************
      */

      // decode the run-length scheme
      RunLengthDecodeBinary(allDataZ, blockNumber, i, compression_dataZ, &runLengthDecoded);

      // undo the zigzag scan
      ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

      // undo the scaling from the quantizer
      // TK: Do it on the raw pointer instead
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataZ, &decodedBlock);
      DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataZ, decodedBlock.data());

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[2] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksZ)[i] = decodedBlock;

      // fetch the V matrix to undo the SVD
      Matrix3d V_i = (*vList)[i];

      // undo the SVD transformation
      V_i = (*vList)[i];
      Vector3d postSVD_i = V_i * preSVD_i;

      // set the column
      submatrix->col(i) = postSVD_i;
    }

    // now that the cache has been filled, set the block number
    data->set_cachedBlockNumber(blockNumber);

  }
}

////////////////////////////////////////////////////////
// given a start row, computes the 3 x numCols submatrix
// given the compression data. assumes start row is
// divisible by 3! assumes no SVD!
////////////////////////////////////////////////////////
void GetSubmatrixNoSVD(int startRow, MATRIX_COMPRESSION_DATA* data, MatrixXd* submatrix)
{
  TIMER functionTimer(__FUNCTION__);

  assert ( startRow % 3 == 0 );

  // start fetching useful parameters
  COMPRESSION_DATA* compression_dataX = data->get_compression_dataX();
  const VEC3I& dims = compression_dataX->get_dims();
  int numCols = compression_dataX->get_numCols();

  // resize if necessary
  submatrix->resize(3, numCols);

  // compute the block number and index
  int blockNumber, blockIndex;
  ComputeBlockNumber(startRow, dims, &blockNumber, &blockIndex);

  // compare the cache to the current
  int cachedBlockNumber = data->get_cachedBlockNumber();
  if (blockNumber == cachedBlockNumber) { // if we've already decoded this block
    //TIMER cacheTimer("cached block");

    // cout << "Used cache!" << endl;

    // load the previously decoded data
    vector<FIELD_3D>* cachedBlocksX = data->get_cachedBlocksX();
    vector<FIELD_3D>* cachedBlocksY = data->get_cachedBlocksY();
    vector<FIELD_3D>* cachedBlocksZ = data->get_cachedBlocksZ();

    // SVD transformation matrix container
    Matrix3d V_i;

    for (int i = 0; i < numCols; i++) {

      // a dummy container for each column of the matrix
      Vector3d col_i;
      col_i[0] = (*cachedBlocksX)[i][blockIndex];
      col_i[1] = (*cachedBlocksY)[i][blockIndex];
      col_i[2] = (*cachedBlocksZ)[i][blockIndex];

      // set the column
      submatrix->col(i) = col_i;
    }
  return;
  }

  else { // no cache; have to compute it from scratch
    // cout << "Didn't use cache!" << endl;
    //TIMER uncachedTimer("uncached block");

    // initialize useful data
    int numBlocks = compression_dataX->get_numBlocks();
    const INTEGER_FIELD_3D& zigzagArray = compression_dataX->get_zigzagArray();

    int* allDataX = data->get_dataX();
    int* allDataY = data->get_dataY();
    int* allDataZ = data->get_dataZ();
    COMPRESSION_DATA* compression_dataY = data->get_compression_dataY();
    COMPRESSION_DATA* compression_dataZ = data->get_compression_dataZ();

    // container for the encoded blocks
    vector<FIELD_3D> blocks(numBlocks);

    // variable to store decoded run-length blocks
    // TK: this is always the same size, so set it once
    VectorXi runLengthDecoded(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

    // vector prior to undoing the svd
    Vector3d preSVD;

    // data for doing the inverse DCT
    const fftw_plan& plan = data->get_plan();
    double* in = data->get_dct_in();
    int direction = -1;


    // pointers to the cache
    vector<FIELD_3D>* cachedBlocksX = data->get_cachedBlocksX();
    vector<FIELD_3D>* cachedBlocksY = data->get_cachedBlocksY();
    vector<FIELD_3D>* cachedBlocksZ = data->get_cachedBlocksZ();

    for (int i = 0; i < numCols; i++) {

      // a container for column i
      Vector3d preSVD_i;
      /*
       *****************************************************************
       * X coordinate
       *****************************************************************
      */

      // decode the run-length scheme
      RunLengthDecodeBinary(allDataX, blockNumber, i, compression_dataX, &runLengthDecoded);

      // undo the zigzag scan
      INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
      ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

      // undo the scaling from the quantizer
      FIELD_3D decodedBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
      // TK: Do it on the raw pointer instead
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataX, &decodedBlock);
      DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataX, decodedBlock.data());

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[0] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksX)[i] = decodedBlock;

      /*
       *****************************************************************
       * Y coordinate
       *****************************************************************
      */

      // decode the run-length scheme
      RunLengthDecodeBinary(allDataY, blockNumber, i, compression_dataY, &runLengthDecoded);

      // undo the zigzag scan
      ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

      // undo the scaling from the quantizer
      // TK: Do it on the raw pointer instead
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataY, &decodedBlock);
      DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataY, decodedBlock.data());

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[1] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksY)[i] = decodedBlock;

      /*
       *****************************************************************
       * Z coordinate
       *****************************************************************
      */

      // decode the run-length scheme
      RunLengthDecodeBinary(allDataZ, blockNumber, i, compression_dataZ, &runLengthDecoded);

      // undo the zigzag scan
      ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

      // undo the scaling from the quantizer
      // TK: Do it on the raw pointer instead
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataZ, &decodedBlock);
      DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataZ, decodedBlock.data());

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[2] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksZ)[i] = decodedBlock;

      // set the column
      submatrix->col(i) = preSVD_i;
    }

    // now that the cache has been filled, set the block number
    data->set_cachedBlockNumber(blockNumber);

  }
}

////////////////////////////////////////////////////////
// given a start row, computes the 3 x numCols submatrix
// given the compression data. assumes start row is
// divisible by 3! assumes no SVD!
////////////////////////////////////////////////////////
void GetSubmatrixNoSVDSparse(int startRow, MATRIX_COMPRESSION_DATA* data, MatrixXd* submatrix)
{
  TIMER functionTimer(__FUNCTION__);
  assert ( startRow % 3 == 0 );

  // start fetching useful parameters
  COMPRESSION_DATA* compression_dataX = data->get_compression_dataX();
  const VEC3I& dims = compression_dataX->get_dims();
  int numCols = compression_dataX->get_numCols();

  // resize if necessary
  submatrix->resize(3, numCols);

  // compute the block number and index
  int blockNumber, blockIndex;
  ComputeBlockNumber(startRow, dims, &blockNumber, &blockIndex);

  // compare the cache to the current
  int cachedBlockNumber = data->get_cachedBlockNumber();
  if (blockNumber == cachedBlockNumber) { // if we've already decoded this block
    //TIMER cacheTimer("cached block");

    // cout << "Used cache!" << endl;

    // load the previously decoded data
    vector<FIELD_3D>* cachedBlocksX = data->get_cachedBlocksX();
    vector<FIELD_3D>* cachedBlocksY = data->get_cachedBlocksY();
    vector<FIELD_3D>* cachedBlocksZ = data->get_cachedBlocksZ();

    // SVD transformation matrix container
    Matrix3d V_i;

    for (int i = 0; i < numCols; i++) {

      // a dummy container for each column of the matrix
      Vector3d col_i;
      col_i[0] = (*cachedBlocksX)[i][blockIndex];
      col_i[1] = (*cachedBlocksY)[i][blockIndex];
      col_i[2] = (*cachedBlocksZ)[i][blockIndex];

      // set the column
      submatrix->col(i) = col_i;
    }
  return;
  }

  else { // no cache; have to compute it from scratch
    // cout << "Didn't use cache!" << endl;
    //TIMER uncachedTimer("uncached block");

    // initialize useful data
    int numBlocks = compression_dataX->get_numBlocks();
    //const INTEGER_FIELD_3D& zigzagArray = compression_dataX->get_zigzagArray();
    const INTEGER_FIELD_3D& reverseZigzag = compression_dataX->get_reverseZigzag();

    int* allDataX = data->get_dataX();
    int* allDataY = data->get_dataY();
    int* allDataZ = data->get_dataZ();
    COMPRESSION_DATA* compression_dataY = data->get_compression_dataY();
    COMPRESSION_DATA* compression_dataZ = data->get_compression_dataZ();

    // container for the encoded blocks
    vector<FIELD_3D> blocks(numBlocks);

    // variable to store decoded run-length blocks
    // TK: this is always the same size, so set it once
    VectorXi runLengthDecoded(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

    // vector prior to undoing the svd
    Vector3d preSVD;

    // data for doing the inverse DCT
    const fftw_plan& plan = data->get_plan();
    double* in = data->get_dct_in();
    int direction = -1;

    // pointers to the cache
    vector<FIELD_3D>* cachedBlocksX = data->get_cachedBlocksX();
    vector<FIELD_3D>* cachedBlocksY = data->get_cachedBlocksY();
    vector<FIELD_3D>* cachedBlocksZ = data->get_cachedBlocksZ();

    INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
    FIELD_3D decodedBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
    NONZERO_ENTRIES nonZeros;
    for (int i = 0; i < numCols; i++)
    {
      // a container for column i
      Vector3d preSVD_i;
      /*
       *****************************************************************
       * X coordinate
       *****************************************************************
      */

      // DEBUG: while the end of the iteration doesn't clear it entirely
      //unflattened.clear();
      decodedBlock.clear();

      nonZeros.clear();
      //vector<int> nonZeros;
      //nonZeros.reserve(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

      // decode the run-length scheme
      RunLengthDecodeBinaryInPlaceSparse(allDataX, blockNumber, i, reverseZigzag, compression_dataX, unflattened, nonZeros);

      // undo the scaling from the quantizer
      //DecodeBlockWithCompressionDataSparse(unflattened, blockNumber, i, compression_dataX, decodedBlock.data(), nonZeros);
      DecodeBlockWithCompressionDataSparseQuantized(unflattened, blockNumber, i, compression_dataX, decodedBlock.data(), nonZeros);

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[0] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksX)[i] = decodedBlock;

      // clear things out for the next coordinate
      unflattened.clear(nonZeros.data());

      // TK: Why doesn't this clear do what it's supposed to?
      //decodedBlock.clear(nonZeros);
      //unflattened.clear();
      decodedBlock.clear();
      nonZeros.clear();
      //nonZeros.reserve(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

      /*
       *****************************************************************
       * Y coordinate
       *****************************************************************
      */

      // decode the run-length scheme
      //RunLengthDecodeBinary(allDataY, blockNumber, i, compression_dataY, &runLengthDecoded);
      //ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);
      RunLengthDecodeBinaryInPlaceSparse(allDataY, blockNumber, i,
          reverseZigzag, compression_dataY, unflattened, nonZeros);

      // undo the scaling from the quantizer
      // TK: Do it on the raw pointer instead
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataY, decodedBlock.data());
      //DecodeBlockWithCompressionDataSparse(unflattened, blockNumber, i, compression_dataY, decodedBlock.data(), nonZeros);
      DecodeBlockWithCompressionDataSparseQuantized(unflattened, blockNumber, i, compression_dataY, decodedBlock.data(), nonZeros);

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[1] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksY)[i] = decodedBlock;

      // clear things out for the next coordinate
      unflattened.clear(nonZeros.data());

      // TK: Why doesn't this clear do what it's supposed to?
      //decodedBlock.clear(nonZeros);
      //unflattened.clear();
      decodedBlock.clear();
      nonZeros.clear();
      //nonZeros.reserve(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

      /*
       *****************************************************************
       * Z coordinate
       *****************************************************************
      */

      // decode the run-length scheme
      //RunLengthDecodeBinary(allDataZ, blockNumber, i, compression_dataZ, &runLengthDecoded);
      //ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);
      RunLengthDecodeBinaryInPlaceSparse(allDataZ, blockNumber, i,
          reverseZigzag, compression_dataZ, unflattened, nonZeros);

      // undo the scaling from the quantizer
      // TK: Do it on the raw pointer instead
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataZ, &decodedBlock);
      //DecodeBlockWithCompressionData(unflattened, blockNumber, i, compression_dataZ, decodedBlock.data());
      //DecodeBlockWithCompressionDataSparse(unflattened, blockNumber, i, compression_dataZ, decodedBlock.data(), nonZeros);
      DecodeBlockWithCompressionDataSparseQuantized(unflattened, blockNumber, i, compression_dataZ, decodedBlock.data(), nonZeros);

      // inverse DCT
      DCT_Smart_Unitary(plan, direction, in, &decodedBlock);

      preSVD_i[2] = decodedBlock[blockIndex];

      // update the cache
      (*cachedBlocksZ)[i] = decodedBlock;


      // clear things out for the next coordinate
      unflattened.clear(nonZeros.data());

      // TK: Why doesn't this clear do what it's supposed to?
      //decodedBlock.clear(nonZeros);
      //unflattened.clear();
      decodedBlock.clear();
      nonZeros.clear();
      //nonZeros.reserve(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

      // set the column
      submatrix->col(i) = preSVD_i;
    }

    // now that the cache has been filled, set the block number
    data->set_cachedBlockNumber(blockNumber);

  }
}

////////////////////////////////////////////////////////
// generates the header information in the binary file
// need to include:
//
////////////////////////////////////////////////////////
void WriteMetaData(const char* filename, COMPRESSION_DATA& compression_data)
{

  TIMER functionTimer(__FUNCTION__);

    FILE* pFile;
    pFile = fopen(filename, "wb");
    if (pFile == NULL) {
      perror ("Error opening file.");
    }
    else {

      // write nBits to the binary file
      int nBits = compression_data.get_nBits();
      fwrite(&nBits, sizeof(int), 1, pFile);

      // write dims, numCols, and numBlocks
      const VEC3I& dims = compression_data.get_dims();
      int xRes = dims[0];
      int yRes = dims[1];
      int zRes = dims[2];
      fwrite(&xRes, sizeof(int), 1, pFile);
      fwrite(&yRes, sizeof(int), 1, pFile);
      fwrite(&zRes, sizeof(int), 1, pFile);
      int numCols = compression_data.get_numCols();
      int numBlocks = compression_data.get_numBlocks();
      int blocksXcols = numBlocks * numCols;
      fwrite(&numCols, sizeof(int), 1, pFile);
      fwrite(&numBlocks, sizeof(int), 1, pFile);

      MatrixXd* sListMatrix = compression_data.get_sListMatrix();
      assert( sListMatrix->rows() * sListMatrix->cols() == blocksXcols );

      // write the matrix data for sList, blockLengths, and blockIndices.
      // note that Eigen uses column-major format!
      fwrite(sListMatrix->data(), sizeof(double), blocksXcols, pFile);

      MatrixXd* gammaListMatrix = compression_data.get_gammaListMatrix();
      assert( gammaListMatrix->rows() * gammaListMatrix->cols() == blocksXcols );

      fwrite(gammaListMatrix->data(), sizeof(double), blocksXcols, pFile);

      MatrixXi* blockLengthsMatrix = compression_data.get_blockLengthsMatrix();
      assert( blockLengthsMatrix->rows() * blockLengthsMatrix->cols() == blocksXcols);

      fwrite(blockLengthsMatrix->data(), sizeof(int), blocksXcols, pFile);

      MatrixXi* blockIndicesMatrix = compression_data.get_blockIndicesMatrix();
      assert( blockIndicesMatrix->rows() * blockIndicesMatrix->cols() == blocksXcols);

      fwrite(blockIndicesMatrix->data(), sizeof(int), blocksXcols, pFile);

      fclose(pFile);
    }
  }


////////////////////////////////////////////////////////
// concatenate two binary files and put them into
// a new binary file
////////////////////////////////////////////////////////
void PrefixBinary(string prefix, string filename, string newFile) {
  TIMER functionTimer(__FUNCTION__);
  string command = "cat " + prefix + ' ' + filename + "> " + newFile;
  const char* command_c = command.c_str();
  system(command_c);
}


////////////////////////////////////////////////////////
// destroy the no longer needed metadata binary file
////////////////////////////////////////////////////////
void CleanUpPrefix(const char* prefix, const char* filename) {
  TIMER functionTimer(__FUNCTION__);
  string prefix_string(prefix);
  string filename_string(filename);
  string command1 = "rm " + prefix_string;
  string command2 = "rm " + filename_string;
  const char* command1_c = command1.c_str();
  const char* command2_c = command2.c_str();

  system(command1_c);
  system(command2_c);
}


////////////////////////////////////////////////////////
// write the singular values and V matrices to a binary file
////////////////////////////////////////////////////////
void WriteSVDData(const char* filename, COMPRESSION_DATA* data)
{
  FILE* pFile;
  pFile = fopen(filename, "wb");
  if (pFile == NULL) {
    perror ("Error opening file.");
  }
  else {
    vector<Vector3d>* singularList = data->get_singularList();
    vector<Matrix3d>* vList        = data->get_vList();

    int numCols = data->get_numCols();

    assert(singularList->size() == numCols && vList->size() == numCols);

    // for each column, first write the singular values, then write
    // the values in the V matrix. there are 3 singular values,
    // and the V matrices are all 3 x 3.
    for (int col = 0; col < numCols; col++) {
      fwrite((*singularList)[col].data(), sizeof(double), 3, pFile);
      fwrite((*vList)[col].data(), sizeof(double), 3 * 3, pFile);
    }
  }
  fclose(pFile);
}

////////////////////////////////////////////////////////
// compress all of the scalar field components
// of a matrix (which represents a vector field) and write them to
// a binary file. applies svd coordinate transform first.
////////////////////////////////////////////////////////
void CompressAndWriteMatrixComponents(const char* filename, const MatrixXd& U,
      COMPRESSION_DATA* data0, COMPRESSION_DATA* data1, COMPRESSION_DATA* data2)
{
  TIMER functionTimer(__FUNCTION__);

  // wipe any pre-existing binary file of the same name, since we will be opening
  // in append mode!
  DeleteIfExists(filename);

  // write to component X, Y, and Z accordingly
  // initialize strings 0, 1, and 2 for the final result
  string filenameX(filename); filenameX += 'X';
  string filename0(filename); filename0 += '0';
  DeleteIfExists(filenameX.c_str());
  DeleteIfExists(filename0.c_str());

  string filenameY(filename); filenameY += 'Y';
  string filename1(filename); filename1 += '1';
  DeleteIfExists(filenameY.c_str());
  DeleteIfExists(filename1.c_str());

  string filenameZ(filename); filenameZ += 'Z';
  string filename2(filename); filename2 += '2';
  DeleteIfExists(filenameZ.c_str());
  DeleteIfExists(filename2.c_str());

  // get useful compression data. it should be the same across all three,
  // so just fetch it from data0
  const VEC3I& dims = data0->get_dims();
  int xRes = dims[0];
  int yRes = dims[1];
  int zRes = dims[2];
  int numCols = U.cols();

  for (int col = 0; col < numCols; col++) {

    // for each column, grab the vector field that it corresponds to
    VECTOR3_FIELD_3D V(U.col(col), xRes, yRes, zRes);

    // do the svd coordinate transform in place and update the data for
    // vList and singularList. only data0 contains these!
    // TransformVectorFieldSVDCompression(&V, data0);


    // write the components to an (appended) binary file

    CompressAndWriteField(filenameX.c_str(), V.scalarField(0), col, data0);
    CompressAndWriteField(filenameY.c_str(), V.scalarField(1), col, data1);
    CompressAndWriteField(filenameZ.c_str(), V.scalarField(2), col, data2);

    // progress printout for the impatient user
    PrintProgress(col, numCols);
  }

  // once we've gone through each column, we can write the full SVD data
  // WriteSVDData("U.preadvect.SVD.data", data0);

  // we can also build the block indices matrix for each component
  BuildBlockIndicesMatrix(data0);
  BuildBlockIndicesMatrix(data1);
  BuildBlockIndicesMatrix(data2);

  // write the metadata for each component one at a time
  const char* metafile = "metadata.bin";
  WriteMetaData(metafile, *data0);

  // appends the metadata as a header to the main binary file and pipes them into final_string
  PrefixBinary(metafile, filenameX, filename0);

  // removes the now-redundant files
  CleanUpPrefix(metafile, filenameX.c_str());

  // do the same for component 1
  WriteMetaData(metafile, *data1);
  PrefixBinary(metafile, filenameY, filename1);
  CleanUpPrefix(metafile, filenameY.c_str());

  // do the same for component 2
  WriteMetaData(metafile, *data2);
  PrefixBinary(metafile, filenameZ, filename2);
  CleanUpPrefix(metafile, filenameZ.c_str());
}

////////////////////////////////////////////////////////
// compress all of the scalar field components
// of a matrix (which represents a vector field) and write them to
// a binary file. applies svd coordinate transform first.
// uses gamma as zero everywhere for debugging!
////////////////////////////////////////////////////////
void CompressAndWriteMatrixComponentsDebug(const char* filename, const MatrixXd& U,
      COMPRESSION_DATA* data0, COMPRESSION_DATA* data1, COMPRESSION_DATA* data2)
{
  TIMER functionTimer(__FUNCTION__);

  // wipe any pre-existing binary file of the same name, since we will be opening
  // in append mode!
  DeleteIfExists(filename);

  // write to component X, Y, and Z accordingly
  // initialize strings 0, 1, and 2 for the final result
  string filenameX(filename); filenameX += 'X';
  string filename0(filename); filename0 += '0';
  DeleteIfExists(filenameX.c_str());
  DeleteIfExists(filename0.c_str());

  string filenameY(filename); filenameY += 'Y';
  string filename1(filename); filename1 += '1';
  DeleteIfExists(filenameY.c_str());
  DeleteIfExists(filename1.c_str());

  string filenameZ(filename); filenameZ += 'Z';
  string filename2(filename); filename2 += '2';
  DeleteIfExists(filenameZ.c_str());
  DeleteIfExists(filename2.c_str());

  // get useful compression data. it should be the same across all three,
  // so just fetch it from data0
  const VEC3I& dims = data0->get_dims();
  int xRes = dims[0];
  int yRes = dims[1];
  int zRes = dims[2];
  int numCols = U.cols();

  for (int col = 0; col < numCols; col++) {

    // for each column, grab the vector field that it corresponds to
    VECTOR3_FIELD_3D V(U.col(col), xRes, yRes, zRes);

    // do the svd coordinate transform in place and update the data for
    // vList and singularList. only data0 contains these!
    // TransformVectorFieldSVDCompression(&V, data0);


    // write the components to an (appended) binary file

    CompressAndWriteFieldDebug(filenameX.c_str(), V.scalarField(0), col, data0);
    CompressAndWriteFieldDebug(filenameY.c_str(), V.scalarField(1), col, data1);
    CompressAndWriteFieldDebug(filenameZ.c_str(), V.scalarField(2), col, data2);

    // progress printout for the impatient user
    PrintProgress(col, numCols);
  }

  // once we've gone through each column, we can write the full SVD data
  // WriteSVDData("U.preadvect.SVD.data", data0);

  // we can also build the block indices matrix for each component
  BuildBlockIndicesMatrix(data0);
  BuildBlockIndicesMatrix(data1);
  BuildBlockIndicesMatrix(data2);

  // write the metadata for each component one at a time
  const char* metafile = "metadata.bin";
  WriteMetaData(metafile, *data0);

  // appends the metadata as a header to the main binary file and pipes them into final_string
  PrefixBinary(metafile, filenameX, filename0);

  // removes the now-redundant files
  CleanUpPrefix(metafile, filenameX.c_str());

  // do the same for component 1
  WriteMetaData(metafile, *data1);
  PrefixBinary(metafile, filenameY, filename1);
  CleanUpPrefix(metafile, filenameY.c_str());

  // do the same for component 2
  WriteMetaData(metafile, *data2);
  PrefixBinary(metafile, filenameZ, filename2);
  CleanUpPrefix(metafile, filenameZ.c_str());
}

////////////////////////////////////////////////////////
// print four different percents for how far along each
// column we are
////////////////////////////////////////////////////////
void PrintProgress(int col, int numCols)
{
  // percent total progress
  double percent = col / (double) numCols;

  // four checkpoints
  int checkPoint1 = (numCols - 1) / 4;
  int checkPoint2 = (numCols - 1) / 2;
  int checkPoint3 = (3 * (numCols - 1)) / 4;
  int checkPoint4 = numCols - 1;

  // if we're at any of the checkpoints, print the progress
  if (col == checkPoint1 || col == checkPoint2 || col == checkPoint3 || col == checkPoint4) {
    cout << "    Percent: " << percent << flush;
    if (col == checkPoint4) {
      cout << endl;
    }
  }
}


////////////////////////////////////////////////////////
// decode an entire scalar field of a particular column from matrix compression data
////////////////////////////////////////////////////////
void DecodeScalarField(COMPRESSION_DATA* compression_data, int* allData,
    int col, FIELD_3D* decoded)
{
  TIMER functionTimer(__FUNCTION__);

  // get the dims from data and construct a field of the appropriate size
  const VEC3I& dims = compression_data->get_dims();

  // fill in the paddings
  VEC3I newDims;
  GetPaddings(dims, &newDims);

  // update the resolutions
  newDims += dims;

  // fetch numBlocks and precomputed zigzag table
  int numBlocks = compression_data->get_numBlocks();
  const INTEGER_FIELD_3D& zigzagArray = compression_data->get_zigzagArray();

  // container for the encoded blocks
  vector<FIELD_3D> blocks(numBlocks);

  // variable to store decoded run-length blocks
  // TK: This is always the same size, so set it once
  VectorXi runLengthDecoded(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  FIELD_3D decodedBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  for (int blockNumber = 0; blockNumber < numBlocks; blockNumber++) {

    // decode the run-length scheme
    RunLengthDecodeBinary(allData, blockNumber, col,
        compression_data, &runLengthDecoded);

    // undo the zigzag scan
    ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

    // undo the scaling from the quantizer and push to the block container
    // TK: Do it on the raw pointer instead
    //DecodeBlockWithCompressionData(unflattened, blockNumber, col, compression_data, &decodedBlock);
    DecodeBlockWithCompressionData(unflattened, blockNumber, col, compression_data, decodedBlock.data());

    blocks[blockNumber] = decodedBlock;
  }

  // perform the IDCT on each block
  // -1 <-- inverse
  UnitaryBlockDCT(-1, &blocks);

  // reassemble the blocks into one large scalar field
  FIELD_3D padded_result(newDims[0], newDims[1], newDims[2]);
  AssimilateBlocks(newDims, blocks, &padded_result);

  // strip the padding
  *decoded = padded_result.subfield(0, dims[0], 0, dims[1], 0, dims[2]);

}

////////////////////////////////////////////////////////
// decode an entire scalar field of a particular column from matrix compression data
// *without* going back to the spatial domain (or the SVD transform). leave them
// in a list of blocks as well.
////////////////////////////////////////////////////////
void DecodeScalarFieldEigen(COMPRESSION_DATA* compression_data, int* allData,
    int col, vector<VectorXd>* decoded)
{
  TIMER functionTimer(__FUNCTION__);

  // get the dims from data and construct a field of the appropriate size
  const VEC3I& dims = compression_data->get_dims();

  // fill in the paddings
  VEC3I newDims;
  GetPaddings(dims, &newDims);

  // update the resolutions
  newDims += dims;

  // fetch numBlocks and precomputed zigzag table
  int numBlocks = compression_data->get_numBlocks();
  //const INTEGER_FIELD_3D& zigzagArray = compression_data->get_zigzagArray();
  const INTEGER_FIELD_3D& reverseZigzag = compression_data->get_reverseZigzag();

  //INTEGER_FIELD_3D reverseZigzag(zigzagArray);
  //for (int x = 0; x < zigzagArray.totalCells(); x++)
  //  reverseZigzag[zigzagArray[x]] = x;

  // resize the container we will return
  decoded->resize(numBlocks);
  for (int x = 0; x < numBlocks; x++)
    (*decoded)[x].resize(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  // variable to store decoded run-length blocks
  // TK: This is always the same size, so set it once
  VectorXi runLengthDecoded(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  //TIMER loopTimer("DecodeScalar loop timer");

  // TK: Trying making these static so it doesn't allocate and deallocate every
  // time. Might cause problems with OpenMP later.
  static FIELD_3D decodedBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  static INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  for (int blockNumber = 0; blockNumber < numBlocks; blockNumber++)
  {
    unflattened.clear();

    // decode the run length scheme
    //RunLengthDecodeBinary(allData, blockNumber, col,
    //    compression_data, &runLengthDecoded);
    RunLengthDecodeBinaryInPlace(allData, blockNumber, col,
        reverseZigzag, compression_data, unflattened);

    // undo the zigzag scan
    // TK: folded this into the in-place decode above
    //ZigzagUnflatten(runLengthDecoded, zigzagArray, &unflattened);

    /*
    // undo the scaling from the quantizer and push to the block container
    DecodeBlockWithCompressionData(unflattened, blockNumber, col,
        compression_data, &decodedBlock);

    TIMER copyTimer("DecodeScalar copy timer");
    // TK: this is a pretty expensive copy, try it the memcpy way
    //(*decoded)[blockNumber] = decodedBlock.flattenedEigen();
    VectorXd& final = (*decoded)[blockNumber];
    final.resize(decodedBlock.totalCells());
    Real* blockData = decodedBlock.data();
    double* finalData = final.data();
    memcpy(finalData, blockData, sizeof(double) * decodedBlock.totalCells());
    */

    // TK: Do the decode in a way that avoids the need to copy into an Eigen
    // VectorXd at the end
    VectorXd& final = (*decoded)[blockNumber];
    //final.resize(decodedBlock.totalCells());
    DecodeBlockWithCompressionData(unflattened, blockNumber, col,
        compression_data, final.data());
  }
}

////////////////////////////////////////////////////////
// decode an entire scalar field of a particular column from matrix compression data
// *without* going back to the spatial domain (or the SVD transform). leave them
// in a list of blocks as well.
////////////////////////////////////////////////////////
void DecodeScalarFieldEigenSparseStackless(COMPRESSION_DATA* compression_data, int* allData,
    int col, vector<VectorXd>* decoded, vector<NONZERO_ENTRIES>& allNonZeros)
{
  TIMER functionTimer(__FUNCTION__);

  // get the dims from data and construct a field of the appropriate size
  const VEC3I& dims = compression_data->get_dims();

  // fill in the paddings
  VEC3I newDims;
  GetPaddings(dims, &newDims);

  // update the resolutions
  newDims += dims;

  // fetch numBlocks and precomputed zigzag table
  const int numBlocks = compression_data->get_numBlocks();
  //const INTEGER_FIELD_3D& reverseZigzag = compression_data->get_reverseZigzag();

  // TK: Trying making these static so it doesn't allocate and deallocate every
  // time. Might cause problems with OpenMP later.
  static INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);

  // store globally once
  _compression_data = compression_data;
  _allData = allData;
  _col = col;
  _unflattened = &unflattened;
  _allNonZeros = &allNonZeros;
  _decoded = decoded;

  NONZERO_ENTRIES* nonZeros = &(allNonZeros[0]);
  for (int blockNumber = 0; blockNumber < numBlocks; blockNumber++)
  {
    _blockNumber = blockNumber;

    //NONZERO_ENTRIES* nonZeros = allNonZeros[blockNumber];
    nonZeros = &(allNonZeros[blockNumber]);
    //nonZeros.clear();
    nonZeros->clear();

    // decode the run length scheme
    //RunLengthDecodeBinaryInPlaceSparse(allData, blockNumber, col,
    //    reverseZigzag, compression_data, unflattened, nonZeros);
    RunLengthDecodeBinaryInPlaceSparseStackless();

    // TK: Do the decode in a way that avoids the need to copy into an Eigen
    // VectorXd at the end
    //VectorXd& final = (*decoded)[blockNumber];
    //DecodeBlockWithCompressionDataSparse(unflattened, blockNumber, col,
    //    compression_data, final.data(), nonZeros);
    DecodeBlockWithCompressionDataSparseStackless();

    //unflattened.clear(nonZeros.data(), nonZeros.size());
    unflattened.clear(nonZeros->data(), nonZeros->size());
  }
}

////////////////////////////////////////////////////////
// decode an entire scalar field of a particular column from matrix compression data
// *without* going back to the spatial domain (or the SVD transform). leave them
// in a list of blocks as well.
////////////////////////////////////////////////////////
void DecodeScalarFieldEigenSparse(COMPRESSION_DATA* compression_data, int* allData,
    int col, vector<VectorXd>* decoded, vector<NONZERO_ENTRIES>& allNonZeros)
{
  TIMER functionTimer(__FUNCTION__);

  // get the dims from data and construct a field of the appropriate size
  const VEC3I& dims = compression_data->get_dims();

  // fill in the paddings
  VEC3I newDims;
  GetPaddings(dims, &newDims);

  // update the resolutions
  newDims += dims;

  // fetch numBlocks and precomputed zigzag table
  const int numBlocks = compression_data->get_numBlocks();
  const INTEGER_FIELD_3D& reverseZigzag = compression_data->get_reverseZigzag();

  // TK: Trying making these static so it doesn't allocate and deallocate every
  // time. Might cause problems with OpenMP later.
  static INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);

  //TIMER functionTimer2("DecodeScalarFieldEigenSparse, block 2");
  //NONZERO_ENTRIES& nonZeros = allNonZeros[0];
  //VectorXd& final = (*decoded)[0];

  //for (int blockNumber = 0; blockNumber < numBlocks; blockNumber++)
  for (int blockNumber = 0; blockNumber < numBlocks; ++blockNumber)
  {
    NONZERO_ENTRIES& nonZeros = allNonZeros[blockNumber];
    nonZeros.clear();

    // decode the run length scheme
    RunLengthDecodeBinaryInPlaceSparse(allData, blockNumber, col,
        reverseZigzag, compression_data, unflattened, nonZeros);

    // TK: Do the decode in a way that avoids the need to copy into an Eigen
    // VectorXd at the end
    VectorXd& final = (*decoded)[blockNumber];
    /*DecodeBlockWithCompressionDataSparse(unflattened, blockNumber, col,
        compression_data, final.data(), nonZeros);*/
     DecodeBlockWithCompressionDataSparseQuantized(unflattened, blockNumber, col,
        compression_data, final.data(), nonZeros);
    unflattened.clear(nonZeros.data(), nonZeros.size());
  }
}

////////////////////////////////////////////////////////
// decode an entire scalar field of a particular column from matrix compression data
// *without* going back to the spatial domain (or the SVD transform). leave them
// in a list of blocks as well.
////////////////////////////////////////////////////////
void DecodeScalarFieldEigenSparse(COMPRESSION_DATA* compression_data, int* allData,
    int col, vector<VectorXd>* decoded)
{
  TIMER functionTimer(__FUNCTION__);

  // get the dims from data and construct a field of the appropriate size
  const VEC3I& dims = compression_data->get_dims();

  // fill in the paddings
  VEC3I newDims;
  GetPaddings(dims, &newDims);

  // update the resolutions
  newDims += dims;

  // fetch numBlocks and precomputed zigzag table
  const int numBlocks = compression_data->get_numBlocks();
  const INTEGER_FIELD_3D& reverseZigzag = compression_data->get_reverseZigzag();

  // resize the container we will return
  //decoded->resize(numBlocks);
  for (int x = 0; x < numBlocks; x++)
  {
    //(*decoded)[x].resize(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);
    (*decoded)[x].setZero();
  }

  // TK: Trying making these static so it doesn't allocate and deallocate every
  // time. Might cause problems with OpenMP later.
  //static FIELD_3D decodedBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  static INTEGER_FIELD_3D unflattened(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  static NONZERO_ENTRIES nonZeros;
  for (int blockNumber = 0; blockNumber < numBlocks; blockNumber++)
  {
    //vector<int> nonZeros;
    //nonZeros.reserve(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);
    nonZeros.clear();

    // decode the run length scheme
    RunLengthDecodeBinaryInPlaceSparse(allData, blockNumber, col,
        reverseZigzag, compression_data, unflattened, nonZeros);

    // TK: Do the decode in a way that avoids the need to copy into an Eigen
    // VectorXd at the end
    VectorXd& final = (*decoded)[blockNumber];
    /*DecodeBlockWithCompressionDataSparse(unflattened, blockNumber, col,
        compression_data, final.data(), nonZeros);*/
   DecodeBlockWithCompressionDataSparseQuantized(unflattened, blockNumber, col,
        compression_data, final.data(), nonZeros);
    //DecodeBlockWithCompressionData(unflattened, blockNumber, col,
    //    compression_data, final.data());
    unflattened.clear(nonZeros.data(), nonZeros.size());
  }
}


////////////////////////////////////////////////////////
// uses DecodeScalarField three times to reconstruct
// a lossy vector field
////////////////////////////////////////////////////////

void DecodeVectorField(MATRIX_COMPRESSION_DATA* data, int col,
    VECTOR3_FIELD_3D* decoded)
{
  TIMER functionTimer(__FUNCTION__);

  COMPRESSION_DATA* compression_dataX = data->get_compression_dataX();
  COMPRESSION_DATA* compression_dataY = data->get_compression_dataY();
  COMPRESSION_DATA* compression_dataZ = data->get_compression_dataZ();

  const VEC3I& dims = compression_dataX->get_dims();
  int xRes = dims[0];
  int yRes = dims[1];
  int zRes = dims[2];

  int* allDataX = data->get_dataX();
  int* allDataY = data->get_dataY();
  int* allDataZ = data->get_dataZ();

  FIELD_3D scalarX, scalarY, scalarZ;
  DecodeScalarField(compression_dataX, allDataX, col, &scalarX);
  DecodeScalarField(compression_dataY, allDataY, col, &scalarY);
  DecodeScalarField(compression_dataZ, allDataZ, col, &scalarZ);

  // copy the data into the resuting vector field data structure
  (*decoded) = VECTOR3_FIELD_3D(scalarX.data(), scalarY.data(), scalarZ.data(),
      xRes, yRes, zRes);

  // vector<Matrix3d>* vList = compression_dataX->get_vList();
  // UntransformVectorFieldSVD((*vList)[col], decoded);

}

////////////////////////////////////////////////////////
// reconstruct an entire (lossy) matrix
////////////////////////////////////////////////////////
void DecodeMatrix(MATRIX_COMPRESSION_DATA* data, MatrixXd* decoded)
{
  TIMER functionTimer(__FUNCTION__);

  COMPRESSION_DATA* compression_dataX = data->get_compression_dataX();
  const VEC3I& dims = compression_dataX->get_dims();
  int numRows = 3 * dims[0] * dims[1] * dims[2];
  int numCols = compression_dataX->get_numCols();
  decoded->resize(numRows, numCols);

  VECTOR3_FIELD_3D decodedVecfield;
  for (int i = 0; i < numCols; i++) {
    DecodeVectorField(data, i, &decodedVecfield);
    decoded->col(i) = decodedVecfield.flattenedEigen();
  }

}



//////////////////////////////////////////////////////////////////////
// unproject the reduced coordinate into the peeled cells in this field
// using compression data
//////////////////////////////////////////////////////////////////////

void PeeledCompressedUnproject(MATRIX_COMPRESSION_DATA* U_data, const VectorXd& q,
    VECTOR3_FIELD_3D* V)
{
  TIMER functionTimer(__FUNCTION__);

  int xRes = V->xRes();
  int yRes = V->yRes();
  int zRes = V->zRes();
  COMPRESSION_DATA* dataX = U_data->get_compression_dataX();
  int totalColumns = dataX->get_numCols();
  const VEC3I& dims = dataX->get_dims();

  // verify that the (peeled) dimensions match
  void GetPaddings(const VEC3I& v, VEC3I* paddings);
  VEC3I paddings;
  GetPaddings(dims, &paddings);
  paddings += dims;
  assert( xRes == paddings[0] && yRes == paddings[1] && zRes == paddings[2] );
  assert ( totalColumns == q.size() );

  // times 3 since it is a vec3 field
  const int numRows = 3 * dims[0] * dims[1] * dims[2];

  VectorXd result(numRows);
  result.setZero();
  VECTOR3_FIELD_3D decodedVecField;

  for (int col = 0; col < totalColumns; col++) {
    DecodeVectorField(U_data, col, &decodedVecField);
    result += ( q[col] * decodedVecField.flattenedEigen() );
  }

  V->setWithPeeled(result);

}

//////////////////////////////////////////////////////////////////////
// scale a vector<VectorXd> each by the same scalar
//////////////////////////////////////////////////////////////////////
void ScaleVectorEigen(double alpha, vector<VectorXd>* V)
{
  for (int i = 0; i < V->size(); i++) {
    (*V)[i] *= alpha;
  }
}

//////////////////////////////////////////////////////////////////////
// increment each element of a vector<VectorXd> correspondingly
// by another vector<VectorXd>
//////////////////////////////////////////////////////////////////////
void AddVectorEigen(const vector<VectorXd> V, vector<VectorXd>* W)
{
  assert( V.size() == W->size() );
  for (int i = 0; i < V.size(); i++) {
    (*W)[i] += V[i];
  }
}

//////////////////////////////////////////////////////////////////////
// unproject the reduced coordinate into the peeled cells in this field
// using compression data. stays in the frequency domain until the end
//////////////////////////////////////////////////////////////////////
void PeeledCompressedUnprojectTransform(MATRIX_COMPRESSION_DATA* U_data, const VectorXd& q,
    VECTOR3_FIELD_3D* V)
{
  TIMER functionTimer(__FUNCTION__);

  // fetch the relevant data
  int xRes = V->xRes();
  int yRes = V->yRes();
  int zRes = V->zRes();
  COMPRESSION_DATA* dataX = U_data->get_compression_dataX();
  int totalColumns = dataX->get_numCols();
  int numBlocks = dataX->get_numBlocks();
  const VEC3I& dims = dataX->get_dims();
  int* allDataX = U_data->get_dataX();
  COMPRESSION_DATA* dataY = U_data->get_compression_dataY();
  int* allDataY = U_data->get_dataY();
  COMPRESSION_DATA* dataZ = U_data->get_compression_dataZ();
  int* allDataZ = U_data->get_dataZ();

  // verify that the (peeled) dimensions match
  //void GetPaddings(const VEC3I& v, VEC3I* paddings);
  VEC3I paddings;
  GetPaddings(dims, &paddings);
  paddings += dims;

  // TK: Trying a different res. This seems quite restrictive?
  //assert( xRes == paddings[0] && yRes == paddings[1] && zRes == paddings[2] );
  assert ( totalColumns == q.size() );

  // a dummy container for each decoded column
  vector<VectorXd> blocks;

  // size it once and for all
  blocks.resize(numBlocks);
  for (int x = 0; x < numBlocks; x++)
  {
    blocks[x].resize(BLOCK_SIZE * BLOCK_SIZE *BLOCK_SIZE);
    blocks[x].setZero();
  }

  // track where all the nonzeros are
  vector<NONZERO_ENTRIES> allNonZeros(numBlocks);

  // three containers for the x, y, and z components
  vector<VectorXd> resultX(numBlocks);
  vector<VectorXd> resultY(numBlocks);
  vector<VectorXd> resultZ(numBlocks);

  // initialize the result blocks to contain all zeros
  for (int i = 0; i < numBlocks; i++) {
    resultX[i].setZero(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);
    resultY[i].setZero(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);
    resultZ[i].setZero(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);
  }

  TIMER axpyTimer("Unproject core");
  // loop through the columns and interpret the matrix-vector multiply as
  // a linear combination of the columns
  for (int col = 0; col < totalColumns; col++) {

    // decode the data but stay in the freq domain
    //DecodeScalarFieldEigenSparse(dataX, allDataX, col, &blocks);
    DecodeScalarFieldEigenSparse(dataX, allDataX, col, &blocks, allNonZeros);
    //DecodeScalarFieldEigenSparseStackless(dataX, allDataX, col, &blocks, allNonZeros);

    // TK: give Eigen the chance to fuse the axpy
    // Most architectures have a special instruction that does a single multiply and add
    // in a single clock cycle, so you should try to apply them simultaneously whenever
    // you can so that the compiler can try to fuse the two ops
#define SPARSE_AXPY 1
#if !SPARSE_AXPY
    for (int x = 0; x < numBlocks; x++)
      resultX[x] += q[col] * blocks[x];
#else
    const Real qCurrent = q[col];
    for (int x = 0; x < numBlocks; x++)
    {
      const NONZERO_ENTRIES& nonZeros = allNonZeros[x];
      const VectorXd& block = blocks[x];
      VectorXd& result = resultX[x];
      for (int y = 0; y < nonZeros.size(); y++)
      {
        const int i = nonZeros[y];
        result[i] += qCurrent * block[i];
      }
    }
#endif
    clearNonZeros(blocks, allNonZeros);

    //DecodeScalarFieldEigenSparse(dataY, allDataY, col, &blocks);
    DecodeScalarFieldEigenSparse(dataY, allDataY, col, &blocks, allNonZeros);
    //DecodeScalarFieldEigenSparseStackless(dataY, allDataY, col, &blocks, allNonZeros);

    // TK:  give Eigen the chance to fuse the axpy
#if !SPARSE_AXPY
    for (int x = 0; x < numBlocks; x++)
      resultY[x] += q[col] * blocks[x];
#else
    for (int x = 0; x < numBlocks; x++)
    {
      const NONZERO_ENTRIES& nonZeros = allNonZeros[x];
      const VectorXd& block = blocks[x];
      VectorXd& result = resultY[x];
      for (int y = 0; y < nonZeros.size(); y++)
      {
        const int i = nonZeros[y];
        result[i] += qCurrent * block[i];
      }
    }
#endif
    clearNonZeros(blocks, allNonZeros);

    //DecodeScalarFieldEigenSparse(dataZ, allDataZ, col, &blocks);
    DecodeScalarFieldEigenSparse(dataZ, allDataZ, col, &blocks, allNonZeros);
    //DecodeScalarFieldEigenSparseStackless(dataZ, allDataZ, col, &blocks, allNonZeros);

    // give Eigen the chance to fuse the axpy
#if !SPARSE_AXPY
    for (int x = 0; x < numBlocks; x++)
      resultZ[x] += q[col] * blocks[x];
#else
    for (int x = 0; x < numBlocks; x++)
    {
      const NONZERO_ENTRIES& nonZeros = allNonZeros[x];
      const VectorXd& block = blocks[x];
      VectorXd& result = resultZ[x];
      for (int y = 0; y < nonZeros.size(); y++)
      {
        const int i = nonZeros[y];
        result[i] += qCurrent * block[i];
      }
    }
#endif
    clearNonZeros(blocks, allNonZeros);
  }

  // TK: lump the rest of the timing into the function's timer
  TIMER functionTimer2(__FUNCTION__);

  // now go back to the spatial domain
  int direction = -1;
  UnitaryBlockDCTEigen(direction, &resultX);
  UnitaryBlockDCTEigen(direction, &resultY);
  UnitaryBlockDCTEigen(direction, &resultZ);

  // assimilate the blocks into the three component scalar fields
  FIELD_3D Xpart(paddings[0], paddings[1], paddings[2]);
  FIELD_3D Ypart(paddings[0], paddings[1], paddings[2]);
  FIELD_3D Zpart(paddings[0], paddings[1], paddings[2]);
  AssimilateBlocksEigen(paddings, &resultX, &Xpart);
  AssimilateBlocksEigen(paddings, &resultY, &Ypart);
  AssimilateBlocksEigen(paddings, &resultZ, &Zpart);

  // form the corresponding vector field from the peeled scalar fields
  VECTOR3_FIELD_3D result(Xpart.subfield(0, dims[0], 0, dims[1], 0, dims[2]).data(),
      Ypart.subfield(0, dims[0], 0, dims[1], 0, dims[2]).data(),
      Zpart.subfield(0, dims[0], 0, dims[1], 0, dims[2]).data(), dims[0], dims[1], dims[2]);

  // set the result
  V->setWithPeeled(result.flattenedEigen());

}

//////////////////////////////////////////////////////////////////////
// compute the block-wise dot product between two lists and sum them into one
// large dot product
//////////////////////////////////////////////////////////////////////
double GetDotProductSumSparse(const vector<VectorXd>& Vlist, const vector<VectorXd>& Wlist, const vector<NONZERO_ENTRIES>& allNonZeros)
{
  assert(Vlist.size() == Wlist.size());

  int size = Vlist.size();
  double dotProductSum = 0.0;

  for (int i = 0; i < size; i++)
  {
    double dotProduct_i = 0;
    const VectorXd& V = Vlist[i];
    const VectorXd& W = Wlist[i];
    const NONZERO_ENTRIES& nonZeros = allNonZeros[i];

    for (int x = 0; x < nonZeros.size(); x++)
    {
      const int y = nonZeros[x];
      dotProduct_i += V[y] * W[y];
    }
    dotProductSum += dotProduct_i;
  }

  return dotProductSum;

}


//////////////////////////////////////////////////////////////////////
// compute the block-wise dot product between two lists and sum them into one
// large dot product
//////////////////////////////////////////////////////////////////////
double GetDotProductSum(const vector<VectorXd>& Vlist, const vector<VectorXd>& Wlist) {

  assert(Vlist.size() == Wlist.size());

  int size = Vlist.size();

  double dotProductSum = 0.0;
  for (int i = 0; i < size; i++) {
    double dotProduct_i = Vlist[i].dot(Wlist[i]);
    dotProductSum += dotProduct_i;
  }

  return dotProductSum;

}


//////////////////////////////////////////////////////////////////////
// do peeled compressed projection naively in the regular spatial domain
//////////////////////////////////////////////////////////////////////
void PeeledCompressedProject(const VECTOR3_FIELD_3D& V, MATRIX_COMPRESSION_DATA* U_data,
    VectorXd* q)
{
  TIMER functionTimer(__FUNCTION__);

  COMPRESSION_DATA* dataX = U_data->get_compression_dataX();
  const int numCols = dataX->get_numCols();

  // preallocate result to the appropriate size
  q->setZero(numCols);

  for (int col = 0; col < numCols; col++) {
    TIMER columnLoopTimer("peeled project column loop");

    VECTOR3_FIELD_3D decoded_i;
    DecodeVectorField(U_data, col, &decoded_i);
    (*q)[col] = decoded_i.flattenedEigen().dot( V.peelBoundary().flattenedEigen() );

  }

}


//////////////////////////////////////////////////////////////////////
// helper function for frequency domain projection. transforms
// V first by the SVD and then DCT. fills up the three vectors
// with the three components.
//////////////////////////////////////////////////////////////////////
void TransformSVDAndDCT(int col, const VECTOR3_FIELD_3D& V,
    MATRIX_COMPRESSION_DATA* U_data,
    vector<VectorXd>* Xpart, vector<VectorXd>* Ypart, vector<VectorXd>* Zpart)
{
  COMPRESSION_DATA* dataX = U_data->get_compression_dataX();
  vector<Matrix3d>* vList = dataX->get_vList();

  // make a copy of V to write over since V is read-only
  VECTOR3_FIELD_3D Vcopy = V;
  TransformVectorFieldSVDCached(&(*vList)[col], &Vcopy);

  FIELD_3D V_X, V_Y, V_Z;
  GetScalarFields(Vcopy.peelBoundary(), &V_X, &V_Y, &V_Z);

  GetBlocksEigen(V_X, Xpart);
  GetBlocksEigen(V_Y, Ypart);
  GetBlocksEigen(V_Z, Zpart);

  UnitaryBlockDCTEigen(1, Xpart);
  UnitaryBlockDCTEigen(1, Ypart);
  UnitaryBlockDCTEigen(1, Zpart);

}

//////////////////////////////////////////////////////////////////////
// helper function for frequency domain projection. transforms
// V by the DCT. fills up the three vectors
// with the three components.
//////////////////////////////////////////////////////////////////////
void TransformDCT(const VECTOR3_FIELD_3D& V,
    MATRIX_COMPRESSION_DATA* U_data, vector<VectorXd>* Xpart, vector<VectorXd>* Ypart,
    vector<VectorXd>* Zpart)
{
  TIMER functionTimer(__FUNCTION__);

  FIELD_3D V_X, V_Y, V_Z;
  GetScalarFields(V.peelBoundary(), &V_X, &V_Y, &V_Z);

  GetBlocksEigen(V_X, Xpart);
  GetBlocksEigen(V_Y, Ypart);
  GetBlocksEigen(V_Z, Zpart);

  UnitaryBlockDCTEigen(1, Xpart);
  UnitaryBlockDCTEigen(1, Ypart);
  UnitaryBlockDCTEigen(1, Zpart);

}
//////////////////////////////////////////////////////////////////////
// projection, implemented in the frequency domain
//////////////////////////////////////////////////////////////////////
void PeeledCompressedProjectTransform(const VECTOR3_FIELD_3D& V,
    MATRIX_COMPRESSION_DATA* U_data, VectorXd* q)
{
  TIMER functionTimer(__FUNCTION__);

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  // fetch the compression data and the full data buffer for each component
  COMPRESSION_DATA* dataX = U_data->get_compression_dataX();
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  int* allDataX = U_data->get_dataX();
  COMPRESSION_DATA* dataY = U_data->get_compression_dataY();
  int* allDataY = U_data->get_dataY();
  COMPRESSION_DATA* dataZ = U_data->get_compression_dataZ();
  int* allDataZ = U_data->get_dataZ();

  // preallocate the resulting vector
  int totalColumns = dataX->get_numCols();
  cout << "total columns: " << totalColumns << endl;
  q->resize(totalColumns);

  vector<VectorXd> Xpart, Ypart, Zpart;
  vector<VectorXd> blocks;


  for (int col = 0; col < totalColumns; col++) {
    // transform V with an SVD and a DCT to mimic compression

    TransformSVDAndDCT(col, V, U_data, &Xpart, &Ypart, &Zpart);
    double totalSum = 0.0;
    DecodeScalarFieldEigen(dataX, allDataX, col, &blocks);
    totalSum += GetDotProductSum(blocks, Xpart);

    DecodeScalarFieldEigen(dataY, allDataY, col, &blocks);
    totalSum += GetDotProductSum(blocks, Ypart);

    DecodeScalarFieldEigen(dataZ, allDataZ, col, &blocks);
    totalSum += GetDotProductSum(blocks, Zpart);

    (*q)[col] = totalSum;

  }

}

//////////////////////////////////////////////////////////////////////
// projection, implemented in the frequency domain. assumes no SVD!
//////////////////////////////////////////////////////////////////////
void PeeledCompressedProjectTransformNoSVD(const VECTOR3_FIELD_3D& V,
    MATRIX_COMPRESSION_DATA* U_data, VectorXd* q)
{
  TIMER functionTimer(__FUNCTION__);

  // fetch the compression data and the full data buffer for each component
  COMPRESSION_DATA* dataX = U_data->get_compression_dataX();
  int* allDataX = U_data->get_dataX();
  COMPRESSION_DATA* dataY = U_data->get_compression_dataY();
  int* allDataY = U_data->get_dataY();
  COMPRESSION_DATA* dataZ = U_data->get_compression_dataZ();
  int* allDataZ = U_data->get_dataZ();

  // preallocate the resulting vector
  int totalColumns = dataX->get_numCols();
  q->resize(totalColumns);

  vector<VectorXd> Xpart, Ypart, Zpart;

  // size the dummy blocks once and for all
  const int numBlocks = dataX->get_numBlocks();
  vector<VectorXd> blocks(numBlocks);
  for (int x = 0; x < numBlocks; x++)
  {
    blocks[x].resize(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);
    blocks[x].setZero();
  }

  // track where all the nonzeros are
  vector<NONZERO_ENTRIES> allNonZeros(numBlocks);

  TransformDCT(V, U_data, &Xpart, &Ypart, &Zpart);

  TIMER functionTimer2("Project core");
  for (int col = 0; col < totalColumns; col++)
  {
    double totalSum = 0.0;
    //DecodeScalarFieldEigen(dataX, allDataX, col, &blocks);
    //DecodeScalarFieldEigenSparse(dataX, allDataX, col, &blocks);
    DecodeScalarFieldEigenSparse(dataX, allDataX, col, &blocks, allNonZeros);
    //DecodeScalarFieldEigenSparseStackless(dataX, allDataX, col, &blocks, allNonZeros);
    //totalSum += GetDotProductSum(blocks, Xpart);
    totalSum += GetDotProductSumSparse(blocks, Xpart, allNonZeros);
    clearNonZeros(blocks, allNonZeros);

    //DecodeScalarFieldEigen(dataY, allDataY, col, &blocks);
    //DecodeScalarFieldEigenSparse(dataY, allDataY, col, &blocks);
    DecodeScalarFieldEigenSparse(dataY, allDataY, col, &blocks, allNonZeros);
    //DecodeScalarFieldEigenSparseStackless(dataY, allDataY, col, &blocks, allNonZeros);
    //totalSum += GetDotProductSum(blocks, Ypart);
    totalSum += GetDotProductSumSparse(blocks, Ypart, allNonZeros);
    clearNonZeros(blocks, allNonZeros);

    //DecodeScalarFieldEigen(dataZ, allDataZ, col, &blocks);
    //DecodeScalarFieldEigenSparse(dataZ, allDataZ, col, &blocks);
    DecodeScalarFieldEigenSparse(dataZ, allDataZ, col, &blocks, allNonZeros);
    //DecodeScalarFieldEigenSparseStackless(dataZ, allDataZ, col, &blocks, allNonZeros);
    //totalSum += GetDotProductSum(blocks, Zpart);
    totalSum += GetDotProductSumSparse(blocks, Zpart, allNonZeros);
    clearNonZeros(blocks, allNonZeros);

    (*q)[col] = totalSum;
  }
}

//////////////////////////////////////////////////////////////////////
// set zeros at the places where we have artificially padded
//////////////////////////////////////////////////////////////////////
void SetZeroPadding(vector<VectorXd>* blocks, COMPRESSION_DATA* data)
{
  TIMER functionTimer(__FUNCTION__);

  const VEC3I& paddedDims = data->get_paddedDims();
  int xPadded = paddedDims[0];
  int yPadded = paddedDims[1];
  int zPadded = paddedDims[2];

  const VEC3I& dims = data->get_dims();
  int xRes = dims[0];
  int yRes = dims[1];
  int zRes = dims[2];

  FIELD_3D assimilated(xPadded, yPadded, zPadded);
  AssimilateBlocksEigen(paddedDims, blocks, &assimilated);
  FIELD_3D sub = assimilated.subfield(0, xRes, 0, yRes, 0, zRes);

  VEC3I paddings = paddedDims - dims;
  FIELD_3D zeroPadded = sub.zeroPad_xyz(paddings);
  GetBlocksEigen(zeroPadded, blocks);

}

//////////////////////////////////////////////////////////////////////
// clear out all the non-zeros from the passed in blocks
//////////////////////////////////////////////////////////////////////
void clearNonZeros(vector<VectorXd>& blocks, const vector<NONZERO_ENTRIES>& allNonZeros)
{
  //TIMER functionTimer(__FUNCTION__);
  assert(blocks.size() == allNonZeros.size());

  for (unsigned int x = 0; x < blocks.size(); x++)
  {
    VectorXd& block = blocks[x];
    const NONZERO_ENTRIES& nonZeros = allNonZeros[x];

    //block.setZero();
    //continue;

    for (int i = 0; i < nonZeros.size(); i++)
    {
      assert(nonZeros[i] < block.size());
      block[nonZeros[i]] = 0.0;
    }
  }
}
