/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
// SPARSE_PETSC_MATRIX.h: interface for the SPARSE_PETSC_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#include "SPARSE_PETSC_MATRIX.h"
#include <TIMER.h>

#if USING_PETSC
#include <slepceps.h>
#include <slepcsvd.h>
#endif

bool SPARSE_PETSC_MATRIX::singleton = true;

//////////////////////////////////////////////////////////////////////
// Constructor for the sparse matrix
//////////////////////////////////////////////////////////////////////
SPARSE_PETSC_MATRIX::SPARSE_PETSC_MATRIX() :
  _rows(0), 
  _cols(0),
  _firstSolve(true), 
  _finalized(false),
  _iterations(1000000),
  _eps(1e-8), 
  _totalTime(0.0),
  _whichPreconditioner(JACOBI),
  _whichSolver(MINRES)
{
#if USING_PETSC
  if (singleton)
  {
    PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
    SlepcInitialize(NULL,NULL,NULL,NULL);
    singleton = false;
  }
#endif
}

SPARSE_PETSC_MATRIX::SPARSE_PETSC_MATRIX(int rows, int cols) :
  _rows(rows), 
  _cols(cols),
  _firstSolve(true), 
  _iterations(1000000), 
  _totalTime(0.0),
  _whichPreconditioner(JACOBI),
  _whichSolver(MINRES)
{
}

SPARSE_PETSC_MATRIX::SPARSE_PETSC_MATRIX(MATRIX& A) :
  _firstSolve(true), 
  _iterations(1000000), 
  _totalTime(0.0),
  _whichPreconditioner(JACOBI),
  _whichSolver(MINRES)
{
#if USING_PETSC
  if (singleton)
  {
    PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
    SlepcInitialize(NULL,NULL,NULL,NULL);
    singleton = false;
  }

  _rows = A.rows();
  _cols = A.cols();

  // set sparsity to dense
  PetscInt* nonZeros = new PetscInt[_rows];
  int maxNonZeros = _cols;
  for (int x = 0; x < _rows; x++)
    nonZeros[x] = _cols;

  // create A
  MatCreateSeqAIJ(PETSC_COMM_WORLD,_rows, _cols, maxNonZeros, nonZeros, &_APetsc);
  MatSetFromOptions(_APetsc);

  // copy values into sparse form
  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _cols; y++)
      set(x,y, A(x,y));

  delete[] nonZeros;
#endif
}

SPARSE_PETSC_MATRIX::~SPARSE_PETSC_MATRIX()
{
#if USING_PETSC
  /*
  if (!singleton)
  {
    PetscFinalize();
    singleton = true;
  }
  */

  if (!_firstSolve)
  {
    // Causes seg fault?
    //KSPDestroy(_ksp);
    VecDestroy(_BPetsc);
    VecDestroy(_XPetsc);
    MatDestroy(_APetsc);
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// set the sparsity pattern based on another sparse matrix
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::setSparsity(SPARSE_MATRIX& A)
{
#if USING_PETSC
  _rows = A.rows();
  _cols = A.cols();

  // get preallocation hints for PetSc
  vector<int> nonZerosVec = A.nonZerosPerRow();
  PetscInt* nonZeros = new PetscInt[nonZerosVec.size()];
  int maxNonZeros = 0;
  for (unsigned int x = 0; x < nonZerosVec.size(); x++)
  {
    if (nonZerosVec[x] > maxNonZeros)
      maxNonZeros = nonZerosVec[x];
    nonZeros[x] = nonZerosVec[x];
  }

  // create A
  MatCreateSeqAIJ(PETSC_COMM_WORLD,_rows,_cols,maxNonZeros,nonZeros, &_APetsc);
  MatSetFromOptions(_APetsc);

  // create P
  MatCreateSeqAIJ(PETSC_COMM_WORLD,_rows,_cols,maxNonZeros,nonZeros, &_PPetsc);
  MatSetFromOptions(_PPetsc);

  // keep a copy of the sparsity structure around
  _sparsity.clear();
  const map<pair<int,int>,Real>& data = A.matrix();
  map<pair<int,int>, Real>::const_iterator i = data.begin();
  while (i != data.end())
  {
    _sparsity[i->first] = true;
    i++;
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// solve an SPD system using PCG
//////////////////////////////////////////////////////////////////////
bool SPARSE_PETSC_MATRIX::solveCG(VECTOR& x, VECTOR& b)
{
#if USING_PETSC
  TIMER total;

  PetscInt N = _rows;
  int        i;

  if (_firstSolve)
  {
    TIMER initial;

    /*
    // get preallocation hints for PetSc
    vector<int> nonZerosVec = nonZerosPerRow();
    PetscInt* nonZeros = new PetscInt[nonZerosVec.size()];
    int maxNonZeros = 0;
    for (int x = 0; x < nonZerosVec.size(); x++)
    {
      if (nonZerosVec[x] > maxNonZeros)
        maxNonZeros = nonZerosVec[x];
      nonZeros[x] = nonZerosVec[x];
    }

    // create A
    MatCreateSeqAIJ(PETSC_COMM_WORLD,N,N,maxNonZeros,nonZeros, &_APetsc);
    //MatCreateMPIAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,N,N,&_APetsc);
    MatSetFromOptions(_APetsc);

    // create P
    MatCreateSeqAIJ(PETSC_COMM_WORLD,N,N,maxNonZeros,nonZeros, &_PPetsc);
    MatSetFromOptions(_PPetsc);

    // cache matrix entries
    entries(_rowsIndices, _colsIndices, _entries);
    */

    // duplicate the current matrix as the preconditioner

    // create X and B
    VecCreateSeq(PETSC_COMM_WORLD,N,&_BPetsc);
    //VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, N,&_BPetsc);
    VecDuplicate(_BPetsc,&_XPetsc);

    /*
    // delete preallocation hints for PetSc
    delete[] nonZeros;
    */

    _timingBreakdown["Initial Setup"] += initial.timing();
  }

  TIMER dataCopy;

  // copy new values into B
  for (i=0; i<N; i++)
  {
    double bValue = b[i];
    double xValue = x[i];

    //VecSetValues(_BPetsc,1,&i,&b[i],INSERT_VALUES);
    //VecSetValues(_XPetsc,1,&i,&x[i],INSERT_VALUES);
    VecSetValues(_BPetsc,1,&i,&bValue,INSERT_VALUES);
    VecSetValues(_XPetsc,1,&i,&xValue,INSERT_VALUES);
  }

  /*
  // copy entries into generic format
  cout << " Setting up system ... "; flush(cout);
  for (i = 0; i < _rowsIndices.size(); i++)
  {
    // format is:
    //
    // MatSetValues(matrix, 
    //  number of rows, row index, 
    //  number of cols, col index, 
    //  actual value, INSERT_VALUES);
    PetscScalar 	 v = *(_entries[i]);
    MatSetValues(_APetsc,1,&(_rowsIndices[i]),1,&(_colsIndices[i]),&v,INSERT_VALUES);
  }

  // the first time, construct the preconditioner as well
  if (_firstSolve)
  {
    for (i = 0; i < _rowsIndices.size(); i++)
    {
      PetscScalar 	 v = *(_entries[i]);
      if (_rowsIndices[i] == _colsIndices[i])
        v *= 1.1;
      MatSetValues(_PPetsc,1,&(_rowsIndices[i]),1,&(_colsIndices[i]),&v,INSERT_VALUES);
    }
  }
  flush(cout);
  cout << " done. " << endl; flush(cout);
  */

  _timingBreakdown["A and b copy"] += dataCopy.timing();

  // finalize matrices
  TIMER finalizeTimer;
  if (_firstSolve)
  {
    MatAssemblyBegin(_APetsc,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_APetsc,MAT_FINAL_ASSEMBLY);
    //_PPetsc = MatDuplicate(_APetsc);
    MatDuplicate(_APetsc, MAT_COPY_VALUES, &_PPetsc);
    MatAssemblyBegin(_PPetsc,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_PPetsc,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(_XPetsc);
    VecAssemblyEnd(_XPetsc);
    VecAssemblyBegin(_BPetsc);
    VecAssemblyEnd(_BPetsc);

    KSPCreate(PETSC_COMM_WORLD,&_ksp);
    //KSPSetOperators(_ksp,_APetsc,_APetsc,SAME_NONZERO_PATTERN);
    KSPSetOperators(_ksp,_APetsc,_PPetsc,SAME_NONZERO_PATTERN);

    switch (_whichSolver)
    {
      case PCG:
        KSPSetType(_ksp, KSPCG);
        cout << " Using PCG solver " << endl;
        break;
      case MINRES:
        KSPSetType(_ksp, KSPMINRES);
        cout << " Using MINRES solver " << endl;
        break;
      case GMRES:
        KSPSetType(_ksp, KSPGMRES);
        cout << " Using GMRES solver " << endl;
        break;
      default:
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Invalid Krylov solver setting! " << endl;
        break;
    }
    KSPSetInitialGuessNonzero(_ksp,PETSC_TRUE);

    // create preconditioner
    KSPGetPC(_ksp,&_preconditioner);
    switch (_whichPreconditioner)
    {
      case NONE:
        PCSetType(_preconditioner, PCNONE);
        cout << " Using no preconditioner " << endl;
        break;
      case JACOBI:
        PCSetType(_preconditioner, PCJACOBI);
        cout << " Using Jacobi preconditioner " << endl;
        break;
      case SOR:
        PCSetType(_preconditioner, PCSOR);
        cout << " Using SOR preconditioner " << endl;
        break;
      case ILU:
        PCSetType(_preconditioner, PCILU);
        cout << " Using ILU preconditioner " << endl;
        break;
      case ICC:
        PCSetType(_preconditioner, PCICC);
        cout << " Using ICC preconditioner " << endl;
        break;
      default:
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Invalid preconditioner setting!" << endl;
        break;
    }

    KSPSetFromOptions(_ksp);
    KSPSetUp(_ksp);

    // set tolerances
    PetscReal rtol;
    PetscReal abstol;
    PetscReal dtol;
    PetscInt maxits;
    KSPGetTolerances(_ksp, &rtol, &abstol, &dtol, &maxits);
    //abstol = 1e-8;
    //rtol = 1e-6;
    //abstol = 1e-7;
    //rtol = 1e-5;
    abstol = _eps;
    rtol = _eps;
    //abstol = 1e-15;
    //rtol = 1e-15;
    dtol = 1e8;

    //cout << __FILE__ << " " << __LINE__ << " : " << endl;
    //cout << " DEBUG!!!! HUGE NUMBER OF PCG ITERATIONS" << endl;
    //maxits = 1000000;
    //maxits = 100000;
    maxits = _iterations;
    
    cout << " rtol: " << rtol << endl;
    cout << " abstol: " << abstol << endl;
    cout << " dtol: " << dtol << endl;
    cout << " maxits: " << maxits << endl;
    KSPSetTolerances(_ksp, rtol, abstol, dtol, maxits);

    _firstSolve = false;
  }
  else
  {
    MatAssemblyBegin(_APetsc,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_APetsc,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(_XPetsc);
    VecAssemblyEnd(_XPetsc);
    VecAssemblyBegin(_BPetsc);
    VecAssemblyEnd(_BPetsc);
    KSPSetOperators(_ksp,_APetsc,_APetsc,SAME_PRECONDITIONER);
  }
  _timingBreakdown["Finalize matrices"] += finalizeTimer.timing();

  // solve the system
  TIMER solve;
  int its;
  KSPSolve(_ksp,_BPetsc,_XPetsc);
  KSPGetIterationNumber(_ksp,&its);
  cout << " Converged in " << its << " iterations " << endl;
  //printf("Convergence in %d iterations.\n",(int)its);
  _timingBreakdown["Solver"] += solve.timing();

  PetscReal rnorm;
  KSPGetResidualNorm(_ksp, &rnorm);
  //cout << " residual: " << rnorm << endl;

  // see if it diverged
  KSPConvergedReason reason;
  KSPGetConvergedReason(_ksp, &reason);

  if (reason < 0)
  {
    if (reason != -4)
      cout << " CG Solve diverged!!!" << endl;

    switch (reason) {
      case -10:
        cout << " Matrix is indefinite!" << endl;
        break;
      case -7:
        cout << " Matrix is not symmetric!" << endl;
        break;
      case -5:
        cout << " Breakdown! " << endl;
        break;
      case -4:
        //cout << " Solution is already within tolerance! " << endl;
        cout << " Solution exceeded divegence tolerance! " << endl;
        break;
      case -9:
        cout << " Nan encountered! " << endl;
        break;
      case -8:
        cout << " Preconditioner construction died! " << endl;
        break;
      case -3:
        cout << " Max iteration count exceeded! " << endl;
        break;
      default:
        cout << " Reason: " << reason << endl;
        break;
    }

    //cout << __FILE__ << " " << __LINE__ << " : " << endl;
    //exit(0);
    return false;
  }

  // get final results
  TIMER finalCopy;
  PetscScalar* data;
  VecGetArray(_XPetsc, &data);
  for (int i = 0; i < N; i++)
    x[i] = data[i];
  VecRestoreArray(_XPetsc, &data);
  _timingBreakdown["Final copy"] += finalCopy.timing();

  _totalTime += total.timing();
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " PetSc not supported under Win32!!!!" << endl;
#endif
  return true;
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::axpy(Real scalar, SPARSE_MATRIX& A)
{
  // iterate through all the entries
  //map<int, Real>::const_iterator i;
  //const map<int, Real>& data = A.matrix();
  map<pair<int,int>, Real>::const_iterator i;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (i = data.begin(); i != data.end(); i++)
  {
    int row = i->first.first;
    int col = i->first.second;
    const Real value = i->second;
    add(row, col, scalar * value);
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::eraseRowColumn(int rowCol)
{
  pair<int,int> entry(rowCol, rowCol);

  map<pair<int,int>, bool>::iterator i = _sparsity.find(entry);

  if (i == _sparsity.end())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Row/Col erasure failed!" << endl;
    return;
  }
  while (i->first.first == rowCol)
  {
    int row = i->first.first;
    int col = i->first.second;
    set(row, col, 0.0);
    set(col, row, 0.0);
    i++;
  }

  map<pair<int,int>, bool>::iterator ri = _sparsity.find(entry);
  while (ri->first.first == rowCol)
  {
    int row = ri->first.first;
    int col = ri->first.second;
    set(row, col, 0.0);
    set(col, row, 0.0);
    ri--;
  }
}

//////////////////////////////////////////////////////////////////////
// get a specified number of eigenvalues of a system
//////////////////////////////////////////////////////////////////////
VECTOR SPARSE_PETSC_MATRIX::eigenvalues(bool smallest, int totalEigenvalues)
{
  // if no limit is specified, solve for all the eigenvalues
  if (totalEigenvalues == -1)
    totalEigenvalues = _rows;

  VECTOR values(totalEigenvalues);

#if USING_PETSC
  Mat& A = _APetsc;
  PetscReal   	 error, tol;
  int         	 maxit, i, its;

  // finalize the current matrix
  finalize();

  // Create eigensolver context
  EPS         	 eps;
  EPSCreate(PETSC_COMM_WORLD,&eps);

  // Set operators. In this case, it is a non-generalized eigenvalue problem
  EPSSetOperators(eps,A,PETSC_NULL);
  EPSSetProblemType(eps,EPS_HEP);

  // set the tolerances
  Real tolerance = 1e-2;
  EPSSetTolerances(eps, tolerance, PETSC_DECIDE);

  // Slepc recommends using at least twice as many eigenvectors
  // for the search space
  EPSSetDimensions(eps, totalEigenvalues, PETSC_DECIDE, PETSC_DECIDE);
  if (smallest)
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
  else
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);

  PetscTruth hermitian = PETSC_TRUE;
  EPSIsHermitian(eps, &hermitian);

  // Send Slepc the solver parameters
  EPSSetFromOptions(eps);

  // Solve the eigensystem
  TIMER eigenTimer;
  cout << " Solving eigensystem ... ";
  EPSSolve(eps);
  EPSGetIterationNumber(eps, &its);
  EPSGetTolerances(eps,&tol,&maxit);
  cout << "done in " << its << " of " << maxit << " max iterations " << endl;
  cout << "took " << eigenTimer.timing() / 60.0 << " minutes" << endl;

  // Display solution and clean up
  int totalFound = 0;
  EPSGetConverged(eps, &totalFound);

  if (totalFound < totalEigenvalues)
    cout << __FILE__ << " " << __LINE__ << " : Not enough eigenvalues were found! " << endl;

  // make sure no imaginary components were found
  bool imaginaryFound = false;

  if (totalFound > 0) {
    int currentMode = 0;
    for (i = 0; i < totalEigenvalues; i++) {
      PetscScalar realValue, imagValue;

      // Get converged eigenpairs: i-th eigenvalue is stored in 
      // kr (real part) and ki (imaginary part)
      EPSGetEigenpair(eps, i, &realValue, &imagValue, PETSC_NULL, PETSC_NULL);

      // Compute the relative error associated to each eigenpair
      EPSComputeRelativeError(eps,i,&error);

      Real re = realValue;
      Real im = imagValue;

      // make sure the eigenvalue is non-zero within the
      // solver tolerance before storing, and also make
      // sure we haven't already stored enough eigenvectors
      if (fabs(re) > tolerance && currentMode < totalEigenvalues)
      {
        // store the eigenvalue
        values(currentMode) = re;

        if (im!=0.0) 
          imaginaryFound = true;

        currentMode++;
      }
    }
  }
  
  // Free work space
  EPSDestroy(eps);
  MatDestroy(A);
  SlepcFinalize();

  if (imaginaryFound)
    cout << __FILE__ << " " << __LINE__ << " : Eigenvalue solver found an imaginary component! There must be a problem with the stiffness matrix. " << endl;
#endif

  return values;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::clear() { 
#if USING_PETSC
  MatZeroEntries(_APetsc); 
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::finalize() { 
  MatAssemblyBegin(_APetsc,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(_APetsc,MAT_FINAL_ASSEMBLY);
}

//////////////////////////////////////////////////////////////////////
// set equal to a GMRES_MATRIX
// assumes setSparsity has been called
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::equals(BLOCK_MATRIX& matrix)
{
  // wipe all previous
  clear();
  
  int blockRows = matrix.blockRows();
  int blockCols = matrix.blockCols();

  int currentRow = 0;
  for (int x = 0; x < blockRows; x++)
  {
    int currentCol = 0;
    for (int y = 0; y < blockCols; y++)
    {
      if (matrix.exists(x,y))
      {
        const MATRIX& toCopy = matrix.constEntry(x,y);
        for (int i = 0; i < toCopy.rows(); i++)
          for (int j = 0; j < toCopy.cols(); j++)
            add(currentRow + i, currentCol + j, toCopy(i,j));
      }
      currentCol += matrix.subCols(y);
    }
    currentRow += matrix.subRows(x);
  }
}
//////////////////////////////////////////////////////////////////////
// set equal to a SPARSE_MATRIX
// assumes setSparsity has been called
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::equals(SPARSE_MATRIX& matrix) {
  // wipe all previous
  clear();

  // copy everything into the PetSc matrix
  const map<pair<int,int>, Real>& m = matrix.matrix();
  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = m.begin(); iter != m.end(); iter++)
  {
    const int row = iter->first.first;
    const int col = iter->first.second;
    const Real entry = iter->second;
    add(row, col, entry);
  }
}

//////////////////////////////////////////////////////////////////////
// set equal to a MATRIX
// assumes setSparsity has been called
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::equals(MATRIX& matrix) {
#if USING_PETSC
  _rows = matrix.rows();
  _cols = matrix.cols();

  // set sparsity to dense
  PetscInt* nonZeros = new PetscInt[_rows];
  int maxNonZeros = _cols;
  for (int x = 0; x < _rows; x++)
    nonZeros[x] = _cols;

  // create A
  MatCreateSeqAIJ(PETSC_COMM_WORLD,_rows, _cols, maxNonZeros, nonZeros, &_APetsc);
  MatSetFromOptions(_APetsc);

  // copy values into sparse form
  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _cols; y++)
      set(x,y, matrix(x,y));

  delete[] nonZeros;
#endif
}

//////////////////////////////////////////////////////////////////////
// perform PCA getting the first 'r' principal components
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::PCA(MATRIX& data, int rank, MATRIX& components, VECTOR& values)
{
  data.subtractRowMeans();
  this->equals(data);
  meanNormalizedPCA(rank, components, values);
}

//////////////////////////////////////////////////////////////////////
// perform PCA getting the first 'r' principal components
//
// perform PCA, but the matrix being passed in should already have
// subtracted its mean!
//
// Assumes each column of 'data' is a new snapshot
//////////////////////////////////////////////////////////////////////
void SPARSE_PETSC_MATRIX::meanNormalizedPCA(int rank, MATRIX& components, VECTOR& values)
{
  SPARSE_PETSC_MATRIX& A = *this;

  // finalize the current matrix
  A.finalize();

  // create the SVD solver
  SVD svd;
  SVDCreate(PETSC_COMM_WORLD, &svd);

  // set it to solve the current matrix
  SVDSetOperator(svd, A.petscData());

  // set tolerances
  int maxIterations = 1000000;
  SVDSetTolerances(svd, 1e-8, maxIterations);

  // set number of singular values to solve for
  //SVDSetDimensions(svd, rank, PETSC_DECIDE);
  // SLEPc 3.0
  SVDSetDimensions(svd, rank, PETSC_DECIDE, PETSC_DECIDE);
  SVDSetFromOptions(svd);
  SVDSolve(svd);

  int its, maxit;
  PetscReal tol;
  SVDGetIterationNumber(svd, &its);
  SVDGetTolerances(svd,&tol,&maxit);

  // see how many triplets converged
  PetscInt numConverged;
  SVDGetConverged(svd, &numConverged);

  // see how many we want to keep
  int totalKept = (numConverged < rank) ? numConverged : rank;

  // resize the results containers
  components.resizeAndWipe(_rows, totalKept);
  values.resizeAndWipe(totalKept);

  // copy the results into the results containers
  Vec u,v;
  VecCreateSeq(PETSC_COMM_SELF, _rows, &u);
  VecCreateSeq(PETSC_COMM_SELF, _cols, &v);
  PetscReal sigma;
  for (int x = 0; x < totalKept; x++)
  {
    PetscInt i = x;
    SVDGetSingularTriplet(svd, i, &sigma, u, v);

    // copy the V into the container
    PetscScalar* uArray;
    VecGetArray(u, &uArray);
    for (int y = 0; y < _rows; y++)
      components(y,x) = uArray[y];
    VecRestoreArray(u, &uArray);
    
    // principal components are the singular values squared
    values[x] = sigma * sigma;
  }
  // clean up
  SVDDestroy(svd);
  VecDestroy(u);
  VecDestroy(v);
}

/*
{
  // subtract out the means
  int rows = data.rows();
  int cols = data.cols();
  for (int x = 0; x < rows; x++)
  {
    Real mean = 0.0;
    for (int y = 0; y < cols; y++)
      mean += data(x,y);
    mean /= cols;
    cout << " mean: " << mean << endl;

    for (int y = 0; y < cols; y++)
      data(x,y) -= mean;
  }

  cout << "data after centering: " << data << endl;
 
  // divide through by sqrt(N - 1)
  data *= 1.0 / sqrt(cols - 1);

  cout << "data after sqrt: " << data << endl;
  data = data.transpose();

  // copy everything to a SPARSE_PETSC_MATRIX
  SPARSE_PETSC_MATRIX A(data);

  // finalize the current matrix
  A.finalize();

  // create the SVD solver
  SVD svd;
  SVDCreate(PETSC_COMM_WORLD,&svd);

  // set it to solve the current matrix
  SVDSetOperator(svd, A._APetsc);

  // set number of singular values to solve for
  SVDSetDimensions(svd, rank, rank * 2);
  SVDSetFromOptions(svd);
  SVDSolve(svd);

  // see how many triplets converged
  PetscInt numConverged;
  SVDGetConverged(svd, &numConverged);

  // resize the results containers
  components.resizeAndWipe(rows, rank);
  values.resizeAndWipe(rank);

  // copy the results into the results containers
  Vec u,v;
  PetscReal sigma;
  for (int x = 0; x < numConverged; x++)
  {
    SVDGetSingularTriplet(svd, x, &sigma, u, v);

    // copy the V into the container
    PetscScalar* vArray;
    VecGetArray(v, &vArray);
    for (int y = 0; y < rows; y++)
      components(y,x) = vArray[y];

    // principal components are the singular values squared
    values[x] = sigma * sigma;
  }
 
  // clean up
  SVDDestroy(svd);
}
*/
/*
{
  Mat         	 A, B;
  EPS         	 eps;
  PetscReal   	 error, tol;
  int         	 maxit, i, its;

  SlepcInitialize(NULL,NULL,NULL,NULL);

  // get preallocation hints for PetSc
  PetscInt* nonZeros = new PetscInt[10];
  int maxNonZeros = 1;
  for (int x = 0; x < 10; x++)
    nonZeros[x] = 1;

  // create A
  MatCreateSeqAIJ(PETSC_COMM_WORLD,10,10,maxNonZeros,nonZeros, &A);
  MatSetFromOptions(A);
 
  // create B
  MatCreateSeqAIJ(PETSC_COMM_WORLD,10,10,1, PETSC_NULL, &B);
  MatSetFromOptions(B);

  // delete preallocation hints for PetSc
  delete[] nonZeros;

  for (int x = 0; x < 10; x++)
  {
    PetscScalar 	 v = 1;
    MatSetValues(A,1,&x,1,&x,&v,INSERT_VALUES);
  }
  cout << " done. " << endl;
  flush(cout);

  // copy entries into mass matrix
  // Do it a generic but memory-hogging way
  cout << " Setting up mass eigensystem ... ";
  for (int x = 0; x < 10; x++)
  {
    PetscScalar 	 v = 10.0;
    MatSetValues(B,1,&x,1,&x,&v,INSERT_VALUES);
  }
  cout << " done. " << endl;

  // finalize matrices
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

  // Create eigensolver context
  EPSCreate(PETSC_COMM_WORLD,&eps);

  // Set operators. In this case, it is a generalized eigenvalue problem
  EPSSetOperators(eps,A,B);
  EPSSetProblemType(eps,EPS_GHEP);
  EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);

  int totalEigenvalues = 1;

  // Slepc recommends using at least twice as many eigenvectors
  // for the search space
  int totalEigenvectors = totalEigenvalues * 3 + 1;
  EPSSetDimensions(eps, totalEigenvalues, totalEigenvectors);

  PetscTruth hermitian = PETSC_TRUE;
  EPSIsHermitian(eps, &hermitian);

  // set the tolerances
  int maxIterations = 1000000;
  EPSSetTolerances(eps, 1e-15, maxIterations);

  // set the solver -- only if you have ARPACK installed!
  //EPSSetType(eps, EPSARPACK);

  // Send Slepc the solver parameters
  EPSSetFromOptions(eps);

  // Solve the eigensystem
  cout << " Solving eigensystem ... ";
  EPSSolve(eps);
  EPSGetIterationNumber(eps, &its);
  EPSGetTolerances(eps,&tol,&maxit);

  // Free work space
  EPSDestroy(eps);
  MatDestroy(A);
  SlepcFinalize();
}
*/

