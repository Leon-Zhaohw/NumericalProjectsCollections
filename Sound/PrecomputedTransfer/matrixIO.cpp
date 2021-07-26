#include "matrixIO.h"

// writes out the m x n matrix onto the stream, in binary format
// LAPACK column-major order
template <class real>
int WriteMatrixToStream(FILE * file, int m, int n, real * matrix)
{
  if ((int)(fwrite(matrix,sizeof(real),m*n,file)) < m*n)
    return 1;

  return 0;
}

// writes out the m x n matrix header onto the stream
int WriteMatrixHeaderToStream(FILE * file, int m, int n)
{
  if (fwrite(&m,sizeof(int),1,file) < 1)
    return 1;

  if (fwrite(&n,sizeof(int),1,file) < 1)
    return 1;

  return 0;
}


template <class real>
int WriteMatrixToDisk(char* filename, int m, int n, real * matrix)
{
  FILE * file;
  file = fopen(filename,"wb");
  if (!file)
  {
    printf ("Can't open output file: %s.\n",filename);
    return 1;
  }

  if (WriteMatrixHeaderToStream(file,m,n) != 0)
  {
    printf ("Error writing the matrix header to disk file: %s.\n",filename);
    return 1;
  }

  if (WriteMatrixToStream(file,m,n,matrix) != 0)
  {
    printf ("Error writing the matrix to disk file: %s.\n",filename);
    return 1;
  }

  fclose(file);

  return 0;
}


// read the m x n matrix from the stream, in binary format
template <class real>
int ReadMatrixFromStream(FILE * file, int M, int N, real * matrix)
{

  if (fread(matrix,sizeof(real),M*N,file) < (unsigned int)M*N)
    return 1;

  return 0;
}

template <class real>
int ReadMatrixFromDisk(char * filename, int * m, int * n, real ** matrix)
{
  FILE * file;
  file = fopen(filename,"rb");
  if (!file)
  {
    printf ("Can't open input matrix file: %s.\n",filename);
    return 1;
  }

  if (ReadMatrixSizeFromStream(file,m,n) != 0)
  {
    printf ("Error reading the matrix from disk file: %s.\n",filename);
    return 1;
  }

  *matrix = (real *) malloc (sizeof(real)*(*m)*(*n));

  if (ReadMatrixFromStream(file,*m,*n,*matrix) != 0)
  {
    printf ("Error reading the matrix from disk file: %s.\n",filename);
    return 1;
  }

  fclose(file);

  return 0;
}

template <class real>
int WriteMatrixToDisk_(char* filename, int m, int n, real * matrix)
{
  if (WriteMatrixToDisk(filename, m, n, matrix) != 0)
    exit(1);

  return 0;
}

template <class real>
int ReadMatrixFromDisk_(char* filename, int * m, int * n, real ** matrix)
{
  if (ReadMatrixFromDisk(filename, m, n, matrix) != 0)
    exit(1);

  return 0;
}


int ReadMatrixSizeFromDisk(char * filename, int * m, int * n)
{
  FILE * file;
  file = fopen(filename,"rb");
  if (!file)
  {
    printf ("Can't open input matrix file: %s.\n",filename);
    return 1;
  }

  if (fread(m,sizeof(int),1,file) < 1)
    return 1;

  if (fread(n,sizeof(int),1,file) < 1)
    return 1;

  fclose(file);

  return 0;
}

void ReadMatrixSizeFromDisk_(char * filename, int * m, int * n)
{
  if (ReadMatrixSizeFromDisk(filename,m,n) != 0)
    exit(1);
}

int ReadMatrixSizeFromStream(FILE * file, int * m, int * n)
{
  if (fread(m,sizeof(int),1,file) < 1)
    return 1;

  if (fread(n,sizeof(int),1,file) < 1)
    return 1;

  return 0;
}

int OpenFile_(char * filename, FILE ** fin, char * mode)
{
  *fin = fopen(filename,mode);
  if (!(*fin))
  {
    printf("Error: could not open file %s.\n",filename);
    exit(1);
  }

  return 0;
}

int Assert_(int a, int b, int positionIdentifier)
{
  if (a != b)
  {
    printf("Error: assertion %d == %d failed (%d).\n",a,b,positionIdentifier);
    exit(1);
  }
  return 0;
}

int Assert_(bool cond, int positionIdentifier)
{
  if (!cond)
  {
    printf("Error: assertion failed (%d).\n",positionIdentifier);
    exit(1);
  }
  return 0;
}

template int WriteMatrixToStream<double>(FILE * file, int m, int n, double * matrix);
template int WriteMatrixToStream<float>(FILE * file, int m, int n, float * matrix);

template int WriteMatrixToDisk<double>(char* filename, int m, int n, double * matrix);
template int WriteMatrixToDisk<float>(char* filename, int m, int n, float * matrix);

template int WriteMatrixToDisk_<double>(char* filename, int m, int n, double * matrix);
template int WriteMatrixToDisk_<float>(char* filename, int m, int n, float * matrix);

template int ReadMatrixFromStream<double>(FILE * file, int m, int n, double * matrix);
template int ReadMatrixFromStream<float>(FILE * file, int m, int n, float * matrix);

template int ReadMatrixFromDisk<double>(char* filename, int * m, int * n, double ** matrix);
template int ReadMatrixFromDisk<float>(char* filename, int * m, int * n, float ** matrix);

template int ReadMatrixFromDisk_<double>(char* filename, int * m, int * n, double ** matrix);
template int ReadMatrixFromDisk_<float>(char* filename, int * m, int * n, float ** matrix);
