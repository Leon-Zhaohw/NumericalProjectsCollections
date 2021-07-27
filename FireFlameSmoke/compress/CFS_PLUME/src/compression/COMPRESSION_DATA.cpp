#include "COMPRESSION_DATA.h"

COMPRESSION_DATA::COMPRESSION_DATA() {}


COMPRESSION_DATA::COMPRESSION_DATA(const VEC3I& dims, int numCols, int nBits, double percent) :
  _dims(dims), _numCols(numCols), _nBits(nBits), _percent(percent) {}

/*
COMPRESSION_DATA::COMPRESSION_DATA(COMPRESSION_DATA& data) :
  _dims(data.get_dims()), _paddedDims(data.get_paddedDims()), _numCols(data.get_numCols()), _nBits(data.get_nBits()), _percent(data.get_percent()),
  _maxIterations(data.get_maxIterations()), _numBlocks(data.get_numBlocks()) 
{
   _blockLengths = data.get_blockLengths();
   _blockIndices = data.get_blockIndices();
   _blockLengthsMatrix = *(data.get_blockLengthsMatrix());
   _blockIndicesMatrix = *(data.get_blockIndicesMatrix());
   _sList = *(data.get_sList());
   _gammaList = *(data.get_gammaList());
   _singularList = *(data.get_singularList());
   _vList = *(data.get_vList());
   _dampingArray = data.get_dampingArray();
   _zigzagArray = data.get_zigzagArray();
   _reverseZigzag = data.get_reverseZigzag();
   _dct_plan = data.get_dct_plan();

   // memory allocations
   if (data.get_dct_in()) {
     _dct_in = (double*)fftw_malloc(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
     }
   else { 
     _dct_in = NULL;
     }
   if (data.get_dct_out()) {
     _dct_out = (double*)fftw_malloc(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
   }
   else {
     _dct_out = NULL;
   }
}
*/

/*  
COMPRESSION_DATA& COMPRESSION_DATA::operator=(COMPRESSION_DATA& data)
{
   _dims = data.get_dims();
   _paddedDims = data.get_paddedDims();
   _numCols = data.get_numCols();
   _nBits = data.get_nBits();
   _percent = data.get_percent();
   _maxIterations = data.get_maxIterations();
   _numBlocks = data.get_numBlocks();
   _blockLengths = data.get_blockLengths();
   _blockIndices = data.get_blockIndices();
   _blockLengthsMatrix = *(data.get_blockLengthsMatrix());
   _blockIndicesMatrix = *(data.get_blockIndicesMatrix());
   _sList = *(data.get_sList());
   _gammaList = *(data.get_gammaList());
   _singularList = *(data.get_singularList());
   _vList = *(data.get_vList());
   _dampingArray = data.get_dampingArray();
   _zigzagArray = data.get_zigzagArray();
   _reverseZigzag = data.get_reverseZigzag();
   _dct_plan = data.get_dct_plan();

   // memory allocations
   if (data.get_dct_in()) {
     _dct_in = (double*)fftw_malloc(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
     }
   else { 
     _dct_in = NULL;
     }
   if (data.get_dct_out()) {
     _dct_out = (double*)fftw_malloc(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
   }
   else {
     _dct_out = NULL;
   }

   return *this;
}
*/

/*
COMPRESSION_DATA::~COMPRESSION_DATA()
{
  if (_dct_in) {
    fftw_free(_dct_in);
    _dct_in = NULL;
  }
  if (_dct_out) {
    fftw_free(_dct_out);
    _dct_out = NULL;
  }
}
*/
