#include "MATRIX_COMPRESSION_DATA.h"

MATRIX_COMPRESSION_DATA::MATRIX_COMPRESSION_DATA()
{
}

MATRIX_COMPRESSION_DATA::MATRIX_COMPRESSION_DATA(int* dataX, int* dataY, int* dataZ,
        COMPRESSION_DATA* compression_dataX, COMPRESSION_DATA* compression_dataY, COMPRESSION_DATA* compression_dataZ) :

  _dataX(dataX), _dataY(dataY), _dataZ(dataZ),
  _compression_dataX(*compression_dataX), _compression_dataY(*compression_dataY), _compression_dataZ(*compression_dataZ)
{
}

/*
MATRIX_COMPRESSION_DATA::MATRIX_COMPRESSION_DATA(MATRIX_COMPRESSION_DATA& m) :
  _dataX(m.get_dataX()), _dataY(m.get_dataY()), _dataZ(m.get_dataZ()),
  _compression_dataX(*(m.get_compression_dataX())), _compression_dataY(*(m.get_compression_dataY())), _compression_dataZ(*m.get_compression_dataZ())

{
  _dct_plan = m.get_plan();
  _cachedBlocksX = *(m.get_cachedBlocksX());
  _cachedBlocksY = *(m.get_cachedBlocksY());
  _cachedBlocksZ = *(m.get_cachedBlocksZ());

  if (m.get_dct_in()) {
    _dct_in = (double*)fftw_malloc(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
  }
  else {
    _dct_in = NULL;
  }

}
*/

/*
MATRIX_COMPRESSION_DATA& MATRIX_COMPRESSION_DATA::operator=(MATRIX_COMPRESSION_DATA& m) 
{
  _dataX = m.get_dataX();
  _dataY = m.get_dataY();
  _dataZ = m.get_dataZ();
  _compression_dataX = *(m.get_compression_dataX());
  _compression_dataY = *(m.get_compression_dataY());
  _compression_dataZ = *(m.get_compression_dataZ());

  _dct_plan = m.get_plan();
  _cachedBlocksX = *(m.get_cachedBlocksX());
  _cachedBlocksY = *(m.get_cachedBlocksY());
  _cachedBlocksZ = *(m.get_cachedBlocksZ());

  if (m.get_dct_in()) {
    _dct_in = (double*)fftw_malloc(BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
  }
  else {
    _dct_in = NULL;
  }
  return *this;
}
*/

/*
MATRIX_COMPRESSION_DATA::~MATRIX_COMPRESSION_DATA()
{
  if (_dct_in) {
    fftw_free(_dct_in);
    _dct_in = NULL;
  }
}
*/
