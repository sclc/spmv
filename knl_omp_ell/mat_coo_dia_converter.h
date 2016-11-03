//
#include <stdio.h>
#include <stdlib.h>
#include "DataTypes.h"
#include <assert.h>
#include <string.h>

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifndef DIA_MAGIC_NUM
#define DIA_MAGIC_NUM 8888
#endif

#ifndef COO2DIA_DEBUG

#define COO2DIA_DEBUG_A
#define COO2DIA_DEBUG_B

#define COO2DIA_DEBUG
#endif

#ifdef	__cplusplus
extern "C" {
#endif

void Converter_Coo2Dia (cooMat src, diaMat * target);

#ifdef	__cplusplus
}
#endif