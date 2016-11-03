#include <stdio.h>
#include <stdlib.h>
#include "DataTypes.h"
#include <assert.h>
#include <omp.h>

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifndef DIA_SPMV_DEBUG
#define DIA_SPMV_DEBUG_A

#define DIA_SPMV_DEBUG
#endif

#ifdef	__cplusplus
extern "C" {
#endif
	void spmv_ell(denseMat *vec_result, ellMat mat, denseMat vec);
	void omp_spmv_ell_v1(denseMat *vec_result, ellMat mat, denseMat vec);


#ifdef	__cplusplus
}
#endif