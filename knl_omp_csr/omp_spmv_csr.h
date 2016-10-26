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

#ifndef SPMV_DEBUG
#define SPMV_DEBUG_A

#endif

#ifdef	__cplusplus
extern "C" {
#endif
	void spmv_csr(denseMat vec_result, csrMat mat, denseMat vec);

#ifdef	__cplusplus
}
#endif