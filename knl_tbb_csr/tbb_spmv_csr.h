#include <stdio.h>
#include <stdlib.h>
#include "DataTypes.h"
#include <assert.h>

#include "tbb/tbb.h"
// #include "tbb/parallel_for.h"
// #include "tbb/blocked_range.h"
using namespace tbb;

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifndef CSR_SPMV_DEBUG
//#define CSR_SPMV_DEBUG_A

#define CSR_SPMV_DEBUG
#endif

void tbb_spmv_csr_v1(denseMat vec_result, csrMat mat, denseMat vec);

#ifdef	__cplusplus
extern "C" {
#endif
	// to add functions

#ifdef	__cplusplus
}
#endif