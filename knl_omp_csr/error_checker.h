#ifndef DATATYPES_H
#include "DataTypes.h"
#endif
#include <iostream>
#include <assert.h>
#include <stdlib.h>

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifndef CHECKER_BOUND
#define CHECKER_BOUND 1e-12
#endif

#ifdef	__cplusplus
extern "C" {
#endif

void spmv_coo(denseMat vec_result, cooMat mat, denseMat vec);
void results_comparsion (denseMat vec_being_checked, denseMat vec_baseline);

#ifdef	__cplusplus
}
#endif
