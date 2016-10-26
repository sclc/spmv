#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifndef DATATYPES_H
#include "DataTypes.h"
#endif

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifdef	__cplusplus
extern "C" {
#endif

void generate_dense_mat(denseMat * mat, IDX_TYPE num_rows, IDX_TYPE num_cols,\
                     VAL_TYPE ranMin, VAL_TYPE ranMax);

#ifdef	__cplusplus
}
#endif