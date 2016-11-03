//
#ifndef MAT_COO_CSR_CONVERTER

#include <stdio.h>
#include <stdlib.h>
#include "DataTypes.h"

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifdef	__cplusplus
extern "C" {
#endif

void Converter_Coo2Csr (cooMat src, csrMat * target);

#ifdef	__cplusplus
}
#endif

#define MAT_COO_CSR_CONVERTER

#endif /*MAT_COO_CSR_CONVERTER*/