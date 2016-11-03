//
#ifndef MAT_CSR_ELL_CONVERTER_

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

void Converter_Csr2Ell (csrMat src, ellMat * target);

#ifdef	__cplusplus
}
#endif

#define MAT_CSR_ELL_CONVERTER_

#endif /*MAT_CSR_ELL_CONVERTER_*/