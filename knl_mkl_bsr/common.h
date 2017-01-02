#ifndef COMMON_

#include "mat_coo_csr_converter.h"
#include "DataTypes.h"
#include <iostream>
#include <assert.h>

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifdef	__cplusplus
extern "C" {
#endif

	void delete_csrMat(csrMat * mat);
	void delete_diaMat(diaMat * mat);
	void delete_cooMat(cooMat * mat);
	void delete_ellMat(ellMat * mat);

	void draw_cooMat_pattern(cooMat mat);

	IDX_TYPE calculate_bsr_num_blocks (csrMat mat, IDX_TYPE block_size); // mat must be square

#ifdef	__cplusplus
}
#endif


#define COMMON_
#endif /*COMMON_*/