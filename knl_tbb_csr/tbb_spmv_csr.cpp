#include "tbb_spmv_csr.h"

class tbb_csr_calculator_1 {

public:
	double *result;
	double *vec;
	csrMat *csrmat;

	void operator() (const blocked_range<IDX_TYPE>& range) const {
		for (IDX_TYPE rowidx=range.begin(); rowidx!=range.end(); rowidx++){

			IDX_TYPE i_start = csrmat->row_start[rowidx];
			IDX_TYPE i_end   = csrmat->row_start[rowidx+1];
			//(*result)[rowidx]=0.0;
			* (result+rowidx) = 0.0;
			for (IDX_TYPE mat_col_idx=i_start; mat_col_idx<i_end; mat_col_idx++)
			{
				// (*result)[rowidx] += (*csrmat).csrdata[rowidx] 
				// 					 * (*vec)[ (*csrmat).col_idx[rowidx] ];
				*(result+rowidx)  += csrmat->csrdata[mat_col_idx]
									 * *(vec + csrmat->col_idx[mat_col_idx]);
			}
		}
	}

};

void tbb_spmv_csr_v1(denseMat vec_result, csrMat mat, denseMat vec)
{

#ifdef CSR_SPMV_DEBUG_A
	assert (vec_result.global_num_col == 1);
	assert (vec.global_num_col == 1);
#endif

	tbb_csr_calculator_1 cal_csr;

	cal_csr.result = vec_result.data;
	cal_csr.vec    = vec.data;
	cal_csr.csrmat = &mat;
	IDX_TYPE grainsize = 256;

	//parallel_for(blocked_range<IDX_TYPE>(0,mat.num_rows, grainsize) ,cal_csr);
	parallel_for(blocked_range<IDX_TYPE>(0,mat.num_rows) ,cal_csr);
}
