#include "omp_spmv_csr.h"

void spmv_csr(denseMat vec_result, csrMat mat, denseMat vec)
{

#ifdef CSR_SPMV_DEBUG_A
	assert (vec_result.global_num_col == 1);
	assert (vec.global_num_col == 1);
#endif

	IDX_TYPE mat_row_idx, mat_col_idx;

	for (mat_row_idx = 0; mat_row_idx<mat.num_rows; mat_row_idx++)
	{
		IDX_TYPE i_start = mat.row_start[mat_row_idx];
		IDX_TYPE i_end   = mat.row_start[mat_row_idx+1];
		vec_result.data[mat_row_idx]=0.0;

		for (mat_col_idx=i_start; mat_col_idx<i_end; mat_col_idx++)
		{
			vec_result.data[mat_row_idx] += mat.csrdata[mat_col_idx] 
									 * vec.data[ mat.col_idx[mat_col_idx] ];
		}
	}

}

void omp_spmv_csr_v1(denseMat vec_result, csrMat mat, denseMat vec)
{

#ifdef CSR_SPMV_DEBUG_A
	assert (vec_result.global_num_col == 1);
	assert (vec.global_num_col == 1);
#endif

	IDX_TYPE mat_row_idx, mat_col_idx;

//#pragma omp parallel num_threads(8)
#pragma omp parallel 
	{
		// int tid, nthreads;
		// tid = omp_get_thread_num();
		// nthreads = omp_get_num_threads();
		// if (tid==0)
		// 	printf("total #threads:%d, myid:%d\n",nthreads,tid);

		#pragma omp for
		for (mat_row_idx = 0; mat_row_idx<mat.num_rows; mat_row_idx++)
		{
			IDX_TYPE i_start = mat.row_start[mat_row_idx];
			IDX_TYPE i_end   = mat.row_start[mat_row_idx+1];
			vec_result.data[mat_row_idx]=0.0;

			for (mat_col_idx=i_start; mat_col_idx<i_end; mat_col_idx++)
			{
				vec_result.data[mat_row_idx] += mat.csrdata[mat_col_idx] 
										 * vec.data[ mat.col_idx[mat_col_idx] ];
			}
		}

	}


}

void omp_spmv_csr_v2(denseMat vec_result, csrMat mat, denseMat vec)
{

#ifdef CSR_SPMV_DEBUG_A
	assert (vec_result.global_num_col == 1);
	assert (vec.global_num_col == 1);
#endif

	IDX_TYPE mat_row_idx, mat_col_idx;

	for (mat_row_idx = 0; mat_row_idx<mat.num_rows; mat_row_idx++)
	{
		IDX_TYPE i_start = mat.row_start[mat_row_idx];
		IDX_TYPE i_end   = mat.row_start[mat_row_idx+1];
		VAL_TYPE temp_val= 0.0;

		#pragma omp parallel for reduction (+:temp_val)
		for (mat_col_idx=i_start; mat_col_idx<i_end; mat_col_idx++)
		{
			temp_val += mat.csrdata[mat_col_idx] 
						* vec.data[ mat.col_idx[mat_col_idx] ];
		}
		vec_result.data[mat_row_idx] = temp_val;
	}


}