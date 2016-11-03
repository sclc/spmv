#include "omp_spmv_ell.h"

void spmv_ell(denseMat * vec_result, ellMat mat, denseMat vec)
{
#ifdef DIA_SPMV_DEBUG_A
	assert (vec_result->global_num_row == mat.ell_num_rows);
	assert (vec_result->global_num_col == 1);
	assert (vec.global_num_col == 1);
#endif

	if(vec_result->data != NULL) free(vec_result->data);
	vec_result->data = (VAL_TYPE*)calloc( vec_result->global_num_row * vec_result->global_num_col, 
											sizeof(VAL_TYPE));

	IDX_TYPE ell_lda = mat.ell_num_rows;
	IDX_TYPE ell_col_idx;

	for (ell_col_idx = 0; ell_col_idx < mat.ell_row_length; ell_col_idx++)
	{
		IDX_TYPE ell_row_idx;
		for(ell_row_idx=0; ell_row_idx<mat.ell_num_rows; ell_row_idx++)
		{
			// position of the element to touch
			IDX_TYPE pos = ell_row_idx + ell_lda * ell_col_idx;
			IDX_TYPE rowid = ell_row_idx;
			IDX_TYPE colid = mat.ell_col_idx[pos];
			VAL_TYPE val   = mat.ell_data[pos];

			vec_result->data[rowid] += val * vec.data[colid];
		}
	}


}

void omp_spmv_ell_v1(denseMat * vec_result, ellMat mat, denseMat vec)
{
#ifdef DIA_SPMV_DEBUG_A
	assert (vec_result->global_num_row == mat.ell_num_rows);
	assert (vec_result->global_num_col == 1);
	assert (vec.global_num_col == 1);
#endif

	if(vec_result->data != NULL) free(vec_result->data);
	vec_result->data = (VAL_TYPE*)calloc( vec_result->global_num_row * vec_result->global_num_col, 
											sizeof(VAL_TYPE));

	IDX_TYPE ell_lda = mat.ell_num_rows;
	IDX_TYPE ell_col_idx;

	for (ell_col_idx = 0; ell_col_idx < mat.ell_row_length; ell_col_idx++)
	{
		IDX_TYPE ell_row_idx;
		#pragma omp parallel for
		for(ell_row_idx=0; ell_row_idx<mat.ell_num_rows; ell_row_idx++)
		{
			// position of the element to touch
			IDX_TYPE pos = ell_row_idx + ell_lda * ell_col_idx;
			IDX_TYPE rowid = ell_row_idx;
			IDX_TYPE colid = mat.ell_col_idx[pos];
			VAL_TYPE val   = mat.ell_data[pos];

			//#pragma omp atomic
			vec_result->data[rowid] += val * vec.data[colid];
		}
	}


}

