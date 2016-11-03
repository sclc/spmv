#include "omp_spmv_dia.h"

void spmv_dia(denseMat * vec_result, diaMat mat, denseMat vec)
{
#ifdef DIA_SPMV_DEBUG_A
	assert (vec_result->global_num_row == mat.num_rows);
	assert (vec_result->global_num_col == 1);
	assert (vec.global_num_col == 1);
#endif

	// for profiling the kernel, new re-allocate vec_result for every call
	if(vec_result->data != NULL) free(vec_result->data);
	vec_result->data = (VAL_TYPE *) calloc ( vec_result->global_num_row * vec_result->global_num_col,
											sizeof(VAL_TYPE) );

	// assume mat is square now
	IDX_TYPE diag_length = mat.num_rows;

	IDX_TYPE idx, offset, on_diag_idx;
	IDX_TYPE rowid, colid, diag_offset_start;

	for (idx=0; idx < mat.num_diags; idx++)
	{
		diag_offset_start = idx * diag_length;
		offset = mat.diag_offset_id[idx];
		if (offset < 0)
		{
			for (on_diag_idx = -offset; on_diag_idx<diag_length; on_diag_idx++)
			{
				rowid = on_diag_idx;
				colid = offset + on_diag_idx;
				vec_result->data[rowid] += mat.diag_zone[diag_offset_start + on_diag_idx] 
									* vec.data[colid];
			}
		}
		else
		{
			for (on_diag_idx = 0; on_diag_idx<diag_length -offset ; on_diag_idx++)
			{
				rowid = on_diag_idx;
				colid = offset + on_diag_idx;
				vec_result->data[rowid] += mat.diag_zone[diag_offset_start + on_diag_idx] 
									* vec.data[colid];
			}
		}
	}
}

void omp_spmv_dia_v1(denseMat * vec_result, diaMat mat, denseMat vec)
{
#ifdef DIA_SPMV_DEBUG_A
	assert (vec_result->global_num_row == mat.num_rows);
	assert (vec_result->global_num_col == 1);
	assert (vec.global_num_col == 1);
#endif

	// for profiling the kernel, new re-allocate vec_result for every call
	if(vec_result->data != NULL) free(vec_result->data);
	vec_result->data = (VAL_TYPE *) calloc ( vec_result->global_num_row * vec_result->global_num_col,
											sizeof(VAL_TYPE) );

	// assume mat is square now
	IDX_TYPE diag_length = mat.num_rows;

	IDX_TYPE idx;

#pragma omp parallel for
	for (idx=0; idx < mat.num_diags; idx++)
	{
		IDX_TYPE offset, on_diag_idx;
		IDX_TYPE rowid, colid, diag_offset_start;

		diag_offset_start = idx * diag_length;
		offset = mat.diag_offset_id[idx];
		if (offset < 0)
		{
			for (on_diag_idx = -offset; on_diag_idx<diag_length; on_diag_idx++)
			{
				rowid = on_diag_idx;
				colid = offset + on_diag_idx;
				#pragma omp atomic
				vec_result->data[rowid] += mat.diag_zone[diag_offset_start + on_diag_idx] 
									* vec.data[colid];
			}
		}
		else
		{
			for (on_diag_idx = 0; on_diag_idx<diag_length -offset ; on_diag_idx++)
			{
				rowid = on_diag_idx;
				colid = offset + on_diag_idx;
				#pragma omp atomic
				vec_result->data[rowid] += mat.diag_zone[diag_offset_start + on_diag_idx] 
									* vec.data[colid];
			}
		}
	}
}
