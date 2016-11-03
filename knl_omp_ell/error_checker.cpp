#include "error_checker.h"

void spmv_coo(denseMat vec_result, cooMat mat, denseMat vec)
{
	IDX_TYPE idx;

	assert(vec.global_num_row == mat.num_rows);
	assert(vec.global_num_col == 1);

	for (idx=0; idx< vec_result.global_num_row; idx++)
	{
		vec_result.data[idx] = 0.0;
	}

	for (idx=0; idx<mat.nnz; idx++)
	{
		vec_result.data[ mat.rowIdx[idx] ] += mat.coodata[ idx ] * vec.data[ mat.colIdx[idx] ];
	} 
}

void results_comparsion (denseMat vec_being_checked, denseMat vec_baseline)
{
	assert(vec_being_checked.global_num_row == vec_baseline.global_num_row);
	assert(vec_being_checked.global_num_col == vec_baseline.global_num_col);

	IDX_TYPE idx;

	for (idx =0; idx< vec_being_checked.global_num_row * vec_being_checked.global_num_col; idx++)
	{
		if ( abs(vec_being_checked.data[idx] - vec_baseline.data[idx] ) / vec_being_checked.data[idx] 
			 > CHECKER_BOUND)
		{
			std::cout<< "idx:"<<idx<<" error happening"<<std::endl;
			std::cout<< "value of result:"<<vec_being_checked.data[idx] 
			         <<", the real value:"<< vec_baseline.data[idx] <<std::endl;

			return;
		}
	}

	std::cout<< "checker_on_coo_spmv passed!!!"<<std::endl;

}
