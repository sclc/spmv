#include <iostream>
#include <omp.h>
#include <string>
#include <assert.h>
#include <iomanip>

#include "DataTypes.h"
#include "readMTX.h"
#include "mat_coo_csr_converter.h"
#include "dense_mat_generator.h"
#include "omp_spmv_csr.h"
#include "error_checker.h"
#include "peformance_profiling.h"
#include "common.h"

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifndef EXP_NUM
#define EXP_NUM 100
#endif

#define CHECK_CSR_WITH_COO

//#define DEBUG_A
//#define DEBUG_B
//#define DEBUG_B_1
//#define DEBUG_B_2
//#define DEBUG_C
//#define DEBUG_D
//#define DEBUG_E

using namespace std;

int main (int argc, char* argv[])
{

	assert(argc == 4);
	string mat_path(argv[1]);
	string mat_filename(argv[2]);
	string size_rhs(argv[3]);

#ifdef DEBUG_A
	cout<< mat_path<<endl;
	cout<< mat_filename<<endl;
	cout<< size_rhs<<endl;
#endif

	matInfo mat_info;
	csrMat csrA;	
	cooMat cooA;

	double t1,t2;


	
	readMtx_info_and_ordered_coo(mat_path, mat_filename, &mat_info, &cooA);

#ifdef DEBUG_B

	cooMat cooA_tmp;
	matInfo mat_info_tmp;
	readMtx_info_and_coo(mat_path, mat_filename, &mat_info_tmp, &cooA_tmp);
	IDX_TYPE db_idx;
	cout<< "#row="<<cooA.num_rows<<", #col="<<cooA.num_cols<<", nnz="<<cooA.nnz<<endl;

	// for (db_idx=0;db_idx<cooA.nnz;db_idx++)
	// {
	// 	if (cooA_tmp.rowIdx[db_idx] - cooA.rowIdx[db_idx]  != 0 
	// 		|| cooA_tmp.colIdx[db_idx] -  cooA.colIdx[db_idx] != 0 
	// 		|| cooA_tmp.coodata[db_idx] - cooA.coodata[db_idx] != 0.0)
	// 	{
	// 		std::cout<< "value differs at:"<<db_idx<<std::endl;
	// 		std::cout<< cooA_tmp.rowIdx[db_idx] - cooA.rowIdx[db_idx] << ", "
	// 				<< cooA_tmp.colIdx[db_idx] -  cooA.colIdx[db_idx] << ", "
	// 				<< cooA_tmp.coodata[db_idx] - cooA.coodata[db_idx] << std::endl;

	// 		//std::cout<< cooA_tmp.coodata[db_idx]<< ", "<<cooA.coodata[db_idx] << std::endl;
	// 		//std::cout<< cooA_tmp.rowIdx[db_idx+5]<< ", "<<cooA.rowIdx[db_idx+5] << std::endl;
	// 		return 1;
	// 	}
	// }


	// std::cout<< "readMtx_info_and_coo == readMtx_info_and_ordered_coo"<<std::endl;

	for (db_idx=0;db_idx<cooA_tmp.nnz -1;db_idx++)
	{
		if (cooA_tmp.colIdx[db_idx] - cooA_tmp.colIdx[db_idx + 1] > 0)
		{
			std::cout<< "cooA_tmp wrong col order at "<<db_idx<<std::endl;
			return 1;
		}
	}

	delete_cooMat(&cooA_tmp);

#endif /*DEBUG_B*/

#ifdef	DEBUG_B_1
	IDX_TYPE db_idx;
	cout<< "#row="<<cooA.num_rows<<", #col="<<cooA.num_cols<<", nnz="<<cooA.nnz<<endl;

	for (db_idx=0;db_idx<cooA.nnz -1;db_idx++)
	{
		if (cooA.colIdx[db_idx] - cooA.colIdx[db_idx + 1] > 0)
		{
			std::cout<< "sorting did not work at "<<db_idx<<std::endl;
			return 1;
		}
	}
	std::cout<< "Qsorting works good "<<std::endl;


#endif /* DEBUG_B_1 */

	Converter_Coo2Csr (cooA, &csrA);
#ifdef DEBUG_C
	
	IDX_TYPE db_c_row_idx;
	cout<< "#row="<<csrA.num_rows<<", #col="<<csrA.num_cols<<", nnz="<<csrA.nnz<<endl;

	for (db_c_row_idx=0; db_c_row_idx<csrA.num_rows;db_c_row_idx++)
	{
		IDX_TYPE i_start = csrA.row_start[db_c_row_idx];
		IDX_TYPE i_end   = csrA.row_start[db_c_row_idx+1];
		IDX_TYPE db_c_col_idx;

		cout<<"row "<<db_c_row_idx<<":"<<endl;

		for(db_c_col_idx = i_start; db_c_col_idx < i_end; db_c_col_idx++)
		{
			cout<<"col: "<<csrA.col_idx[db_c_col_idx]<<", val:"<<csrA.csrdata[db_c_col_idx]<<" ;	";

		} 
		cout<<endl;
	}

#endif /* DEBUG_C */

	denseMat vec_result, vec, vec_checker;
	VAL_TYPE val_min, val_max;
	val_min = 10.;
	val_max = 100.;
	IDX_TYPE vec_ncol = 1;

	generate_dense_mat(&vec, csrA.num_rows, vec_ncol, val_min, val_max);
	generate_dense_mat_uniform_val(&vec_result, csrA.num_rows, vec_ncol, 66.66);

#ifdef CHECK_CSR_WITH_COO
	generate_dense_mat_uniform_val(&vec_checker, csrA.num_rows, vec_ncol, 66.66);
#endif /*CHECK_CSR_WITH_COO*/

#ifdef	DEBUG_B_2
	cooMat cooA_tmp;
	matInfo mat_info_tmp;
	denseMat vec_checker_2;
	readMtx_info_and_coo(mat_path, mat_filename, &mat_info_tmp, &cooA_tmp);
	generate_dense_mat_uniform_val(&vec_checker_2, csrA.num_rows, vec_ncol, 66.66);
	spmv_coo(vec_checker_2, cooA_tmp, vec);
#endif /* DEBUG_B_2 */

#ifdef DEBUG_D

	//cout<< "#vec_row:"<<vec.global_num_row<<" , #vec_col:"<<vec.global_num_col<<endl;
	cout<< "#vec_row:"<<vec_result.global_num_row<<" , #vec_col:"<<vec_result.global_num_col<<endl;

	IDX_TYPE db_d_idx;
	for (db_d_idx = 0; db_d_idx<vec.global_num_row*vec.global_num_col; db_d_idx++)
//		cout<<vec.data[db_d_idx]<<", ";
		cout<<vec_result.data[db_d_idx]<<", ";
	cout<<endl;

#endif/*DEBUG_D*/

#ifdef CHECK_CSR_WITH_COO
	spmv_coo(vec_checker, cooA, vec);
#endif /*CHECK_CSR_WITH_COO*/

#ifdef	DEBUG_B_2

	results_comparsion (vec_checker, vec_checker_2);

	std::cout<<" results_comparsion (vec_checker, vec_checker_2); pass"<<std::endl;

	delete_cooMat(&cooA_tmp);
#endif /* DEBUG_B_2 */
	//del coo
	delete_cooMat(&cooA);

	t1 = mysecond();
	
	for(int exp_idx=0; exp_idx<EXP_NUM; exp_idx++)
	{
		spmv_csr(vec_result, csrA, vec);
		//omp_spmv_csr_v1(vec_result, csrA, vec);
		//omp_spmv_csr_v2(vec_result, csrA, vec);
	}
	t2 = mysecond();

	VAL_TYPE per_loop_time_cost = (t2-t1)/(double)EXP_NUM;
	cout<< setprecision(12)<<per_loop_time_cost<<"sec pased"<<endl;
	cout<< "performance:"<<setprecision(12)<<( (double)(2*csrA.nnz)/(double)(1024*1024*1024) )/per_loop_time_cost <<"GFlops"<<endl;

#ifdef DEBUG_E

	cout<< "#vec_row:"<<vec_result.global_num_row<<" , #vec_col:"<<vec_result.global_num_col<<endl;

	IDX_TYPE db_e_idx;
	for (db_e_idx = 0; db_e_idx<vec.global_num_row*vec.global_num_col; db_e_idx++)
		cout<<vec_result.data[db_e_idx]<<", ";
	cout<<endl;
#endif


#ifdef CHECK_CSR_WITH_COO
	results_comparsion (vec_result, vec_checker);
#endif

	return 0;
}
