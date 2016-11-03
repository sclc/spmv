#include <iostream>
#include <omp.h>
#include <string>
#include <assert.h>
#include <iomanip>

#include "DataTypes.h"
#include "readMTX.h"
#include "mat_coo_csr_converter.h"
#include "mat_csr_ell_converter.h"
#include "dense_mat_generator.h"
#include "error_checker.h"
#include "peformance_profiling.h"
#include "common.h"
#include "omp_spmv_ell.h"

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifndef EXP_NUM
#define EXP_NUM 50
#endif

#define CHECK_RES_WITH_COO

//#define DEBUG_A
//#define DEBUG_B
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
	csrMat csrA_temp;	
	cooMat cooA;
	ellMat ellA;

	double t1,t2;


	readMtx_info_and_coo(mat_path, mat_filename, &mat_info, &cooA);
#ifdef DEBUG_B
	
	IDX_TYPE db_idx;
	cout<< "#row="<<cooA.num_rows<<", #col="<<cooA.num_cols<<", nnz="<<cooA.nnz<<endl;

	for (db_idx=0;db_idx<cooA.nnz;db_idx++)
	{
		cout<<cooA.rowIdx[db_idx]<<", "<<cooA.colIdx[db_idx]<<", "<<setprecision(10)<<cooA.coodata[db_idx]<<endl;
	}
#endif

	Converter_Coo2Csr (cooA, &csrA_temp);
	Converter_Csr2Ell (csrA_temp, &ellA);
	delete_csrMat (&csrA_temp);

#ifdef DEBUG_C
	
	//IDX_TYPE db_c_idx, db_c_ondiag_idx;
	cout<< "#row="<<ellA.ell_num_rows<<", #col="<<ellA.ell_num_cols<<", nnz="<<ellA.ell_nnz
	<< ", ell_row_length="<<ellA.ell_row_length<<endl;

	//draw_cooMat_pattern(cooA);

	return 0;

#endif



	denseMat vec_result, vec, vec_checker;
	VAL_TYPE val_min, val_max;
	val_min = 10.;
	val_max = 100.;
	IDX_TYPE vec_ncol = 1;

	generate_dense_mat(&vec, ellA.ell_num_rows, vec_ncol, val_min, val_max);
	generate_dense_mat_uniform_val(&vec_result, ellA.ell_num_rows, vec_ncol, 66.66);
	// generate_dense_mat_uniform_val(&vec_checker, diaA.num_rows, vec_ncol, 66.66);

#ifdef CHECK_RES_WITH_COO
	generate_dense_mat_uniform_val(&vec_checker, cooA.num_rows, vec_ncol, 66.66);
#endif /*CHECK_RES_WITH_COO*/

#ifdef DEBUG_D

	//cout<< "#vec_row:"<<vec.global_num_row<<" , #vec_col:"<<vec.global_num_col<<endl;
	cout<< "#vec_row:"<<vec_result.global_num_row<<" , #vec_col:"<<vec_result.global_num_col<<endl;

	IDX_TYPE db_d_idx;
	for (db_d_idx = 0; db_d_idx<vec.global_num_row*vec.global_num_col; db_d_idx++)
//		cout<<vec.data[db_d_idx]<<", ";
		cout<<vec_result.data[db_d_idx]<<", ";
	cout<<endl;

#endif /*DEBUG_D*/

#ifdef CHECK_RES_WITH_COO
	spmv_coo(vec_checker, cooA, vec);
#endif/*CHECK_RES_WITH_COO*/

	// delete COO
	delete_cooMat(&cooA);

	t1 = mysecond();
	
	for(int exp_idx=0; exp_idx<EXP_NUM; exp_idx++)
	{
		//spmv_ell(&vec_result, ellA, vec);
		omp_spmv_ell_v1(&vec_result, ellA, vec);
	}
	t2 = mysecond();

	VAL_TYPE per_loop_time_cost = (t2-t1)/(double)EXP_NUM;
	cout<< setprecision(12)<<per_loop_time_cost<<"sec pased"<<endl;
	cout<< "performance:"<<setprecision(12)<<( (double)(2*ellA.ell_nnz)/(double)(1024*1024*1024) )/per_loop_time_cost
	    <<"GFlops"<<endl;

#ifdef DEBUG_E

	cout<< "#vec_row:"<<vec_result.global_num_row<<" , #vec_col:"<<vec_result.global_num_col<<endl;

	IDX_TYPE db_e_idx;
	for (db_e_idx = 0; db_e_idx<vec.global_num_row*vec.global_num_col; db_e_idx++)
		cout<<vec_result.data[db_e_idx]<<", ";
	cout<<endl;
#endif /*DEBUG_E*/

	// spmv_coo(vec_checker, cooA, vec);
	// results_comparsion (vec_result, vec_checker);
#ifdef CHECK_RES_WITH_COO
	results_comparsion (vec_result, vec_checker);
#endif /*CHECK_RES_WITH_COO*/

	return 0;
}
