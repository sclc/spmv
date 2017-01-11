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
#include "tbb_spmv_csr.h"
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
#define EXP_NUM 1
#endif

#define CHECK_CSR_WITH_COO

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


	
//	readMtx_info_and_ordered_coo(mat_path, mat_filename, &mat_info, &cooA);
	readMtx_info_and_coo(mat_path, mat_filename, &mat_info, &cooA);


	Converter_Coo2Csr (cooA, &csrA);

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


#ifdef CHECK_CSR_WITH_COO
	spmv_coo(vec_checker, cooA, vec);
#endif /*CHECK_CSR_WITH_COO*/

	//del coo
	delete_cooMat(&cooA);

	t1 = mysecond();
	
	for(int exp_idx=0; exp_idx<EXP_NUM; exp_idx++)
	{
		//spmv_csr(vec_result, csrA, vec);
		//omp_spmv_csr_v1(vec_result, csrA, vec);
		//omp_spmv_csr_v2(vec_result, csrA, vec);
		tbb_spmv_csr_v1(vec_result, csrA, vec);
	}
	t2 = mysecond();

	VAL_TYPE per_loop_time_cost = (t2-t1)/(double)EXP_NUM;
	cout<< setprecision(12)<<per_loop_time_cost<<"sec pased"<<endl;
	cout<< "performance:"<<setprecision(12)<<( (double)(2*csrA.nnz)/(double)(1024*1024*1024) )/per_loop_time_cost <<"GFlops"<<endl;


#ifdef CHECK_CSR_WITH_COO
	results_comparsion (vec_result, vec_checker);
#endif

	return 0;
}
