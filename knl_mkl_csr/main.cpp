#include <iostream>
#include <omp.h>
#include <string>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include "mkl_spblas.h"

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

#ifndef MKL_EXP_NUM
#define MKL_EXP_NUM 100
#endif

#define CHECK_CSR_WITH_COO

// #define DEBUG_MKL_SPMV

//#define DEBUG_A

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

/* MKL SpMV start */

	// declare space for mkl_sparse_d_create_csr
	double * mkl_csrdata = (double*)calloc( csrA.nnz,sizeof(double));
	MKL_INT * mkl_csr_colInd = (MKL_INT *) calloc ( csrA.nnz,sizeof(MKL_INT));
	MKL_INT * mkl_csr_rowInd = (MKL_INT *) calloc ( csrA.num_rows+1, sizeof(MKL_INT) );

	struct matrix_descr mkl_descrA;
	sparse_matrix_t mkl_csrA;

	assert( sizeof(double) == sizeof(VAL_TYPE) );
	assert( sizeof(IDX_TYPE) == sizeof(MKL_INT) );

	memcpy( mkl_csrdata, csrA.csrdata, csrA.nnz * sizeof(VAL_TYPE) );
	memcpy( mkl_csr_colInd, csrA.col_idx, csrA.nnz * sizeof(IDX_TYPE) );
	memcpy( mkl_csr_rowInd, csrA.row_start, (csrA.num_rows+1) * sizeof(IDX_TYPE) );

	// Create handle with matrix stored in CSR format
    mkl_sparse_d_create_csr ( &mkl_csrA, SPARSE_INDEX_BASE_ZERO,
                                    csrA.num_rows,  // number of rows
                                    csrA.num_cols,  // number of cols
                                    mkl_csr_rowInd,
                                    mkl_csr_rowInd+1,
                                    mkl_csr_colInd,
                                    mkl_csrdata );
    // Create matrix descriptor
    mkl_descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_optimize ( mkl_csrA );

    double mkl_t1, mkl_t2;
    double mkl_spmv_overhead_accumulator = 0.0;

    double mkl_alpha= 1.0, mkl_beta = 0.0;
    int mkl_exp_counter;
    double * mkl_spmv_result = NULL;

    for ( mkl_exp_counter = 0; mkl_exp_counter < MKL_EXP_NUM; mkl_exp_counter++)
    {

    	// free and reallocate
    	if (mkl_spmv_result != NULL) free(mkl_spmv_result);
    	mkl_spmv_result = (double *) calloc ( csrA.num_rows, sizeof(double));
    	    // y := alpha*op(A)*x + beta*y
    	// Compute y = alpha * A * x + beta * y

    	mkl_t1 = mysecond();
    	mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE,
          	            mkl_alpha,
          	            mkl_csrA,
          	            mkl_descrA,
           	            vec.data, // x
           	            mkl_beta,
           	            mkl_spmv_result ); // y

    	mkl_t2 = mysecond();
    	mkl_spmv_overhead_accumulator += mkl_t2-mkl_t1;
    }
    mkl_spmv_overhead_accumulator /= (double)MKL_EXP_NUM;
    cout<< setprecision(12)<<mkl_spmv_overhead_accumulator<<"sec passed(per loop)"<<endl;
	cout<< "performance:"<<setprecision(12)<<( (double)(2*csrA.nnz + 3*csrA.num_rows)/(double)(1024*1024*1024) )/mkl_spmv_overhead_accumulator <<"GFlops"<<endl;


	// Release matrix handle and deallocate matrix
    mkl_sparse_destroy ( mkl_csrA );

/* MKL SpMV end */

#ifdef DEBUG_MKL_SPMV

    double error_ratio = 1e-2;
    int db_mkl_idx;
    for (db_mkl_idx = 0; db_mkl_idx< csrA.num_rows; db_mkl_idx++)
    {
    	double db_mkl_diff_ratio = (vec_checker.data[db_mkl_idx] - mkl_spmv_result[db_mkl_idx]) / 
    								vec_checker.data[db_mkl_idx];
    	//printf ("checker:%lf, mkl_res:%lf \n", vec_checker.data[db_mkl_idx], mkl_spmv_result[db_mkl_idx]);

    	if (db_mkl_diff_ratio >  error_ratio )
    	{
    		printf ("idx=%d, diff_ratio = %lf , checker = %lf, mkl_res = %lf \n", 
    				db_mkl_idx, db_mkl_diff_ratio, vec_checker.data[db_mkl_idx], mkl_spmv_result[db_mkl_idx]);
    		return 1;
    	} 
    }
#endif /* DEBUG_MKL_SPMV*/


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
	cout<< setprecision(12)<<per_loop_time_cost<<"sec pased(per loop)"<<endl;
	cout<< "performance:"<<setprecision(12)<<( (double)(2*csrA.nnz)/(double)(1024*1024*1024) )/per_loop_time_cost <<"GFlops"<<endl;


#ifdef CHECK_CSR_WITH_COO
	results_comparsion (vec_result, vec_checker);
#endif

	return 0;
}
