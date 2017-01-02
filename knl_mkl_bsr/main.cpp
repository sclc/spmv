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
#define EXP_NUM 1
#endif

#ifndef MKL_EXP_NUM
#define MKL_EXP_NUM 1
#endif

// #define MKL_SPARSE_BLAS
#define MKL_INS_EXE_SPARSE_BLAS

#define CHECK_CSR_WITH_COO

#define DEBUG_MKL_SPMV

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
#ifdef MKL_SPARSE_BLAS /*MKL_SPARSE_BLAS*/

	assert( sizeof(double) == sizeof(VAL_TYPE) );
	assert( sizeof(IDX_TYPE) == sizeof(MKL_INT) );

    //convert MKL CSR to MKL BSR
    MKL_INT job[8], info=0, mblk=2, mkl_row_nums = (MKL_INT)csrA.num_rows;
    MKL_INT ldAbsr= mblk * mblk;
    
    job[0] = 0; //the matrix in the CSR format is converted to the BSR format;
    job[1] = 0; //zero-based indexing for the matrix in CSR format is used;
    job[2] = 0; //zero-based indexing for the matrix in the BSR format is used
    job[5] = 1; //all output arrays absr, jab, and iab are filled in for the output storage.

	// declare space for mkl_bsr
	// mat must be square
	MKL_INT num_blocks = (MKL_INT)calculate_bsr_num_blocks (csrA, (IDX_TYPE)mblk);

	printf ("mblk = %d, num_blocks = %d\n",  mblk, num_blocks);


	MKL_INT num_block_row = csrA.num_rows / mblk;
	assert (num_block_row * mblk == csrA.num_rows);

	double * mkl_bsrdata = (double*)calloc( num_blocks*ldAbsr, sizeof(double));
	MKL_INT * mkl_bsr_colInd = (MKL_INT *) calloc ( num_blocks, sizeof(MKL_INT));
	MKL_INT * mkl_csr_rowInd = (MKL_INT *) calloc ( num_block_row+1, sizeof(MKL_INT) );


/*
    void mkl_dcsrbsr (const MKL_INT *job , const MKL_INT *m , const MKL_INT *mblk , const MKL_INT *ldabsr 
					, double *acsr , MKL_INT *ja , MKL_INT *ia 
					, double *absr , MKL_INT *jab , MKL_INT *iab , MKL_INT *info );
*/
	mkl_dcsrbsr(job, &mkl_row_nums, &mblk, &ldAbsr
				, csrA.csrdata, csrA.col_idx, csrA.row_start
				, mkl_bsrdata, mkl_bsr_colInd, mkl_csr_rowInd, &info);
	printf ("info=%d \n", info);
	
	// for (int kk=0; kk<csrA.nnz; kk++)
	// {
	// 	printf("diff: %d \n", mkl_bsr_colInd[kk] -  csrA.col_idx[kk]);
	// }

	// return 1;

    char transa = 'N';
    double mkl_t1, mkl_t2;
    double mkl_spmv_overhead_accumulator = 0.0;

    int mkl_exp_counter;
    double * mkl_spmv_result = NULL;

    for ( mkl_exp_counter = 0; mkl_exp_counter < MKL_EXP_NUM; mkl_exp_counter++)
    {

    	// free and reallocate
    	if (mkl_spmv_result != NULL) free(mkl_spmv_result);
    	mkl_spmv_result = (double *) calloc ( csrA.num_rows, sizeof(double));

    	mkl_t1 = mysecond();
/*
void mkl_cspblas_dbsrgemv (const char *transa , const MKL_INT *m , const MKL_INT *lb 
						, const double *a , const MKL_INT *ia , const MKL_INT *ja 
						, const double *x , double *y );
*/
    	mkl_cspblas_dbsrgemv (&transa, &num_block_row , &mblk
						, mkl_bsrdata , mkl_csr_rowInd , mkl_bsr_colInd 
						, vec.data , mkl_spmv_result);

    	mkl_t2 = mysecond();
    	mkl_spmv_overhead_accumulator += mkl_t2-mkl_t1;
    }

    mkl_spmv_overhead_accumulator /= (double)MKL_EXP_NUM;
    cout<< setprecision(12)<<mkl_spmv_overhead_accumulator<<"sec passed(per loop)"<<endl;
	cout<< "performance:"<<setprecision(12)<<( (double)(2*csrA.nnz)/(double)(1024*1024*1024) )/mkl_spmv_overhead_accumulator <<"GFlops"<<endl;

#endif /* MKL_SPARSE_BLAS */

#ifdef MKL_INS_EXE_SPARSE_BLAS /* MKL_INS_EXE_SPARSE_BLAS */

	assert( sizeof(double) == sizeof(VAL_TYPE) );
	assert( sizeof(IDX_TYPE) == sizeof(MKL_INT) );

    //convert MKL CSR to MKL BSR
    MKL_INT job[8], info=0, mblk=2, mkl_row_nums = (MKL_INT)csrA.num_rows;
    MKL_INT ldAbsr= mblk * mblk;
    
    job[0] = 0; //the matrix in the CSR format is converted to the BSR format;
    job[1] = 0; //zero-based indexing for the matrix in CSR format is used;
    job[2] = 0; //zero-based indexing for the matrix in the BSR format is used
    job[5] = 1; //all output arrays absr, jab, and iab are filled in for the output storage.

	// declare space for mkl_bsr
	// mat must be square
	MKL_INT num_blocks = (MKL_INT)calculate_bsr_num_blocks (csrA, (IDX_TYPE)mblk);

	printf ("mblk = %d, num_blocks = %d\n",  mblk, num_blocks);


	MKL_INT num_block_row = csrA.num_rows / mblk;
	assert (num_block_row * mblk == csrA.num_rows);

	double * mkl_bsrdata = (double*)calloc( num_blocks*ldAbsr, sizeof(double));
	MKL_INT * mkl_bsr_colInd = (MKL_INT *) calloc ( num_blocks, sizeof(MKL_INT));
	MKL_INT * mkl_csr_rowInd = (MKL_INT *) calloc ( num_block_row+1, sizeof(MKL_INT) );


/*
    void mkl_dcsrbsr (const MKL_INT *job , const MKL_INT *m , const MKL_INT *mblk , const MKL_INT *ldabsr 
					, double *acsr , MKL_INT *ja , MKL_INT *ia 
					, double *absr , MKL_INT *jab , MKL_INT *iab , MKL_INT *info );
*/
	mkl_dcsrbsr(job, &mkl_row_nums, &mblk, &ldAbsr
				, csrA.csrdata, csrA.col_idx, csrA.row_start
				, mkl_bsrdata, mkl_bsr_colInd, mkl_csr_rowInd, &info);
	printf ("info=%d \n", info);


	struct matrix_descr mkl_descrA;
	// Create matrix descriptor
    mkl_descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

	sparse_matrix_t mkl_bsrA;

/*
sparse_status_t mkl_sparse_d_create_bsr (sparse_matrix_t *A, sparse_index_base_t indexing, sparse_layout_t block_layout
										, MKL_INT rows, MKL_INT cols, MKL_INT block_size
										, MKL_INT *rows_start, MKL_INT *rows_end, MKL_INT *col_indx, double *values);
*/
	mkl_sparse_d_create_bsr(&mkl_bsrA, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR
							, num_block_row, num_block_row, mblk
							, mkl_csr_rowInd, mkl_csr_rowInd+1, mkl_bsr_colInd, mkl_bsrdata);

	// Analyze sparse matrix; choose proper kernels and workload balancing strategy
    mkl_sparse_optimize ( mkl_bsrA );

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

    	mkl_t1 = mysecond();
    	
    	// y := alpha*op(A)*x + beta*y
    	mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE,
          	            mkl_alpha,
          	            mkl_bsrA,
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
    mkl_sparse_destroy ( mkl_bsrA );

#endif /* MKL_INS_EXE_SPARSE_BLAS */

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
    printf("DEBUG_MKL_SPMV passed \n");
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
