#include <iostream>
#include <omp.h>
#include <string>
#include <assert.h>
#include <iomanip>

#include "DataTypes.h"
#include "readMTX.h"
#include "mat_coo_csr_converter.h"

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

//#define DEBUG_A
//#define DEBUG_B
#define DEBUG_C

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


	readMtx_info_and_coo(mat_path, mat_filename, &mat_info, &cooA);
#ifdef DEBUG_B
	
	IDX_TYPE db_idx;
	cout<< "#row="<<cooA.num_rows<<", #col="<<cooA.num_cols<<", nnz="<<cooA.nnz<<endl;

	for (db_idx=0;db_idx<cooA.nnz;db_idx++)
	{
		cout<<cooA.rowIdx[db_idx]<<", "<<cooA.colIdx[db_idx]<<", "<<setprecision(10)<<cooA.coodata[db_idx]<<endl;
	}
#endif

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

#endif

	return 0;
}
