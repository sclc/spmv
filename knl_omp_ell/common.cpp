#include "common.h"

void delete_csrMat(csrMat * mat)
{
	free(mat->col_idx);
	free(mat->row_start);
	free(mat->csrdata);
}
	
void delete_diaMat(diaMat * mat)
{
    free(mat->diag_zone);
    free(mat->diag_offset_id);
}


void delete_cooMat(cooMat * mat)
{
	free(mat->colIdx);
	free(mat->rowIdx);
	free(mat->coodata);
}

void delete_ellMat (ellMat * mat)
{
    free(mat->ell_col_idx);
    free(mat->ell_data);
}

void draw_cooMat_pattern(cooMat mat)
{
	csrMat csrA;

	Converter_Coo2Csr (mat, &csrA);

	IDX_TYPE rowidx, colidx, idx;

	for (rowidx =0; rowidx<csrA.num_rows; rowidx++)
	{
		IDX_TYPE i_start = csrA.row_start[rowidx];
		IDX_TYPE i_end   = csrA.row_start[rowidx+1];

		IDX_TYPE * row_pattern = (IDX_TYPE*)calloc( csrA.num_cols, sizeof(IDX_TYPE));
		for (colidx = i_start; colidx<i_end; colidx++)
		{
			if (csrA.col_idx[colidx] == rowidx)
				// set main diagonal to magic number 66
				row_pattern[ csrA.col_idx[colidx] ] = 66;
			else
				row_pattern[ csrA.col_idx[colidx] ] = csrA.col_idx[colidx] - rowidx;
		}

		for (idx=0; idx < csrA.num_cols; idx++)
		{
			std::cout<< row_pattern[idx] <<"	";
		}
		std::cout<<std::endl;
		free (row_pattern);
	}

	delete_csrMat(&csrA);
}