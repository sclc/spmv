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

IDX_TYPE calculate_bsr_num_blocks (csrMat mat, IDX_TYPE block_size) // mat must be square
{
	IDX_TYPE num_blocks=0;
	IDX_TYPE num_row_col = mat.num_rows;
	assert (mat.num_cols == mat.num_rows);

	long max_num_block_1D = (long)num_row_col / block_size;
	assert (max_num_block_1D * block_size == num_row_col); // test if block_size divides num_row_col 

	printf ("max_num_block_1D = %d\n", max_num_block_1D);

	int * histgram_blocks = (int*) calloc ( max_num_block_1D * max_num_block_1D, sizeof(int));
	assert (histgram_blocks != NULL);

	IDX_TYPE idx_row, idx_col;
	IDX_TYPE mat_col_idx;
	IDX_TYPE histgram_blocks_row_idx, histgram_blocks_col_idx;

	for (idx_row = 0; idx_row < num_row_col; idx_row++)
	{
		IDX_TYPE start_row = mat.row_start[idx_row];
		IDX_TYPE end_row   = mat.row_start[idx_row+1];

		for (idx_col = start_row; idx_col<end_row; idx_col++)
		{
			mat_col_idx = mat.col_idx[idx_col];

			histgram_blocks_row_idx = idx_row / block_size;
			histgram_blocks_col_idx = mat_col_idx / block_size;

			if (histgram_blocks[histgram_blocks_row_idx * max_num_block_1D + histgram_blocks_col_idx]++ == 0 )
			{
				num_blocks++;
			}
		}
	}

	free (histgram_blocks);

	return num_blocks;
}
