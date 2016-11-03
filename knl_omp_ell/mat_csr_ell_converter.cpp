#include "mat_csr_ell_converter.h"

void Converter_Csr2Ell (csrMat src, ellMat * target)
{
    IDX_TYPE row_max_width=0;

    IDX_TYPE idx;
    for (idx=0; idx<src.num_rows; idx++)
    {
        IDX_TYPE temp;
        temp = src.row_start[idx+1] - src.row_start[idx];
        if (temp > row_max_width)
        {
            row_max_width = temp;
        }
    }
    target->ell_row_length = row_max_width;
    target->ell_num_rows = src.num_rows;
    target->ell_num_cols = src.num_cols;
    target->ell_nnz      = src.nnz;

    target->ell_col_idx = (IDX_TYPE*)calloc(target->ell_num_rows * target->ell_row_length
                                        , sizeof(IDX_TYPE));
    target->ell_data    = (VAL_TYPE*)calloc(target->ell_num_rows * target->ell_row_length
                                        , sizeof(VAL_TYPE));

    IDX_TYPE idx_row, idx_col;
    for (idx_row=0; idx_row<src.num_rows; idx_row++)
    {
        // declare this inside loop, for the convience of parallization with OpenMP
        IDX_TYPE id_start = src.row_start[idx_row];
        IDX_TYPE id_end   = src.row_start[idx_row+1];

        IDX_TYPE col_counter=0;
        for (idx_col=id_start; idx_col<id_end; idx_col++)
        {
            IDX_TYPE rowid = idx_row;
            IDX_TYPE colid = src.col_idx[idx_col];
            VAL_TYPE val   = src.csrdata[idx_col];

            IDX_TYPE pos = rowid + col_counter * target->ell_num_rows;
            target->ell_col_idx[pos] = colid;
            target->ell_data[pos]    = val;
            col_counter++;
        }
    }

}