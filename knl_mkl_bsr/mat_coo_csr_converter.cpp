#include "mat_coo_csr_converter.h"

void Converter_Coo2Csr (cooMat src, csrMat * target)
{
    IDX_TYPE idx, cumsum, last;
    
    target->num_rows = src.num_rows;
    target->num_cols = src.num_cols;
    target->nnz      = src.nnz;
    
    target->row_start = (IDX_TYPE *) calloc (target->num_rows + 1, sizeof(IDX_TYPE) );
    target->col_idx   = (IDX_TYPE *) calloc (target->nnz, sizeof (IDX_TYPE));
    target->csrdata   = (VAL_TYPE *) calloc (target->nnz, sizeof(VAL_TYPE));
    
    for (idx = 0; idx < src.num_rows; idx++)
        target->row_start[idx] = 0;

    for (idx = 0; idx < src.nnz; idx++)
        target->row_start[src.rowIdx[idx]]++;


    //cumsum the nnz per row to get Bp[]
    for(idx = 0, cumsum = 0; idx < src.num_rows; idx++){     
        IDX_TYPE temp = target->row_start[idx];
        target->row_start[idx] = cumsum;
        cumsum += temp;
    }
    target->row_start[src.num_rows] = src.nnz;

    // final csr, in row section unsorted by column or essentially sorted for 
    // coo sorted by col
    for(idx = 0; idx < src.nnz; idx++){
        IDX_TYPE row  = src.rowIdx[idx];
        IDX_TYPE dest = target->row_start[row];

        target->col_idx[dest] = src.colIdx [idx];
        target->csrdata[dest] = src.coodata[idx];
        // printf ("%lf\n", target->csrdata[dest]);

        target->row_start[row]++;
    }
    
    last = 0;
    // set back row_start values
    for( idx = 0; idx <= src.num_rows; idx++){
        IDX_TYPE temp = target->row_start[idx];
        target->row_start[idx]  = last;
        last   = temp;
    }
    
}