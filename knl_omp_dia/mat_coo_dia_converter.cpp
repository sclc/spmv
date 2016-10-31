#include "mat_coo_dia_converter.h"

void Converter_Coo2Dia (cooMat src, diaMat * target)
{
    target->num_rows = src.num_rows;
    target->num_cols = src.num_cols;
    target->nnz      = src.nnz;
    
    // scan COO matrix 
    IDX_TYPE pos_offset_zero = src.num_rows - 1;
    IDX_TYPE length_diags_lable = 2*src.num_rows - 1;

    IDX_TYPE * diag_lable = (IDX_TYPE *)calloc(length_diags_lable, sizeof(IDX_TYPE) );
    VAL_TYPE ** diag_list = (VAL_TYPE **)calloc(length_diags_lable, sizeof(VAL_TYPE *) );
    IDX_TYPE diag_list_idx;
    for (diag_list_idx = 0; diag_list_idx < length_diags_lable; diag_list_idx++)
        diag_list[diag_list_idx] = NULL;



    IDX_TYPE rowid, colid, pos_diag_label;
    VAL_TYPE val;
    IDX_TYPE id;

    for(id=0; id<src.nnz;id++)
    {
        rowid = src.rowIdx[id];
        colid = src.colIdx[id];
        val   = src.coodata[id];
        pos_diag_label = colid - rowid + pos_offset_zero;

        if (diag_lable[pos_diag_label] != DIA_MAGIC_NUM)
        {
#ifdef COO2DIA_DEBUG_A
            assert(diag_list[pos_diag_label] == NULL);
#endif /*COO2DIA_DEBUG_A*/
            // assume this is a square matrix
            // calloc() zero-initializes the buffer
            diag_list[pos_diag_label] = (VAL_TYPE *)calloc( src.num_rows, sizeof(VAL_TYPE));
            diag_lable[pos_diag_label] = DIA_MAGIC_NUM;
        }
        
        diag_list[pos_diag_label][rowid] = val;

    }

    // fill DIA matrix
    IDX_TYPE num_diags = 0;
    for (id = 0; id< length_diags_lable; id++)
    {
        if (diag_lable[id] == DIA_MAGIC_NUM)
            num_diags++;
    }

    target->num_diags = num_diags;

    IDX_TYPE * diag_offset_id = (IDX_TYPE *)calloc( target->num_diags, sizeof(IDX_TYPE) );
    // you must use calloc for diag_zone, otherwise you need to explicitly zero the elements
    VAL_TYPE * diag_zone      = (VAL_TYPE *)calloc( target->num_diags * target->num_rows, 
                                sizeof(VAL_TYPE) );

    IDX_TYPE pos_counter=0;
    for(id=0; id<length_diags_lable;id++)
    {
        if (diag_lable[id] == DIA_MAGIC_NUM)
        {
            diag_offset_id[pos_counter] = id - pos_offset_zero;
            memcpy ( diag_zone + target->num_rows*pos_counter, 
                        diag_list[id], 
                        sizeof(VAL_TYPE) * target->num_rows);
            pos_counter++;
        }
    }
#ifdef COO2DIA_DEBUG_B
    assert( pos_counter == target->num_diags);
#endif /*COO2DIA_DEBUG_B*/

    target->diag_offset_id = diag_offset_id;
    target->diag_zone      = diag_zone;

    // delete pos_offset_zero and diag_list
    free(diag_lable);
    for (diag_list_idx = 0; diag_list_idx<length_diags_lable;diag_list_idx++)
    {
        if (diag_list[diag_list_idx] != NULL)
        {
            free( diag_list[diag_list_idx] );
            diag_list[diag_list_idx] = NULL;
        }
    }
    free(diag_list);

}
