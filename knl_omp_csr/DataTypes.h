/* 
 * File:   DataTypes.h
 * Author: scl
 *
 * Created on November 27, 2013, 5:29 PM
 */

#ifndef DATATYPES_H
#define DATATYPES_H

#include <map>
#include <vector>

#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

// #ifdef  __cplusplus
// extern "C" {
// #endif

typedef struct {

    IDX_TYPE num_rows;
    IDX_TYPE num_cols;
    IDX_TYPE nnz;

} matInfo;

typedef struct {

    IDX_TYPE rowid;
    IDX_TYPE colid;
    VAL_TYPE val;

} matElement;

typedef struct
{
                 
    /* Dimensions, number of nonzeros 
     * (m == n for square, m < n for a local "slice") */
    // local matrix property
    IDX_TYPE num_rows, num_cols, nnz;
                     
    /* Starts of rows owned by local processor
     * row_start[0] == 0 == offset of first nz in row start[MY_THREAD] 
     */
    IDX_TYPE *row_start;
                                          
    /* Column indices and values of matrix elements at local processor */
    IDX_TYPE *col_idx;

    VAL_TYPE *csrdata;
                                          
} csrMat;


typedef struct 
{
    IDX_TYPE num_rows, num_cols, nnz;
    IDX_TYPE * rowIdx;
    IDX_TYPE * colIdx;
    VAL_TYPE * coodata;
} cooMat;

typedef struct
{

    IDX_TYPE global_num_row;
    IDX_TYPE global_num_col;
        
    VAL_TYPE * data;
} denseMat;

typedef struct
{
    IDX_TYPE num_rows, num_cols, nnz, num_diags;
    IDX_TYPE * diag_offset_id;
    VAL_TYPE * diag_zone;
    
}diaMat;

typedef struct 
{
    IDX_TYPE ell_num_rows, ell_num_cols, ell_nnz;
    IDX_TYPE ell_row_length;
    IDX_TYPE * ell_col_idx;
    VAL_TYPE * ell_data;
    
}ellMat;


// #ifdef  __cplusplus
// }
// #endif

#endif  /* DATATYPES_H */
