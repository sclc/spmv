#include "dense_mat_generator.h"

void generate_dense_mat(denseMat * mat, IDX_TYPE num_rows, IDX_TYPE num_cols,\
                     VAL_TYPE ranMin, VAL_TYPE ranMax)
{
	mat->global_num_row = num_rows;
	mat->global_num_col = num_cols;

    srand (time(NULL));
    VAL_TYPE ranRange = ranMax - ranMin;
    IDX_TYPE totalEle = num_rows * num_cols;
    mat->data = (VAL_TYPE*) calloc (totalEle,sizeof(VAL_TYPE));

    IDX_TYPE idx;
    for (idx=0; idx<totalEle; idx++)
    {
        mat->data[idx] = ranMin + ( (VAL_TYPE)rand() / (VAL_TYPE)RAND_MAX ) * ranRange;
    }

}
