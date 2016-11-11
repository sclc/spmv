//
#include "readMTX.h"

void readMtx_coo(char* path, char* name, cooMat mtr, matInfo info) {

}

int comparator ( const void *p, const void *q)
{
    IDX_TYPE l = ((matElement *) p) -> colid;
    IDX_TYPE r = ((matElement *) q) -> colid;

    return (l-r);
}

void readMtx_info_and_coo(string path, string name, matInfo * info, cooMat * mat) {
    int ret_code;
    int m, n, nnz;
    MM_typecode matcode;
    FILE *fid;
    IDX_TYPE idx;
    string full_path = path+name;

    fid = fopen(full_path.c_str(), "r");

    if (fid == NULL) {
        printf("Unable to open file %s\n", full_path.c_str());
        exit(1);
    }

    if (mm_read_banner(fid, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (!mm_is_valid(matcode)) {
        printf("Invalid Matrix Market file.\n");
        exit(1);
    }

    if (!((mm_is_real(matcode) || mm_is_integer(matcode) || mm_is_pattern(matcode)) && mm_is_coordinate(matcode) && mm_is_sparse(matcode))) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        printf("Only sparse real-valued or pattern coordinate matrices are supported\n");
        exit(1);
    }
    ////////////////////////////////////////

    if ((ret_code = mm_read_mtx_crd_size(fid, &m, &n, &nnz)) != 0)
        exit(1);

    info->num_rows = m;
    info->num_cols = n;
    info->nnz = nnz;
    // printf ("%ld, %ld, %ld", info->num_rows, info->num_cols, info->nnz);
    // exit(0);
    
    mat->rowIdx = (IDX_TYPE *) calloc(info->nnz, sizeof (IDX_TYPE));
    mat->colIdx = (IDX_TYPE *) calloc(info->nnz, sizeof (IDX_TYPE));
    mat->coodata = (VAL_TYPE *) calloc(info->nnz, sizeof (VAL_TYPE));
    // printf ("so far so good K\n");
    // exit(0);


//    printf("Reading sparse matrix from file ( %s ):\n", concatStr(path, name));
    printf("Reading sparse matrix from file ( %s ):\n", full_path.c_str());

    if (mm_is_pattern(matcode)) {
        // pattern matrix defines sparsity pattern, but not values
        for (idx = 0; idx < info->nnz; idx++) {
//            assert(fscanf(fid, " %ld %ld \n", &(mat->rowIdx[idx]), &(mat->colIdx[idx])) == 2);
            assert(fscanf(fid, " %d %d \n", &(mat->rowIdx[idx]), &(mat->colIdx[idx])) == 2);
            // possible bug place ,%d or %ld
            // printf ("%ld, %ld\n", mat->rowIdx[idx], mat->colIdx[idx]);
            // exit(0);

            mat->rowIdx[idx]--; //adjust from 1-based to 0-based indexing
            mat->colIdx[idx]--;
            mat->coodata[idx] = 1.0; //use value 1.0 for all nonzero entries 
        }
    } else if (mm_is_real(matcode) || mm_is_integer(matcode)) {
        for (idx = 0; idx < info->nnz; idx++) {
            IDX_TYPE I, J;
            VAL_TYPE V; // always read in a VAL_TYPE and convert later if necessary

//            assert(fscanf(fid, " %ld %ld %lf \n", &I, &J, &V) == 3);
            assert(fscanf(fid, " %d %d %lf \n", &I, &J, &V) == 3);

            // printf("%ld, %ld, %lf\n", I, J, V);
        

            mat->rowIdx[idx] = I - 1;
            mat->colIdx[idx] = J - 1;
            mat->coodata[idx] = V;
        }
    } else {
        printf("Unrecognized data type\n");
        exit(1);
    }


    fclose(fid);
    printf(" done reading\n");

    if (mm_is_symmetric(matcode)) { //duplicate off diagonal entries
        IDX_TYPE off_diagonals = 0;
        // for spmv feature
        for (idx = 0; idx < info->nnz; idx++) {
            if (mat->rowIdx[idx] != mat->colIdx[idx])
                off_diagonals++;
        }

        IDX_TYPE true_nonzeros = 2 * off_diagonals + (info->nnz - off_diagonals);

        IDX_TYPE* new_I = (IDX_TYPE *) calloc(true_nonzeros, sizeof (IDX_TYPE)); 
        IDX_TYPE* new_J = (IDX_TYPE *) calloc(true_nonzeros, sizeof (IDX_TYPE)); 
        VAL_TYPE* new_V = (VAL_TYPE *) calloc(true_nonzeros, sizeof (VAL_TYPE)); 

        IDX_TYPE ptr = 0;
        for (idx = 0; idx < info->nnz; idx++) {
            if (mat->rowIdx[idx] != mat->colIdx[idx]) {
                new_I[ptr] = mat->rowIdx[idx];
                new_J[ptr] = mat->colIdx[idx];
                new_V[ptr] = mat->coodata[idx];
                ptr++;
                new_J[ptr] = mat->rowIdx[idx];
                new_I[ptr] = mat->colIdx[idx];
                new_V[ptr] = mat->coodata[idx];
                ptr++;
            } else {
                new_I[ptr] = mat->rowIdx[idx];
                new_J[ptr] = mat->colIdx[idx];
                new_V[ptr] = mat->coodata[idx];
                ptr++;
            }
        }
        assert (ptr == true_nonzeros);
        free(mat->rowIdx);
        free(mat->colIdx);
        free(mat->coodata);

        mat->rowIdx = new_I;
        mat->colIdx = new_J;
        mat->coodata = new_V;
        mat->nnz = true_nonzeros;
        mat->num_rows = info->num_rows;
        mat->num_cols = info->num_cols;
        
        info->nnz = true_nonzeros;

    } //end symmetric case

}

void readMtx_info_and_ordered_coo(string path, string name, matInfo * info, cooMat * mat) 
{
    int ret_code;
    int m, n, nnz;
    MM_typecode matcode;
    FILE *fid;
    IDX_TYPE idx;
    string full_path = path+name;

    fid = fopen(full_path.c_str(), "r");

    if (fid == NULL) {
        printf("Unable to open file %s\n", full_path.c_str());
        exit(1);
    }

    if (mm_read_banner(fid, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (!mm_is_valid(matcode)) {
        printf("Invalid Matrix Market file.\n");
        exit(1);
    }

    if (!((mm_is_real(matcode) || mm_is_integer(matcode) || mm_is_pattern(matcode)) && mm_is_coordinate(matcode) && mm_is_sparse(matcode))) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        printf("Only sparse real-valued or pattern coordinate matrices are supported\n");
        exit(1);
    }
    ////////////////////////////////////////

    if ((ret_code = mm_read_mtx_crd_size(fid, &m, &n, &nnz)) != 0)
        exit(1);

    info->num_rows = m;
    info->num_cols = n;
    info->nnz = nnz;
    // printf ("%ld, %ld, %ld", info->num_rows, info->num_cols, info->nnz);
    // exit(0);
    
    // mat->rowIdx = (IDX_TYPE *) calloc(info->nnz, sizeof (IDX_TYPE));
    // mat->colIdx = (IDX_TYPE *) calloc(info->nnz, sizeof (IDX_TYPE));
    // mat->coodata = (VAL_TYPE *) calloc(info->nnz, sizeof (VAL_TYPE));
    matElement * matList_1 = (matElement *)malloc( sizeof(matElement) * info->nnz);

    // printf ("so far so good K\n");
    // exit(0);


//    printf("Reading sparse matrix from file ( %s ):\n", concatStr(path, name));
    printf("Reading sparse matrix from file ( %s ):\n", full_path.c_str());

    if (mm_is_pattern(matcode)) {
        // pattern matrix defines sparsity pattern, but not values
        for (idx = 0; idx < info->nnz; idx++) {
//            assert(fscanf(fid, " %ld %ld \n", &(mat->rowIdx[idx]), &(mat->colIdx[idx])) == 2);
//            assert(fscanf(fid, " %d %d \n", &(mat->rowIdx[idx]), &(mat->colIdx[idx])) == 2);
            assert(fscanf(fid, " %d %d \n", &(matList_1+idx)->rowid, &(matList_1+idx)->colid)  == 2);
            // possible bug place ,%d or %ld
            // printf ("%ld, %ld\n", mat->rowIdx[idx], mat->colIdx[idx]);
            // exit(0);

            // mat->rowIdx[idx]--; //adjust from 1-based to 0-based indexing
            // mat->colIdx[idx]--;
            // mat->coodata[idx] = 1.0; //use value 1.0 for all nonzero entries 
            (matList_1+idx)->rowid --;
            (matList_1+idx)->colid --;
            (matList_1+idx)->val = 1.0;
        }
    } else if (mm_is_real(matcode) || mm_is_integer(matcode)) {
        for (idx = 0; idx < info->nnz; idx++) {
            IDX_TYPE I, J;
            VAL_TYPE V; // always read in a VAL_TYPE and convert later if necessary

//            assert(fscanf(fid, " %ld %ld %lf \n", &I, &J, &V) == 3);
           assert(fscanf(fid, " %d %d %lf \n", &I, &J, &V) == 3);
            
            // printf("%ld, %ld, %lf\n", I, J, V);
        

            // mat->rowIdx[idx] = I - 1;
            // mat->colIdx[idx] = J - 1;
            // mat->coodata[idx] = V;
           (matList_1+idx)->rowid = I - 1;
           (matList_1+idx)->colid = J - 1;
           (matList_1+idx)->val   = V;

        }
    } else {
        printf("Unrecognized data type\n");
        exit(1);
    }


    fclose(fid);
    printf(" done reading\n");

    if (mm_is_symmetric(matcode)) { //duplicate off diagonal entries
        IDX_TYPE off_diagonals = 0;
        // for spmv feature
        for (idx = 0; idx < info->nnz; idx++) {
            // if (mat->rowIdx[idx] != mat->colIdx[idx])
            if ((matList_1+idx)->rowid != (matList_1+idx)->colid)
                off_diagonals++;
        }

        IDX_TYPE true_nonzeros = 2 * off_diagonals + (info->nnz - off_diagonals);

        // IDX_TYPE* new_I = (IDX_TYPE *) calloc(true_nonzeros, sizeof (IDX_TYPE)); 
        // IDX_TYPE* new_J = (IDX_TYPE *) calloc(true_nonzeros, sizeof (IDX_TYPE)); 
        // VAL_TYPE* new_V = (VAL_TYPE *) calloc(true_nonzeros, sizeof (VAL_TYPE)); 
        
        matElement * matList_final = (matElement *)malloc( sizeof(matElement) * true_nonzeros);

        IDX_TYPE ptr = 0;
        for (idx = 0; idx < info->nnz; idx++) {
            if ((matList_1+idx)->rowid != (matList_1+idx)->colid) {
                // new_I[ptr] = mat->rowIdx[idx];
                // new_J[ptr] = mat->colIdx[idx];
                // new_V[ptr] = mat->coodata[idx];
                (matList_final + ptr)->rowid = (matList_1 + idx)->rowid;
                (matList_final + ptr)->colid = (matList_1 + idx)->colid;
                (matList_final + ptr)->val = (matList_1 + idx)->val;

                ptr++;
                // new_J[ptr] = mat->rowIdx[idx];
                // new_I[ptr] = mat->colIdx[idx];
                // new_V[ptr] = mat->coodata[idx];

                (matList_final + ptr)->rowid = (matList_1 + idx)->colid;
                (matList_final + ptr)->colid = (matList_1 + idx)->rowid;
                (matList_final + ptr)->val = (matList_1 + idx)->val;

                ptr++;
            } else {
                // new_I[ptr] = mat->rowIdx[idx];
                // new_J[ptr] = mat->colIdx[idx];
                // new_V[ptr] = mat->coodata[idx];
                (matList_final + ptr)->rowid = (matList_1 + idx)->rowid;
                (matList_final + ptr)->colid = (matList_1 + idx)->rowid;
                (matList_final + ptr)->val = (matList_1 + idx)->val;

                ptr++;
            }
        }
        assert (ptr == true_nonzeros);
        // free(mat->rowIdx);
        // free(mat->colIdx);
        // free(mat->coodata);
        free(matList_1);


        // re-ordering here
        // TO-DO
        qsort (matList_final,  true_nonzeros, sizeof(matElement), comparator);

        IDX_TYPE* new_I = (IDX_TYPE *) calloc(true_nonzeros, sizeof (IDX_TYPE)); 
        IDX_TYPE* new_J = (IDX_TYPE *) calloc(true_nonzeros, sizeof (IDX_TYPE)); 
        VAL_TYPE* new_V = (VAL_TYPE *) calloc(true_nonzeros, sizeof (VAL_TYPE)); 

        for (idx=0; idx<true_nonzeros; idx++)
        {
            new_I[idx] = (matList_final + idx)->rowid;
            new_J[idx] = (matList_final + idx)->colid;
            new_V[idx] = (matList_final + idx)->val;
        }

        mat->rowIdx = new_I;
        mat->colIdx = new_J;
        mat->coodata = new_V;
        mat->nnz = true_nonzeros;
        mat->num_rows = info->num_rows;
        mat->num_cols = info->num_cols;
        
        info->nnz = true_nonzeros;

        free(matList_final);

    } //end symmetric case

}

