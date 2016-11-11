#include "DataTypes.h"
#include "mmio.h"
#include <string>
#include <assert.h>

using namespace std;


#ifndef VAL_TYPE
#define VAL_TYPE double
#endif

#ifndef IDX_TYPE
#define IDX_TYPE int
#endif

#ifdef	__cplusplus
extern "C" {
#endif

void readMtx_coo(char* path, char* name, cooMat mtr, matInfo info);

#ifdef	__cplusplus
}
#endif

void readMtx_info_and_coo(string path, string name, matInfo* info, cooMat* mat);
void readMtx_info_and_ordered_coo(string path, string name, matInfo * info, cooMat * mat);
