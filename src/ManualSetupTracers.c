//  Copyright (C) 2018, Ludwig Jens Papenfort
//                      <papenfort@th.physik.uni-frankfurt.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define length(X) (sizeof(X)/sizeof(*X))
#define MYH5CHECK(val) assert(val >= 0)

// Stores all the tracers info
typedef struct {
    CCTK_REAL * x_coords;
    CCTK_REAL * y_coords;
    CCTK_REAL * z_coords;
    const int npoints;
} myTracers;

static myTracers * trcrs = NULL;

static void trcrs_free() {
    free(trcrs->x_coords);
    free(trcrs->y_coords);
    free(trcrs->z_coords);
    free(trcrs);
    trcrs = NULL;
}

int TracerFlow_Read_Positions(const char *filename) {

    double * coords = NULL;
    int rank;

    int const npoints;

    /* Decl. HDF vars */
    hid_t fid, dataid, fapl;
    hsize_t *dim = NULL, size;
    hid_t datatype, dataspace;
    herr_t hErrVal;

    printf("Reading from file: %s\n",filename);

    /* Load the library -- not required for most platforms. */
    hErrVal = H5open(); MYH5CHECK(hErrVal);

		/*
    * Open file, thereby enforcing proper file close
    * semantics
    */
    fapl = H5Pcreate(H5P_FILE_ACCESS); MYH5CHECK(fapl);
    H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);
    fid = H5Fopen(filename, H5F_ACC_RDONLY,fapl); MYH5CHECK(fid);
    H5Pclose(fapl);

		/*
    * open and read the coords dataset with the initial positions
    */
    dataid = H5Dopen(fid,"/initial_pos/coords",H5P_DEFAULT); MYH5CHECK(dataid);
    dataspace = H5Dget_space(dataid); MYH5CHECK(dataspace);
    rank = H5Sget_simple_extent_ndims(dataspace); MYH5CHECK(rank);
    dim = malloc(rank*sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dim, NULL);
    size = dim[0]*dim[1];
    coords = malloc(size*sizeof(double));

    /*
    * Coordinates are stored as double
    */
    H5Dread(dataid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords);
    H5Dclose(dataid);
    H5Sclose(dataspace);
    H5Fclose(fid);

		printf("Successfully read data from file: %s\n",filename);

    /*
    * Allocate memory for the coordinate vectors
    */
    trcrs->x_coords = malloc(dim[0]*sizeof(CCTK_REAL)); assert(trcrs->x_coords);
    trcrs->y_coords = malloc(dim[0]*sizeof(CCTK_REAL)); assert(trcrs->y_coords);
    trcrs->z_coords = malloc(dim[0]*sizeof(CCTK_REAL)); assert(trcrs->z_coords);
    *(int *)&trcrs->npoints  = dim[0];

    int index = 0;

		printf("Begin setting up tracers\n");
    for(int i = 0; i < size; i = i+3)
    {
        //printf("Index %i \n", index);
        trcrs->x_coords[index] = coords[i];
        trcrs->y_coords[index] = coords[i+1];
        trcrs->z_coords[index] = coords[i+2];
        /*
         * DEBUG
         *
        if (index < 50){
            printf("X %8.2f \n",trcrs->x_coords[index]);
            printf("Y %8.2f \n",trcrs->y_coords[index]);
            printf("Z %8.2f \n",trcrs->z_coords[index]);
        }
        */
        index++;
    }
    // Free coords and dim
    free(coords);
    free(dim);
    coords = NULL;
    dim = NULL;

    //printf("X %8.2f Y %8.2f Z %8.2f\n",trcrs->x_coords[0],trcrs->y_coords[0],trcrs->z_coords[0]);
    printf("TracerFlow:: Number of tracers %i \n", index);
    printf("TracerFlow:: Done reading the file.\n");

    return 0;
}

void TracerFlow_Manual_SetupTracers(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int  Err = 0;

    CCTK_VInfo(CCTK_THORNSTRING, "Enter Manual_SetupTracers");

    assert(NULL == trcrs);
    trcrs = malloc(sizeof(myTracers)); assert(trcrs);

    Err =  TracerFlow_Read_Positions(tracers_pos_file); assert(Err==0);

    int group = CCTK_GroupIndex("TracerFlow::tracer_evol");
    cGroupDynamicData data;
    (void)CCTK_GroupDynamicData (cctkGH, group, &data);
    int myntracers = data.lsh[0];

    int myoffset = data.lbnd[0];
    const int siz = myntracers;

    int npoints = trcrs->npoints;

    /* For now we assume that the tracers are "massless".
     * We assign to them only an initial position.
     */

    for(int i = 0 ; i < siz; ++i) {
                tracer_x[i] = trcrs->x_coords[i+myoffset];
                tracer_y[i] = trcrs->y_coords[i+myoffset];
                tracer_z[i] = trcrs->z_coords[i+myoffset];
    };

    trcrs_free();

    CCTK_VInfo(CCTK_THORNSTRING, "Exit Manual_SetupTracers");
}
