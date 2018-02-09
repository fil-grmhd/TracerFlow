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

#define H5Acreate_vers 2
#define H5Gcreate_vers 2
#define H5Gopen_vers   2

#include <hdf5.h>
#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "OutputUtils.h"

#include <util_String.h>

static void WriteHeader(hid_t file_id, Header * header)
{
  hid_t group_id, dspace_id, attr_id;

  hsize_t const size_one = 1;
  hsize_t const size_six = 6;

  group_id = H5Gcreate(file_id, "Header", H5P_DEFAULT, H5P_DEFAULT,
          H5P_DEFAULT); MYH5CHECK(group_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "BoxSize", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->box_size);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_Cooling", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_cooling);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_DoublePrecision", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_double_precision);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_Feedback", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_feedback);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_IC_Info", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_ic_info);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_Metals", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_metals);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_Sfr", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_sfr);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_StellarAge", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_stellar_age);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "HubbleParam", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->hubble_param);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_six, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "MassTable", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->mass_table);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "NumFilesPerSnapshot", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->num_files_per_snapshot);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_six, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "NumPart_ThisFile", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->num_part_this_file);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_six, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "NumPart_Total", H5T_NATIVE_UINT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_UINT, &header->num_part_total);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_six, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "NumPart_Total_HighWord", H5T_NATIVE_UINT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_UINT, &header->num_part_total_high_word);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Omega0", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->Omega0);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "OmegaLambda", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->OmegaLambda);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Redshift", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->redshift);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Time", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->time);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  H5Gclose(group_id);
}

static void WriteIDs(hid_t file_id, unsigned int const * IDs, int const siz)
{
  hid_t group_id, dspace_id, dset_id;
  herr_t h5err;

  hsize_t hsiz = siz;
  group_id = H5Gcreate(file_id, "PartType0", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT); MYH5CHECK(group_id);
  dspace_id = H5Screate_simple(1, &hsiz, NULL); MYH5CHECK(dspace_id);
  dset_id = H5Dcreate2(group_id, "ParticleIDs", H5T_NATIVE_UINT, dspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(dset_id);
  h5err = H5Dwrite(dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, IDs); MYH5CHECK(h5err);

  H5Dclose(dset_id);
  H5Sclose(dspace_id);
  H5Gclose(group_id);
}

static void WriteField(hid_t file_id, float const * var, char const * name,
    int ncomp, int const siz)
{
  hid_t group_id, dspace_id, dset_id;
  herr_t h5err;

  hsize_t hsiz = siz;
  group_id = H5Gopen(file_id, "/PartType0", H5P_DEFAULT); MYH5CHECK(group_id);
  if(ncomp == 1)
  {
    dspace_id = H5Screate_simple(1, &hsiz, NULL); MYH5CHECK(dspace_id);
  }
  else if(ncomp == 3)
  {
    hsize_t dspace_siz[2] = {hsiz, 3};
    dspace_id = H5Screate_simple(2, dspace_siz, NULL); MYH5CHECK(dspace_id);
  }
  else
  {
    char msg[STRLEN];
    snprintf(msg, STRLEN, "Trying to write a field of unsupported ncomp: %d",
        ncomp);
    CCTK_ERROR(msg);
  }
  dset_id = H5Dcreate2(group_id, name, H5T_NATIVE_FLOAT, dspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(dset_id);
  h5err = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, var); MYH5CHECK(h5err);

  H5Dclose(dset_id);
  H5Sclose(dspace_id);
  H5Gclose(group_id);
}

static void WriteVisItFile(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char *fn;
  FILE *fh;
  static int first_time = 1;
  int truncate, nprocs;

  nprocs = CCTK_nProcs(cctkGH);

  fn = NULL;
  Util_asprintf(&fn, "%s/%s.visit", out_dir, tracersfile);
  assert(fn != NULL);

  // open file and write header if needed
  truncate = first_time && IO_TruncateOutputFiles(cctkGH);
  fh = fopen(fn, truncate ? "w" : "a");
  if(fh == NULL)
  {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not open file '%s' for writing: %s\n", fn,
                strerror(errno));
  }
  if(truncate)
    fprintf(fh, "!NBLOCKS %d\n", one_file_per_timestep ? 1 : nprocs);

  // one line per output file
  if(one_file_per_timestep)
  {
    fprintf(fh, "%s.%d.hdf5\n", tracersfile, cctk_iteration);
  }
  else
  {
    for(int i = 0 ; i < nprocs ; i++)
    {
      fprintf(fh, "%s.%d.%d.hdf5\n", tracersfile, cctk_iteration, i);
    }
  }

  fclose(fh);
  free(fn);

  first_time = 0;
}

static void WriteSingleVariableChunked(cGH const * cctkGH, int ntracers_on_this_proc,
    int components, CCTK_REAL * const restrict * restrict data, char const * label,
    hid_t file_id)
{
  float * local_data = malloc(ntracers_on_this_proc * components * sizeof(float));
  for(int i = 0; i < ntracers_on_this_proc; ++i)
    for(int c = 0; c < components; ++c)
      local_data[i*components + c] = (float)data[c][i];
  WriteField(file_id, local_data, label, components, ntracers_on_this_proc);
  free(local_data);
}

static void CollectAndWriteSingleVariable(const cGH *cctkGH, int ntracers_on_this_proc,
                                          int components, CCTK_REAL * const restrict * restrict data,
                                          const char *label,
                                          hid_t file_id)
{
  DECLARE_CCTK_PARAMETERS;

  const int nProcs = CCTK_nProcs(cctkGH);
  const int myProc = CCTK_MyProc(cctkGH);

  // offsets and counts stay zero on the non-Root procs
  int *offsets = calloc(nProcs, sizeof(*offsets));
  int *counts = calloc(nProcs, sizeof(*counts));
  float *local_data, *global_data;

  int count = ntracers_on_this_proc * components;
  int count_total;

  const MPI_Comm comm = CCTK_IsFunctionAliased("GetMPICommWorld") ?
                         *(const MPI_Comm*)GetMPICommWorld(cctkGH) :
                         MPI_COMM_WORLD;

  assert(counts);
  assert(offsets);

  // local size on each process
  MPI_Gather(&count, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);

  // turns size into offsets
  count_total = 0;
  for(int i = 0 ; i < nProcs ; i++)
  {
    offsets[i] = count_total;
    count_total += counts[i];
  }
  assert(count_total == ntracers*components || myProc != 0);

  // data buffers to collect data in
  local_data = malloc(count * sizeof(*local_data));
  global_data = malloc(count_total * sizeof(*global_data));
  assert(local_data);
  assert(global_data);
  for(int i = 0 ; i < ntracers_on_this_proc ; i++)
    for(int c = 0 ; c < components ; c++)
      local_data[i*components+c] = (float)data[c][i];

  // get all data
  MPI_Gatherv(local_data, count, MPI_FLOAT, global_data,
              counts, offsets, MPI_FLOAT, 0, comm);

  // write to disk
  if(myProc == 0)
    WriteField(file_id, global_data, label, components, ntracers);

  // clean up
  free(global_data);
  free(local_data);
  free(counts);
  free(offsets);
}

