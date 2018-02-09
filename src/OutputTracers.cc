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
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <util_String.h>

#include "Extras.hh"

extern "C" {
#include "OutputUtils.h"
}

typedef struct {
  char const * name;
  CCTK_REAL * data;
} Field;

void TracerFlow_OutputTracers(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  if(outTracers_every == 0) {
    if(cctk_iteration != outTracers_iteration) {
      return;
    }
  }
  else {
    if(cctk_iteration % outTracers_every != 0 && cctk_iteration != outTracers_iteration) {
      return;
    }
  }

  const int nProcs = CCTK_nProcs(cctkGH);
  const int myProc = CCTK_MyProc(cctkGH);
  int ierr;

  int ntracers_on_this_proc;
  int group;
  cGroupDynamicData groupdata;

  // collect information about tracers get data
  group = CCTK_GroupIndex("TracerFlow::tracervars");
  ierr = CCTK_GroupDynamicData(cctkGH, group, &groupdata);
  if(ierr != 0)
  {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to obtain information about tracervars group: %d",
                ierr);
  }
  ntracers_on_this_proc = groupdata.lsh[0];
  int const ntracers_to_write = one_file_per_timestep ? ntracers :
    ntracers_on_this_proc;

  // open file for the output
  char * fn = NULL;
  hid_t file_id = -1;
  if(one_file_per_timestep)
  {
    if(myProc == 0)
    {
      Util_asprintf(&fn, "%s/%s.%d.hdf5", out_dir, tracersfile,
          cctk_iteration);
      assert(fn);
      file_id = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        MYH5CHECK(file_id);
    }
  }
  else
  {
    Util_asprintf(&fn, "%s/%s.%d.%d.hdf5", out_dir, tracersfile,
        cctk_iteration, myProc);
    assert(fn);
      file_id = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        MYH5CHECK(file_id);
  }

  // Construct and write header
  Header header;
  memset(&header, 0, sizeof(Header));
  header.num_files_per_snapshot = one_file_per_timestep ? 1 : nProcs;
  header.num_part_this_file[0] = ntracers_to_write;
  header.num_part_total[0] = ntracers;
  header.time = cctk_time;
  if(file_id >= 0)
    WriteHeader(file_id, &header);

  // Construct the index
  std::vector<unsigned int> ID(ntracers_to_write);
  if(one_file_per_timestep)
  {
    for(int i = 0; i < ntracers; ++i)
     ID[i] = i;
  }
  else
  {
    for(int i = 0; i < ntracers_to_write; ++i)
      ID[i] = i + groupdata.lbnd[0];
  }
  if(file_id >= 0)
    WriteIDs(file_id, ID.data(), ntracers_to_write);

  // All of the scalar fields to output
  // CHECK: any leakage quantities to output?
  Field fields_base[] = {
    {"Frozen",                  tracer_frozen},
    {"LorentzFactor",           tracer_wlorentz},
    {"KineticEnergyAtInfinity", tracer_eninf},
    {"Density",                 tracer_rho},
    {"Temperature",             tracer_temp},
    {"Ye",                      tracer_ye}
  };
  int const n_fields_base = length(fields_base);

  std::vector<Field> fields_extra;
  for(int i = 0; i < extra_var_names.size(); ++i) {
    Field extra_field = {extra_var_names[i].c_str(), tracer_extras[i].data()};
    fields_extra.push_back(extra_field);
  }

  // Choose which function to use for writing the data
  WriteVar writer;
  if(one_file_per_timestep)
    writer = CollectAndWriteSingleVariable;
  else
    writer = WriteSingleVariableChunked;

  // Output vector vars
  CCTK_REAL * Pos[3] = {tracer_x, tracer_y, tracer_z};
  writer(cctkGH, ntracers_on_this_proc, 3, Pos, "Coordinates", file_id);
  CCTK_REAL * Vel[3] = {tracer_velx_out, tracer_vely_out, tracer_velz_out};
  writer(cctkGH, ntracers_on_this_proc, 3, Vel, "Velocities", file_id);

  // Output scalar vars
  for(int f = 0; f < n_fields_base; ++f)
  {
    writer(cctkGH, ntracers_on_this_proc, 1, &fields_base[f].data,
        fields_base[f].name, file_id);
  }

  // Output extra vars
  for(auto it = fields_extra.begin(); it < fields_extra.end(); ++it) {
    writer(cctkGH, ntracers_on_this_proc, 1, &(it->data),
        it->name, file_id);
  }

  // Write VisIt file
  if(myProc == 0 && write_VisIt_file)
    WriteVisItFile(CCTK_PASS_CTOC);

  // Cleanup
  if(file_id >= 0)
  {
    H5Fclose(file_id);
  }

  free(fn);
}
