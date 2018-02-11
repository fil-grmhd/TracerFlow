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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <vector>
#include <sstream>

std::vector<std::string> extra_var_names;
std::vector<std::vector<CCTK_REAL>> tracer_extras;

extern "C" void TracerFlow_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  // initialize RK4 counter
  *RK4_counter = 1;

  // figure out how many tracers each process has
  int group = CCTK_GroupIndex("TracerFlow::tracer_evol");
  cGroupDynamicData data;
  int retval = CCTK_GroupDynamicData(cctkGH, group, &data);
  // (processor) local size of the array, i.e. local number of tracers
  int num_tracers = data.lsh[0];

  // at the beginning none of the tracers are frozen (at the boundary)
  for(int i = 0; i < num_tracers; ++i) {
    tracer_frozen[i] = 0;
  }

  // parse extra vars
  if(!CCTK_Equals(extra_vars, "")) {
    // read vars line by line
    std::istringstream extras_stream(extra_vars);
    std::string line;
    while(std::getline(extras_stream, line)) {
      if(line != "" && line.at(0) != '#')
        extra_var_names.push_back(line);
    }

    // now check variable indices and create storage for them
    for(auto it = extra_var_names.begin(); it < extra_var_names.end(); ++it) {
      CCTK_INT var_idx = CCTK_VarIndex(it->c_str());
      if(var_idx < 0) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Couldn't find extra variable '%s'",it->c_str());
      }
    }
    tracer_extras = std::vector<std::vector<CCTK_REAL>>(extra_var_names.size(), std::vector<CCTK_REAL>(num_tracers));
  }
}
