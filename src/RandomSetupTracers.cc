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

#include <random>

extern "C" void TracerFlow_Random_SetupTracers(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_VInfo(CCTK_THORNSTRING, "Enter Random_SetupTracers");

    int group = CCTK_GroupIndex("TracerFlow::tracer_evol");
    cGroupDynamicData data;
    (void) CCTK_GroupDynamicData (cctkGH, group, &data);
    int num_tracers = data.lsh[0];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist_x(xmin,xmax);
    std::uniform_real_distribution<> dist_y(ymin,ymax);

    for(int i = 0; i < num_tracers; ++i) {
      tracer_x[i] = dist_x(gen);
      tracer_y[i] = dist_y(gen);
      tracer_z[i] = 0;
    };

    CCTK_VInfo(CCTK_THORNSTRING, "Exit Random_SetupTracers");
}
