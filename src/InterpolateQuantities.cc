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
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <util_Table.h>

#include "InterpolateQuantities.hh"

#define length(X) (sizeof(X)/sizeof(*X))

void TracerFlow_InterpAdvectionVelocity(cGH *cctkGH,
                                        int const num_tracers,
                                        double * t_x,
                                        double * t_y,
                                        double * t_z,
                                        double * t_velx,
                                        double * t_vely,
                                        double * t_velz,
                                        int * t_frozen)
{
  DECLARE_CCTK_PARAMETERS;

  // CHECK: move this? add frozen check?
  for(int i = 0; i < num_tracers; ++i) {
    if(fabs(t_x[i]) > coord_bound ||
       fabs(t_y[i]) > coord_bound ||
       fabs(t_z[i]) > coord_bound) {

      if(verbose > 0)
        CCTK_VInfo(CCTK_THORNSTRING, "InterpAdvectionVelocity: Tracer %i moved out of bounds: (%e,%e,%e)",i,t_x[i],t_y[i],t_z[i]);

      t_velx[i] = 0;
      t_vely[i] = 0;
      t_velz[i] = 0;

      t_frozen[i] = 1
    }
  }

  // get interpolator
  int const interp_handle = CCTK_InterpHandle(interpolator);
  if(interp_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't find interpolator \"%s\"!",interpolator);

  // get options into a table and coordinate handle
  int const options_handle =
    Util_TableCreateFromString(interpolator_options);
  assert (options_handle >= 0);
  if(options_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't create interp option table with \"%s\"!",interpolator_options);

  int const coords_handle = CCTK_CoordSystemHandle(coordinate_system);
  if(coords_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't get handle for coordinate system \"%s\"!",coordinate_system);

  void const *const interp_coords[] = {
      reinterpret_cast<void *>(t_x),
      reinterpret_cast<void *>(t_y),
      reinterpret_cast<void *>(t_z)
  };

  CCTK_INT const input_array_indices[] = {
      CCTK_VarIndex("TracerFlow::adv_velx"),
      CCTK_VarIndex("TracerFlow::adv_vely"),
      CCTK_VarIndex("TracerFlow::adv_velz"),
  };
  int const ninputs = length(input_array_indices);

  CCTK_INT const output_array_types[] = {
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL
  };
  assert(ninputs == length(output_array_types));

  void * output_arrays[] = {
      reinterpret_cast<void *>(t_velx),
      reinterpret_cast<void *>(t_vely),
      reinterpret_cast<void *>(t_velz)
  };
  assert(ninputs == length(output_arrays));

  // Actual interpolation
  int const ierr =
    CCTK_InterpGridArrays (cctkGH, 3,
                           interp_handle, options_handle, coords_handle,
                           num_tracers, CCTK_VARIABLE_REAL, interp_coords,
                           ninputs, input_array_indices,
                           ninputs, output_array_types, output_arrays);
  assert (ierr == 0);

  Util_TableDestroy (options_handle);

  return;
}

extern "C" void TracerFlow_InterpAllQuantities(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /*
    Check if tracers are about to be written to file in this iteration.
    Interpolate all important quantities to the current position.
  */

  // if one of these is true, just skip the interpolation (not needed for the evolution)
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

  // get interpolator
  int const interp_handle = CCTK_InterpHandle(interpolator);
  if(interp_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't find interpolator \"%s\"!",interpolator);

  // get options into a table and coordinate handle
  int const options_handle =
    Util_TableCreateFromString(interpolator_options);
  assert (options_handle >= 0);
  if(options_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't create interp option table with \"%s\"!",interpolator_options);

  int const coords_handle = CCTK_CoordSystemHandle(coordinate_system);
  if(coords_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Couldn't get handle for coordinate system \"%s\"!",coordinate_system);

  void const *const interp_coords[] = {
      reinterpret_cast<void *>(tracer_x),
      reinterpret_cast<void *>(tracer_y),
      reinterpret_cast<void *>(tracer_z)
  };

  // CHECK: not interpolating the adv velocity here again
  //        could be slightly off due to RK4 substep
  CCTK_INT const input_array_indices[] = {
      CCTK_VarIndex("HydroBase::rho"),
      CCTK_VarIndex("HydroBase::temperature"),
      CCTK_VarIndex("HydroBase::Y_e"),
      CCTK_VarIndex("HydroBase::w_lorentz"),
      CCTK_VarIndex("TracerFlow::eninf")
  };
  int const ninputs = length(input_array_indices);

  CCTK_INT const output_array_types[] = {
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL
  };
  assert(ninputs == length(output_array_types));

  void * output_arrays[] = {
      reinterpret_cast<void *>(tracer_rho),
      reinterpret_cast<void *>(tracer_temperature),
      reinterpret_cast<void *>(tracer_ye),
      reinterpret_cast<void *>(tracer_wlorentz),
      reinterpret_cast<void *>(tracer_einf)
  };
  assert(ninputs == length(output_arrays));

  // Actual interpolation
  int const ierr =
    CCTK_InterpGridArrays (cctkGH, 3,
                           interp_handle, options_handle, coords_handle,
                           num_tracers, CCTK_VARIABLE_REAL, interp_coords,
                           ninputs, input_array_indices,
                           ninputs, output_array_types, output_arrays);
  assert (ierr == 0);

  Util_TableDestroy (options_handle);

  return;
}
