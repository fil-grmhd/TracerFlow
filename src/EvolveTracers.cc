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
#include "util_Table.h"

#include <assert.h>
#include <memory>

#include "InterpolateQuantities.hh"

/*
  RK4 steps for tracer advection evolution
*/

#define DEBUG_A (-1.31)
void debug_analytic_test(int num_tracers,double *x_posn,double *y_posn,double *z_posn,  double *tracer_velx,double *tracer_vely,double *tracer_velz) {
  for(int i = 0; i<num_tracers; i++) {
    double xc = x_posn[i];
    double yc = y_posn[i];
    double zc = z_posn[i];
    tracer_velx[i] = DEBUG_A*xc;
    tracer_vely[i] = DEBUG_A*yc;
    tracer_velz[i] = DEBUG_A*zc;
  }
}

extern "C" void TracerFlow_Advect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // skip any evolution
  if(step_freq <= 0) return;

  // figure out how many tracers each process has
  int group = CCTK_GroupIndex("TracerFlow::tracer_evol");
  cGroupDynamicData data;
  int retval = CCTK_GroupDynamicData(cctkGH, group, &data);
  (void)retval;
  // (processor) local size of the array, i.e. local number of tracers
  int num_tracers = data.lsh[0];

  if(cctk_iteration % step_freq == 0 && cctk_iteration >= start_iteration) {
    if(verbose > 1)
      CCTK_VInfo(CCTK_THORNSTRING, "Advect: Got a dt of %e.",CCTK_DELTA_TIME);

    // The factor of 2 comes from RK4 going from t0 to t0+dt/2 to t0+dt.
    double dt = CCTK_DELTA_TIME * step_freq * 2.0;

    /* We must solve the following 3 ODEs to get the next particle position
       dx/dt = vx
       dy/dt = vy
       dz/dt = vz

       For one variable:
       y' = f(y_n,t_n)
    */

    if(*RK4_counter == 1) {
      /*
        RK4 step 1:
        k_1 = dt f(y_n,t_n) = dt f(y_n) <- no explicit time dependence in f
      */

      // Interpolate advection velocity to the right locations
      TracerFlow_InterpAdvectionVelocity(cctkGH,num_tracers,tracer_x,
                                                            tracer_y,
                                                            tracer_z,
                                                            tracer_velx,
                                                            tracer_vely,
                                                            tracer_velz,
                                                            tracer_frozen);


      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_analytic_test(num_tracers,tracer_x,tracer_y,tracer_z, tracer_velx,tracer_vely,tracer_velz);
      /***************************************************************/

      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
        tracer_x_k1[i] = dt*tracer_velx[i];
        tracer_y_k1[i] = dt*tracer_vely[i];
        tracer_z_k1[i] = dt*tracer_velz[i];
      }

      *RK4_counter = 2;

      return;

    } else if(*RK4_counter == 2) {
      // make some temporary arrays for shifted positions
//      auto shifted_tracer_x = std::make_unique<double[]>(num_tracers);
//      auto shifted_tracer_y = std::make_unique<double[]>(num_tracers);
//      auto shifted_tracer_z = std::make_unique<double[]>(num_tracers);

      /*
        RK4 step 2:
        k_2 = dt f(y_n + 0.5*k_1,t_n + 0.5 dt) = dt f(y_n + 0.5*k_1) <- no explicit time dependence in f
      */

      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
        shifted_tracer_x[i] = tracer_x[i] + 0.5*tracer_x_k1[i];
        shifted_tracer_y[i] = tracer_y[i] + 0.5*tracer_y_k1[i];
        shifted_tracer_z[i] = tracer_z[i] + 0.5*tracer_z_k1[i];
      }

      TracerFlow_InterpAdvectionVelocity(cctkGH,num_tracers,shifted_tracer_x,
                                                            shifted_tracer_y,
                                                            shifted_tracer_z,
                                                            tracer_velx,
                                                            tracer_vely,
                                                            tracer_velz,
                                                            tracer_frozen);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_analytic_test(num_tracers,shifted_tracer_x,shifted_tracer_y,shifted_tracer_z, tracer_velx,tracer_vely,tracer_velz);
      /***************************************************************/

      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
        tracer_x_k2[i] = dt*tracer_velx[i];
        tracer_y_k2[i] = dt*tracer_vely[i];
        tracer_z_k2[i] = dt*tracer_velz[i];
      }

      /*
        RK4 step 3:
        k_3 = dt f(y_n + 0.5*k_2,t_n + 0.5 dt) = dt f(y_n + 0.5 k_2) <- no explicit time dependence in f
      */

      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
        shifted_tracer_x[i] = tracer_x[i] + 0.5*tracer_x_k2[i];
        shifted_tracer_y[i] = tracer_y[i] + 0.5*tracer_y_k2[i];
        shifted_tracer_z[i] = tracer_z[i] + 0.5*tracer_z_k2[i];
      }

      TracerFlow_InterpAdvectionVelocity(cctkGH,num_tracers,shifted_tracer_x,
                                                            shifted_tracer_y,
                                                            shifted_tracer_z,
                                                            tracer_velx,
                                                            tracer_vely,
                                                            tracer_velz,
                                                            tracer_frozen);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_analytic_test(num_tracers,shifted_tracer_x,shifted_tracer_y,shifted_tracer_z, tracer_velx,tracer_vely,tracer_velz);
      /***************************************************************/

      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
        tracer_x_k3[i] = dt*tracer_velx[i];
        tracer_y_k3[i] = dt*tracer_vely[i];
        tracer_z_k3[i] = dt*tracer_velz[i];
      }

      *RK4_counter = 4;
      return;
    } else if(*RK4_counter == 4) {
      // make some temporary arrays for shifted positions
//      auto shifted_tracer_x = std::make_unique<double[]>(num_tracers);
//      auto shifted_tracer_y = std::make_unique<double[]>(num_tracers);
//      auto shifted_tracer_z = std::make_unique<double[]>(num_tracers);

      /*
        RK4 step 4:
        k_4 = dt f(y_n + k_3,t_n + dt) = dt f(y_n + k_3) <- no explicit time dependence in f
      */
      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
        shifted_tracer_x[i] = tracer_x[i] + tracer_x_k3[i];
        shifted_tracer_y[i] = tracer_y[i] + tracer_y_k3[i];
        shifted_tracer_z[i] = tracer_z[i] + tracer_z_k3[i];
      }

      TracerFlow_InterpAdvectionVelocity(cctkGH,num_tracers,shifted_tracer_x,
                                                            shifted_tracer_y,
                                                            shifted_tracer_z,
                                                            tracer_velx,
                                                            tracer_vely,
                                                            tracer_velz,
                                                            tracer_frozen);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_analytic_test(num_tracers,shifted_tracer_x,shifted_tracer_y,shifted_tracer_z, tracer_velx,tracer_vely,tracer_velz);
      /***************************************************************/

      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
        tracer_x_k4[i] = dt*tracer_velx[i];
        tracer_y_k4[i] = dt*tracer_vely[i];
        tracer_z_k4[i] = dt*tracer_velz[i];
      }

      /******************************/
      /* DEBUG MODE: compute errors */
      if(debug) {
        double sumx,sumy,sumz,sumxe,sumye,sumze;
        sumx=sumy=sumz=sumxe=sumye=sumze=0;
        for(int i = 0; i<num_tracers; i++) {
          double orig_x = tracer_x[i];
          double orig_y = tracer_y[i];
          double orig_z = tracer_z[i];

          //dx/dt = DEBUG_A*x -> x = x_0 e^{DEBUG_A t} = 1 + (DEBUG_A) t + (DEBUG_A t)^2/2 + (DEBUG_A t)^3/6 + (DEBUG_A t)^4/24 + O(t^5). RK4 should get all but the O(t^5) term exact.
          double a_deltat=DEBUG_A*dt;
          double x_exact_to_4th_order = orig_x * (1. + a_deltat + pow(a_deltat,2)/2. + pow(a_deltat,3)/6. + pow(a_deltat,4)/24.);
          double y_exact_to_4th_order = orig_y * (1. + a_deltat + pow(a_deltat,2)/2. + pow(a_deltat,3)/6. + pow(a_deltat,4)/24.);
          double z_exact_to_4th_order = orig_z * (1. + a_deltat + pow(a_deltat,2)/2. + pow(a_deltat,3)/6. + pow(a_deltat,4)/24.);

          double x_num = (orig_x + (1.0/6.0)*(tracer_x_k1[i] + 2.0*tracer_x_k2[i] + 2.0*tracer_x_k3[i] + tracer_x_k4[i]));
          double y_num = (orig_y + (1.0/6.0)*(tracer_y_k1[i] + 2.0*tracer_y_k2[i] + 2.0*tracer_y_k3[i] + tracer_y_k4[i]));
          double z_num = (orig_z + (1.0/6.0)*(tracer_z_k1[i] + 2.0*tracer_z_k2[i] + 2.0*tracer_z_k3[i] + tracer_z_k4[i]));
          sumx += x_exact_to_4th_order - x_num;
          sumy += y_exact_to_4th_order - y_num;
          sumz += z_exact_to_4th_order - z_num;

          sumxe+= fabs(x_exact_to_4th_order);
          sumye+= fabs(y_exact_to_4th_order);
          sumze+= fabs(z_exact_to_4th_order);
        }
        printf("PARTICLE_TRACER DEBUG OUTPUT: Error in x=%e. Error in y=%e. Error in z=%e.\n",sumx/sumxe,sumy/sumye,sumz/sumze);
        printf(" (Errors within ~8x machine precision -> particle tracer RK4 routine is coded correctly.)\n");
      }
      /******************************/

      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
      /*
        Update particle positions now. For RK4:
        x_{n+1} = x_n + 1/6 ( k_1 + 2 k_2 + 2 k_3 + k_4 )
      */
        tracer_x[i] += (1.0/6.0)*(tracer_x_k1[i] + 2.0*tracer_x_k2[i] + 2.0*tracer_x_k3[i] + tracer_x_k4[i]);
        tracer_y[i] += (1.0/6.0)*(tracer_y_k1[i] + 2.0*tracer_y_k2[i] + 2.0*tracer_y_k3[i] + tracer_y_k4[i]);
        tracer_z[i] += (1.0/6.0)*(tracer_z_k1[i] + 2.0*tracer_z_k2[i] + 2.0*tracer_z_k3[i] + tracer_z_k4[i]);
      }

      /*
        Ready for RK4 step 1 again!
        k_1 = dt f(y_n,t_n)
      */
      TracerFlow_InterpAdvectionVelocity(cctkGH,num_tracers,tracer_x,
                                                            tracer_y,
                                                            tracer_z,
                                                            tracer_velx,
                                                            tracer_vely,
                                                            tracer_velz,
                                                            tracer_frozen);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_analytic_test(num_tracers,tracer_x,tracer_y,tracer_z, tracer_velx,tracer_vely,tracer_velz);
      /***************************************************************/

      #pragma omp simd
      for(int i = 0; i<num_tracers; i++) {
        tracer_x_k1[i] = dt*tracer_velx[i];
        tracer_y_k1[i] = dt*tracer_vely[i];
        tracer_z_k1[i] = dt*tracer_velz[i];

        tracer_velx_out[i] = tracer_velx[i];
        tracer_vely_out[i] = tracer_vely[i];
        tracer_velz_out[i] = tracer_velz[i];
      }

      *RK4_counter = 2;
    }
  }
}
