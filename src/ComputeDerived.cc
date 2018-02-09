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

extern "C" void TracerFlow_ComputeDerivedQuantities(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // skip any computation
  if(step_freq <= 0)
    return;

  int const gsiz = cctkGH->cctk_ash[0]*cctkGH->cctk_ash[1]*cctkGH->cctk_ash[2];
  CCTK_REAL const * velx = &vel[0*gsiz];
  CCTK_REAL const * vely = &vel[1*gsiz];
  CCTK_REAL const * velz = &vel[2*gsiz];

    #pragma omp parallel for schedule(static)
    for(int k = 0; k < cctk_lsh[2]; ++k) {
      for(int j = 0; j < cctk_lsh[1]; ++j) {
        for(int i = 0; i < cctk_lsh[0]; ++i) {
          const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          adv_velx[ijk] = alp[ijk] * velx[ijk] - betax[ijk];
          adv_vely[ijk] = alp[ijk] * vely[ijk] - betay[ijk];
          adv_velz[ijk] = alp[ijk] * velz[ijk] - betaz[ijk];

          // compute covariant three velocity
          CCTK_REAL const v_x = gxx[i]*velx[i] + gxy[i]*vely[i] + gxz[i]*velz[i];
          CCTK_REAL const v_y = gxy[i]*velx[i] + gyy[i]*vely[i] + gyz[i]*velz[i];
          CCTK_REAL const v_z = gxz[i]*velx[i] + gyz[i]*vely[i] + gzz[i]*velz[i];

          // compute approximation of E_inf, i.e. - u_t - 1
          eninf[ijk] = - w_lorentz[ijk]
                     * (v_x * betax[ijk]
                      + v_y * betay[ijk]
                      + v_z * betaz[ijk]
                      - alp[ijk])
                     - 1.0;
        }
      }
    }

}
