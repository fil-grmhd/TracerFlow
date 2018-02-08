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



#ifndef INTERP_QUANT_HH
#define INTERP_QUANT_HH


// Interpolates advection velocity to given tracer positions
void TracerFlow_InterpAdvectionVelocity(cGH *cctkGH,
                                        int const num_tracers,
                                        double * t_x,
                                        double * t_y,
                                        double * t_z,
                                        double * t_velx,
                                        double * t_vely,
                                        double * t_velz,
                                        int * t_frozen);

#endif
