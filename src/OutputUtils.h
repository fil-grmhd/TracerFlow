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

#ifndef TRACERFLOW_OUTPUT_H
#define TRACERFLOW_OUTPUT_H

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <hdf5.h>

#define STRLEN 1024
#define MYH5CHECK(val) assert(val >= 0)

#define length(X) (sizeof(X)/sizeof(*X))

typedef struct {
  double box_size;
  int flag_cooling;
  int flag_double_precision;
  int flag_feedback;
  int flag_ic_info;
  int flag_metals;
  int flag_sfr;
  int flag_stellar_age;
  double hubble_param;
  double mass_table[6];
  int num_files_per_snapshot;
  int num_part_this_file[6];
  unsigned int num_part_total[6];
  unsigned int num_part_total_high_word[6];
  double Omega0;
  double OmegaLambda;
  double redshift;
  double time;
} Header;

static void WriteVisItFile(CCTK_ARGUMENTS);

static void WriteHeader(hid_t file_id, Header * header);

static void WriteIDs(hid_t file_id, unsigned int const * IDs, int const siz);

static void WriteField(hid_t file_id, float const * var, char const * name, int ncomp, int const siz);

typedef void (*WriteVar)(cGH const *, int, int, CCTK_REAL * const restrict *
        restrict, char const *, hid_t);

static void WriteSingleVariableChunked(cGH const * cctkGH, int ntracers_on_this_proc,
                                       int components, CCTK_REAL * const restrict * restrict data,
                                       char const * label,hid_t file_id);

static void CollectAndWriteSingleVariable(const cGH *cctkGH, int ntracers_on_this_proc,
                                          int components, CCTK_REAL * const restrict * restrict data,
                                          const char *label,
                                          hid_t file_id);

#endif
