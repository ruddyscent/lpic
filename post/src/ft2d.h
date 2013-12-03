/*
   This file is part of LPIC++, a particle-in-cell code for
   simulating the interaction of laser light with plasma.

   Copyright (C) 1994-1997 Roland Lichters

   LPIC++ is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#ifndef FT2_H
#define FT2_H

#include <common.h>
#include <utilities.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string.h>

class FFT2D {
 public:

// input
  int   periods_input_1, periods_input_2;
  int   steps_pp_input_1, steps_pp_input_2;
  int   screen_input_1, screen_input_2;

// internal
  int   steps_1, steps_2, steps_half_1, steps_half_2, steps_input_1, steps_input_2;
  float dt_1, dt_input_1, dt_2, dt_input_2;
  float df_1, df_2;
  float **local;
  float *data;
  float *frequency_1, *frequency_2;
  float **si;
  float **co;
  float **power;
  int   *nn;

  char errname[filename_size];

      FFT2D            ( int Q_kw, int periods_1, int steps_pp_1, int screen_1,
			           int periods_2, int steps_pp_2, int screen_2 );
  int   steps_ft       ( int steps_input );
  float window         ( float t, int screen_input, int periods_input );
  void  fftn           ( float* data, int *nn, int ndim, int isign );
  void  RealFt         ( float** input );
};

#endif

