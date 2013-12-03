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

#ifndef FT_H
#define FT_H

#include <common.h>
#include <utilities.h>
#include <math.h>
#include <fstream.h>
#include <iomanip.h>
#include <string.h>

class FFT {
 public:

// input
  int   periods_input;
  int   steps_pp_input;
  int   periods_screen;
  int   padding;

// internal
  int   steps, steps_half, steps_input;
  double dt, dt_input;
  double df;
  double *local;
  double *data;
  double *frequency;
  double *si;
  double *co;
  double *phase;
  double *power;
  double *corr;

  char errname[filename_size];

                      FFT( int periods, int steps_pp, int screen );
  int            steps_ft( void );
  void                fft( double* data, int nn, int isign );
  int              RealFt( float* input );
  double           window( double );
  int         correlation( float *input, double mid, double width );
  double frequency_filter( double freq, double mid, double width );
  void                out( float* input, char *path, char *appendix );
};

#endif

