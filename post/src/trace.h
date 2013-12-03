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

#ifndef TRACE_H
#define TRACE_H

#include <common.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <parameter.h>
#include <utilities.h>
#include <ft.h>
#include <math.h>


class input_trace {
private:
  char errname[filename_size];

public:

  // input from input.post:

  int region;
  int period_start;
  int period_stop;
  int screen;

  int Q_fp, Q_fm, Q_gp, Q_gm;
  int Q_ex, Q_ey, Q_ez, Q_by, Q_bz;
  int Q_Sr, Q_Pr, Q_Si, Q_Pi;
  int Q_de, Q_di, Q_jx, Q_jy, Q_jz;

  // input from arbitrary trace file and derived values
  int traces;                        // to be read from a trace file
  int steps_pp;                      // to be read from a trace file
  int periods;

  readfile rf;
  void save( parameter &p );

  input_trace( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class trace {

 private:
  input_trace input;
  FFT ft;                      // class 'Fourier Transforms'

  int region;
  int period;
  int period_start;
  int period_stop;
  int periods;
  int traces;
  float *position;            // position of traces
  int steps_pp;

  int Q;
  float **fp, **fm, **gp, **gm, **ex, **ey, **ez, **by, **bz;
  float **de, **di;
  float **jx, **jy, **jz;
  float *vector_read;

  double **power_spectrum;
  char *power_name, *trace_name;
  std::ofstream powerfile, tracefile;
  FILE* read_file;
  char *read_name;
  char *path;
  char errname[filename_size];

 public:
       trace          ( parameter &p );
  void read_traces    ( void );
  void transform      ( parameter &p );
  void write_transform( parameter &p, char* appendix );
  void write_traces   ( parameter &p, float** input, char* appendix );

};


#endif
