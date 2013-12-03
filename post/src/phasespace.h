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

#ifndef PHASESPACE_H
#define PHASESPACE_H

#include <common.h>
#include <fstream.h>
#include <iomanip.h>
#include <string.h>
#include <parameter.h>
#include <utilities.h>
#include <error.h>
#include <math.h>


class input_phasespace {
private:
  char errname[filename_size];

public:

  int    Q_vx, Q_vy, Q_vz;
  int    Q_el, Q_ion;
  double period_start;
  double period_stop;
  double period_step;
  double xmax;
  double xoffset;

  readfile rf;
  void save( parameter &p );

  input_phasespace( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class phasespace {

 private:

  input_phasespace input;

  int    dim;
  int    Q_vx, Q_vy, Q_vz;
  int    Q_el, Q_ion;
  double period_start;
  double period_stop;
  double period_step;
  double xmax, xoffset;
  unsigned char **matrix_read, **matrix_write;
  int    **matrix_inter;
  char   *input_path;
  char   *output_path;
  char   errname[filename_size];

 public:
   phasespace( parameter &p );
  void           concat( void );
  int              read( char *unit, char *spec, double time );
  void            write( char *unit, char *spec, double time );
  void write_idl_header( char *spec, char *direct );
};

#endif

