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

#ifndef DIAGNOSTIC_TRACE_H
#define DIAGNOSTIC_TRACE_H

#include <diagnostic_stepper.h>
#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>
#include <readfile.h>

class input_trace {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );

public:
  stepper_param stepper;
  int           domain_number;
  int           traces;
  int           *tracepos;
  int           Q_restart;
  char          restart_file[filename_size];
  int           Q_restart_save;
  char          restart_file_save[filename_size];
  char          *path;

  input_trace( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class trace {
private:
  char        errname[filename_size];
  readfile    rf;
  input_trace input;

public:
  diagnostic_stepper stepper;

  int         traces;
  int         *cell_number;
  struct cell **cell;
  float       **fp, **fm, **gp, **gm, **ex;
  float       **dens_e, **dens_i;
  float       **jx, **jy, **jz;
  char        *name;
  FILE        *file;

  trace             ( parameter &p );
  void store_traces ( domain* grid );
  void write_traces ( double time, parameter &p );
  void restart_save ( void );
};

//////////////////////////////////////////////////////////////////////////////////////////

#endif
