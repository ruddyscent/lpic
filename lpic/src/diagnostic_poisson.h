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

#ifndef DIAGNOSTIC_POISSON_H
#define DIAGNOSTIC_POISSON_H

#include <diagnostic_stepper.h>
#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>
#include <readfile.h>

class input_poisson {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );

public:
  stepper_param stepper;
  int           n_domains;
  double        time_start, time_stop;
  int           Q_restart;
  char          restart_file[filename_size];

  input_poisson( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class poisson {

private:

  readfile      rf;
  input_poisson input;
  int           spp;
  int           domain_number;
  char          errname[filename_size];
  char          *output_path;

  void fft            ( double* data, int nn, int isign );
  double dif          ( double in );
  double smooth       ( double in );

public:
  diagnostic_stepper stepper;

  double *ex, *rhok, *phik;
  char *name;
  std::ofstream file;

       poisson ( parameter &p, domain *grid );
  void solve   ( domain* grid );
  void write   ( double time, domain* grid );
};

//////////////////////////////////////////////////////////////////////////////////////////

#endif
