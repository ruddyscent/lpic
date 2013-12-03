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

#ifndef DIAGNOSTIC_VELOCITY_H
#define DIAGNOSTIC_VELOCITY_H

#include <diagnostic_stepper.h>
#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>


class input_velocity {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );
  char          species_name[filename_size];

public:
  stepper_param stepper;
  int           Q_restart;
  char          restart_file[filename_size];

  input_velocity( parameter &p, char *species_name_input );
};


//////////////////////////////////////////////////////////////////////////////////////////


class velocity {
private:

  readfile       rf;
  input_velocity input;
  char           errname[filename_size];
  int            dim;
  int            *x, *y, *z, *a;
  double         Beta, Gamma, vcut;
  int            species;
  char           species_name[filename_size];

public:
  diagnostic_stepper stepper;

  char *name;

  velocity ( parameter &p,
	     int species_input, char *species_name_input );
  void write_velocity( double time, parameter &p, domain *grid );
};


class el_velocity : public velocity {
public:

  el_velocity( parameter &p,
	       int species_input=0, char *species_name_input="el" );
};


class ion_velocity : public velocity {
public:

  ion_velocity( parameter &p,
		int species_input=1, char *species_name_input="ion" );
};

//////////////////////////////////////////////////////////////////////////////////////////

#endif
