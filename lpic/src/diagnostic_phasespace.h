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

#ifndef DIAGNOSTIC_PHASESPACE_H
#define DIAGNOSTIC_PHASESPACE_H

#include <diagnostic_stepper.h>
#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>
#include <readfile.h>


class input_phasespace {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );
  char          species_name[filename_size];

public:
  stepper_param stepper;
  int           n_domains;
  int           cells, cells_per_wl;
  int           Q_restart;
  char          restart_file[filename_size];

  input_phasespace( parameter &p, char *species_name_input );
};


//////////////////////////////////////////////////////////////////////////////////////////


class phasespace {
private:
  char             errname[filename_size];
  readfile         rf;
  input_phasespace input;
  int              species;
  char             species_name[filename_size];

public:
  diagnostic_stepper stepper;

  int    dim;
  int    box_cells;
  double box_length;
  double Beta, Gamma;
  double vcut;
  unsigned char **x, **y, **z;
  double vx, vy, vz;
  int bx, bvx, bvy, bvz;
  char *name;
  FILE *file_x, *file_y, *file_z;

  phasespace           ( parameter &p,
			 int species_input, char *species_name_input );
  void write_phasespace( double time, parameter &p, domain *grid );
};


class el_phasespace : public phasespace {
public:

  el_phasespace( parameter &p,
		 int species_input=0, char *species_name_input="el" );
};


class ion_phasespace : public phasespace {
public:

  ion_phasespace( parameter &p,
		  int species_input=1, char *species_name_input="ion" );
};

//////////////////////////////////////////////////////////////////////////////////////////

#endif
