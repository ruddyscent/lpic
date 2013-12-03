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

#ifndef DIAGNOSTIC_ENERGY_H
#define DIAGNOSTIC_ENERGY_H

#include <diagnostic_stepper.h>
#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>

class input_energy {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );

public:
  stepper_param stepper;
  int           Q_restart;
  char          restart_file[filename_size];

  input_energy( parameter &p );
};

//////////////////////////////////////////////////////////////////////////////////////////

class energy {
private:
  char         errname[filename_size];
  readfile     rf;
  input_energy input;

public:
  diagnostic_stepper stepper;

  double flux;
  double field, field_0;
  double field_t, field_t_0;
  double field_l, field_l_0;
  double kinetic, kinetic_0;
  double total, total_0;
  char *name;
  ofstream file;

  energy              ( parameter &p, domain* grid );
  void get_energies   ( domain* grid );
  void write_energies ( double time );
  void average_reflex ( domain *grid );
};

//////////////////////////////////////////////////////////////////////////////////////////

#endif
