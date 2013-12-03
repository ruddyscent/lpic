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

#ifndef DIAGNOSTIC_REFLEX_H
#define DIAGNOSTIC_REFLEX_H

#include <diagnostic_stepper.h>
#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>


class input_reflex {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );

public:
  stepper_param stepper;
  int           Q_restart;
  char          restart_file[filename_size];

  input_reflex( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class reflex {
private:
  char         errname[filename_size];
  readfile     rf;
  input_reflex input;

public:
  diagnostic_stepper stepper;

  double   *buf;
  char     *name;
  std::ofstream file;

  reflex              ( parameter &p );
  void average_reflex ( domain* grid );
  void write_reflex   ( double time );
};

//////////////////////////////////////////////////////////////////////////////////////////

#endif
