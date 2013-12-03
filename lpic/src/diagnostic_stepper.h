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

#ifndef DIAGNOSTIC_STEPPER_H
#define DIAGNOSTIC_STEPPER_H

#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>


typedef struct stepper_param {
  int Q;
  double t_start;                // in periods
  double t_stop;                 // in periods
  double t_step;                 // in periods
  int x_start;                   // in cells
  int x_stop;                    // in cells
  int x_step;                    // in cells
} stepper_param;


//////////////////////////////////////////////////////////////////////////////////////////


class diagnostic_stepper {
private:
  char     errname[filename_size];

public:
  int Q;

  int   t_start;
  int   t_stop;
  int   t_step;
  int   t_count;

  int x_start;
  int x_stop;
  int x_step;

  diagnostic_stepper( stepper_param &t, parameter &p );
};

//////////////////////////////////////////////////////////////////////////////////////////
#endif




