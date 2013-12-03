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

#ifndef DIAGNOSTIC_SPACETIME_H
#define DIAGNOSTIC_SPACETIME_H

#include <diagnostic_stepper.h>
#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>


class input_spacetime {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );

public:
  stepper_param stepper_de, stepper_di,
                stepper_jx, stepper_jy, stepper_jz,
                stepper_ex, stepper_ey, stepper_ez,
                stepper_bx, stepper_by, stepper_bz,
                stepper_edens;
  int           Q_restart;
  char          restart_file[filename_size];

  input_spacetime( parameter &p );
};

//////////////////////////////////////////////////////////////////////////////////////////

class spacetime {

private:
  char            errname[filename_size];
  readfile        rf;
  input_spacetime input;

public:
  int             output_period_de, output_period_di,
                  output_period_jx, output_period_jy, output_period_jz,
                  output_period_ex, output_period_ey, output_period_ez,
                  output_period_bx, output_period_by, output_period_bz,
                  output_period_edens;

  diagnostic_stepper stepper_de, stepper_di,
                     stepper_jx, stepper_jy, stepper_jz,
                     stepper_ex, stepper_ey, stepper_ez,
                     stepper_bx, stepper_by, stepper_bz,
                     stepper_edens;

  char *name_de, *name_di,
       *name_jx, *name_jy, *name_jz,
       *name_ex, *name_ey, *name_ez,
       *name_bx, *name_by, *name_bz,
       *name_edens;

  spacetime            ( parameter &p );
  void boundaries      ( float *x_start, float *x_stop, int *x_steps,
			 diagnostic_stepper *stepper, domain *grid );
  void write_de        ( domain *grid, int time_out_count, parameter &p );
  void write_di        ( domain *grid, int time_out_count, parameter &p );
  void write_jx        ( domain *grid, int time_out_count, parameter &p );
  void write_jy        ( domain *grid, int time_out_count, parameter &p );
  void write_jz        ( domain *grid, int time_out_count, parameter &p );
  void write_ex        ( domain *grid, int time_out_count, parameter &p );
  void write_ey        ( domain *grid, int time_out_count, parameter &p );
  void write_ez        ( domain *grid, int time_out_count, parameter &p );
  void write_bx        ( domain *grid, int time_out_count, parameter &p );
  void write_by        ( domain *grid, int time_out_count, parameter &p );
  void write_bz        ( domain *grid, int time_out_count, parameter &p );
  void write_edens     ( domain *grid, int time_out_count, parameter &p );
};

//////////////////////////////////////////////////////////////////////////////////////////

#endif



