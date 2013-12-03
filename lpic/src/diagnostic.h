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

#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include <common.h>
#include <error.h>
#include <parameter.h>
#include <domain.h>
#include <matrix.h>
#include <math.h>
#include <readfile.h>

#include <diagnostic_stepper.h>
#include <diagnostic_trace.h>
#include <diagnostic_spacetime.h>
#include <diagnostic_energy.h>
#include <diagnostic_reflex.h>
#include <diagnostic_flux.h>
#include <diagnostic_snapshot.h>
#include <diagnostic_velocity.h>
#include <diagnostic_phasespace.h>
#include <diagnostic_poisson.h>

class input_diagnostic {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );

public:
  int Q_restart;
  char restart_file[filename_size];
  int Q_restart_save;
  char restart_file_save[filename_size];

  input_diagnostic( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class diagnostic {

private:

  readfile         rf;
  input_diagnostic input;

  int              time_steps, time_out;

  int              domain_number;
  char             errname[filename_size];
  char             output_path[filename_size];

public:

  diagnostic ( parameter &p, domain* grid );
  void             out( double time, domain* grid, parameter &p );
  void           count( void );
  int     write_window( int time_steps, diagnostic_stepper *stepper );

  int     public_time_steps;
  int     time_out_count;

  poisson        poi;
  snapshot       sna;
  el_velocity    vel_el;
  ion_velocity   vel_ion;
  flux           flu;
  reflex         ref;
  spacetime      spa;
  energy         ene;
  trace          tra;
  el_phasespace  pha_el;
  ion_phasespace pha_ion;

};

//////////////////////////////////////////////////////////////////////////////////////////

#endif



