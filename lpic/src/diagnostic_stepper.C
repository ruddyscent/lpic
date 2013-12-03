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

#include <diagnostic_stepper.h>


diagnostic_stepper::diagnostic_stepper( stepper_param &t, parameter &p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("stepper_diagnostic::Constructor",errname);

  Q = t.Q;

  t_start = (int) floor( t.t_start * p.spp + 0.5 );
  t_stop  = (int) floor( t.t_stop * p.spp + 0.5 );
  if (t.t_step<TINY) t_step = 1;
  else               t_step  = (int) floor( t.t_step * p.spp + 0.5 );
  t_count = 0;

  x_start = t.x_start;
  x_stop  = t.x_stop;
  x_step  = t.x_step;
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof
