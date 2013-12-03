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

#ifndef MAIN_H
#define MAIN_H

#include <common.h>
#include <iostream.h>
#include <error.h>
#include <parameter.h>
#include <box.h>
#include <diagnostic.h>
#include <uhr.h>
#include <propagate.h>

int main(int argc, char **argv);

int main_exit( parameter &p, box &sim )
{
  char                 errname[filename_size];
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("main_exit",errname);

#ifdef LPIC_PARALLEL
  sim.talk.end_task();
#endif

  bob.message("done");

  return(0);
}

#endif
