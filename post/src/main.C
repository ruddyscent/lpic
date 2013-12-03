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

//////////////////////////////////////////////////////////////////////////////////////////
//
// postprocessor for lpic++
//
//////////////////////////////////////////////////////////////////////////////////////////

#include <main.h>

//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  parameter            p(argc, argv);

  char errname[filename_size];
  sprintf( errname, "%s/error", p.output_path );
  static error_handler bob("main",errname);

  //  trace                tr(p);
  //  tr.transform(p);           // calculate Fourier transforms from trace data

  spacetime            sp(p);
  sp.select();               // convert spacetime files into plot format


  phasespace           ph(p);
  ph.concat();               // concatenate phasespace files of different domains

  bob.message("done");

  exit(0);
}

//////////////////////////////////////////////////////////////////////////////////////////
//EOF


