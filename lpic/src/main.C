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

#include <main.h>


//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    // initialize classes ////////////////////////////////////////////////////////////////

    parameter p(argc,argv);                         // read parameters


    char                 errname[filename_size];
    sprintf( errname, "%s/error-%d", p.path, p.domain_number );
    static error_handler bob("main",errname);       // error handler for main.C


    box        sim(p);                              // init domain, cells, particles
                                                    // spawn task for the following domain


    pulse      laser_front(p,"&pulse_front");       // init laser pulses
    pulse      laser_rear(p,"&pulse_rear");

    diagnostic diag(p,&(sim.grid));                 // init diagnostics
    propagate  prop(p,sim.grid);                    // init propagator

    // main loop /////////////////////////////////////////////////////////////////////////

    prop.loop(p,sim,laser_front,laser_rear,diag);

    // exit //////////////////////////////////////////////////////////////////////////////

    return main_exit(p,sim);                        // stop parallel task

}

//////////////////////////////////////////////////////////////////////////////////////////
//EOF


