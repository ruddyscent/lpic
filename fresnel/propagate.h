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

#ifndef PROPAGATE_H
#define PROPAGATE_H

#include <common.h>
#include <error.h>
#include <parameter.h>
#include <uhr.h>

class propagate {
public:
                  propagate( parameter &p );
    void               loop( parameter &p );
    double     reflectivity( double w, double wp, double angle, int polarization );
    double            phase( double w, double wp, double angle, int polarization );
    double              dif( double x );
    void                out( parameter &p );
    void            out_avg( parameter &p );
    void            out_int( parameter &p );

private:

    int        t_avg;
    double     t, w, wmax;      // in periods and laser frequency, respectively
    double     dt, dw;
    double     r, rs, rp;       // plane wave reflectivity
    double     phi, phis, phip; // reflected plane wave's phase
    double     reflex;          // as a function of time for the rectangular pulse
    double     reflex_avg;
    double     integrand;
    double     T;               // pulse duration in periods

    uhr        zeit;

    char errname[filename_size];
    ofstream reflex_file, reflex_file_avg, reflex_file_int;
};
#endif

