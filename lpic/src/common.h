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
// physical constants in SI units
// mathematical constants
// common definitions
// common procedures
//
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef COMMON_H
#define COMMON_H

#include <config.h>
#include <debug.h>

#ifdef LPIC_PARALLEL
//#define SLOW 1
#undef SLOW
#endif

#define C    2.9979246e+8    // m/s    velocity of light in vacuum
#define E    1.6021773e-19   // C      electron charge
#define M    9.1093897e-31   // kg     electron mass
#define EPS  8.8541878e-12   // C/Vm   dielectric permability
#define EV   1.6021773e-19   // J      electron volt
#define HB   1.0545727e-34   // Js     Planck's h/2pi
#define AB   5.2917726e-11   // m      Bor radius
#define EN   4.3597483e-18   // J      atomic energy unit = 2 * 13.6 eV
#define TI   2.4188843e-17   // s      atomic time unit = HB / EN
#define PI   M_PI

#define TINY 1e-10
#define MASK 10              // -> propagate::mask_current()

#define filename_size 100

inline double sqr(double x) { return (x*x); }

#endif





