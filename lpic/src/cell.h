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

#ifndef CELL_H
#define CELL_H

#include <particle.h>

struct cell {

  int    number;                 // number of this cell
  int    domain;                 // domain number it belongs to
  struct cell *prev;             // pointer to the previous (left) cell
  struct cell *next;             // pointer to the next (right) cell

  double x;                      // position of the left cell boundary in wavelengths
  double charge;                 // charge density in units e*n_c
  double jx, jy, jz;             // current density in units e*n_c*c
  double ex, ey, ez;             // electric fields in units m*omega*c/e
  double bx, by, bz;             // magnetic fields in units m*omega/e
  double fp, fm, gp, gm;         //
  double dens[2];                // densities for each species in units n_c

  int             np[2];         // # of electrons [0] and ions [1]
  int             npart;         // # particles

  struct particle *first;        // pointer to the first particle
  struct particle *last;         // pointer to the last particle
  struct particle *insert;       // pointer to particle, in front of which new particles
                                 // have to be inserted
};

#endif


