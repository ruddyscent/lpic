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

#ifndef PARTICLE_H
#define PARTICLE_H

#include <cell.h>

struct particle {

  int             number;        // number of this particle
  int             species;       // particle species, 0=electron, 1=ion
  struct cell     *cell;         // pointer to the cell this particle belongs to
  struct particle *prev;         // pointer to the previous particle in this cell
  struct particle *next;         // pointer to the next particle in this cell

  int    fix;                    // fixed species? 0->no, 1->yes
  double z;                      // charge of the micro particle in units of e
  double m;                      // mass of the micro particle in units of m_e
  double zm;                     // specific charge, z/m
  double x, dx;                  // position and shift within one timestep
  double igamma;                 // inverse gamma factor
  double ux, uy, uz;             // gamma * velocity
  double n;                      // particle density in units of n_c
  double zn;                     // contribution of the particle to the charge density
                                 // in units of n_c ( = z * n )
};

#endif
