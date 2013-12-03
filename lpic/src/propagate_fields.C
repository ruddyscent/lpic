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

#include <propagate.h>

void propagate::fields( domain &grid, pulse &laser_front, pulse &laser_rear )
{
  static error_handler bob("propagate::fields", errname);

  struct cell *cell;
  double fp, fp_old, gm, gm_old;
  double pidt = PI * dt;

  cell=grid.left;

  fp_old   = cell->fp;
  gm_old   = cell->gm;

  if (domain_number==1) {
    cell->fp = laser_front.Qy * laser_front.field( time );
    cell->gm = laser_front.Qz * laser_front.field( time + laser_front.shift );
  }
  else {
    cell->fp = grid.lbuf->fp - pidt * grid.lbuf->jy;
    cell->gm = grid.lbuf->gm - pidt * grid.lbuf->jz;
  }

  cell->fm = cell->next->fm - pidt * cell->jy;
  cell->gp = cell->next->gp - pidt * cell->jz;

  cell->ey = cell->fp + cell->fm;   cell->ez = cell->gp + cell->gm;
  cell->bz = cell->fp - cell->fm;   cell->by = cell->gp - cell->gm;

  cell->ex -= 2.0 * pidt * cell->jx;

  for( cell=grid.left->next; cell->next!=grid.rbuf; cell=cell->next )
    {
      fp = cell->fp;                                   // could be done in a more elegant
      cell->fp = fp_old - pidt * cell->prev->jy;       // way: i.e. do the loop twice,
      fp_old = fp;                                     // forward for fm and gp,
                                                       // backward for fp and gm
      gm = cell->gm;
      cell->gm = gm_old - pidt * cell->prev->jz;
      gm_old = gm;

      cell->fm = cell->next->fm - pidt * cell->jy;
      cell->gp = cell->next->gp - pidt * cell->jz;

      cell->ey = cell->fp + cell->fm;   cell->bz = cell->fp - cell->fm;
      cell->ez = cell->gp + cell->gm;   cell->by = cell->gp - cell->gm;

      cell->ex -= 2.0 * pidt * cell->jx;
    }

  cell->fp = fp_old - pidt * cell->prev->jy;
  cell->gm = gm_old - pidt * cell->prev->jz;

  if (domain_number==n_domains) {
    cell->fm = laser_rear.Qy * laser_rear.field( time );
    cell->gp = laser_rear.Qz * laser_rear.field( time + laser_rear.shift);
  }
  else {
    cell->fm = grid.rbuf->fm - pidt * cell->jy;
    cell->gp = grid.rbuf->gp - pidt * cell->jz;
  }

  cell->ey = cell->fp + cell->fm;   cell->bz = cell->fp - cell->fm;
  cell->ez = cell->gp + cell->gm;   cell->by = cell->gp - cell->gm;

  cell->ex -= 2*pidt * cell->jx;
}


//////////////////////////////////////////////////////////////////////////////////////////
//EOF
