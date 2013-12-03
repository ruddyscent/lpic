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

#ifndef DOMAIN_H
#define DOMAIN_H

#include <common.h>
#include <fstream.h>
#include <stdio.h>
#include <iomanip.h>
#include <math.h>
#include <error.h>
#include <cell.h>
#include <particle.h>
#include <parameter.h>
#include <readfile.h>

class input_domain {
private:
  char errname[filename_size];

public:
  int      Q_restart;           // start from t>0, using restart files
  char     restart_file[filename_size];
  int      Q_restart_save;

  int      n_domains;           // grid
  int      cells;
  int      cells_per_wl;
  int      cells_left;
  int      cells_ramp;
  int      cells_plasma;
  double   dx;

  double   n_ion_over_nc;       // plasma density
  double   n_el_over_nc;

  int      nsp;                 // particles
  int      *ppc;
  int      *fix;
  double   *z, *zmax, *m;
  double   *vtherm;

  double   angle;                    // angles of incidence
  double   Gamma, Beta;              // Lorentz transformation
  int      spp;                      // steps per laser period

  readfile rf;
  input_domain( parameter &p );
  void save( parameter &p );
};

//////////////////////////////////////////////////////////////////////////////////////////

class domain {
 private:

  char         errname[filename_size];
  int          domain_number;
  int          n_domains;
  char         path[filename_size];
  input_domain input;

  void restart_configuration( void );
  void        set_boundaries( void );
  void           chain_cells( void );
  void            init_cells( void );
  void       chain_particles( void );
  void        init_particles( void );
  double        gauss_rand48( void );
  double    exponential_rand( double ); // ## exponential velocity distribution

public:

  int         n_left;     // cell number at the left boundary
  int         n_right;    // cell number at the right boundary
  int         n_cells;    // number of cells in this domain

  double dx;              // cell width

  struct cell *Lbuf;      // lhs: left buffer cell  left of left
  struct cell *lbuf;      //      right buffer      left of left
  struct cell *left;      // pointer to the first occupied cell
  struct cell *right;     // pointer to the last occupied cell
  struct cell *rbuf;      // rhs: left buffer cell  right of right
  struct cell *Rbuf;      //      right buffer      right of right
  struct cell *dummy;     // definitely the last one

  int n_el;               // # of electrons
  int n_ion;              // # of ions
  int n_part;             // total # particles

                    domain( parameter &p );
  void     count_particles( void );
  void               check( void );

  void         reo_to_prev( int request_to_prev, int *cells_to_prev, int *parts_to_prev );
  void         reo_to_next( int request_to_next, int *cells_to_next, int *parts_to_next );
  void  reo_delete_to_prev( int cells_to_prev, int parts_to_prev );
  void  reo_delete_to_next( int cells_to_next, int parts_to_next );
  void reo_alloc_from_prev( int cells_from_prev, int parts_from_prev );
  void reo_alloc_from_next( int cells_from_next, int parts_from_next );
  void reo_update_n_el_n_ion( int el_count, int ion_count );
};


#endif



