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

#ifndef BOX_H
#define BOX_H

#include <common.h>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <error.h>
#include <cell.h>
#include <particle.h>
#include <parameter.h>
#include <domain.h>
#include <diagnostic.h>
#include <uhr.h>
#include <readfile.h>

#ifdef LPIC_PARALLEL
#include <network.h>
#endif


class input_box {
private:
  char errname[filename_size];

public:
  int      Q_restart;           // start from t>0, using restart files
  char     restart_file[filename_size];
  int      Q_restart_save;
  char     restart_file_save[filename_size];

  int      n_domains;
  int      Q_reorganize;
  int      delta_reo;

  int      nsp;

  readfile rf;
  void save( parameter &p );

  input_box( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class box {
private:
  char errname[filename_size];
  input_box input;
#ifdef LPIC_PARALLEL
  readfile rf;
#endif

public:
  box( parameter &p );

#ifdef LPIC_PARALLEL
  void new_global_particle_numbers( domain &grid, network &talk );
  void  com_total_particle_numbers( domain &grid, network &talk );
  void                reorganize_f( domain &grid, network &talk );
  void                  reorganize( domain &grid, network &talk, double time );
  void             init_reorganize( parameter &p );
  void            count_reorganize( void );
  void               particle_load( domain &grid );

  network  talk;

  struct reo_struct {
  int Q_reorganize;
  int delta_reo;
  int count_reo;
  char file[filename_size];
  } reo;
#endif

  struct rest_struct {
  int Q_restart_save;
  int delta_rest;
  int count_rest;
  } rest;

  void    init_restart( parameter &p );
  void   count_restart( void );
  void    restart_save( diagnostic &diag, double time, parameter &p,
			uhr &zeit, uhr &zeit_particles, uhr &zeit_fields,
			uhr &zeit_diagnostic);

  domain   grid;

  int n_domains;          // # of domains

  int n_el;               // total # of electrons
  int n_ion;              // total # of ions
  int n_part;             // total # particles
};

#endif



