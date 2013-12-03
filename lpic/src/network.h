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

#include <config.h>

#ifdef LPIC_PARALLEL
#ifdef LPIC_PVM

#ifndef NETWORK_H
#define NETWORK_H

#include <pvm3.h>

#include <common.h>
#include <parameter.h>
#include <error.h>
#include <cell.h>
#include <particle.h>
#include <domain.h>

class network {

 private:

  int domain_number;
  int n_domains;

  int tid;        // my_tid()
  int tid_prev;   // this task was spawned by 'prev'
  int tid_next;   // this task will have spawned 'next'

  char errname[filename_size];

 public:

  network( parameter &p );

  void      start_next_task( parameter &p );

  void                field( int time_step, domain* grid );
  void        field_get_cpy( struct cell*, int ptid, int time_step );
  void       field_send_cpy( struct cell*, int ptid, int time_step );
  void            particles( int time_step, domain* grid );
  void        particles_get( struct cell* cell, int ptid, int time_step,
			     int *el_count, int *ion_count );
  void       particles_send( struct cell* cell, int ptid, int time_step,
			     int *el_count, int *ion_count );
  void              current( int time_step, domain* grid );
  void          current_get( struct cell* cell, int ptid, int time_step );
  void         current_send( struct cell* cell, int ptid, int time_step );
  void      current_get_cpy( struct cell* cell, int ptid, int time_step );
  void     current_send_cpy( struct cell* cell, int ptid, int time_step );
  void              density( int time_step, domain* grid );
  void          density_get( struct cell* cell, int ptid, int time_step );
  void         density_send( struct cell* cell, int ptid, int time_step );

  void            current_1( int time_step, domain* grid );
  void            current_2( int time_step, domain* grid );
  void       current_get_12( struct cell* cell, int ptid, int time_step );
  void      current_send_12( struct cell* cell, int ptid, int time_step );

  void   get_part_numbers_from_prev( int* number, int n );
  void    send_part_numbers_to_next( int* number, int n );
  void  get_total_numbers_from_next( int* number, int n );
  void   send_total_numbers_to_prev( int* number, int n );

  void       reo_get_mesg_from_prev( int* mesg );
  void        reo_send_mesg_to_next( int* mesg );
  void                reo_from_prev( int *cells_from_prev, int *parts_from_prev );
  void                reo_from_next( int *cells_from_next, int *parts_from_next );
  void reo_recieve_from_prev_and_unpack( int cells_from_prev, int parts_from_prev,
                                     struct cell* firstcell,
				     int *el_count, int *ion_count );
  void reo_recieve_from_next_and_unpack( int cells_from_next, int parts_from_next,
                                     struct cell* lastcell,
				     int *el_count, int *ion_count );
  void                  reo_to_prev( int cells_from_prev, int parts_from_prev );
  void                  reo_to_next( int cells_to_next, int parts_to_next );
  void    reo_pack_and_send_to_prev( int cells_to_prev, int parts_to_prev,
                                     struct cell* firstcell );
  void    reo_pack_and_send_to_next( int cells_to_next, int parts_to_next,
                                     struct cell* lastcell );

  void                pack_particle( struct particle *part );
  void              unpack_particle( struct particle *part );
  void                    pack_cell( struct cell *cell );
  void                  unpack_cell( struct cell *cell );
  void          pack_cell_as_buffer( struct cell *cell );

  void             end_task( void );
};

#endif
#endif
#endif
