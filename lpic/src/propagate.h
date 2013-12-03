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
#include <pulse.h>
#include <box.h>
#include <particle.h>
#include <stack.h>
#include <diagnostic.h>
#include <uhr.h>
#include <readfile.h>

class input_propagate {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );

public:
  double start_time, stop_time;

  int    n_domains;

  int    Q_restart;
  char   restart_file[filename_size];

  input_propagate( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class propagate {
public:
                  propagate( parameter &p, domain &grid );
    void               loop( parameter &p, box &sim,
			     pulse &laser_front, pulse &laser_rear,
			     diagnostic &diag );

private:
    input_propagate input;

    stack      stk;
    readfile   rf;
    double     time, start_time, stop_time;
    double     dt, dx, idx;                  // timestep and grid spacing
    double     Gamma;                        // gamma factor due to Lorentz Transformation
    int        domain_number;                // domain number
    int        n_domains;                    // # of domains

    ofstream grid_file;

    char errname[filename_size];

    void                clear_grid( domain &grid );
    void                    fields( domain &grid, pulse &laser_front, pulse &laser_rear );
    void                 particles( domain &grid );
    void         reflect_particles( domain &grid );
    inline void	        accelerate( struct cell *cell, struct particle *part );
    inline void	      accelerate_1( struct cell *cell, struct particle *part );
    inline void	      accelerate_2( struct cell *cell, struct particle *part );
    inline void               move( struct particle *part );
    inline void has_to_change_cell( struct cell *cell, struct particle *part );
    inline void     do_change_cell( domain &grid );
    inline void     deposit_charge( struct cell *cell, struct particle *part );
    inline void    deposit_current( struct cell *cell, struct particle *part );
    inline void       mask_current( domain &grid );
    inline double             mask( int i );
    inline void           left_one( struct cell *cell, struct particle *part );
    inline void      left_two_left( struct cell *cell, struct particle *part );
    inline void     left_two_right( struct cell *cell, struct particle *part );
    inline void          right_one( struct cell *cell, struct particle *part );
    inline void    right_two_right( struct cell *cell, struct particle *part );
    inline void     right_two_left( struct cell *cell, struct particle *part );

    inline double        weighting( struct cell *cell, struct particle *part );
    inline double      weighting_0( struct cell *cell, struct particle *part );

};
#endif

