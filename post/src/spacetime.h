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

#ifndef SPACETIME_H
#define SPACETIME_H

#include <common.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <parameter.h>
#include <utilities.h>
#include <error.h>
#include <math.h>
#include <readfile.h>
#include <ft.h>
#include <ft2d.h>

class input_spacetime {
private:
  char errname[filename_size];

public:

  int   t_start, t_stop;
  float x_start, x_stop, x_offset;
  int   periods_x, periods_t, cells_per_wl, steps_per_period;
  int   average;
  int   size;
  float contour_1, contour_2, contour_3;
  int   Q_de, Q_di, Q_jx, Q_jy, Q_jz, Q_ex, Q_ey, Q_ez, Q_bx, Q_by, Q_bz, Q_edens,
        Q_de_fi, Q_de_ii;
  float C_de, C_di, C_jx, C_jy, C_jz, C_ex, C_ey, C_ez, C_bx, C_by, C_bz, C_edens,
        C_de_fi, C_de_ii;
  int   Q_kw, Q_kt;
  float C_kw, C_kt;
  float K_cut, W_cut;

  readfile rf;
  void save( parameter &p );

  input_spacetime( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class spacetime {

 private:
  input_spacetime input;
  FFT             ft;
  FFT2D           ft2d;

  int           t_start_in,  t_stop_in, t_steps_in;
  float         x_start_in, x_stop_in;
  int           x_steps_in;
  int           spp, spl;
  int           fnumber;
  float         **matrix_read;
  unsigned char **matrix_write;
  unsigned char *vector_write;
  float         **power_spectrum, **kspace;
  float         **kw;
  char          *input_path;
  char          *output_path;
  char          errname[filename_size];

 public:
                   spacetime( parameter &p );
  void                select( void );
  void              xt_kt_kw( char *unit, float cut, int sign, int scale_write );

  void read_input_array_size( char *input );
  void                  read( char *input );
  void                smooth( float **field );
  void                 scale( float cut, int sign,
			      float **m, unsigned char **mw, int nt, int nx );
  void                 write( char *unit, int scale_write );
  void   write_idl_header_xt( float cut, int sign, char *idlname, char *unit,
		   	      char *axislabel );
  void           transform_k( float **matrix_read );
  void     write_transform_k( char* unit, unsigned char **m );
  void   write_idl_header_kt( float cut, int sign, char *idlname, char *unit,
		  	      char *axislabel );
  void          transform_kw( float **matrix_read );
  void    write_transform_kw( char* unit, unsigned char **m );
  void   write_idl_header_kw( float cut, int sign, char *idlname, char *unit,
			      char *axislabel );
};

#endif

