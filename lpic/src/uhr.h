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

#ifndef UHR_H
#define UHR_H

#include <common.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <fstream.h>
#include <iomanip.h>
#include <error.h>
#include <parameter.h>

class input_uhr {
private:
  char          errname[filename_size];
  readfile      rf;
  void          save( parameter &p );

public:
  int Q_restart;
  char restart_file[filename_size];
  int Q_restart_save;
  char restart_file_save[filename_size];

  input_uhr( parameter &p );
};


//////////////////////////////////////////////////////////////////////////////////////////


class uhr {
 private:

  double    h_cpu, m_cpu, s_cpu, sec_cpu;       // cpu time
  clock_t   start_tics, stop_tics, tics;        // processor clocks
  double    h_sys, m_sys, s_sys, sec_sys;       // system time
  time_t    start_time, stop_time;

  char      errname[filename_size];
  char      uhrname[filename_size];
  char      path[filename_size];

  int       domain_number;
  readfile  rf;
  input_uhr input;

 public:

  uhr( parameter &p, char *name );
  void        reset( void );
  void        start( void );
  void stop_and_add( void );
  void          add( void );
  void          sys( void );
  void  seconds_cpu( void );
  void  seconds_sys( void );
  void      restart( void );
  void restart_save( void );
};


#endif



