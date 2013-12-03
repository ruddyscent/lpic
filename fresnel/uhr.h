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

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <fstream.h>
#include <iomanip.h>
#include <error.h>
#include <parameter.h>

class uhr {
 private:

  double h_cpu, m_cpu, s_cpu;                // total time
  double h_sys, m_sys, s_sys;                // intermediate time
  clock_t start_clock, stop_clock;           // clock ticks
  time_t  start_time, stop_time;             // system time
  char errname[filename_size];

 public:

  uhr( parameter &p );
  void       proc( void );
  void        sys( void );
  void       init( void );
  void       exit( void );
};


#endif
