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

//////////////////////////////////////////////////////////////////////////////////////////
//
// laser pulses : linear, sin, sin^2
//
//////////////////////////////////////////////////////////////////////////////////////////


#ifndef PULSE_H
#define PULSE_H

#include <common.h>
#include <math.h>
#include <error.h>
#include <parameter.h>
#include <fstream.h>
#include <iomanip.h>
#include <readfile.h>

//////////////////////////////////////////////////////////////////////////////////////////

class input_pulse {
private:
  char errname[filename_size];

public:
  int    Q;

  int    shape;         // 1=linear, 2=sin, 3=sinsqr
  int    polarization;  // 1=s, 2=p, 3=circular
  double a0;            // dimensionless amplitude of 1omega light
  double raise;         // raise time in periods
  double duration;      // pulse duration in periods
  double a2, a3;        // amplitudes of 2omega and 3omega
  double p2, p3;        // phases of 2omega and 3omega

  int    Q_save;
  double save_step;

  double time_start, time_stop;  // simulation time interval in periods

  int    Q_restart;

  readfile rf;
  input_pulse( parameter &p, char *side, int pulse_number );
  void save( parameter &p, int pulse_number );
};

//////////////////////////////////////////////////////////////////////////////////////////

class pulse {

private:

  static int pulse_number; // counts elements of type pulse

  input_pulse input;

  char   errname[filename_size];
  char   *path;

public:

  double Qy;            // polarization dependent
  double Qz;            // "
  double shift;         // relative shift between y and z for circular polarization


            pulse( parameter &p, char *side );
  void       save( double t_start, double t_stop, double t_step );
  double    field( double t );
  double envelope( double t );
};


//////////////////////////////////////////////////////////////////////////////////////////

#endif


