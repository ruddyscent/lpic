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

/////////////////////////////////////////////
//
// header file, input parameters to fresnel.C
//
/////////////////////////////////////////////


#ifndef PARAMETER_H
#define PARAMETER_H

#include <common.h>
#include <error.h>
#include <fstream.h>
#include <iomanip.h>
#include <string.h>
#include <utilities.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////////////////////////

struct parameter  {

  //--------------------------------------------------------------------------------------

  struct plasma_struct {
    double wp;                   // plasma frequency in units of laser frequency
  } plasma;

  //--------------------------------------------------------------------------------------

  struct pulse_struct {
    double amplitude;            // dimensionless laser amplitude (lpic++)
    double angle;                // angle of incidence in degree
    double angle_rad;            // angle of incidence in radiant
    int    polarization;         // s=1, p=2
    int    duration;             // pulse duration in periods
  } pulse;

  //--------------------------------------------------------------------------------------

  struct propagate_struct {
    double time_start;              // start time in periods
    double time_stop;               // stop time in periods
    int spp;                     // steps per period
  } prop;

  //--------------------------------------------------------------------------------------

  char   *read_filename;
  char   *my_name;
  char   *path;
  char   errname[filename_size];

  parameter(int argc, char **argv);
  void read(char *fname);
  void save(char *path, char *errname);
};

//////////////////////////////////////////////////////////////////////////////////////////
#endif





