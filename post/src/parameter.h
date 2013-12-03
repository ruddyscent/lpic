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
//
//////////////////////////////////////////////////////////////////////////////////////////


#ifndef PARAMETER_H
#define PARAMETER_H

#include <common.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <utilities.h>
#include <error.h>
#include <readfile.h>

//////////////////////////////////////////////////////////////////////////////////////////

class parameter  {

private:

  readfile rf;

public:

  char   *read_filename;
  char   *file_path;
  char   *output_path;
  char   *save_filename;
  char   *save_path_name;
  char   *my_name;
  char   *errname;

       parameter(int argc, char **argv);
  void      save(char *fname);
//  void      read(char *fname);
};

//////////////////////////////////////////////////////////////////////////////////////////
#endif





