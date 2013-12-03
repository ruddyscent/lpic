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
// tools for reading ascii files
// e.g. 'namelist' input
//
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef READFILE_H
#define READFILE_H

#define MAX_LINE_LENGTH 1000
#define MAX_COL 50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class readfile {
 public:

  readfile();

  void           openinput( char* );
  void          closeinput( void );
  int             setinput( char* );
  char*           getinput( char* );
  char*             setget( char*, char* );
  int        read_one_line( void );
  void      write_one_line( void );

  int            read_line( int *narg, char **arg );
  int             read_col( double *data, int col, int max_rows );

  void           copy_file( char *input, char *output );
  int        compare_files( char *name1, char *nsma2 );

  char**           cmatrix( long nrh, long nch );
  void        free_cmatrix( char **m );
  char*            cvector( long nch );
  void        free_cvector( char *m );

 private:

  int  already_open;
  FILE *fd;
  char *buffer;
  char *result;
};

#endif

