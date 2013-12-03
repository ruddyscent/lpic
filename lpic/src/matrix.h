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

#ifndef MATRIX_H
#define MATRIX_H

#include <common.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <fstream.h>
#include <string.h>
#include <iomanip.h>
double** dmatrix( long nrl, long nrh, long ncl, long nch );
void     delete_dmatrix( double **m, long nrl, long nrh, long ncl, long nch );
float**  matrix( long nrl, long nrh, long ncl, long nch );
void     delete_matrix( float **m, long nrl, long nrh, long ncl, long nch );
float**  fmatrix( long nrl, long nrh, long ncl, long nch );
void     delete_fmatrix( float **m, long nrl, long nrh, long ncl, long nch );
int**    imatrix( long nrl, long nrh, long ncl, long nch );
void     delete_imatrix( int **m, long nrl, long nrh, long ncl, long nch );
unsigned char **ucmatrix(long nrl, long nrh, long ncl, long nch);
void     delete_ucmatrix(unsigned char **m, long nrl, long nrh, long ncl, long nch);

void error(char* s1, char*  s2="", char* s3="", char*  s4="");
void error(char* s1, double d2,    char* s3="", char*  s4="");

#endif
