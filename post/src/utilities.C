/*
   This file is part of LPIC++, a particle-in-cell code for
   simulating the interaction of laser light with plasma.

   Copyright (C) 2001, 2002 Andreas Kemp
   Copyright (C) 1994-1997  Roland Lichters

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

#include <utilities.h>

// original version by R.Lichters taken from Numerical Recipies.
// changes by A.Kemp indicated by ## and date of change

#define NR_END 1
#define FREE_ARG char*

// iostream.h added by R.L. on 06.08.05, needed by cout and endl
// to compile with gcc 4 on Mac OS X 10.4
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////
void error(char* s1, char* s2, char *s3, char *s4)
{
  cout << "FAILURE: " << s1 << ' ' << s2 << s3 << s4 << endl;

  exit(1);
}

void error(char* s1, double d2, char *s3, char *s4)
{
  cout << "FAILURE: " << s1 << ' ' << d2 << s3 << s4 << endl;

  exit(1);
}

//////////////////////////////////////////////////////////////////////////////////////////

float **fmatrix(long nrl, long nrh, long ncl, long nch)
// allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  float **m;

  // allocate pointers to rows
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) error("allocation failure 1 in fmatrix()");

  m += NR_END;  // ##31.10.01
  m -= nrl;    // ##

  // allocate rows and set pointers to them
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) error("allocation failure 2 in fmatrix()");
  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  // return pointer to array of pointers to rows
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////

void delete_fmatrix( float **m, long nrl, long nrh, long ncl, long nch)
// free a float matrix allocated by fmatrix()
{
  nrh = nch = 0; // ##
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

//////////////////////////////////////////////////////////////////////////////////////////

double **dmatrix(long nrl, long nrh, long ncl, long nch)
// allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

  // allocate pointers to rows
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) error("allocation failure 1 in dmatrix()");

  m += NR_END;  // ##31.10.01
  m -= nrl;     // ##

  // allocate rows and set pointers to them
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) error("allocation failure 2 in dmatrix()");
  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  // return pointer to array of pointers to rows
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////

void delete_dmatrix( double **m, long nrl, long nrh, long ncl, long nch)
// free a double matrix allocated by dmatrix()
{
  nrh = nch = 0; // ##31.10.01
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/////////////////////////////////////////////////////////////////////////////////////////

unsigned char **ucmatrix(long nrl, long nrh, long ncl, long nch)
// allocate a unsigned char matrix with subscript range m[nrl..nrh][ncl..nch]
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  unsigned char **m;

  // allocate pointers to rows
  m=(unsigned char **) malloc((size_t)((nrow+NR_END)*sizeof(unsigned char*)));
  if (!m) error("allocation failure 1 in ucmatrix()");

  m += NR_END;  // ## 31.10.01
  m -= nrl;     // ##

  // allocate rows and set pointers to them
  m[nrl]=(unsigned char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned char)));
  if (!m[nrl]) error("allocation failure 2 in ucmatrix()");
  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  // return pointer to array of pointers to rows
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////

void delete_ucmatrix( unsigned char **m, long nrl, long nrh, long ncl, long nch)
// free a unsigned char matrix allocated by ucmatrix()
{
  nrh = nch = 0; // ## 31.10.01
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/////////////////////////////////////////////////////////////////////////////////////////

int **imatrix(long nrl, long nrh, long ncl, long nch)
// allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  int **m;

  // allocate pointers to rows
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) error("allocation failure 1 in imatrix()");

  m += NR_END;  // ##31.10.01
  m -= nrl;     // ##

  // allocate rows and set pointers to them
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) error("allocation failure 2 in ucmatrix()");
  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  // return pointer to array of pointers to rows
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////

void delete_imatrix( int **m, long nrl, long nrh, long ncl, long nch)
// free a unsigned char matrix allocated by imatrix()
{
  nrh = nch = 0; //##
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/////////////////////////////////////////////////////////////////////////////////////////

ifstream& operator>>(ifstream& input, Trash& trash)
// reads a file to the next occuring ':' by using  'file >> trash'
{
  do input.getline(trash.string,2);
  while( strstr(trash.string,":")==0 );

  return input;
}

/////////////////////////////////////////////////////////////////////////////////////////
//eof





