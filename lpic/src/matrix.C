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

#include <iostream>
#include <matrix.h>

using namespace std;

#define NR_END 1
#define FREE_ARG char*

//////////////////////////////////////////////////////////////////////////////////////////

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

double **dmatrix(long nrl, long nrh, long ncl, long nch)
// allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

  // allocate pointers to rows
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) error("allocation failure 1 in dmatrix()");

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
  nrh = nch = 0;
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/////////////////////////////////////////////////////////////////////////////////////////

float **matrix(long nrl, long nrh, long ncl, long nch)
// allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  float **m;

  // allocate pointers to rows
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) error("allocation failure 1 in dmatrix()");

  // allocate rows and set pointers to them
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) error("allocation failure 2 in dmatrix()");
  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  // return pointer to array of pointers to rows
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////

void delete_matrix( float **m, long nrl, long nrh, long ncl, long nch)
// free a float matrix allocated by matrix()
{
  nrh = nch = 0;
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
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
  nrh = nch = 0;
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

//////////////////////////////////////////////////////////////////////////////////////////

int **imatrix(long nrl, long nrh, long ncl, long nch)
// allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  int **m;

  // allocate pointers to rows
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) error("allocation failure 1 in imatrix()");

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
  nrh = nch = 0;
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/////////////////////////////////////////////////////////////////////////////////////////

unsigned char **ucmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a unsigned char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  unsigned char **m;

  /* allocate pointers to rows */
  m=(unsigned char **) malloc((size_t)((nrow+NR_END)*sizeof(unsigned char*)));
  if (!m) error("allocation failure 1 in ucmatrix()");

  /* allocate rows and set pointers to them */
  m[nrl]=(unsigned char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned char)));
  if (!m[nrl]) error("allocation failure 2 in ucmatrix()");
  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/////////////////////////////////////////////////////////////////////////////////////////

void delete_ucmatrix( unsigned char **m, long nrl, long nrh, long ncl, long nch)
// free a unsigned char matrix allocated by ucmatrix()
{
  nrh = nch = 0;
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/////////////////////////////////////////////////////////////////////////////////////////
//eof






