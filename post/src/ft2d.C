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

#include <ft2d.h>

//////////////////////////////////////////////////////////////////////////////////////////


FFT2D::FFT2D( int Q_kw, int periods_1, int steps_pp_1, int screen_1,
	                int periods_2, int steps_pp_2, int screen_2 )
{
  if (Q_kw){
    periods_input_1  = periods_1;
    periods_input_2  = periods_2;
    steps_pp_input_1 = steps_pp_1;
    steps_pp_input_2 = steps_pp_2;
    screen_input_1   = screen_1;
    screen_input_2   = screen_2;

    steps_input_1  = periods_1 * steps_pp_1;
    steps_input_2  = periods_2 * steps_pp_2;
    dt_input_1     = (float) periods_input_1 / ( steps_input_1 - 1 );
    dt_input_2     = (float) periods_input_2 / ( steps_input_2 - 1 );
    steps_1        = steps_ft( steps_input_1 );
    steps_2        = steps_ft( steps_input_1 );
    steps_half_1   = (int) floor( 0.5*steps_1 + 0.5 );
    steps_half_2   = (int) floor( 0.5*steps_2 + 0.5 );
    dt_1           = (float) periods_input_1 / ( steps_1 - 1 );     // time step in periods
    dt_2           = (float) periods_input_2 / ( steps_2 - 1 );     // time step in periods
    df_1           = (float) 1.0 / periods_input_1;                 // frequency step
    df_2           = (float) 1.0 / periods_input_2;                 // frequency step

    nn             = new int [ 3 ];
    nn[1]          = steps_1;
    nn[2]          = steps_2;
    local          = fmatrix( 0, steps_1, 0, steps_2 );
    data           = new float [ 2 * steps_1 * steps_2 + 1 ];
    frequency_1    = new float [ steps_1 ];
    frequency_2    = new float [ steps_2 ];
    //  co  = fmatrix( 0, steps_half_1+1, 0, steps_2+1 );// pos and neg k, pos frequency!
    //  si  = fmatrix( 0, steps_half_1+1, 0, steps_2+1 );
    power          = fmatrix( 0, steps_half_1+1, 0, steps_2+1 );

    if (!nn || !local || !data || !frequency_1 || !frequency_2 || !power) {
      printf( "\n allocation failure in FFT2D::constructor" );
      exit(-1);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


int FFT2D::steps_ft( int steps_input )
{
  int MAX_EXPO = 14;
  int expo, steps_ft;

  expo = (int) ceil( log(steps_input)/log(2) );  // steps_ft = next power of 2
                                                 //            larger than steps_input
  if (expo > MAX_EXPO) {
    expo = MAX_EXPO;
    printf( "\n FFT2D: number of steps is 2^%d < %d\n", MAX_EXPO, steps_input );
  }

  steps_ft = (int) pow( 2, expo );

  return steps_ft;
}


//////////////////////////////////////////////////////////////////////////////////////////


float FFT2D::window( float t, int screen_input, int periods_input )
// input: time t in periods
// window cuts off smoothly first and last period
{
  float tscr = (float) screen_input;
  float toff = (float) periods_input - tscr;

  if ( t < tscr ) return pow( sin(0.5*PI*t/tscr), 2 );
  else {
    if ( t > toff ) return pow( sin(0.5*PI*(1.0-(t-toff)/tscr)), 2 );
    else            return 1.0;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void FFT2D::RealFt( float **input )
/*
  takes      real input matrix [0, ..., steps_input_1 - 1] [0, ..., steps_input_2 - 1 ]
  interpolates to local matrix [0, ..., steps_1 - 1] [0, ..., steps_2 - 1],
              where steps_1 and steps_2 are powers of 2, usually larger than steps_input

                                      i.e.  t_1=0 .... (steps_input_1-1) dt_1
                                      i.e.  t_2=0 .... (steps_input_2-1) dt_2

  returns   real output arrays [0, ... , steps_half_1] [0, ..., steps_2]

                                      i.e.  f_1= 0           .... 1/(2*dt_1)
                                      i.e.  f_2= -1/(2*dt_2) .... 1/(2*dt_2)

  returns cosine-transform  co  = re
          sine-transform    si  = im
          power spectrum    pow = 2*(re*re+im*im)
*/
{
  int   i, index, j;
  float t, x;
  int   th, tl, xh, xl;
  float re, im;

  for( i=0; i<steps_1; i++ )                   // interpolate
    {
      for( j=0; j<steps_2; j++ )
	{
	  t = (float) i * (steps_input_1-1) / (steps_1-1);
	  x = (float) j * (steps_input_2-1) / (steps_2-1);
	  tl = (int) floor( t );
	  xl = (int) floor( x );
	  if ( i<steps_1-1) th = tl+1;
	  else              th = tl;
	  if ( j<steps_2-1) xh = xl+1;
	  else              xh = xl;

	  t  = th - t;
	  x  = xh - x;

	  local[i][j] =        t*x*input[tl][xl] +       t*(1.0-x)*input[tl][xh];
	  local[i][j] += (1.0-t)*x*input[th][xl] + (1.0-t)*(1.0-x)*input[th][xh];
	  local[i][j] *= window( dt_1*i, screen_input_1, periods_input_1 );
	  local[i][j] *= window( dt_2*j, screen_input_2, periods_input_2 );

	  data[2*i*steps_2 + 2*j+1] = local[i][j];
	  data[2*i*steps_2 + 2*j+2] = 0;
	}
    }

  fftn(data,nn,2,1);

  for( i=0; i<=steps_half_1; i++ )
    {
      if (i==0) index = 0;                  // only negative frequencies
      else      index = steps_1 - i;        // only negative frequencies
      //      index = i;                            // only positive frequencies

      frequency_1[i] = df_1 * i;

      for( j=0; j<steps_2; j++ )            // positive and negative k-values
	{
	  if (j<steps_half_2) frequency_2[j+steps_half_2] = df_2 * j;
	  else                frequency_2[j-steps_half_2] = df_2 * (j-steps_2);

	  re = data[2*index*steps_2 + 2*j+1] / (steps_1 * steps_2);     // cos-part
	  im = data[2*index*steps_2 + 2*j+2] / (steps_1 * steps_2);     // sin-part

	  if (j<steps_half_2) {
	    //	       co[ i ][ j + steps_half_1 ] = re;
	    //	       si[ i ][ j + steps_half_1 ] = im;
	    power[ i ][ j + steps_half_1 ] = ( re*re+im*im );
	  }
	  else {
	    //	       co[ i ][ j - steps_half_1 ] = re;
	    //	       si[ i ][ j - steps_half_1 ] = im;
	    power[ i ][ j - steps_half_1 ] = ( re*re+im*im );
	  }
      }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void FFT2D::fftn(float *data, int *nn, int ndim, int isign)
// numerical recipies routine "fourn.c"
{
	int i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	int ibit,idim,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	float theta,wi,wpi,wpr,wr,wtemp;

	ntot=1;
	for (idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=wr*data[k2]-wi*data[k2+1];
						tempi=wr*data[k2+1]+wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}

#undef SWAP


//////////////////////////////////////////////////////////////////////////////////////////
//eof

