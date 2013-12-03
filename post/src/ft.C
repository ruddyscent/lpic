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

#include <ft.h>

//////////////////////////////////////////////////////////////////////////////////////////


FFT::FFT( int periods, int steps_pp, int screen )
{
  periods_input  = periods;
  steps_pp_input = steps_pp;
  periods_screen = screen;

  steps_input    = periods * steps_pp;
  dt_input       = (double) periods_input / ( steps_input - 1 );
  steps          = steps_ft();
  steps_half     = (int) floor( 0.5*steps + 0.5 );
  dt             = (double) periods_input / ( steps - 1 );         // time step in periods
  df             = (double) 1.0 / periods_input;                   // frequency step

  local          = new double [steps];
  data           = new double [2*steps+1];
  frequency      = new double [steps];
  co             = new double [steps];
  si             = new double [steps];
  power          = new double [steps];
  corr           = new double [steps];
}


//////////////////////////////////////////////////////////////////////////////////////////


int FFT::steps_ft( void )
{
  int max_expo = 14;
  int expo, steps_ft;

  expo = (int) ceil( log(steps_input)/log(2) );  // steps_ft = next power of 2
                                                 //            larger than steps_input
  if (expo > max_expo) expo = max_expo;

  steps_ft = (int) pow( 2, expo );

  return steps_ft;
}


//////////////////////////////////////////////////////////////////////////////////////////


double FFT::window( double t )
// input: time t in periods
// window cuts off smoothly first and last period
{
  double tscr = (double) periods_screen;
  double toff = (double) periods_input - tscr;

  if ( t < tscr ) return pow( sin(0.5*PI*t/tscr), 2 );
  else {
    if ( t > toff ) return pow( sin(0.5*PI*(1.0-(t-toff)/tscr)), 2 );
    else            return 1.0;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


int FFT::RealFt( float* input )
/*
  takes      real input 0, ... , m-1
  interpolates to local 0, ... , n-1, where n is a power of 2, usually larger than m
                                      i.e.  t=0 .... t=(n-1)dt
  returns   real output 0, ... , n/2  i.e.  f=0 .... f=1/(2*dt)

  returns cosine-transform  co  = re
          sine-transform    si  = im
          power spectrum    pow = 2*(re*re+im*im)

*/
{
  int i, low, high;
  double t, th, tl;
  double re, im;

  for( i=0; i<steps; i++ )                     // interpolates vector 1, ... , steps_input
    {                                          //           to vector 0, ... , steps-1
      low = (int) floor( (double) i * (steps_input-1) / (steps-1) );
      if ( i<steps-1) high = low+1;
      else            high = low;

      tl = dt_input * low;
      th = dt_input * high;
      t  = dt * i;

      local[i]    = (th-t)/dt_input * input[low] + (t-tl)/dt_input * input[high];
      data[2*i+1] = window(t) * local[i];
      data[2*i+2] = 0;
    }

  fft(data,steps,1);

  for( i=0; i<=0.5*steps; i++ )     // positive frequency part, since input real
    {
      frequency[i] = df * i;

      re = data[2*i+1]/steps;     // cos-part
      im = data[2*i+2]/steps;     // sin-part

      if ( i==0 ) {
	co[i]    = re;
	si[i]    = im;
	power[i] = (re*re+im*im);
      }
      else {
	co[i]    = 2*re;
	si[i]    = 2*im;
	power[i] = 2*( re*re+im*im );
      }
    }

  return steps;
}


//////////////////////////////////////////////////////////////////////////////////////////


double FFT::frequency_filter( double freq, double mid, double width )
{
  if ( fabs(freq)-mid < -0.5*width )     return 0;
  else if ( fabs(freq)-mid < 0.5*width ) return cos(M_PI*(fabs(freq)-mid)/width);
  else                                   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////


int FFT::correlation( float* input, double mid, double width )
/*
  the correlation of 'input' is calculated and returned in vector 'corr'

  this is done modifying function FFT::RealFt():
  - input is Fourier transformed
  - the power spectrum is calculated
  - the inverse transform of the power spectrum is taken
  - the result is normalized, divided by the result at delay t=0
  'mid' and 'width' denote a frequency window that is used to filter the
  power spectrum
*/
{
  int i, low, high;
  double t, th, tl;
  double re, im;

  for( i=0; i<steps; i++ )                     // interpolates vector 1, ... , steps_input
    {                                          //           to vector 0, ... , steps-1
      low = (int) floor( (double) i * (steps_input-1) / (steps-1) );
      if ( i<steps-1) high = low+1;
      else            high = low;

      tl = dt_input * low;
      th = dt_input * high;
      t  = dt * i;

      local[i]    = (th-t)/dt_input * input[low] + (t-tl)/dt_input * input[high];
      data[2*i+1] = window(t) * local[i];
      data[2*i+2] = 0;
    }

  fft(data,steps,1);

  for( i=0; i<steps; i++ )         // for all frequencies
    {
      if (i<=steps_half) frequency[i] = df * i;
      else               frequency[i] = df * ( i - steps );

      re = data[2*i+1]/steps;      // cos-part
      im = data[2*i+2]/steps;      // sin-part

      re *= frequency_filter( frequency[i], mid, width );  // filter in frequency space
      im *= frequency_filter( frequency[i], mid, width );

      data[2*i+1] = re*re + im*im; // power spectrum
      data[2*i+2] = 0;             // real!
    }

  fft(data,steps,-1);              // inverse transform should be real

  for( i=0; i<steps; i++ )         // for all time delays
    {
      corr[i] = data[2*i+1] / data[1]; // normalize, correlation is real
    }

  return steps;
}


//////////////////////////////////////////////////////////////////////////////////////////


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void FFT::fft( double* data, int nn, int isign)
// numerical recipies routine "four1.c"
{
	int n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=2*mmax;
		theta=6.28318530717959/(isign*mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

#undef SWAP


//////////////////////////////////////////////////////////////////////////////////////////


void FFT::out( float* input, char* path, char* appendix )
{

  int MAX_STEPS = 5000;
  int sample    = 1 + (int) floor( (double) steps_input / MAX_STEPS );
  FILE* f;
  int i;
  double time;
  char file[filename_size];

  sprintf( file, "%s/ft.%s.trace_1", path, appendix );
  f=fopen( file, "w" );

  for( i=0; i<steps_input; i+=sample )
    {
      time=dt_input*i;
      fprintf( f, "\n %.5e %.5e", time, input[i] );
    }

  fclose( f );

  sprintf( file, "%s/ft.%s.trace_2", path, appendix );
  f=fopen( file, "w" );

  for( i=0; i<steps; i++ )
    {
      time=dt*i;
      fprintf( f, "\n %.5e %.5e", time, local[i] );
    }

  fclose( f );

  sprintf( file, "%s/ft.%s", path, appendix );
  f=fopen( file, "w" );

  for( i=0; i<=0.5*steps; i++ )
    {
      fprintf( f, "\n %.5e", frequency[i] );
      fprintf( f, " %.5e",   co[i] );
      fprintf( f, " %.5e",   si[i] );
      fprintf( f, " %.5e",   power[i] );
    }

  fclose( f );

}


//////////////////////////////////////////////////////////////////////////////////////////
//eof

