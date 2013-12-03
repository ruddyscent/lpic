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

#include <propagate.h>


//////////////////////////////////////////////////////////////////////////////////////////

propagate::propagate(parameter &p)
    : zeit(p)
{
  sprintf( errname, "%s/error", p.path );
  static error_handler bob("propagate::Constructor",errname);

  dt    = 1.0 / p.prop.spp;
  dw    = 1.0 / ( p.prop.time_stop - p.prop.time_start ) / 50;
  wmax  = 5;
  t     = p.prop.time_start;
  T     = p.pulse.duration;

  reflex_avg = 0;

  strcpy(p.errname,errname);
}

//////////////////////////////////////////////////////////////////////////////////////////

void propagate::loop(parameter &p)
{
  static error_handler bob("propagate::loop",errname);

  out_int(p);      // plot the integrand

  for( t = p.prop.time_start, t_avg=0; t <= p.prop.time_stop+dt; t += dt, t_avg++ )
    {
      w = 0.5*dw; reflex = 0;
      do
	{
	  r = reflectivity( w, p.plasma.wp, p.pulse.angle_rad, p.pulse.polarization );
	  phi = phase( w, p.plasma.wp, p.pulse.angle_rad, p.pulse.polarization );
	  integrand = r * ( dif((w-1)*PI*T) - dif((w+1)*PI*T) ) * sin( w*PI*(2*t-T)-phi );
	  reflex += dw * T * integrand;

	  w += dw;
	}
      while( w<wmax );

      out(p);    // plot the flux (f(t)/f0)^2

      reflex_avg += 2 * pow( reflex, 2 ) / p.prop.spp;

                 // plot the ratio of time averaged fluxes once per period
      if ( t_avg==p.prop.spp ) { out_avg(p); t_avg=0; reflex_avg=0; }

      zeit.proc();              // update clock
    }

  zeit.exit();
}


//////////////////////////////////////////////////////////////////////////////////////////

double propagate::reflectivity( double w, double wp, double angle, int polarization )
{
  static error_handler bob("propagate::reflectivity",errname);

  if (w*cos(angle)>wp) {

    double n = sqrt( 1 - pow(wp/w,2) );

    rs = cos(angle) - sqrt( n*n-sin(angle)*sin(angle) );
    rs /= ( cos(angle) + sqrt( n*n-sin(angle)*sin(angle) ) );

    rp = sqrt( 1-pow(sin(angle)/n,2) ) - n*cos(angle);
    rp /= ( sqrt( 1-pow(sin(angle)/n,2) ) + n*cos(angle) );
  }
  else {
    rs = rp = 1;
  }

  if (polarization==1) return rs;
  else                 return rp;
}

//////////////////////////////////////////////////////////////////////////////////////////

double propagate::phase( double w, double wp, double angle, int polarization )
{
  static error_handler bob("propagate::phase",errname);
  double brev = w*cos(angle)/wp;


  if (brev>=1) {

    double n = sqrt( 1 - pow(wp/w,2) );

    phis = phip = 0;
  }
  else {

    // from Born and Wolf, page 49
    phis = - 2 * atan( sqrt( 1-brev*brev ) / brev );

    // from Born and Wolf, page 49
    phip = + 2 * atan( wp/w * cos(angle) * ( pow(w/wp,2)-1 ) / sqrt( 1-brev*brev ) );
  }

  if (polarization==1) return phis;
  else                 return phip;
}

//////////////////////////////////////////////////////////////////////////////////////////

double propagate::dif( double x )
{
  static error_handler bob("propagate::dif",errname);

  if ( x==0 ) return 1;
  else return sin(x)/x;
}

//////////////////////////////////////////////////////////////////////////////////////////

void propagate::out( parameter &p )
{
  static error_handler bob("propagate::out",errname);

  char filename[filename_size];

  sprintf(filename,"%s/reflex", p.path );

  reflex_file.open(filename,ios::app);
  if (!reflex_file) bob.error( "cannot open file", filename );

  reflex_file.precision( 3 );
  reflex_file.setf( ios::showpoint | ios::scientific );

  reflex_file << setw(10)  << t  << setw(10) << pow( p.pulse.amplitude*reflex,2 ) << endl;

  reflex_file.close();
}

//////////////////////////////////////////////////////////////////////////////////////////

void propagate::out_avg( parameter &p )
{
  static error_handler bob("propagate::out_avg",errname);

  char filename[filename_size];

  sprintf(filename,"%s/reflex-avg", p.path );

  reflex_file_avg.open(filename,ios::app);
  if (!reflex_file_avg) bob.error( "cannot open file", filename );

  reflex_file_avg.precision( 3 );
  reflex_file_avg.setf( ios::showpoint | ios::scientific );

  reflex_file_avg << setw(10)  << t  << setw(10) << reflex_avg << endl;

  reflex_file_avg.close();
}

//////////////////////////////////////////////////////////////////////////////////////////

void propagate::out_int( parameter &p )
{
  static error_handler bob("propagate::out_int",errname);

  char filename[filename_size];

  sprintf(filename,"%s/reflex-int", p.path );

  reflex_file_int.open(filename,ios::app);
  if (!reflex_file_int) bob.error( "cannot open file", filename );

  reflex_file_int.precision( 3 );
  reflex_file_int.setf( ios::showpoint | ios::scientific );

  reflex_file_int << "#"
                  << setw(12)  << "w"     << setw(13) << "r(w)"
		  << setw(13)  << "phase"
                  << setw(13)  << "dif()" << setw(13) << "r(w)*dif()" << endl;

  w = 0; reflex = 0;
  do
    {
      r = reflectivity( w, p.plasma.wp, p.pulse.angle_rad, p.pulse.polarization );
      phi = phase( w, p.plasma.wp, p.pulse.angle_rad, p.pulse.polarization );
      integrand = r * ( dif( (w-1)*PI*T ) - dif( (w+1)*PI*T) );

      reflex_file_int << setw(13) << w
                      << setw(13) << r
		      << setw(13) << phi
		      << setw(13) << dif( (w-1)*PI*T ) - dif( (w+1)*PI*T)
		      << setw(13) << integrand << endl;

      w += 5*dw;
    }
  while( w<2 );

  reflex_file_int.close();
}

//////////////////////////////////////////////////////////////////////////////////////////
//EOF



