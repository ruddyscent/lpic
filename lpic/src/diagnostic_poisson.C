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

#include <diagnostic_poisson.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////


poisson::poisson( parameter &p, domain* grid )
  : rf(),
    input(p),
    stepper( input.stepper, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("poisson::Constructor",errname);

  double nach, vor;
  int    n = grid->n_cells;                  // number of grid points

  spp           = p.spp;
  domain_number = p.domain_number;
  output_path   = new char [filename_size];

  strcpy(output_path,p.path);

  if (input.n_domains>1) {
    stepper.Q = 0;
    bob.message( "parallel processing -> no fft" );
  }
  else {
    if ( (nach=modf(log(n)/log(2),&vor)) > TINY ) {
      stepper.Q   = 0;
      bob.message( "# cells not a power of 2 -> no fft" );
    }
    else {
      stepper.Q   = 1;

      ex   = new( double [grid->n_cells] );
      rhok = new( double [2*grid->n_cells + 1] );
      phik = new( double [2*grid->n_cells + 1] );

      name  = new( char [filename_size] );
    }
  }

  if( input.Q_restart == 1 ) {
    char fname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    stepper.t_count = atoi( rf.getinput( "poi.stepper.t_count" ) );
    rf.closeinput();
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_poisson::input_poisson( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_poisson::Constructor",errname);

  rf.openinput( p.input_file_name );

  n_domains         = p.n_domains;

  stepper.Q         = 1;                // not needed
  stepper.t_start   = atof( rf.setget( "&propagate", "prop_start" ) );
  stepper.t_stop    = atof( rf.setget( "&propagate", "prop_stop" ) );
  stepper.t_step    = 1;                // once per cycle
  stepper.x_start   = -1;               // not needed
  stepper.x_stop    = -1;               // not needed
  stepper.x_step    = -1;               // not needed

  Q_restart       = atoi( rf.setget( "&restart", "Q"     ) );
  strcpy( restart_file, rf.setget( "&restart", "file"    ) );

  rf.closeinput();

  bob.message("parameter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_poisson::save( parameter &p )
{
  static error_handler bob("input_poisson::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic poisson" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Q                : " << stepper.Q       << endl;
  outfile << "t_start          : " << stepper.t_start << endl;
  outfile << "t_stop           : " << stepper.t_stop  << endl;
  outfile << "t_step           : " << stepper.t_step  << endl;
  outfile << "Q_restart        : " << Q_restart       << endl;
  outfile << "restart_file     : " << restart_file    << endl << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void poisson::write( double time, domain *grid )
{
  static error_handler bob("diagnostic::write_poisson",errname);

  struct cell *cell;
  int i;

  sprintf(name,"%s/poisson-%d-%.3f", output_path, domain_number, time);

  file.open(name);
  if (!file) bob.error("cannot open snapshot file", name );

  file.precision( 3 );
  file.setf( ios::showpoint | ios::scientific );

  file << "#"  << setw(11) << "x"
               << setw(12) << "Ex-Current"
               << setw(12) << "Ex-Poisson" << endl;

  for( i=0, cell=grid->left; cell!=grid->rbuf; cell=cell->next, i++ )
    {
      file << setw(12) << cell->x
	           << setw(12) << cell->ex
		   << setw(12) << ex[i] << endl;
    }

  file.close();
}


//////////////////////////////////////////////////////////////////////////////////////////


void poisson::solve( domain* grid )
// n_cells = power of 2 is assumed
// solves Poisson equation to obtain electric field Ex'
// compares Ex' with Ex obtained from Jx during simulation
// can be used to initialize Ex from an initial charge distribution
{
  static error_handler bob("diagnostic::solve_poisson",errname);

  int     n  = grid->n_cells;
  double  dx = grid->dx;                       // grid constant
  double  dk = 2*PI/(dx*n);                    // grid constant
  double  kn, Kn;                              // slitfunction
  double  b;                                   // homogeneous solution
  int i;
  struct cell *cell;

  // fourier transform of the charge density

  for( i=1, cell=grid->left; cell!=grid->rbuf; cell=cell->next, i++ )
    {
      rhok[2*i-1] = cell->charge;                   // real part
      rhok[2*i]   = 0;                              // imaginary part
    }
  fft(rhok,n,1);                                    // transform

  for( i=1; i<n; i++ )                                      // potential in k-space
    {                                                       // smooth at large k
      if (i<=0.5*n) kn=dk*i;
      else          kn=dk*(i-n);

      Kn=kn*dif(kn*dx/2);                                   // LOCAL differences
      // new units: 1/eps -> (2pi)^2 in Poisson's equation!
      phik[2*i+1] = sqr(2*PI/Kn) * rhok[2*i+1] * smooth(kn*dx/2);
      phik[2*i+2] = sqr(2*PI/Kn) * rhok[2*i+2] * smooth(kn*dx/2);
    }
  phik[1] = phik[2] = 0;                    // FT_phi(k=0):=0

  fft(phik,n,-1);

  // take real part and add homogeneous solution --> potential
  // calculate electric field and electrostatic energy density

  // new units: E = - div phi --> E = - 1/2pi div phi

  b = ( phik[1] - phik[3] ) / (2*PI*dx*n);
  ex[0] = - (double) ( phik[3] - phik[1] ) / (2*PI*dx*n) - b;

  for(i=1;i<=n-2;i++)
    ex[i] = - ( phik[2*i+3] - phik[2*i-1] ) / (4*PI*dx*n) - b;

  ex[n-1] = - ( phik[2*n-3] - phik[2*n-5] ) / (2*PI*dx*n) - b;
}


//////////////////////////////////////////////////////////////////////////////////////////


double poisson::dif( double in )
{
  double out;

  if (in==0) out=1;
  else out=sin(in)/in;
  return out;
}


//////////////////////////////////////////////////////////////////////////////////////////


double poisson::smooth( double in )
/*
  cutoff at short wavelengths : a2
  corrects omega_p(k) to order (kdx)^4 : a1
*/
{
  double out;
  double a1=0.5;
  double a2=0.1; /* vorher a2=0.01 */

  if (in==0) out=1;
  else
    {
      out = a1 * pow( sin(in),2 ) - a2 * pow( sin(in)/cos(in), 4 );
      out = exp( 2 * out );
    }
  return out;
}


//////////////////////////////////////////////////////////////////////////////////////////


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void poisson::fft( double* data, int nn, int isign)
// numerical recipies routine "four1.c"
// changed data[], tempr, tempi from float to double
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
//eof
