/*
   This file is part of LPIC++, a particle-in-cell code for
   simulating the interaction of laser light with plasma.

   Copyright (C) 2002      Andreas Kemp
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

#include <domain.h>

using namespace std;

// exponential energy distribution introduced by A.Kemp, 04/02
domain::domain( parameter &p )
  : input(p)
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("domain::Constructor",errname);

  domain_number = p.domain_number;
  n_domains     = input.n_domains;
  dx            = input.dx;

  strcpy( path, p.path );

  n_el    = 0;                           // will be set in domain::chain_particles()
  n_ion   = 0;                           //  ''
  n_part  = 0;                           //  ''

  if(input.Q_restart == 0){
    // simulation box -----------------------

    set_boundaries();                      // determines boundaries from domain number

    // init data structure-------------------

    chain_cells();                         // create chained list of cells
                                           // set cell numbers
                                           // set domain pointers to left and right cells

    init_cells();                          // set fields equal to zero
                                           // set normalized densities

    chain_particles();                     // allocate and link particles to cells
                                           // set number, species, cell, position

  // physical part ------------------------

    init_particles();                      // set charge, mass, velocities, gamma

    check();                               // check and save
                                           // particel numbers, positions, total charge
                                           // plot cell, x, densities, part. per cell
  }
  else {
    restart_configuration();               // set restart configuration
  }

  bob.message( "sizeof(struct cell)     =", sizeof(struct cell), "Byte" );
  bob.message( "sizeof(struct particle) =", sizeof(struct particle), "Byte" );
}


//////////////////////////////////////////////////////////////////////////////////////////


input_domain::input_domain( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_domain::Constructor",errname);

  rf.openinput(p.input_file_name);

  Q_restart      = atoi( rf.setget( "&restart", "Q" ) );
  strcpy( restart_file, rf.setget( "&restart", "file" ) );
  Q_restart_save = atoi( rf.setget( "&restart", "Q_save" ) );

  // read grid ////////////////////////////////////////////////////

  n_domains    = p.n_domains;
  cells        = atoi( rf.setget( "&box", "cells" ) );
  cells_per_wl = atoi( rf.setget( "&box", "cells_per_wl" ) );
  cells_left   = atoi( rf.setget( "&box", "cells_left" ) );
  cells_ramp   = atoi( rf.setget( "&box", "cells_ramp" ) );
  cells_plasma = atoi( rf.setget( "&box", "cells_plasma" ) );
  n_ion_over_nc= atof( rf.setget( "&box", "n_ion_over_nc" ) );

  dx           = 1.0 / cells_per_wl;

  // read particles /////////////////////////////////////////////////

  nsp                 = p.nsp;

  fix                 = new(int[nsp]);
  z                   = new(double[nsp]);
  m                   = new(double[nsp]);
  ppc                 = new(int[nsp]);
  vtherm              = new(double[nsp]);

  // electrons -----------------------------------------------------

  fix[0]         = atoi( rf.setget( "&electrons", "fix" ) );
  z[0]           = -1;        // DEFAULT, SHOULD NOT BE CHANGED
  m[0]           = +1;        // DEFAULT, SHOULD NOT BE CHANGED
  ppc[0]         = atoi( rf.setget( "&electrons", "ppc" ) );
  vtherm[0]      = atof( rf.setget( "&electrons", "vtherm" ) );

  // ions ----------------------------------------------------------

  fix[1]         = atoi( rf.setget( "&ions", "fix" ) );
  z[1]           = atoi( rf.setget( "&ions", "z" ) );
  m[1]           = atof( rf.setget( "&ions", "m" ) );
  ppc[1]         = atoi( rf.setget( "&ions", "ppc" ) );
  vtherm[1]      = atof( rf.setget( "&ions", "vtherm" ) );

  if ( z[1] == 0 && ppc[0] > 0 ) {
    cout << " WARNING     : You selected neutral atoms and free electrons" << endl;
    cout << "               Setting # MacroElectrons := 0" << endl;
    ppc[0] = 0;
  }

  n_el_over_nc = 1.0 * z[1] * n_ion_over_nc;

  // read spp, adjusted angle, Beta and Gamma

  spp         = p.spp;
  angle       = p.angle;
  Beta        = p.Beta;
  Gamma       = p.Gamma;

  rf.closeinput();

  bob.message("parameter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_domain::save( parameter &p )
{
  static error_handler bob("input_domain::save",errname);
  ofstream outfile;
  int i;

  outfile.open(p.outname,ios::app);

  outfile << "domain: plasma" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Q_restart          : " << Q_restart      << endl;
  outfile << "restart_file       : " << restart_file   << endl;
  outfile << "Q_restart_save     : " << Q_restart_save << endl;
  outfile << "N_domains          : " << n_domains      << endl;
  outfile << "cells              : " << cells          << endl;
  outfile << "cells_per_wl       : " << cells_per_wl   << endl;
  outfile << "cells_left         : " << cells_left     << endl;
  outfile << "cells_ramp         : " << cells_ramp     << endl;
  outfile << "cells_plasma       : " << cells_plasma   << endl;
  outfile << "dx                 : " << dx             << endl << endl;
  outfile << "n_ion_over_nc      : " << n_ion_over_nc  << endl;
  outfile << "n_el_over_nc       : " << n_el_over_nc   << endl << endl << endl;

  outfile << "domain: particles" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "nsp                : " << setw(8) << nsp << endl;
  outfile << "fix                : ";
  for(i=0;i<nsp;i++) outfile << setw(8) << fix[i];
  outfile << endl << "charge             : ";
  for(i=0;i<nsp;i++) outfile << setw(8) << z[i];
  outfile << endl << "mass               : ";
  for(i=0;i<nsp;i++) outfile << setw(8) << m[i];
  outfile << endl << "ppc                : ";
  for(i=0;i<nsp;i++) outfile << setw(8) << ppc[i];
  outfile << endl << "vtherm             : ";
  for(i=0;i<nsp;i++) outfile << setw(8) << vtherm[i];
  outfile << endl << endl << endl;

  outfile << "domain: Lorentz transformation" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "spp                : " << spp             << endl;
  outfile << "adjusted angle     : " << angle           << endl;
  outfile << "Beta               : " << Beta            << endl;
  outfile << "Gamma              : " << Gamma           << endl << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void domain::set_boundaries()
// p.box.N domains of equal size
{
  error_handler bob("domain::set_boundaries",errname);

  int cells_per_domain = (int) floor( (double) input.cells / input.n_domains );

  n_left  = 1 + ( domain_number - 1 ) * cells_per_domain;
  n_right = n_left + cells_per_domain - 1;
  n_cells = cells_per_domain;

  if (domain_number==n_domains)
    {
      n_right = input.cells;
      n_cells = n_right - n_left + 1;
    }

  bob.message("n_left  = ", n_left );
  bob.message("n_right = ", n_right );
  bob.message("n_cells = ", n_cells );
}


//////////////////////////////////////////////////////////////////////////////////////////


void domain::chain_cells( void )
{
  error_handler bob("domain::chain_cells",errname);
  struct cell *cell_old, *cell_new;

  Lbuf         = new ( struct cell );
  if (!Lbuf) bob.error("allocation error: Lbuf");
  Lbuf->number = n_left - 2;

  lbuf         = new ( struct cell );
  if (!lbuf) bob.error("allocation error: lbuf");
  lbuf->number = n_left - 1;
  lbuf->prev   = Lbuf;

  left         = new ( struct cell );
  if (!left) bob.error("allocation error: left");
  left->number = n_left;
  left->prev   = lbuf;
  cell_old     = left;

  for( int i=n_left+1; i<=n_right; i++ )
    {
      cell_new         = new ( struct cell );
      if (!cell_new) bob.error("allocation error: cell_new");
      cell_new->prev   = cell_old;
      cell_old->next   = cell_new;
      cell_old         = cell_new;
      cell_new->number = i;
    }

  right         = cell_old;

  rbuf         = new ( struct cell );
  if (!rbuf) bob.error("allocation error: rbuf");
  rbuf->prev   = right;
  rbuf->number = n_right + 1;

  Rbuf         = new ( struct cell );
  if (!Rbuf) bob.error("allocation error: Rbuf");
  Rbuf->prev   = rbuf;
  Rbuf->number = n_right + 2;

  dummy         = new ( struct cell );
  if (!dummy) bob.error("allocation error: dummy");
  dummy->prev   = Rbuf;
  dummy->number = n_right + 3;

  right->next  = rbuf;
  rbuf->next   = Rbuf;
  Rbuf->next   = dummy;
  dummy->next  = Lbuf;
  Lbuf->next   = lbuf;
  Lbuf->prev   = dummy;
  lbuf->next   = left;

  // closed ring of cells with two buffer cells left and right and dummy cell
  // connecting Rbuf and Lbuf
}


//////////////////////////////////////////////////////////////////////////////////////////


void domain::init_cells( void )
{
  error_handler bob("domain::init_cells",errname);
  struct cell *cell;

  for( cell=Lbuf; cell!=dummy; cell=cell->next )
    {
      if (!cell) bob.error("allocation error");

      cell->domain = domain_number;            // domain number
      cell->x = dx * ( cell->number - 1 );     // cell coordinate, left boundary

      // set up the normalized particle densities for a ramp profile

      if ( cell->number < input.cells_left ) cell->dens[0] = cell->dens[1] = 0;
      else
	{
	  if (cell->number < input.cells_left + input.cells_ramp)
	    {
	      cell->dens[0] = (double)(cell->number - input.cells_left) /
                               input.cells_ramp;
	      cell->dens[1] = cell->dens[0];
	    }
	  else
	    {
	      if (cell->number < input.cells_left + input.cells_plasma )
		{
		  cell->dens[0] = cell->dens[1] = 1.0;
		}
	      else cell->dens[0] = cell->dens[1] = 0;
	    }
	}

      if (cell->number < n_left || cell->number > n_right)  // initially empty buffers !
	cell->dens[0] = cell->dens[1] = 0;

      cell->charge = 0;
      cell->jx = cell->jy = cell->jz = 0;
      cell->ex = cell->ey = cell->ez = 0;
      cell->bx = cell->by = cell->bz = 0;
      cell->fp = cell->fm = cell->gp = cell->gm = 0;

      for( int j=0;j<input.nsp; j++ ) cell->np[j] = 0;
      cell->npart = 0;
      // will be set in domain::chain_particles()
    }
}


//////////////////////////////////////////////////////////////////////////////////////////


void domain::chain_particles( void )
  // particle numbers are set in each domain separately, beginning with 1
  // this is corrected in the box::Constructor using communicate_particle_numbers(...)
  // in order to keep domain free of networking
{
  error_handler bob("domain::chain_particles",errname);

  int             i;
  int             *number;         // count particles for each species sperately
  double          delta;
  struct cell     *cell;
  struct particle *pn, *po;

  number = new( int [input.nsp] );
  if (!number) bob.error("allocation error");

  for ( i=0; i<input.nsp; i++ ) number[i] = 0;

  for( cell=Lbuf; cell!=dummy; cell=cell->next )   // for all cells including all buffers
    {
      cell->first = NULL;
      cell->last  = NULL;

      po = cell->first;

      for( int j=0; j<input.nsp; j++ )             // for all species
	{
	  cell->np[j] = (int) floor( cell->dens[j] * input.ppc[j] + 0.5 );

	  if (j==0)  n_el += cell->np[j];          // count particles by species
	  else      n_ion += cell->np[j];
	  cell->npart     += cell->np[j];          // particles per cell
	  n_part          += cell->np[j];          // particles per domain

	  if (cell->np[j]!=0)                      // for occupied cells
	    {
	      delta  = dx / cell->np[j];

	      for( i=1; i<=cell->np[j]; i++ )      // for all particles of
		{                                  // kind j in this cell
		  number[j]++;

		  pn          = new ( struct particle );
		  if (!pn) bob.error("allocation error");

		  pn->prev    = po;
		  if (po==NULL) cell->first = pn;
		  else          po->next    = pn;
		  pn->number  = number[j];
		  pn->species = j;
		  pn->cell    = cell;
		  pn->x       = cell->x + ((double)i-0.50000001) * delta;

		  po          = pn;
		}
	    }
	}
      if (po!=NULL) {
	po->next = NULL;
	cell->last = po;
      }
    }

  if ( n_el != number[0] ) bob.error("# allocated electrons incorrect");

  for( i=2; i<input.nsp; i++ ) number[1]+=number[i];
  if ( n_ion != number[1] ) bob.error("# allocated ions incorrect");

  delete number;
}

//////////////////////////////////////////////////////////////////////////////////////////


void domain::init_particles( void )
{
  error_handler bob("domain::init_particles",errname);

  struct cell     *cell;
  struct particle *part;
  double Gamma   = input.Gamma;                   // gamma factor due to Lorentz transformation
  double Beta    = input.Beta;
  double vx, vy, vz;

  for( cell=left; cell!=rbuf; cell=cell->next )   // for all cells
    {
      if (cell->npart!=0)                         // for occupied cells
	{
	  for( part=cell->first; part!=NULL; part=part->next )
	    {                                      // for all particles in this cell
	      if (!part) bob.error("allocation: part");

	      part->fix     = input.fix[part->species];
	      part->z       = input.z[part->species];
	      part->m       = input.m[part->species];
	      part->zm      = part->z / part->m;

	      // part->n  = density / critical density
	      // part->zn = charge state * density / critical density
	      // ---------------------------------------------------------
	      // Lorentz-Transformation: part->zn is scaled up with Gamma^3!
	      // L-Contraction in y-direction leads to n_M = Gamma n_L
	      // and Doppler shift leads to n_c_M = 1/Gamma^2 * n_c_L
	      // ---------------------------------------------------------
	      // part->zn is designed such that the sum of MacroParticle charges
	      // (electrons and ions) is zero in each cell initially


	      if (part->z == 0) {     // neutral atoms
		part->n  = pow(Gamma,3) * input.n_ion_over_nc / input.ppc[part->species];
		part->zn = 0;
	      }
	      else {                  // electrons or ions
		part->n  = pow(Gamma,3) * fabs( 1.0 * input.z[1] / part->z )
		           * input.n_ion_over_nc / input.ppc[part->species];
		part->zn = part->z * part->n;
	      }
	                                                    // thermal velocities
	      do
		{
		  vx            = input.vtherm[part->species] * gauss_rand48();
		  vy            = input.vtherm[part->species] * gauss_rand48();
		  vz            = input.vtherm[part->species] * gauss_rand48();
		  //		  vx = exponential_rand( input.vtherm[part->species] ); vy = vz = 0.0;
		}
	      while( vx*vx + vy*vy + vz*vz >= 1.0);         // make sure that |v| < c

	                                                    // L-transform to the M frame

	      vx            = vx * sqrt(1.0-Beta*Beta) / ( 1 - vy*Beta );
	      vz            = vz * sqrt(1.0-Beta*Beta) / ( 1 - vy*Beta );
	      vy            = ( vy - Beta ) / ( 1 - vy*Beta );

                                                            // determine gamma*v

	      part->igamma  = sqrt( 1.0 - vx*vx - vy*vy - vz*vz );
	      part->ux      = vx / part->igamma;
	      part->uy      = vy / part->igamma;
	      part->uz      = vz / part->igamma;
	    }
	}
    }
}


//////////////////////////////////////////////////////////////////////////////////////////


double domain::exponential_rand( double tm )
{
  // one dimensional exponential energy distribution ##
  // tm == kT / Me
  static error_handler bob( "domain::exponential_rand", errname );
  double r1;

  r1 = drand48();
  return sqrt( 1.0 - 1.0/ sqr( (1.0 - tm * log( 1.0 - r1)) ));
}
///////////////////////////////////////////////////////////////////////////

double domain::gauss_rand48( void )
{
  double r1, r2;

  r1 = drand48();
  r2 = drand48();

  return sqrt( -2.0 * log( 1.0 - r1 ) ) * sin( 2*PI*r2 );
}


//////////////////////////////////////////////////////////////////////////////////////////


void domain::check( void )
{
  error_handler bob("domain::check_and_save",errname);
  struct cell *cell;
  struct particle *part;
  int count[2];
  double charge=0;

  for( cell=Lbuf; cell!=dummy; cell=cell->next )
    {
      count[0] = count[1] = 0;

      if (cell->npart!=0)
	{
	  for( part=cell->first; part!=NULL; part=part->next )
	    {
	      if (!part) bob.error("allocation: part");
	      if ( part->x < cell->x || part->x > cell->x+dx )
		bob.error("particle position");
	      count[part->species]++;
	      charge += part->zn;
	    }
	}
      if (cell->np[0] != count[0]) bob.error("number of electrons");
      if (cell->np[1] != count[1]) bob.error("number of ions");
    }

  if ( fabs(charge)> 1e-6 ) bob.error("domain is charged");

  char fname[filename_size];
  sprintf(fname,"%s/domain-%d", path, domain_number);

  ofstream domain_file(fname);
  if (!domain_file) bob.error("cannot open output file: ", fname );

  domain_file.precision( 3 );
  domain_file.setf( ios::showpoint );
  domain_file.setf( ios::scientific );

  domain_file << setw(6)  << "# cell"
	      << setw(12) << "x"
	      << setw(12) << "rho_el."
	      << setw(12) << "rho_ion"
	      << setw(7)  << "n_el"
	      << setw(7)  << "n_ion" << endl;

  for( cell=Lbuf; cell!=dummy; cell=cell->next )
    {
      domain_file << setw(6)  << cell->number
		  << setw(12) << cell->x
		  << setw(12) << cell->dens[0]
		  << setw(12) << cell->dens[1]
		  << setw(7)  << cell->np[0]
		  << setw(7)  << cell->np[1] << endl;
    }

  domain_file.close();
}


//////////////////////////////////////////////////////////////////////////////////////////


void domain::count_particles()
  // counts particles in domain -> n_el, n_ion, n_part
{
  static error_handler bob("domain::count_particles",errname);

  struct cell *cell;
  struct particle *part;

  int nparts_1, nparts_2;
  int n_el_old   = n_el;       // keep in mind the old numbers
  int n_ion_old  = n_ion;
  int n_part_old = n_part;

  n_el   = 0;
  n_ion  = 0;
  n_part = 0;

  for( cell=Lbuf; cell!=dummy; cell=cell->next )
    {
      nparts_1 = nparts_2 = 0;

      if (cell->npart != 0) {

	nparts_1 += cell->npart;

	for( part=cell->first; part!=NULL; part=part->next )
	  {
	    nparts_2++;

	    if (part->species==0) n_el  ++;
	    else                  n_ion ++;
	    n_part ++;
	  }

	if (nparts_1 != nparts_2) {
	  bob.message( "particle numbers incorrect" );
	  bob.message( "                 in cell:", cell->number );
	  bob.message( "             cell_parts =", nparts_1 );
	  bob.message( "                  parts =", nparts_2 );

	  bob.error("");
	}

      }
    }

  bob.message( "particle numbers", n_el, n_ion, n_part );

  if (n_el!=n_el_old || n_ion!=n_ion_old || n_part!=n_part_old)
    {
      bob.message( "old particle numbers", n_el_old, n_ion_old, n_part_old );
      bob.error("");
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

void domain::reo_to_prev( int request_prev, int *cells_to_prev, int *parts_to_prev )
{
  static error_handler bob("domain::reo_to_prev",errname);

  struct cell *cell;

  cell = left;
  *cells_to_prev = 1;
  *parts_to_prev = cell->npart;

  while( *parts_to_prev < -request_prev && cell != right->prev ) // to make sure that at
    {                                           // at least one cell stays in the domain
      cell = cell->next;
      *cells_to_prev = *cells_to_prev + 1;
      *parts_to_prev = *parts_to_prev + cell->npart;
    }

  if ( ( -request_prev - (*parts_to_prev - cell->npart) ) / cell->npart < 0.5 )
    {
      *cells_to_prev = *cells_to_prev - 1;
      *parts_to_prev = *parts_to_prev - cell->npart;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

void domain::reo_to_next( int request_next, int *cells_to_next, int *parts_to_next )
{
  static error_handler bob("domain::reo_to_next",errname);

  struct cell *cell;

  cell = right;
  *cells_to_next = 1;
  *parts_to_next = cell->npart;

  while( *parts_to_next < request_next && cell != left->next ) // to make sure that at
    {                                         // at least one cell stays in the domain
      cell = cell->prev;
      *cells_to_next = *cells_to_next + 1;
      *parts_to_next = *parts_to_next + cell->npart;
    }

   if ( ( request_next - (*parts_to_next - cell->npart) ) / cell->npart < 0.5 )
    {
      *cells_to_next = *cells_to_next - 1;
      *parts_to_next = *parts_to_next - cell->npart;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

void domain::reo_delete_to_prev( int cells_to_prev, int parts_to_prev )
{
  static error_handler bob("domain::reo_delete_to_prev",errname);

  int i;
  int partcount = 0;
  int el_count  = 0;
  int ion_count = 0;
  struct cell *cell,*dcell;
  struct particle *part,*dpart;

  cell = left;

  for(i=0;i<cells_to_prev;i++) {

    if(i == cells_to_prev - 2){
      Lbuf->number = cell->number;
      Lbuf->x      = cell->x;
      Lbuf->charge = cell->charge;
      Lbuf->jx     = cell->jx;
      Lbuf->jy     = cell->jy;
      Lbuf->jz     = cell->jz;
      Lbuf->ex     = cell->ex;
      Lbuf->ey     = cell->ey;
      Lbuf->ez     = cell->ez;
      Lbuf->bx     = cell->bx;
      Lbuf->by     = cell->by;
      Lbuf->bz     = cell->bz;
      Lbuf->fp     = cell->fp;
      Lbuf->fm     = cell->fm;
      Lbuf->gp     = cell->gp;
      Lbuf->gm     = cell->gm;
      Lbuf->dens[0]= cell->dens[0];
      Lbuf->dens[1]= cell->dens[1];
      Lbuf->np[0]  = 0;
      Lbuf->np[1]  = 0;
      Lbuf->npart  = 0;
    }

    if(i == cells_to_prev - 1){
      lbuf->number = cell->number;
      lbuf->x      = cell->x;
      lbuf->charge = cell->charge;
      lbuf->jx     = cell->jx;
      lbuf->jy     = cell->jy;
      lbuf->jz     = cell->jz;
      lbuf->ex     = cell->ex;
      lbuf->ey     = cell->ey;
      lbuf->ez     = cell->ez;
      lbuf->bx     = cell->bx;
      lbuf->by     = cell->by;
      lbuf->bz     = cell->bz;
      lbuf->fp     = cell->fp;
      lbuf->fm     = cell->fm;
      lbuf->gp     = cell->gp;
      lbuf->gm     = cell->gm;
      lbuf->dens[0]= cell->dens[0];
      lbuf->dens[1]= cell->dens[1];
      lbuf->np[0]  = 0;
      lbuf->np[1]  = 0;
      lbuf->npart  = 0;
    }

       part      = cell->first;
       while(part != NULL)
        {
	 switch (part->species){
	 case 0:
          el_count ++;
          partcount ++;
          break;
	 case 1:
          ion_count ++;
          partcount ++;
	  break;
         }
         dpart     = part;
         part      = part->next;
         delete dpart;
        }

    dcell = cell;
    cell  = cell->next;
    delete dcell;
  }

  cell->prev = lbuf;
  lbuf->next = cell;

  left       = cell;
  n_left    += cells_to_prev;
  n_cells   -= cells_to_prev;
  n_el      -= el_count;
  n_ion     -= ion_count;
  n_part    -= partcount;

  if (partcount != parts_to_prev) {
   bob.error( "number of particles deleted does NOT match intended number to delete" );}

}

//////////////////////////////////////////////////////////////////////////////////////////

void domain::reo_delete_to_next( int cells_to_next, int parts_to_next )
{
  static error_handler bob("domain::reo_delete_to_next",errname);

  int i;
  int partcount = 0;
  int el_count  = 0;
  int ion_count = 0;
  struct cell *cell,*dcell;
  struct particle *part,*dpart;

  cell = right;

  for(i=0;i<cells_to_next;i++) {

    if(i == cells_to_next - 2){
      Rbuf->number = cell->number;
      Rbuf->x      = cell->x;
      Rbuf->charge = cell->charge;
      Rbuf->jx     = cell->jx;
      Rbuf->jy     = cell->jy;
      Rbuf->jz     = cell->jz;
      Rbuf->ex     = cell->ex;
      Rbuf->ey     = cell->ey;
      Rbuf->ez     = cell->ez;
      Rbuf->bx     = cell->bx;
      Rbuf->by     = cell->by;
      Rbuf->bz     = cell->bz;
      Rbuf->fp     = cell->fp;
      Rbuf->fm     = cell->fm;
      Rbuf->gp     = cell->gp;
      Rbuf->gm     = cell->gm;
      Rbuf->dens[0]= cell->dens[0];
      Rbuf->dens[1]= cell->dens[1];
      Rbuf->np[0]  = 0;
      Rbuf->np[1]  = 0;
      Rbuf->npart  = 0;
    }
    if(i == cells_to_next - 1){
      rbuf->number = cell->number;
      rbuf->x      = cell->x;
      rbuf->charge = cell->charge;
      rbuf->jx     = cell->jx;
      rbuf->jy     = cell->jy;
      rbuf->jz     = cell->jz;
      rbuf->ex     = cell->ex;
      rbuf->ey     = cell->ey;
      rbuf->ez     = cell->ez;
      rbuf->bx     = cell->bx;
      rbuf->by     = cell->by;
      rbuf->bz     = cell->bz;
      rbuf->fp     = cell->fp;
      rbuf->fm     = cell->fm;
      rbuf->gp     = cell->gp;
      rbuf->gm     = cell->gm;
      rbuf->dens[0]= cell->dens[0];
      rbuf->dens[1]= cell->dens[1];
      rbuf->np[0]  = 0;
      rbuf->np[1]  = 0;
      rbuf->npart  = 0;
    }

       part      = cell->first;
       while(part != NULL)
        {
	 switch (part->species){
	 case 0:
          el_count ++;
          partcount ++;
          break;
	 case 1:
          ion_count ++;
          partcount ++;
	  break;
         }
         dpart     = part;
         part      = part->next;
         delete dpart;
        }

    dcell = cell;
    cell  = cell->prev;
    delete dcell;
  }

  cell->next = rbuf;
  rbuf->prev = cell;

  right      = cell;
  n_right   -= cells_to_next;
  n_cells   -= cells_to_next;
  n_el      -= el_count;
  n_ion     -= ion_count;
  n_part    -= partcount;

  if (partcount != parts_to_next) {
   bob.error( "number of particles deleted does NOT match intended number to delete" );}

}

//////////////////////////////////////////////////////////////////////////////////////////

void domain::reo_alloc_from_prev( int cells_from_prev, int parts_from_prev )
{
  static error_handler bob("domain::reo_alloc_from_prev",errname);

  int i,cell_number;
  double cell_x;
  struct cell *cell_new;
  struct particle *part_new;

  cell_new = left;
  cell_number = left->number;
  cell_x      = left->x;

  for(i=0;i<cells_from_prev;i++) {
      cell_new         = new ( struct cell );
      if (!cell_new) bob.error("allocation error: cell_new");
      cell_new->prev   = lbuf;
      cell_new->next   = lbuf->next;
      lbuf->next->prev = cell_new;
      lbuf->next       = cell_new;
      cell_number --;
      cell_new->number = cell_number;
      cell_x -= dx;
      cell_new->x      = cell_x;

      cell_new->first  = NULL;
      cell_new->last   = NULL;
  }

  left     = cell_new;

  n_left  -= cells_from_prev;
  n_cells += cells_from_prev;
  n_part  += parts_from_prev;

  for(i=0;i<parts_from_prev;i++) {
      part_new = new ( struct particle );
      if (!(part_new)) bob.error("allocation error");
      part_new->prev = NULL;
      part_new->next = left->first;
      if (part_new->next!=NULL) left->first->prev = part_new;
      else                      left->last        = part_new;
      left->first    = part_new;

      part_new->cell = left;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void domain::reo_alloc_from_next( int cells_from_next, int parts_from_next )
{
  static error_handler bob("domain::reo_alloc_from_next",errname);

  int i,cell_number;
  double cell_x;
  struct cell *cell_new;
  struct particle *part_new;

  cell_new = right;
  cell_number = right->number;
  cell_x      = right->x;

  for(i=0;i<cells_from_next;i++) {
      cell_new         = new ( struct cell );
      if (!cell_new) bob.error("allocation error: cell_new");
      cell_new->next   = rbuf;
      cell_new->prev   = rbuf->prev;
      rbuf->prev->next = cell_new;
      rbuf->prev       = cell_new;
      cell_number ++;
      cell_new->number = cell_number;
      cell_x += dx;
      cell_new->x      = cell_x;

      cell_new->first = NULL;
      cell_new->last  = NULL;
  }

  right    = cell_new;

  n_right += cells_from_next;
  n_cells += cells_from_next;
  n_part  += parts_from_next;

  for(i=0;i<parts_from_next;i++) {
      part_new = new ( struct particle );
      if (!(part_new)) bob.error("allocation error");
      part_new->prev = NULL;
      part_new->next = right->first;
      if (part_new->next!=NULL) right->first->prev = part_new;
      else                      right->last        = part_new;
      right->first   = part_new;

      part_new->cell = right;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void domain::reo_update_n_el_n_ion( int el_count, int ion_count )
{
  static error_handler bob("domain::reo_update_n_el_n_ion",errname);

  n_el  += el_count;
  n_ion += ion_count;

}

//////////////////////////////////////////////////////////////////////////////////////////

void domain::restart_configuration( void )
{
  error_handler bob("domain::restart_configuration",errname);

  FILE *file;
  char fname[ filename_size ];
  struct cell *cell_old, *cell_new, *cell;
  struct particle *part;
  int i,k;
  int n_el_check, n_ion_check, n_part_check;

  sprintf( fname, "%s/%s-%d-data2", path, input.restart_file, domain_number );
  file = fopen( fname, "rb" );
  if (!file) bob.error( "cannot open file", fname );

  fread( &n_cells, sizeof(int), 1, file );

  // create chained list of cells:

  Lbuf         = new ( struct cell );
  if (!Lbuf) bob.error("allocation error: Lbuf");
  Lbuf->first  = NULL;
  Lbuf->last   = NULL;

  lbuf         = new ( struct cell );
  if (!lbuf) bob.error("allocation error: lbuf");
  lbuf->prev   = Lbuf;
  lbuf->first  = NULL;
  lbuf->last   = NULL;

  left         = new ( struct cell );
  if (!left) bob.error("allocation error: left");
  left->prev   = lbuf;
  left->first  = NULL;
  left->last   = NULL;

  cell_old     = left;

  for( i=0; i<n_cells-1; i++ )
    {
      cell_new         = new ( struct cell );
      if (!cell_new) bob.error("allocation error: cell_new");
      cell_new->prev   = cell_old;
      cell_new->first  = NULL;
      cell_new->last   = NULL;
      cell_old->next   = cell_new;
      cell_old         = cell_new;
    }

  right         = cell_old;
  right->first  = NULL;
  right->last   = NULL;

  rbuf         = new ( struct cell );
  if (!rbuf) bob.error("allocation error: rbuf");
  rbuf->prev   = right;
  rbuf->first  = NULL;
  rbuf->last   = NULL;

  Rbuf         = new ( struct cell );
  if (!Rbuf) bob.error("allocation error: Rbuf");
  Rbuf->prev   = rbuf;
  Rbuf->first  = NULL;
  Rbuf->last   = NULL;

  dummy         = new ( struct cell );
  if (!dummy) bob.error("allocation error: dummy");
  dummy->prev   = Rbuf;

  right->next  = rbuf;
  rbuf->next   = Rbuf;
  Rbuf->next   = dummy;
  dummy->next  = Lbuf;
  Lbuf->next   = lbuf;
  Lbuf->prev   = dummy;
  lbuf->next   = left;

  // read cell information from restart file
  // including particles:

  for( cell=Lbuf; cell!=dummy; cell=cell->next )
    {
      fread( &cell->number , sizeof(int), 1, file );
      fread( &cell->x      , sizeof(double), 1, file );
      fread( &cell->charge , sizeof(double), 1, file );
      fread( &cell->jx     , sizeof(double), 1, file );
      fread( &cell->jy     , sizeof(double), 1, file );
      fread( &cell->jz     , sizeof(double), 1, file );
      fread( &cell->ex     , sizeof(double), 1, file );
      fread( &cell->ey     , sizeof(double), 1, file );
      fread( &cell->ez     , sizeof(double), 1, file );
      fread( &cell->bx     , sizeof(double), 1, file );
      fread( &cell->by     , sizeof(double), 1, file );
      fread( &cell->bz     , sizeof(double), 1, file );
      fread( &cell->fp     , sizeof(double), 1, file );
      fread( &cell->fm     , sizeof(double), 1, file );
      fread( &cell->gp     , sizeof(double), 1, file );
      fread( &cell->gm     , sizeof(double), 1, file );
      fread( &(cell->dens[0]), sizeof(double), 1, file );
      fread( &(cell->dens[1]), sizeof(double), 1, file );
      fread( &(cell->np[0])  , sizeof(int), 1, file );
      fread( &(cell->np[1])  , sizeof(int), 1, file );
      fread( &cell->npart    , sizeof(int), 1, file );

      cell->domain = domain_number;

      for(k=0;k<cell->npart;k++) {

	part          = new ( struct particle );
	if (!part) bob.error("allocation error");

	fread( &part->number , sizeof(int), 1, file );
	fread( &part->species, sizeof(int), 1, file );
	fread( &part->fix    , sizeof(int), 1, file );
	fread( &part->z      , sizeof(double), 1, file );
	fread( &part->m      , sizeof(double), 1, file );
	fread( &part->zm     , sizeof(double), 1, file );
	fread( &part->x      , sizeof(double), 1, file );
	fread( &part->dx     , sizeof(double), 1, file );
	fread( &part->igamma , sizeof(double), 1, file );
	fread( &part->ux     , sizeof(double), 1, file );
	fread( &part->uy     , sizeof(double), 1, file );
	fread( &part->uz     , sizeof(double), 1, file );
	fread( &part->zn     , sizeof(double), 1, file );

	part->cell = cell;
	part->next = NULL;
	part->prev = cell->last;
	if (part->prev==NULL) cell->first = part;
	if (cell->last!=NULL) cell->last->next = part;
	cell->last = part;

	switch (part->species){
	case 0:
	  n_el   ++;
	  n_part ++;
	  break;
	case 1:
	  n_ion  ++;
	  n_part ++;
	  break;
	}

      }
    }
  n_left  = left->number;
  n_right = right->number;

  fread( &n_el_check, sizeof(int), 1, file );
  fread( &n_ion_check, sizeof(int), 1, file );
  fread( &n_part_check, sizeof(int), 1, file );
  fclose( file );

  bob.message("n_cells = ", n_cells);
  bob.message("n_el    = ", n_el);
  bob.message("n_ion   = ", n_ion);
  bob.message("n_part  = ", n_part);

  if( n_el_check!=n_el){
    bob.error("n_el incorrect", n_el);
  }
  if( n_ion_check!=n_ion){
    bob.error("n_ion incorrect");
  }
  if( n_part_check!=n_part){
    bob.error("n_part incorrect");
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
//eof
