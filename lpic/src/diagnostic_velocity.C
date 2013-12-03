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

#include <diagnostic_velocity.h>


//////////////////////////////////////////////////////////////////////////////////////////

el_velocity::el_velocity( parameter &p,
			  int species_input, char *species_name_input )
  : velocity(p,species_input,species_name_input)
{
}

ion_velocity::ion_velocity( parameter &p,
			    int species_input, char *species_name_input )
  : velocity(p,species_input,species_name_input)
{
}

//////////////////////////////////////////////////////////////////////////////////////////


velocity::velocity( parameter &p,
		    int species_input, char *species_name_input )
  : rf(),
    input(p,species_name_input),
    stepper( input.stepper, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("velocity::Constructor",errname);

  species = species_input;
  sprintf( species_name, "%s", species_name_input );

  dim        = 399;       // 400 velocity bins: 0...399

  x          = new int [ dim ];
  y          = new int [ dim ];
  z          = new int [ dim ];
  a          = new int [ dim ];

  Beta       = p.Beta;
  Gamma      = p.Gamma;

  vcut       = 1.0;

  name  = new( char [filename_size] );

  if( input.Q_restart == 1 ){
    char fname[ filename_size ];
    char varname[filename_size];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    sprintf( varname, "vel_%s.stepper.t_count", species_name );
    stepper.t_count = atoi( rf.getinput( varname ) );
    rf.closeinput();
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_velocity::input_velocity( parameter &p, char *species_name_input )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_velocity::Constructor",errname);
  char input_name[filename_size];

  sprintf( species_name, "%s", species_name_input );

  rf.openinput( p.input_file_name );

  sprintf( input_name, "&%s_velocity", species_name );
  stepper.Q         = atoi( rf.setget( input_name, "Q" ) );
  stepper.t_start   = atof( rf.setget( input_name, "t_start" ) );
  stepper.t_stop    = atof( rf.setget( input_name, "t_stop" ) );
  stepper.t_step    = atof( rf.setget( input_name, "t_step" ) );

  stepper.x_start   = 0;
  stepper.x_stop    = atoi( rf.setget( "&box", "cells" ) );
  stepper.x_step    = -1; // not used

  Q_restart       = atoi( rf.setget( "&restart", "Q"     ) );
  strcpy( restart_file, rf.setget( "&restart", "file"    ) );

  rf.closeinput();

  bob.message("parameter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_velocity::save( parameter &p )
{
  static error_handler bob("input_velocity::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic " << species_name << "_velocity" << endl;
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


void velocity::write_velocity( double time, parameter &p, domain *grid )
{
  static error_handler bob("velocity::write_velocity",errname);

  struct cell *cell;
  struct particle *part;
  int i;
  double vx, vy, vz, v, absolut;
  int bvx, bvy, bvz, bv;
  FILE *file;

  for( i=0; i<=dim; i++ )
    x[i]=y[i]=z[i]=a[i]=0;

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {

      if (cell->npart != 0 && cell->number >= stepper.x_start
                           && cell->number < stepper.x_stop) {

	for( part=cell->first; part!=NULL; part=part->next ) {

	    if (part->species == species) {
	      vx = part->ux * part->igamma;
	      vy = part->uy * part->igamma;
	      vz = part->uz * part->igamma;

	      vx = 1.0/Gamma * vx / ( 1 + vy * Beta );
	      vz = 1.0/Gamma * vz / ( 1 + vy * Beta );
	      vy = ( vy + Beta ) / ( 1 + vy * Beta );
	      absolut = sqrt( sqr(vx) + sqr(vy) + sqr(vz) );

	      bvx = (int) floor( 0.5 * dim * (1.0 + vx/vcut) + 0.5 );
	      bvy = (int) floor( 0.5 * dim * (1.0 + vy/vcut) + 0.5 );
	      bvz = (int) floor( 0.5 * dim * (1.0 + vz/vcut) + 0.5 );
	      bv  = (int) floor( 0.5 * dim * (1.0 + absolut/vcut) + 0.5 );

	      if (bvx>=0 && bvx<=dim) x[bvx]++;
	      else bob.error( "velocity bin out of range" );
	      if (bvy>=0 && bvy<=dim) y[bvy]++;
	      else bob.error( "velocity bin out of range" );
	      if (bvz>=0 && bvz<=dim) z[bvz]++;
	      else bob.error( "velocity bin out of range" );
	      if (bv>=0 && bv<=dim)   a[bv]++;
	      else bob.error( "velocity bin out of range" );
	    }
	  }
      }
    }

  sprintf(name,"%s/velocity-%d-sp%d-%.3f", p.path, p.domain_number, species, time);
  file = fopen( name, "w" );

  for( i=0; i<=dim; i++ ) {
    v = -vcut + 2.0*vcut*i/dim;
    fprintf( file, "\n %.4e  %d %d %d %d", v, x[i], y[i], z[i], a[i] );
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof







