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

#include <diagnostic_phasespace.h>


//////////////////////////////////////////////////////////////////////////////////////////

el_phasespace::el_phasespace( parameter &p,
			      int species_input, char *species_name_input )
  : phasespace(p,species_input,species_name_input)
{
}

ion_phasespace::ion_phasespace( parameter &p,
				int species_input, char *species_name_input )
  : phasespace(p,species_input,species_name_input)
{
}

//////////////////////////////////////////////////////////////////////////////////////////


phasespace::phasespace( parameter &p,
			int species_input, char *species_name_input )
  : rf(),
    input(p,species_name_input),
    stepper( input.stepper, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("phasespace::Constructor",errname);

  species = species_input;
  sprintf( species_name, "%s", species_name_input );

  dim        = 399;       // 400 velocity bins: 0...399

  x          = ucmatrix( 0, dim, 0, dim );
  y          = ucmatrix( 0, dim, 0, dim );
  z          = ucmatrix( 0, dim, 0, dim );

  Beta       = p.Beta;
  Gamma      = p.Gamma;

  vcut       = 1.0;

  box_cells  = input.cells;
  box_length = (double) input.cells / input.cells_per_wl;

  name       = new( char [filename_size] );

  if( input.Q_restart == 1 ){
    char fname[ filename_size ];
    char varname[filename_size];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    sprintf( varname, "pha_%s.stepper.t_count", species_name );
    stepper.t_count = atoi( rf.getinput( varname ) );
    rf.closeinput();
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_phasespace::input_phasespace( parameter &p, char *species_name_input )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_phasespace::Constructor",errname);
  char input_name[filename_size];

  sprintf( species_name, "%s", species_name_input );

  rf.openinput( p.input_file_name );

  cells             = atoi( rf.setget( "&box", "cells" ) );
  cells_per_wl      = atoi( rf.setget( "&box", "cells_per_wl" ) );

  sprintf( input_name, "&%s_phasespace", species_name );
  stepper.Q         = atoi( rf.setget( input_name, "Q" ) );
  stepper.t_start   = atof( rf.setget( input_name, "t_start" ) );
  stepper.t_stop    = atof( rf.setget( input_name, "t_stop" ) );
  stepper.t_step    = atof( rf.setget( input_name, "t_step" ) );

  stepper.x_start   = -1;   // not used
  stepper.x_stop    = -1;   // not used
  stepper.x_step    = -1;   // not used

  Q_restart       = atoi( rf.setget( "&restart", "Q"     ) );
  strcpy( restart_file, rf.setget( "&restart", "file"    ) );

  rf.closeinput();

  bob.message("parameter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_phasespace::save( parameter &p )
{
  static error_handler bob("input_phasespace::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic " << species_name << "_phasespace" << endl;
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


void phasespace::write_phasespace( double time, parameter &p, domain *grid )
{
  static error_handler bob("phasespace::write_phasespace",errname);

  struct cell *cell;
  struct particle *part;
  int i, j;

  int dim1 = (int) floor( 1.0 * dim * grid->left->number / box_cells );
  int dim2 = (int) floor( 1.0 * dim * grid->right->number / box_cells );

  for( i=0; i<=dim; i++ )
    for( j=0; j<=dim; j++ )
      x[i][j]=y[i][j]=z[i][j]=0;

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {

      if (cell->npart != 0) {

	for( part=cell->first; part!=NULL; part=part->next ) {

	    if (part->species == species) {
	      vx = part->ux * part->igamma;
	      vy = part->uy * part->igamma;
	      vz = part->uz * part->igamma;

	      vx = 1.0/Gamma * vx / ( 1 + vy * Beta );
	      vz = 1.0/Gamma * vz / ( 1 + vy * Beta );
	      vy = ( vy + Beta ) / ( 1 + vy * Beta );

	      bx  = (int) floor( part->x/box_length * dim + 0.5 );
	      bvx = (int) floor( 0.5 * dim * (1 + vx/vcut) + 0.5 );
	      bvy = (int) floor( 0.5 * dim * (1 + vy/vcut) + 0.5 );
	      bvz = (int) floor( 0.5 * dim * (1 + vz/vcut) + 0.5 );

	      if (bvx>=0 && bvx<=dim) x[bvx][bx]++;
	      else bob.error( "velocity bin out of range" );
	      if (bvy>=0 && bvy<=dim) y[bvy][bx]++;
	      else bob.error( "velocity bin out of range" );
	      if (bvz>=0 && bvz<=dim) z[bvz][bx]++;
	      else bob.error( "velocity bin out of range" );
	    }
	  }
      }
    }

  sprintf(name,"%s/phasex-%d-sp%d-%.3f", p.path, p.domain_number, species, time);
  file_x = fopen( name, "wb" );
  fwrite( &dim1, sizeof(int), 1, file_x );
  fwrite( &dim2, sizeof(int), 1, file_x );

  sprintf(name,"%s/phasey-%d-sp%d-%.3f", p.path, p.domain_number, species, time);
  file_y = fopen( name, "wb" );
  fwrite( &dim1, sizeof(int), 1, file_y );
  fwrite( &dim2, sizeof(int), 1, file_y );

  sprintf(name,"%s/phasez-%d-sp%d-%.3f", p.path, p.domain_number, species,time);
  file_z = fopen( name, "wb" );
  fwrite( &dim1, sizeof(int), 1, file_z );
  fwrite( &dim2, sizeof(int), 1, file_z );

  for( i=0; i<=dim; i++ ) {
      fwrite( x[i]+dim1, sizeof(unsigned char), dim2-dim1+1, file_x );
      fwrite( y[i]+dim1, sizeof(unsigned char), dim2-dim1+1, file_y );
      fwrite( z[i]+dim1, sizeof(unsigned char), dim2-dim1+1, file_z );
    }

  fclose( file_x );
  fclose( file_y );
  fclose( file_z );
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof




