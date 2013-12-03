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

#include <diagnostic_energy.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////


energy::energy( parameter &p, domain* grid )
  : rf(),
    input(p),
    stepper( input.stepper, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("energy::Constructor",errname);

  if( input.Q_restart == 0 ){
    if(stepper.Q){
      name  = new( char [filename_size] );
      sprintf(name, "%s/energy-%d", p.path, p.domain_number);

      file.open(name,ios::out);
      if (!file) bob.error( "cannot open file", name );

      file.precision( 3 );
      file.setf( ios::showpoint | ios::scientific );

      file << "#" << setw(11) << "time"
	   << setw(12) << "flux"
	   << setw(12) << "total"
	   << setw(12) << "field"
	   << setw(12) << "field_tr"
	   << setw(12) << "field_lo"
	   << setw(12) << "kinetic" << endl;

      file.close();
    }

    flux      = 0;

    get_energies( grid );

    field_0   = field;
    field_t_0 = field_t;
    field_l_0 = field_l;
    kinetic_0 = kinetic;
    total_0   = total;
  }
  else{
    if(stepper.Q){
      name  = new( char [filename_size] );
      sprintf(name, "%s/energy-%d", p.path, p.domain_number);

      file.open(name,ios::app);
      if (!file) bob.error( "cannot open file", name );

      file.precision( 3 );
      file.setf( ios::showpoint | ios::scientific );

      file << endl;

      file.close();
    }
    char fname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    flux            = atof( rf.getinput( "ene.flux" ) );
    field_0         = atof( rf.getinput( "ene.field_0" ) );
    field_t_0       = atof( rf.getinput( "ene.field_t_0" ) );
    field_l_0       = atof( rf.getinput( "ene.field_l_0" ) );
    kinetic_0       = atof( rf.getinput( "ene.kinetic_0" ) );
    total_0         = atof( rf.getinput( "ene.total_0" ) );
    stepper.t_count = atoi( rf.getinput( "ene.stepper.t_count" ) );
    rf.closeinput();
  }
}

//////////////////////////////////////////////////////////////////////////////////////////


input_energy::input_energy( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_energy::Constructor",errname);

  rf.openinput( p.input_file_name );

  stepper.Q         = atoi( rf.setget( "&energy", "Q" ) );
  stepper.t_start   = atof( rf.setget( "&energy", "t_start" ) );
  stepper.t_stop    = atof( rf.setget( "&energy", "t_stop" ) );
  stepper.t_step    = atof( rf.setget( "&energy", "t_step" ) );
  stepper.x_start   = -1;   // not used
  stepper.x_stop    = -1;   // not used
  stepper.x_step    = -1;   // not used

  Q_restart       = atoi( rf.setget( "&restart", "Q"     ) );
  strcpy( restart_file, rf.setget( "&restart", "file"    ) );

  rf.closeinput();

  bob.message("paramter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_energy::save( parameter &p )
{
  static error_handler bob("input_energy::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic energy" << endl;
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


void energy::get_energies( domain* grid )
{
  static error_handler bob("energy::get_energies",errname);

  struct cell *cell;
  struct particle *part;

  field   = 0;
  field_l = 0;
  field_t = 0;
  kinetic = 0;
  total   = 0;

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next )
    {
      field_l += 0.5 * sqr(cell->ex);
      field_t += 0.5 * ( sqr(cell->ey)+sqr(cell->ez)+sqr(cell->bz)+sqr(cell->by) );

      if (cell->npart != 0) {

	for( part=cell->first; part!=NULL; part=part->next )
	  {
	    kinetic += part->n * part->m * ( 1.0/part->igamma - 1.0 );
	  }
      }
    }
  field = field_l + field_t;
  total = field + kinetic;
}


//////////////////////////////////////////////////////////////////////////////////////////


void energy::write_energies( double time )
{
  static error_handler bob("energy::write_energies",errname);

  file.open(name,ios::app);
  if (!file) bob.error( "cannot open file", name );

  file.precision( 3 );
  file.setf( ios::showpoint | ios::scientific );

  file << setw(11) << time
              << setw(12) << flux
              << setw(12) << total - total_0
	      << setw(12) << field - field_0
	      << setw(12) << field_t - field_t_0
	      << setw(12) << field_l - field_l_0
	      << setw(12) << kinetic - kinetic_0 << endl;

  file.close();
}


//////////////////////////////////////////////////////////////////////////////////////////


void energy::average_reflex( domain *grid )
{
  static error_handler bob("energy::average_reflex",errname);

  flux += ( sqr(grid->left->fp) + sqr(grid->left->gm) );
  flux -= ( sqr(grid->left->fm) + sqr(grid->left->gp) );
  flux += ( sqr(grid->right->fm) + sqr(grid->right->gp) );
  flux -= ( sqr(grid->right->fp) + sqr(grid->right->gm) );
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof

