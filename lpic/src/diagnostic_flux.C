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

#include <diagnostic_flux.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////


flux::flux( parameter &p )
  : rf(),
    input(p),
    stepper( input.stepper, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("flux::Constructor",errname);

  if( input.Q_restart == 0 ){
    if(stepper.Q){
      name  = new( char [filename_size] );
      sprintf(name, "%s/flux-%d", p.path, p.domain_number);

      file.open(name,ios::out);
      if (!file) bob.error( "cannot open file", name );

      file.precision( 3 );
      file.setf( ios::showpoint | ios::scientific );

      file << "#" << setw(11) << "time"
	   << setw(14) << "S+ left"  << setw(12) << "S- left"
	   << setw(14) << "S+ right" << setw(12) << "S- right" << endl;

      file.close();
    }
  }
  else{
    if(stepper.Q){
      name  = new( char [filename_size] );
      sprintf(name, "%s/flux-%d", p.path, p.domain_number);

      file.open(name,ios::app);
      if (!file) bob.error( "cannot open file", name );

      file.precision( 3 );
      file.setf( ios::showpoint | ios::scientific );
      file << endl;
      file.close();

      char fname[ filename_size ];
      sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
      rf.openinput(fname);
      stepper.t_count = atoi( rf.getinput( "flu.stepper.t_count" ) );
      rf.closeinput();
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_flux::input_flux( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_flux::Constructor",errname);

  rf.openinput( p.input_file_name );

  stepper.Q         = atoi( rf.setget( "&flux", "Q" ) );
  stepper.t_start   = atof( rf.setget( "&flux", "t_start" ) );
  stepper.t_stop    = atof( rf.setget( "&flux", "t_stop" ) );
  stepper.t_step    = atof( rf.setget( "&flux", "t_step" ) );
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


void input_flux::save( parameter &p )
{
  static error_handler bob("input_flux::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic flux" << endl;
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


void flux::write_flux( double time, domain *grid )
{
  static error_handler bob("flux::write_flux",errname);

  file.open(name,ios::app);
  if (!file) bob.error( "cannot open file", name );

  file.precision( 3 );
  file.setf( ios::showpoint | ios::scientific );

  file << setw(12) << time
              << setw(14) << sqr(grid->left->fp) + sqr(grid->left->gm)
              << setw(12) << sqr(grid->left->fm) + sqr(grid->left->gp)
              << setw(14) << sqr(grid->right->fp) + sqr(grid->right->gm)
	      << setw(12) << sqr(grid->right->fm) + sqr(grid->right->gp) << endl;

  file.close();
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof
