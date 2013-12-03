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

#include <diagnostic_snapshot.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////


snapshot::snapshot( parameter &p )
  : rf(),
    input(p),
    stepper( input.stepper, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("snapshot::Constructor",errname);

  name  = new( char [filename_size] );

  if( input.Q_restart == 1 ){
    char fname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    stepper.t_count = atoi( rf.getinput( "sna.stepper.t_count" ) );
    rf.closeinput();
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_snapshot::input_snapshot( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_snapshot::Constructor",errname);

  rf.openinput( p.input_file_name );

  stepper.Q         = atoi( rf.setget( "&snapshot", "Q" ) );
  stepper.t_start   = atof( rf.setget( "&snapshot", "t_start" ) );
  stepper.t_stop    = atof( rf.setget( "&snapshot", "t_stop" ) );
  stepper.t_step    = atof( rf.setget( "&snapshot", "t_step" ) );
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


void input_snapshot::save( parameter &p )
{
  static error_handler bob("input_snapshot::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic snapshot" << endl;
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


void snapshot::write_snap( double time, domain* grid, parameter &p )
{
  static error_handler bob("snapshot::out_snap",errname);
  struct cell *cell;

  sprintf(name,"%s/snap-%d-%.3f", p.path, p.domain_number, time);

  file.open(name);
  if (!file) bob.error("cannot open snapshot file", name );

  file.precision( 3 );
  file.setf( ios::showpoint | ios::scientific );

  file << "#"
            << setw(11) << "x"
	    << setw(12) << "Ex"
            << setw(12) << "Ey"
            << setw(12) << "Ez"
            << setw(12) << "By"
            << setw(12) << "Bz"
	    << setw(12) << "rho_el"
            << setw(12) << "rho_ion"
	    << setw(12) << "jx"
            << setw(12) << "jy"
            << setw(12) << "jz"
            << setw(12) << "#el"
            << setw(12) << "#ion" << endl;

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next )
    {
      file << setw(12) << cell->x
		<< setw(12) << cell->ex
		<< setw(12) << cell->ey
		<< setw(12) << cell->ez
		<< setw(12) << cell->by
		<< setw(12) << cell->bz
		<< setw(12) << cell->dens[0]
		<< setw(12) << cell->dens[1]
                << setw(12) << cell->jx
		<< setw(12) << cell->jy
		<< setw(12) << cell->jz
	        << setw(12) << cell->np[0]
   	        << setw(12) << cell->np[1] << endl;
    }

  file.close();
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof

