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

#include <diagnostic_reflex.h>


//////////////////////////////////////////////////////////////////////////////////////////


reflex::reflex( parameter &p )
  : rf(),
    input(p),
    stepper( input.stepper, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("reflex::Constructor",errname);

  if( input.Q_restart == 0 ){
    buf  = new( double [4] );
    buf[0] = 0;
    buf[1] = 0;
    buf[2] = 0;
    buf[3] = 0;

    if(stepper.Q){
      name  = new( char [filename_size] );
      sprintf(name, "%s/reflex-%d", p.path, p.domain_number);

      file.open(name,ios::app);
      if (!file) bob.error( "cannot open file", name );

      file.precision( 3 );
      file.setf( ios::showpoint | ios::scientific );

      file << "#" << setw(11) << "time"
	   << setw(14) << "R left"  << setw(12) << "R right" << endl;

      file.close();
    }
  }
  else{
    buf  = new( double [4] );

    char fname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    buf[0] = atof( rf.getinput( "ref.buf[0]" ) );
    buf[1] = atof( rf.getinput( "ref.buf[1]" ) );
    buf[2] = atof( rf.getinput( "ref.buf[2]" ) );
    buf[3] = atof( rf.getinput( "ref.buf[3]" ) );
    stepper.t_count = atoi( rf.getinput( "ref.stepper.t_count" ) );
    rf.closeinput();

    if(stepper.Q){
      name  = new( char [filename_size] );
      sprintf(name, "%s/reflex-%d", p.path, p.domain_number);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_reflex::input_reflex( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_reflex::Constructor",errname);

  rf.openinput( p.input_file_name );

  stepper.Q         = atoi( rf.setget( "&reflex", "Q" ) );
  stepper.t_start   = atof( rf.setget( "&reflex", "t_start" ) );
  stepper.t_stop    = atof( rf.setget( "&reflex", "t_stop" ) );
  stepper.t_step    = atof( rf.setget( "&reflex", "t_step" ) );
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


void input_reflex::save( parameter &p )
{
  static error_handler bob("input_reflex::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic reflex" << endl;
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


void reflex::average_reflex( domain *grid )
{
  static error_handler bob("reflex::average_reflex",errname);

  buf[0] += sqr(grid->left->fp) + sqr(grid->left->gm);
  buf[1] += sqr(grid->left->fm) + sqr(grid->left->gp);
  buf[2] += sqr(grid->right->fm) + sqr(grid->right->gp);
  buf[3] += sqr(grid->right->fp) + sqr(grid->right->gm);
}


//////////////////////////////////////////////////////////////////////////////////////////


void reflex::write_reflex( double time )
{
  static error_handler bob("reflex::write_reflex",errname);

  double rl, rr;

  if (buf[0]>0) rl = buf[1] / buf[0];
  else                 rl = 0;
  if (rl>10) rl = 0;


  if (buf[2]>0) rr = buf[3] / buf[2];
  else                 rr = 0;
  if (rr>10) rr = 0;

  buf[0] = 0;
  buf[1] = 0;
  buf[2] = 0;
  buf[3] = 0;

  file.open(name,ios::app);
  if (!file) bob.error( "cannot open file", name );

  file.precision( 3 );
  file.setf( ios::showpoint | ios::scientific );

  file << setw(12) << time
           << setw(14) << rl << setw(12) << rr << endl;

  file.close();
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof
