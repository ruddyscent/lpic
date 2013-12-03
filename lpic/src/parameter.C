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

#include <parameter.h>

parameter::parameter(int argc, char **argv)
  : rf()
{
  my_name         = new char [filename_size];
  input_file_name = new char [filename_size] ;
  errname         = new char [filename_size] ;
  outname         = new char [filename_size] ;
  path            = new char [filename_size] ;

  strcpy(my_name,argv[0]);

  //////////////// check commandline parameters //////////////////////////////////////////

  if (argc==1) {           // no argument: assume domain=1 and input from input/input.list
    domain_number = 1;
    sprintf( input_file_name, "input/input.list" );
    cout << "\n no arguments \n";
  }
  else if (argc==2) {           // one argument: interprete as domain number or input-file
    cout << "\n one argument: ";//               choose default for the other one
    domain_number = atoi( argv[1] );
    if (domain_number==0) {
      cout << "input file \n";
      domain_number = 1;
      strcpy( input_file_name, argv[1] );
    }
    else {
      cout << "domain number \n";
      sprintf( input_file_name, "input/input.list" );
    }
  }
  else {                                       // two arguments: domain number, input file
    cout << "\n two arguments \n";             //         additional arguments are ignored
    domain_number = atoi( argv[1] );
    strcpy( input_file_name, argv[2] );
    if (argc>3) cout << " arguments " << 3 << "-" << argc-1 << " ignored \n";
  }

  if (domain_number<1) {
      cerr << "\n domain number should be larger than zero!" << endl << endl;
      exit(-1);
  }
  else if (domain_number>2048) {
    cerr << "\n this machine does not exist so far!" << endl << endl;
    exit(-1);
  }

  //// read path and number of domains ///////////////////////////////////////////////////

  rf.openinput(input_file_name);

  strcpy( path, rf.setget( "&output", "path" ) );
  n_domains = atoi( rf.setget("&parallel","N_domains") );
  Q_restart = atoi( rf.setget("&restart","Q") );
  rf.closeinput();

  if (Q_restart) cout << " RESTART" << endl;
  cout << " domain      : " << domain_number << endl;
  cout << " input file  : " << input_file_name << endl;
  cout << " output path : " << path << endl << endl;

  //// set initial number of particle species ////////////////////////////////////////////

  nsp = 2;

  //// set in-out file name //////////////////////////////////////////////////////////////

  sprintf( outname, "%s/output.lpi", path );

  //// start first error handler /////////////////////////////////////////////////////////

  sprintf( errname, "%s/error-%d", path, domain_number );
  static error_handler bob("parameter::Constructor", errname);

#ifdef DEBUG
  bob.message("DEBUG is defined");
#else
  bob.message("DEBUG is undefined");
#endif
#ifdef LPIC_PARALLEL
  bob.message("LPIC_PARALLEL is defined");
#ifdef SLOW
  bob.message("SLOW is defined");
#else
  bob.message("SLOW is undefined");
#endif
#else
  bob.message("LPIC_PARALLEL is undefined");
#endif

  if ( !strcmp(my_name,"lpic_plain") && n_domains > 1) {
    n_domains = 1;
    bob.message( "# domains changed for plain version" );
  }

  bob.message( "program            =", my_name );
  bob.message( "domain number      =", domain_number );
  bob.message( "# domains          =", n_domains );
  bob.message( "# species          =", nsp );

  //// adjust angle such that # of steps per period is integer ///////////////////////////
  //// write spp and spl to file 'lpic.steps' for later use in lpic's postprocessor //////

  adjust_angle_write_steps();

  //// check contents of output directory and write in-out file for the first time ///////

  if (domain_number==1) {

    char input_copy[filename_size];
    sprintf( input_copy, "%s/input.lpi", path ); // write copy of input to path/input.lpi

    int diff = rf.compare_files( input_file_name, input_copy );

    if (diff<0) {            // file input_copy does not exist
      rf.copy_file( input_file_name, input_copy );
    }
    else if (diff==0) {      // file input_copy exists and is identical to input_file_name
      cout << " you try to run identical parameters again, exit." << endl << endl;
      bob.message( "you try to run identical parameters again, exit." );
      exit(-1);
    }
    else if (Q_restart==1) { // restart:
      if (diff==1) {         //          only one character changed, ok
	rf.copy_file( input_file_name, input_copy );
      }
      else {                 //          more than one character changed, not ok, exit!
	cout << " " << input_file_name << " and " << input_copy
	     << " differ by more than one character" << endl;
	cout << " if you really want to continue, remove " << input_copy
	     << " first!" << endl;
	cout << " exit." << endl << endl;
	bob.message( input_file_name, "and", input_copy );
	bob.message( "differ by more than one character" );
	bob.message( "if you really want to continue, remove", input_copy );
	bob.message( "first!" );
	bob.message( "exit." );
	exit(-1);
      }
    }
    else {                   // input_copy exists, no restart, exit!
      cout << " " << input_file_name << " and " << input_copy << " differ!" << endl;
      cout << " you try to overwrite old data!" << endl;
      cout << " exit." << endl << endl;
      bob.message( input_file_name, "and", input_copy, "differ!" );
      bob.message( "you try to overwrite old data!" );
      bob.message( "exit." );
      exit(-1);
    }

    save();                                      // save parameters to path/output.lpi
  }
};


//////////////////////////////////////////////////////////////////////////////////////////


void parameter::adjust_angle_write_steps( void )
// read and adjust angle of incidence, calculate steps per period
// write spp and spl to file 'lpic.steps' for later use in lpic's postprocessor
{
  static error_handler bob("parameter::adjust_angle",errname);

  rf.openinput( input_file_name );

  int    front_Q      = atoi( rf.setget( "&pulse_front", "Q" ) );
  int    rear_Q       = atoi( rf.setget( "&pulse_rear", "Q"  ) );
  double front_angle  = atof( rf.setget( "&pulse_front", "angle" ) );
  double rear_angle   = atof( rf.setget( "&pulse_rear", "angle" ) );

  if ( front_Q > 0 ) rear_angle = front_angle;
  else {
    if ( rear_Q > 0 ) front_angle = rear_angle;
    else              rear_angle  = front_angle;
  }

  spl         = atoi( rf.setget( "&box", "cells_per_wl" ) );
  spp         = (int) floor( 1.0/cos(PI/180*front_angle) * spl + 0.5 );
  angle       = 180.0/PI * acos( (double) spl / spp );
  Beta        = sin( PI/180 * angle );
  Gamma       = 1.0 / cos( PI/180 * angle );

  bob.message( "angles adjusted to ", angle, "degree" );
  bob.message( "# spl             =", spl );
  bob.message( "# spp             =", spp );

  rf.closeinput();

  FILE *f;
  char fname[filename_size];
  sprintf( fname, "%s/lpic.steps", path );

  f = fopen( fname, "w" );
  fprintf( f, "spl = %d\n", spl );
  fprintf( f, "spp = %d\n", spp );
  fclose( f );
}


//////////////////////////////////////////////////////////////////////////////////////////


void parameter::save( void )
{
  static error_handler bob("parameter::save",errname);

  ofstream outfile(outname);
  if (!outfile) bob.error("Cannot open outfile: ", outname);

  outfile << "parameter" << endl;
  outfile << "------------------------------------------------------------------" << endl;

  outfile << "program name       : " << my_name         << endl;
  outfile << "input file         : " << input_file_name << endl;
  outfile << "output path        : " << path            << endl;
  outfile << "domain number      : " << domain_number   << endl;
  outfile << "# domains          : " << n_domains       << endl;
  outfile << "# species          : " << nsp             << endl;
  outfile << "# steps per cycle  : " << spp             << endl;
  outfile << "adjusted angle     : " << angle           << endl;
  outfile << "LT-Beta            : " << Beta            << endl;
  outfile << "LT-Gamma           : " << Gamma           << endl << endl << endl;

  outfile.close();
};

//////////////////////////////////////////////////////////////////////////////////////////
//EOF

