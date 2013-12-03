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

//////////////////////////////////////////////////////////////////////////////////////////
//
// postprocessor for lpic++
//
//////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <parameter.h>

using namespace std;

parameter::parameter(int argc, char **argv)
  : rf()
{
  my_name = argv[0];

  errname             = new char [filename_size];
  read_filename       = new char [filename_size];
  save_filename       = new char [filename_size];
  output_path         = new char [filename_size];
  save_path_name      = new char [filename_size];
  file_path           = new char [filename_size];

  if (argc<3) {              // check commandline parameters
    printf( "\n two arguments required: input and output path\n\n");
    exit(0);
  }
  else {
    strcpy( file_path, argv[1] );
    strcpy( output_path, argv[2] );
  }

  sprintf( errname, "%s/error", output_path );
  static error_handler bob("parameter::Constructor", errname );

  bob.message( "sizeof(unsigned char)=",sizeof(unsigned char));
  bob.message( "sizeof(int)=",sizeof(int));

  bob.message("reading lpi data files from ", file_path);
  bob.message("writing post files to       ", output_path);
  bob.message("reading further input from   input.post");

  cout << endl;
  cout << "reading lpi data files from " << file_path << endl;
  cout << "writing post-lpi files to   " << output_path << endl;
  cout << "reading further input from  input.post" << endl << endl;

  strcpy(read_filename,"input.post");
  strcpy(save_filename,"output.post");
  sprintf(save_path_name, "%s/%s", output_path, save_filename);

  //  read(read_filename);

  save(save_path_name);
};

//////////////////////////////////////////////////////////////////////////////////////////

/*
void parameter::read(char *fname)
{
  static error_handler bob("parameter::read", errname);

  rf.openinput(fname);

  rf.closeinput();
}
*/

//////////////////////////////////////////////////////////////////////////////////////////


void parameter::save(char *fname)
{
  static error_handler bob("parameter::save", errname);

  ofstream outfile(fname);
  if (!outfile)
    bob.error("Cannot open outfile: ", fname);

  outfile << "postprocessor parameters:" << endl << endl;

  outfile << "Input File Path" << endl;
  outfile << "--------------------------------------------------" << endl;
  outfile << "input file path: " << file_path << endl << endl;

  outfile << "Output File Path" << endl;
  outfile << "--------------------------------------------------" << endl;
  outfile << "output path:     " << output_path << endl << endl;

  outfile.close();
}

//////////////////////////////////////////////////////////////////////////////////////////
//EOF

