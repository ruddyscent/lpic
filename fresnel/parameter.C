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
// fresnel
//
// initializes and reads parameters
// for documentation see 'parameter.h'
//
//////////////////////////////////////////////////////////////////////////////////////////

#include <parameter.h>

parameter::parameter(int argc, char **argv)
{
  plasma.wp             = 0.8;  // default parameters ////////////////////////////////////

  pulse.amplitude       = 0.01;
  pulse.angle           = 0.0;
  pulse.polarization    = 1;
  pulse.duration        = 60;

  prop.time_start       = 0.0;
  prop.time_stop        = 30.0;
  prop.spp              = 50;

  my_name = new(char[filename_size]);
  strcpy(my_name,argv[0]);

  path    = new char [filename_size];
  strcpy(path,"data");

  if (argc==1) { //////////////// check commandline parameters ///////////////////////////
    cerr << "\n\n no arguments " << endl;
  }
  else {
    cerr << "\n\n " << argc - 1 << " arguments" << endl << endl;
  }

  read_filename = new(char[filename_size]);

  if (argc==1) { //////////// read parameters from default filename //////////////////////
    char *fname = "input.fresnel";
    cout << " reading from " << fname << endl;
    read(fname);
  }
  else { //////////////////// read parameters from specified filename ////////////////////
    cout << " reading from " << argv[1] << endl;
    read(argv[1]);
  }

  /////////////////////////////////////////
  pulse.angle_rad = PI/180 * pulse.angle;//
  /////////////////////////////////////////

  cout << " output in    " << path << endl;

  sprintf( errname, "%s/error", path );
  static error_handler bob("parameter::Constructor", errname);

  save(path,errname);

};

//////////////////////////////////////////////////////////////////////////////////////////

void parameter::read(char *fname)
{
    Trash trash;
    int i;

    ifstream infile(fname);
    if (!infile) {
      cerr << "Cannot open infile: " << fname << endl;
      cerr << "Using default parameters" << endl;
    }
    else {
      cout << "reading from parameter file does not work." << endl;
      cout << "set parameters in parameter.C" << endl;
      exit( 0 );
      //
      infile >> trash >> plasma.wp;

      infile >> trash >> pulse.amplitude;
      infile >> trash >> pulse.angle;
      infile >> trash >> pulse.polarization;
      infile >> trash >> pulse.duration;

      infile >> trash >> prop.time_start;
      infile >> trash >> prop.time_stop;
      infile >> trash >> prop.spp;

      infile >> trash >> path;

      strcpy(read_filename,fname);

      infile.close();
    }
};


//////////////////////////////////////////////////////////////////////////////////////////


void parameter::save(char *path, char *errname)
{
    int i;
    static error_handler bob("parameter::save",errname);
    char *fname;

    fname = new char [filename_size];
    sprintf( fname, "%s/output.fresnel", path );
    ofstream outfile(fname);
    if (!outfile)
	bob.error("Cannot open outfile: ", fname);

    outfile << "fresnel parameters:" << endl << endl;

    outfile << "plasma" << endl;
    outfile << "-----------------------------------" << endl;
    outfile << "plasma frequency   : " << plasma.wp << endl << endl;

    outfile << "laser pulse" << endl;
    outfile << "-----------------------------------" << endl;
    outfile << "incident amplitude : " << pulse.amplitude     << endl;
    outfile << "angle of incidence : " << pulse.angle         << endl;
    outfile << "polarization       : " << pulse.polarization  << endl;
    outfile << "duration           : " << pulse.duration      << endl << endl;

    outfile << "propagation" << endl;
    outfile << "-----------------------------------" << endl;
    outfile << "time_start         : " << prop.time_start << endl;
    outfile << "time_stop          : " << prop.time_stop  << endl;
    outfile << "steps per period   : " << prop.spp        << endl << endl;

    outfile.close();
};

//////////////////////////////////////////////////////////////////////////////////////////
//EOF

