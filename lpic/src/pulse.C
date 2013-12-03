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

#include <pulse.h>

int pulse::pulse_number = 0;

pulse::pulse( parameter &p, char *side )
  : input(p,side,pulse_number+1)
{
  pulse_number ++;

  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("pulse::Constructor",errname);

  path = new char [filename_size];
  strcpy(path,p.path);

  if (input.Q>0) {                  // pulse switched ON

    if (input.polarization==1) { Qz = 1.0; Qy = 0.0; shift = 0.0; }
    if (input.polarization==2) { Qz = 0.0; Qy = 1.0; shift = 0.0; }
    if (input.polarization==3) { Qz = 1.0; Qy = 1.0; shift = 0.25; }

    bob.message("pulse #", pulse_number, "generated");
    input.shape == 1 ? bob.message("shape       linear"):
    input.shape == 2 ? bob.message("shape       sin"):
    bob.message("shape       sinsqr");

    bob.message("amplitude    ", input.a0 );
    bob.message("raise_time   ", input.raise, " periods");
    bob.message("duration     ", input.duration, " periods");
    bob.message("polarization ", input.polarization );
    bob.message("color        ", input.a2, input.a3 );
    bob.message("phase        ", input.p2, input.p3 );

    if (input.Q_save==1 && p.domain_number==1 && input.Q_restart == 0)
      {
	save(input.time_start, input.time_stop, input.save_step);
      }
  }
  else {                              // pulse switched OFF
    Qz = 0.0; Qy = 0.0; shift = 0.0;

    bob.message("pulse #", pulse_number, "is switched off");
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_pulse::input_pulse( parameter &p, char *side, int pulse_number )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_pulse::Constructor",errname);

  rf.openinput(p.input_file_name);

  Q               = atoi( rf.setget( side, "Q"              ) );
  a0              = atof( rf.setget( side, "amplitude"      ) );
  a2              = atof( rf.setget( side, "amplitude2"     ) );
  a3              = atof( rf.setget( side, "amplitude3"     ) );
  p2              = atof( rf.setget( side, "phase2"         ) );
  p3              = atof( rf.setget( side, "phase3"         ) );
  polarization    = atoi( rf.setget( side, "polarization"   ) );
  shape           = atoi( rf.setget( side, "shape"          ) );
  raise           = atof( rf.setget( side, "raise"          ) );
  duration        = atof( rf.setget( side, "duration"       ) );
  Q_save          = atoi( rf.setget( side, "pulse_save"     ) );
  save_step       = atof( rf.setget( side, "pulse_save_step") );

  time_start      = atof( rf.setget( "&propagate", "prop_start") );
  time_stop       = atof( rf.setget( "&propagate", "prop_stop") );

  Q_restart       = atoi( rf.setget( "&restart", "Q"     ) );

  rf.closeinput();

  bob.message("parameter read");

  if (p.domain_number==1) save(p,pulse_number);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_pulse::save( parameter &p, int pulse_number )
{
  static error_handler bob("input_pulse::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "pulse # " << pulse_number << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Q                  : " << Q              << endl;
  outfile << "a0                 : " << a0             << endl;
  outfile << "a2                 : " << a2             << endl;
  outfile << "a3                 : " << a3             << endl;
  outfile << "phase2             : " << p2             << endl;
  outfile << "phase3             : " << p3             << endl;
  outfile << "polarization       : " << polarization   << endl;
  outfile << "shape              : " << shape          << endl;
  outfile << "raise              : " << raise          << endl;
  outfile << "duration           : " << duration       << endl;
  outfile << "Q_save             : " << Q_save         << endl;
  outfile << "save_step          : " << save_step      << endl;
  outfile << "time_start         : " << time_start     << endl;
  outfile << "time_stop          : " << time_stop      << endl;
  outfile << "Q_restart          : " << Q_restart      << endl << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void pulse::save(double t_start, double t_stop, double t_step)
{
    static error_handler bob("pulse::save",errname);
    char filename[filename_size];

    if (t_stop < t_start)
	bob.error("stop_time smaller than start_time");
    if (t_start < 0.)
	bob.error("start_time smaller than 0 not allowed for current pulse");
    if (t_step <0.)
	bob.error("step is negative !");

    sprintf(filename,"%s/pulse#%d", path, pulse_number );

    ofstream pulse_file(filename);
    if (!pulse_file)
	bob.error("cannot open output file: ", filename);

    pulse_file.precision( 4 );
    pulse_file.setf( ios::scientific | ios::showpoint );

    for (double t=t_start; t<t_stop; t+=t_step)
	pulse_file << setw(15) << t << " " << setw(15) << field(t) << endl;

    pulse_file.close();
}


//////////////////////////////////////////////////////////////////////////////////////////

double pulse::envelope( double t )
{
  return (
	  input.shape == 1 ?
	  (
	   t < input.raise ?            ( t / input.raise ) :
	   t < input.duration - input.raise ? 1.0 :
	   t < input.duration ?         ( input.duration - t ) / input.raise :
	   0
	   ) :
	   input.shape == 2 ?
	   (
	    t < input.raise ? 	   sin( 0.5 * PI * t / input.raise ) :
	    t < input.duration - input.raise ? 1.0 :
	    t < input.duration ?         sin( 0.5 * PI * (input.duration - t) / input.raise ) :
	    0
	    ) :
	   // assume shape == 3
	     t < input.raise ?         sqr( sin( 0.5 * PI * t / input.raise ) ) :
	     t < input.duration - input.raise ? 1.0 :
	     t < input.duration ?      sqr( sin( 0.5 * PI * (input.duration - t) / input.raise ) ) :
	     0
	     );
}

//////////////////////////////////////////////////////////////////////////////////////////


double pulse::field( double t )
{
  double env = envelope( t );
  double amp;

  amp =  input.a0 * env * sin( 2.0*PI*t );
  amp += input.a2 * env * sin( 4.0*PI*t + input.p2 );
  // amp += input.a2 * env * sin( 1.4*2.0*PI*t + input.p2 );
  amp += input.a3 * env * sin( 6.0*PI*t + input.p3 );

  return amp;
}

//////////////////////////////////////////////////////////////////////////////////////////
//eof
