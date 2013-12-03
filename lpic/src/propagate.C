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

#include <propagate.h>

//////////////////////////////////////////////////////////////////////////////////////////

propagate::propagate(parameter &p, domain &grid)
    : input(p),
      stk(p),
      rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("propagate::Constructor",errname);

  Gamma         = p.Gamma;
  dx            = grid.dx;
  idx           = 1.0/dx;
  dt            = dx / Gamma;         // <--------------------- LT -----------------------
  domain_number = p.domain_number;

  n_domains   = input.n_domains;

  start_time  = input.start_time;
  stop_time   = input.stop_time;

  if( input.Q_restart == 0 ) start_time = input.start_time;
  else{
    char fname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    start_time = atof( rf.getinput( "time" ) );
    rf.closeinput();
  }

  ofstream outfile;
  outfile.open(p.outname,ios::app);

  outfile << "propagate constructor" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Gamma              : " << Gamma         << endl;
  outfile << "dx                 : " << dx            << endl;
  outfile << "idx                : " << idx           << endl;
  outfile << "dt                 : " << dt            << endl;
  outfile << "domain             : " << domain_number << endl << endl << endl;

  outfile.close();
}

//////////////////////////////////////////////////////////////////////////////////////////


input_propagate::input_propagate( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_propagate::Constructor",errname);

  rf.openinput( p.input_file_name );

  start_time  = atoi( rf.setget( "&propagate", "prop_start" ) );
  stop_time   = atoi( rf.setget( "&propagate", "prop_stop"  ) );

  n_domains   = atoi( rf.setget( "&parallel", "N_domains" ) );

  Q_restart   = atoi( rf.setget( "&restart", "Q" ) );
  strcpy( restart_file, rf.setget( "&restart", "file" ) );

  rf.closeinput();

  bob.message( "parameter read" );

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_propagate::save( parameter &p )
{
  static error_handler bob("input_propagate::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "propagate" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Q_restart          : " << Q_restart      << endl;
  outfile << "restart_file       : " << restart_file   << endl;
  outfile << "prop_start         : " << start_time     << endl;
  outfile << "prop_stop          : " << stop_time      << endl;
  outfile << "N_domains          : " << n_domains      << endl << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void propagate::loop(parameter &p, box &sim, pulse &laser_front, pulse &laser_rear,
                     diagnostic &diag )
  // In case of several domains:
  // >> Usage of "sim.talk.current_1" and "sim.talk.current_2" leads to EXACTLY
  //    the same ordering of operations during adding the current contributions
  //    for any cell as in the case of only one domain. This leads to identical
  //    results for one or several domain, respectively, even in chaotic situations.
  //    This version may be used to check new network routines, but does NOT
  //    increase the job's speed, since a specified task propagates its particles
  //    after its parent task has finished.
  //    For this version, #DEFINE SLOW in header file common.h!
  // >> Usage of "sim.talk.current" instead leads to a slightly DIFFERENT ordering
  //    of operations during adding the current contributions for the buffer cells.
  //    Although PHYSICALLY COMPLETELY EQUIVALENT to the case of only one domain,
  //    slightly different roundoff errors in subtracting large numbers can lead
  //    to DIFFERENT simulation results in chaotic situations compared to the
  //    case of one domain.
  //    For this version #UNDEF SLOW in header file common.h!
{
  static error_handler bob("propagate::loop",errname);

  //////////////////////////////////// clocks ////////////////////////////////////////////
  uhr            zeit(p,"total"     );                                                  //
  uhr  zeit_particles(p,"particles" );                                                  //
  uhr     zeit_fields(p,"fields"    );                                                  //
  uhr zeit_diagnostic(p,"diagnostic");                                                  //
  ////////////////////////////////////////////////////////////////////////////////////////

  zeit.start();

  for( time = start_time; time <= stop_time + dt; time += dt )
    {
      sim.restart_save( diag, time, p, zeit, zeit_particles,
			zeit_fields, zeit_diagnostic );
      sim.count_restart();

      clear_grid( sim.grid );

#ifdef LPIC_PARALLEL
#ifdef SLOW                             // send/recieve current contributions and copies
      sim.talk.current_1(diag.public_time_steps, &(sim.grid) ); // SLOW
#endif
#endif

      zeit_particles.start();
      particles( sim.grid );            // accelerate and move
      reflect_particles( sim.grid );    // reflect particles at box boundaries
      zeit_particles.stop_and_add();

#ifdef LPIC_PARALLEL
#ifdef SLOW
      sim.talk.current_2(diag.public_time_steps, &(sim.grid) ); // SLOW
#else                                   // send/recieve current contributions and copies
      sim.talk.current( diag.public_time_steps, &(sim.grid) );  // FAST
#endif
      sim.talk.particles( diag.public_time_steps, &(sim.grid) );
                                        // send/recieve particles to/from
                                        // neighbour domains
      sim.talk.density( diag.public_time_steps, &(sim.grid) );
                                        // send/recieve density contributions
#endif

      zeit_fields.start();
      fields( sim.grid, laser_front, laser_rear );
                                        // propagate fields
      zeit_fields.stop_and_add();

#ifdef LPIC_PARALLEL
      sim.talk.field( diag.public_time_steps, &(sim.grid) );
                                        // send/recieve field copies to/from
	                                // neighbour domains
      sim.reorganize( sim.grid, sim.talk, time );
	                                // reorganize box
      sim.count_reorganize();           // reorganize counter
#endif

      zeit_diagnostic.start();
      diag.out( time, &(sim.grid), p ); // diagnostics
      diag.count();                     // diagnostic counter
      zeit_diagnostic.stop_and_add();

      zeit.add();                       // update clock
    }

  zeit.stop_and_add();

  zeit_particles.seconds_cpu();
  zeit_fields.seconds_cpu();
  zeit_diagnostic.seconds_cpu();
  zeit.seconds_cpu();
  if( input.Q_restart == 0 ) zeit.seconds_sys();
}


//////////////////////////////////////////////////////////////////////////////////////////


void propagate::clear_grid( domain &grid )
{
  static error_handler bob("propagate::clear_grid",errname);

  struct cell *cell;

  for( cell=grid.Lbuf; cell!=grid.dummy; cell=cell->next )
    {
      cell->charge  = 0;
      cell->dens[0] = 0;
      cell->dens[1] = 0;
      cell->jx      = 0;
      cell->jy      = 0;
      cell->jz      = 0;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////
//EOF



