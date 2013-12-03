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

//
// this file should be always under construction
// any diagnostic can be added, replaced or removed
//
#include <diagnostic.h>

using namespace std;

diagnostic::diagnostic( parameter &p, domain* grid )
  : rf(),
    input(p),
    poi(p,grid),
    sna(p),
    vel_el(p),
    vel_ion(p),
    flu(p),
    ref(p),
    spa(p),
    ene(p,grid),
    tra(p),
    pha_el(p),
    pha_ion(p)
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("diagnostic::Constructor",errname);

  time_out          = p.spp;
  domain_number     = p.domain_number;

  strcpy(output_path,p.path);

  if(input.Q_restart == 0){
    time_steps        = 0;
    time_out_count    = time_out;
  }
  else{
    char fname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, domain_number );
    rf.openinput(fname);
    time_steps     = atoi( rf.getinput( "public_time_steps" ) );
    time_out_count = atoi( rf.getinput( "time_out_count" ) );
    rf.closeinput();
  }
  public_time_steps = time_steps;
}


//////////////////////////////////////////////////////////////////////////////////////////


input_diagnostic::input_diagnostic( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_diagnostic::Constructor",errname);

  rf.openinput( p.input_file_name );

  Q_restart         = atoi( rf.setget( "&restart", "Q" ) );
  strcpy( restart_file, rf.setget( "&restart", "file" ) );
  Q_restart_save    = atoi( rf.setget( "&restart", "Q_save" ) );
  strcpy( restart_file_save, rf.setget( "&restart", "file_save" ) );

  rf.closeinput();

  bob.message("parameter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_diagnostic::save( parameter &p )
{
  static error_handler bob("input_diagnostic::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Q_restart        : " << Q_restart       << endl;
  outfile << "restart_file     : " << restart_file    << endl;
  outfile << "Q_restart_save   : " << Q_restart_save  << endl;
  outfile << "restart_file_save: " << restart_file_save << endl << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void diagnostic::count( void )
{
  static error_handler bob("diagnostic::count",errname);

  time_steps     ++;
  time_out_count ++;
  public_time_steps = time_steps;

  ene.stepper.t_count  ++;
  flu.stepper.t_count  ++;
  sna.stepper.t_count  ++;
  pha_el.stepper.t_count  ++;
  pha_ion.stepper.t_count  ++;
  vel_el.stepper.t_count  ++;
  vel_ion.stepper.t_count  ++;
  ref.stepper.t_count  ++;
  poi.stepper.t_count  ++;
  tra.stepper.t_count  ++;
  spa.stepper_de.t_count    ++;
  spa.stepper_di.t_count    ++;
  spa.stepper_jx.t_count    ++;
  spa.stepper_jy.t_count    ++;
  spa.stepper_jz.t_count    ++;
  spa.stepper_ex.t_count    ++;
  spa.stepper_ey.t_count    ++;
  spa.stepper_ez.t_count    ++;
  spa.stepper_bx.t_count    ++;
  spa.stepper_by.t_count    ++;
  spa.stepper_bz.t_count    ++;
  spa.stepper_edens.t_count ++;
}


//////////////////////////////////////////////////////////////////////////////////////////


int diagnostic::write_window( int time_steps, diagnostic_stepper *stepper )
{
  int ok;

  if ( stepper->Q ) {                                // this diagnostic switched on ?
    if ( time_steps == stepper->t_start )
      stepper->t_count = stepper->t_step;            // start writing now !
    if ( time_steps >= stepper->t_start      &&      // inside time window and at
         time_steps <  stepper->t_stop       &&      //        point of time for writing ?
	 stepper->t_count == stepper->t_step    ) {  // yes:
      stepper->t_count = 0;                          //      reset write counter
      ok = 1;                                        //      return ok
    }
    else ok = 0;                                     // no:  do not write
  }
  else ok = 0;                                       // do not write

  return ok;
}



//////////////////////////////////////////////////////////////////////////////////////////


void diagnostic::out( double time, domain* grid, parameter &p )
{
  static error_handler bob("diagnostic::out",errname);

  if ( time_out_count == time_out ) {
    bob.message( "---------- TIME =", time, "----------" );
    time_out_count = 0;
  }

  // ---- traces -------------------------------------------------------------------------

  if ( tra.stepper.Q ) {
    if ( time_steps == tra.stepper.t_start )
         tra.stepper.t_count = 0;
    if (    time_steps >= tra.stepper.t_start
	 && time_steps <= tra.stepper.t_stop
         && tra.stepper.t_count == tra.stepper.t_step ) {
      tra.write_traces(time,p);
      tra.stepper.t_count = 0;
    }
  }

  if ( tra.stepper.Q ) {                          // store in memory
    if (    time_steps >= tra.stepper.t_start
	 && time_steps <= tra.stepper.t_stop )
	    tra.store_traces(grid);
  }

  // ---- energy, flux, reflectivity -----------------------------------------------------

  ene.average_reflex(grid);

  if ( write_window(time_steps,&(ene.stepper)) ) {
    ene.get_energies(grid);
    ene.write_energies(time);
  }

  if ( write_window(time_steps,&(flu.stepper)) )
    flu.write_flux(time,grid);

  ref.average_reflex(grid);

  if ( write_window(time_steps,&(ref.stepper)) )
    ref.write_reflex(time);

  // ---- snapshot -----------------------------------------------------------------------

  if ( write_window(time_steps,&(sna.stepper)) )
    sna.write_snap(time,grid,p);

  // ---- phasespace ---------------------------------------------------------------------

  if ( write_window(time_steps,&(pha_ion.stepper)) )
    pha_ion.write_phasespace(time,p,grid);

  if ( write_window(time_steps,&(pha_el.stepper)) )
    pha_el.write_phasespace(time,p,grid);

  // ---- velocity -----------------------------------------------------------------------

  if ( write_window(time_steps,&(vel_ion.stepper)) )
    vel_ion.write_velocity(time,p,grid);

  if ( write_window(time_steps,&(vel_el.stepper)) )
    vel_el.write_velocity(time,p,grid);

  // ---- poisson ------------------------------------------------------------------------

  if ( write_window(time_steps,&(poi.stepper)) ) {
    poi.solve(grid);
    poi.write(time,grid);
  }

  //---------------------- de ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_de)) ) spa.write_de(grid,time_out_count,p);

  //---------------------- di ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_di)) ) spa.write_di(grid,time_out_count,p);

  //---------------------- jx ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_jx)) ) spa.write_jx(grid,time_out_count,p);

  //---------------------- jy ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_jy)) ) spa.write_jy(grid,time_out_count,p);

  //---------------------- jz ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_jz)) ) spa.write_jz(grid,time_out_count,p);

  //---------------------- ex ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_ex)) ) spa.write_ex(grid,time_out_count,p);

  //---------------------- ey ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_ey)) ) spa.write_ey(grid,time_out_count,p);

  //---------------------- ez ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_ez)) ) spa.write_ez(grid,time_out_count,p);

  //---------------------- bx ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_bx)) ) spa.write_bx(grid,time_out_count,p);

  //---------------------- by ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_by)) ) spa.write_by(grid,time_out_count,p);

  //---------------------- bz ------------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_bz)) ) spa.write_bz(grid,time_out_count,p);

  //---------------------- edens ---------------------------------------------------------

  if ( write_window(time_steps,&(spa.stepper_edens)) )
    spa.write_edens(grid,time_out_count,p);
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof



