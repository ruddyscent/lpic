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

#include <diagnostic_trace.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////


trace::trace( parameter &p )
  : rf(),
    input(p),
    stepper( input.stepper, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("trace::Constructor",errname);

  int i, j;

  if(stepper.Q){
    name  = new( char [filename_size] );

    traces = input.traces;

    if (traces > 0) {
      cell_number = new( int [traces+1] );
      fp     = matrix( 1, traces, 0, p.spp - 1 );
      fm     = matrix( 1, traces, 0, p.spp - 1 );
      gp     = matrix( 1, traces, 0, p.spp - 1 );
      gm     = matrix( 1, traces, 0, p.spp - 1 );
      ex     = matrix( 1, traces, 0, p.spp - 1 );
      dens_e = matrix( 1, traces, 0, p.spp - 1 );
      dens_i = matrix( 1, traces, 0, p.spp - 1 );
      jx     = matrix( 1, traces, 0, p.spp - 1 );
      jy     = matrix( 1, traces, 0, p.spp - 1 );
      jz     = matrix( 1, traces, 0, p.spp - 1 );
    }

    for( i=1; i<=traces; i++ ) {
      for( j=0; j<p.spp; j++ ) {
	fp[i][j]=0;
	fm[i][j]=0;
	gp[i][j]=0;
	gm[i][j]=0;
	ex[i][j]=0;
	dens_e[i][j]=0;
	dens_i[i][j]=0;
	jx[i][j]=0;
	jy[i][j]=0;
	jz[i][j]=0;
      }
    }

    for( i=0, j=1; i<input.traces; i++, j++ ) cell_number[j] = input.tracepos[i];
  }

  if( input.Q_restart == 1 ){
    char fname[ filename_size ];
    char dataname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    stepper.t_count = atoi( rf.getinput( "tra.stepper.t_count" ) );

    for(i=1; i<=traces; i++){
      sprintf( dataname, "fp[%d][0]", i);
      fp[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "fm[%d][0]", i);
      fm[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "gp[%d][0]", i);
      gp[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "gm[%d][0]", i);
      gm[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "ex[%d][0]", i);
      ex[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "dens_e[%d][0]", i);
      dens_e[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "dens_i[%d][0]", i);
      dens_i[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "jx[%d][0]", i);
      jx[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "jy[%d][0]", i);
      jy[i][0] = atof( rf.getinput( dataname ) );
      sprintf( dataname, "jz[%d][0]", i);
      jz[i][0] = atof( rf.getinput( dataname ) );
    }
    rf.closeinput();
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_trace::input_trace( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_trace::Constructor",errname);
  int i;
  char name[filename_size];

  domain_number     = p.domain_number;

  rf.openinput( p.input_file_name );

  path = new char [ filename_size ];
  strcpy( path, rf.setget( "&output", "path" ) );

  stepper.Q         = atoi( rf.setget( "&traces", "Q" ) );
  stepper.t_start   = atof( rf.setget( "&traces", "t_start" ) );
  stepper.t_stop    = atof( rf.setget( "&traces", "t_stop" ) );
  stepper.t_step    = 1;         // write output once per period, but store each time step
  stepper.x_start   = 0;
  stepper.x_stop    = atoi( rf.setget( "&box", "cells" ) );
  stepper.x_step    = 0;

  traces            = atoi( rf.setget( "&traces", "traces" ) );

  tracepos = new int [traces];
  for(i=0;i<traces;i++) {
    sprintf(name,"t%d",i);        // trace positions expected in variables t0, t1, t2, ...
    tracepos[i] = atoi( rf.setget( "&traces", name ) );
  }

  Q_restart       = atoi( rf.setget( "&restart", "Q"     ) );
  strcpy( restart_file, rf.setget( "&restart", "file"    ) );
  Q_restart_save  = atoi( rf.setget( "&restart", "Q_save"     ) );
  strcpy( restart_file_save, rf.setget( "&restart", "file_save"    ) );

  rf.closeinput();

  bob.message("parameter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_trace::save( parameter &p )
{
  static error_handler bob("input_trace::save",errname);
  ofstream outfile;
  int i;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic trace" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Q                : " << stepper.Q       << endl;
  outfile << "t_start          : " << stepper.t_start << endl;
  outfile << "t_stop           : " << stepper.t_stop  << endl;
  outfile << "t_step           : " << stepper.t_step  << endl;
  outfile << "traces           : " << traces          << endl;
  outfile << "                 : ";
  for(i=0;i<traces;i++) outfile << tracepos[i] << " ";
  outfile << endl;
  outfile << "Q_restart        : " << Q_restart       << endl;
  outfile << "restart_file     : " << restart_file    << endl << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void trace::store_traces( domain* grid )
{
  static error_handler bob("trace::store_traces",errname);
  struct cell *cell;
  int i;

  for( cell=grid->left; cell->next!=grid->rbuf; cell=cell->next )
    {
      for( i=1; i<=traces; i++ ) {

	if ( cell_number[i]==cell->number ) {

	  fp[i][stepper.t_count]     = (float) cell->fp;
	  fm[i][stepper.t_count]     = (float) cell->fm;
	  gp[i][stepper.t_count]     = (float) cell->gp;
	  gm[i][stepper.t_count]     = (float) cell->gm;
	  ex[i][stepper.t_count]     = (float) cell->ex;
	  dens_e[i][stepper.t_count] = (float) cell->dens[0];
	  dens_i[i][stepper.t_count] = (float) cell->dens[1];
	  jx[i][stepper.t_count]     = (float) cell->jx;
	  jy[i][stepper.t_count]     = (float) cell->jy;
	  jz[i][stepper.t_count]     = (float) cell->jz;
	}
      }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////


void trace::write_traces( double time, parameter &p )
{
  static error_handler bob("trace::store_traces",errname);

  int   period = (int) floor( time + 0.5 );
  float position;
  int   i, j;

  sprintf(name,"%s/trace-%d-%d", p.path, p.domain_number, period);

  file = fopen( name, "wb+" );
  if (!file) bob.error("cannot open", name);

  // header contains: time, # traces, # time steps per trace

  fwrite( &period,         sizeof(int), 1, file );
  fwrite( &(traces), sizeof(int), 1, file );
  fwrite( &(stepper.t_step),    sizeof(int), 1, file );

  // main body of the trace file

  for( i=1; i<=traces; i++ ) {

    position = (float) cell_number[i];

    fwrite( &position,   sizeof(float),           1, file );
    fwrite( fp[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( fm[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( gp[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( gm[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( ex[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( dens_e[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( dens_i[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( jx[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( jy[i], sizeof(float)*stepper.t_step, 1, file );
    fwrite( jz[i], sizeof(float)*stepper.t_step, 1, file );

  }

  fclose( file );

  for( i=1; i<=traces; i++ ) {
    for( j=0; j<p.spp; j++ ) {
      fp[i][j]=0;
      fm[i][j]=0;
      gp[i][j]=0;
      gm[i][j]=0;
      ex[i][j]=0;
      dens_e[i][j]=0;
      dens_i[i][j]=0;
      jx[i][j]=0;
      jy[i][j]=0;
      jz[i][j]=0;
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void trace::restart_save( void )
{
  static error_handler bob("trace::restart_save",errname);
  ofstream file1;
  char fname[ filename_size ];
  char  dataname[filename_size];
  int i;

  sprintf( fname, "%s/%s-%d-data1", input.path, input.restart_file_save,
	   input.domain_number);
  file1.open(fname,ios::app);
  if (!file1) bob.error( "cannot open file", fname );

  file1.precision( 20 );
  file1.setf( ios::showpoint | ios::scientific );

  for(i=1; i<=traces; i++){

    sprintf( dataname, "fp[%d][0]     = ", i);
    file1 << dataname <<  fp[i][0] << endl;
    sprintf( dataname, "fm[%d][0]     = ", i);
    file1 << dataname <<  fm[i][0] << endl;
    sprintf( dataname, "gp[%d][0]     = ", i);
    file1 << dataname <<  gp[i][0] << endl;
    sprintf( dataname, "gm[%d][0]     = ", i);
    file1 << dataname <<  gm[i][0] << endl;
    sprintf( dataname, "ex[%d][0]     = ", i);
    file1 << dataname <<  ex[i][0] << endl;
    sprintf( dataname, "dens_e[%d][0] = ", i);
    file1 << dataname <<  dens_e[i][0] << endl;
    sprintf( dataname, "dens_i[%d][0] = ", i);
    file1 << dataname <<  dens_i[i][0] << endl;
    sprintf( dataname, "jx[%d][0]     = ", i);
    file1 << dataname <<  jx[i][0] << endl;
    sprintf( dataname, "jy[%d][0]     = ", i);
    file1 << dataname <<  jy[i][0] << endl;
    sprintf( dataname, "jz[%d][0]     = ", i);
    file1 << dataname <<  jz[i][0] << endl;
  }
  file1.close();
}

//////////////////////////////////////////////////////////////////////////////////////////
//eof

