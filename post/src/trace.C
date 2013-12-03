/*
   This file is part of LPIC++, a particle-in-cell code for
   simulating the interaction of laser light with plasma.

   Copyright (C) 2002      Andreas Kemp
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
// changes by A.Kemp, 2002 denoted by ##
#include <trace.h>


trace::trace( parameter &p )
        : input(p),
	  ft( input.periods, input.steps_pp, input.screen )
{
  sprintf( errname, "%s/error", p.output_path );
  static error_handler bob("trace::Constructor",errname);

  period_start = input.period_start;
  period_stop  = input.period_stop;
  periods      = input.periods;
  traces       = input.traces;
  steps_pp     = input.steps_pp;

  bob.message("reached steps:0. # traces = ", traces);

  read_name    = new char[filename_size];
  power_name   = new char[filename_size];
  trace_name   = new char[filename_size];
  path         = new char[filename_size];

  strcpy(path,p.file_path);

  position     = new float [traces+1];

  Q = 0;

  vector_read = new float [ steps_pp ];
  if ( input.Q_ey || input.Q_bz || input.Q_Pi || input.Q_Pr )
    Q = input.Q_fp = input.Q_fm = 1;
  if ( input.Q_ez || input.Q_by || input.Q_Si || input.Q_Sr )
    Q = input.Q_gp = input.Q_gm = 1;
  if ( input.Q_fp ) { fp = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_fm ) { fm = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_gm ) { gm = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_gp ) { gp = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_ex ) { ex = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_ey ) { ey = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_ez ) { ez = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_by ) { by = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_bz ) { bz = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_de ) { de = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_di ) { di = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_jx ) { jx = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_jy ) { jy = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }
  if ( input.Q_jz ) { jz = fmatrix( 1, traces, 0, periods*steps_pp - 1 ); Q = 1; }


  if (Q) {
    printf( "reading traces ... " ); fflush(stdout);
    read_traces();                                  // read from trace files
    printf( "done\n" );
  }

  bob.message("Q              =", Q );
  bob.message("traces         =", traces );
  bob.message("periods        =", ft.periods_input );
  bob.message("steps_pp_input =", ft.steps_pp_input );
  bob.message("steps_input    =", ft.steps_input );
  bob.message("steps          =", ft.steps );
}


//////////////////////////////////////////////////////////////////////////////////////////


input_trace::input_trace( parameter &p )
  : rf()
{
  strcpy( errname, p.errname );
  static error_handler bob("input_trace::Constructor",errname);
  char filename[filename_size];
  FILE *file;

  rf.openinput( p.read_filename );

  period_start = atoi( rf.setget( "&traces", "period_start"  ) );
  period_stop  = atoi( rf.setget( "&traces", "period_stop"   ) );
  screen       = atoi( rf.setget( "&traces", "period_screen" ) );

  Q_fp = atoi( rf.setget( "&traces", "fp" ) );
  Q_fm = atoi( rf.setget( "&traces", "fm" ) );
  Q_gp = atoi( rf.setget( "&traces", "gp" ) );
  Q_gm = atoi( rf.setget( "&traces", "gm" ) );
  Q_ex = atoi( rf.setget( "&traces", "ex" ) );
  Q_ey = atoi( rf.setget( "&traces", "ey" ) );
  Q_ez = atoi( rf.setget( "&traces", "ez" ) );
  Q_by = atoi( rf.setget( "&traces", "by" ) );
  Q_bz = atoi( rf.setget( "&traces", "bz" ) );
  Q_Pi = atoi( rf.setget( "&traces", "Pi" ) );
  Q_Pr = atoi( rf.setget( "&traces", "Pr" ) );
  Q_Sr = atoi( rf.setget( "&traces", "Sr" ) );
  Q_Si = atoi( rf.setget( "&traces", "Si" ) );
  Q_de = atoi( rf.setget( "&traces", "de" ) );
  Q_di = atoi( rf.setget( "&traces", "di" ) );
  Q_jx = atoi( rf.setget( "&traces", "jx" ) );
  Q_jy = atoi( rf.setget( "&traces", "jy" ) );
  Q_jz = atoi( rf.setget( "&traces", "jz" ) );

  rf.closeinput();

  // read trace-file headers to determine number of traces -------------------------------

  int period=1; // ##
  int region=1, traces_read, traces_min=10000, traces_max=0;

  if ( Q_fp || Q_fm || Q_gp || Q_gm || Q_ex || Q_ey || Q_ez || Q_by || Q_bz ||
       Q_Pi || Q_Pr || Q_Sr || Q_Si || Q_de || Q_di || Q_jx || Q_jy || Q_jz ) {

    sprintf( filename, "%s/trace-%d-%d", p.file_path, region, period );

    while( (file = fopen( filename, "rb" )) )
      {
	fread( &period,      sizeof(int), 1, file );
	fread( &traces_read, sizeof(int), 1, file );              // read number of traces
	fread( &steps_pp,    sizeof(int), 1, file );              // read steps_pp
	rewind( file );
	fclose( file );

	if (traces_read>traces_max) traces_max=traces_read;
	if (traces_read<traces_min) traces_min=traces_read;

	sprintf( filename, "%s/trace-%d-%d", p.file_path, ++region, period );
      }

    if( traces_max!=traces_min) bob.error( "number of traces not unique" );

  }
  bob.message("traces.min = ",traces_min, "traces_max = ",traces_max); //##
  traces = traces_max;                       // input parameter input.traces ! // ##
  periods = period_stop - period_start + 1;  // input parameter input.periods  // ##

  // -------------------------------------------------------------------------------------

  bob.message("parameter read");

  save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_trace::save( parameter &p )
{
  static error_handler bob("input_trace::save",errname);
  ofstream outfile;

  outfile.open(p.save_path_name,ios::app);

  outfile << "Traces and Fourier Transforms with respect to time" << endl;
  outfile << "--------------------------------------------------" << endl;
  outfile << "region         : " << region       << endl;
  outfile << "period_start   : " << period_start << endl;
  outfile << "period_stop    : " << period_stop  << endl;
  outfile << "periods_screen : " << screen       << endl;
  outfile << "traces         : " << traces       << endl;
  outfile << "steps_pp       : " << steps_pp     << endl << endl;

  outfile << "fp fm gp gm   ex ey ez by bz   Pi Pr Sr Si   de di jx jy jz" << endl;
  outfile << ":" << Q_fp << " ";
  outfile << ":" << Q_fm << " ";
  outfile << ":" << Q_gp << " ";
  outfile << ":" << Q_gm << "   ";
  outfile << ":" << Q_ex << " ";
  outfile << ":" << Q_ey << " ";
  outfile << ":" << Q_ez << " ";
  outfile << ":" << Q_by << " ";
  outfile << ":" << Q_bz << "   ";
  outfile << ":" << Q_Pi << " ";
  outfile << ":" << Q_Pr << " ";
  outfile << ":" << Q_Sr << " ";
  outfile << ":" << Q_Si << "   ";
  outfile << ":" << Q_de << " ";
  outfile << ":" << Q_di << " ";
  outfile << ":" << Q_jx << " ";
  outfile << ":" << Q_jy << " ";
  outfile << ":" << Q_jz << " " << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void trace::read_traces( void )
{
  static error_handler bob("trace::read_traces",errname);

  int i, j, k, n_curr, domain;
  int period_read, traces_read, traces_previous, steps_pp_read;

  for( i=period_start, n_curr=0; i<=period_stop; i++, n_curr+=steps_pp ){// forall periods

    traces_read     = 0;
    traces_previous = 0;
    domain          = 1;

    sprintf( read_name, "%s/trace-%d-%d", path, domain, i );

    while( (read_file = fopen( read_name, "rb" )) )
      {
	//	bob.message( "opening file:", read_name );

	fread( &period_read, sizeof(int), 1, read_file );            // consistency checks
	fread( &traces_read, sizeof(int), 1, read_file );
	fread( &steps_pp_read, sizeof(int), 1, read_file );

	bob.message( "file   ", read_name);
	bob.message( "period ", period_read );
	bob.message( "traces ", traces_read );
	bob.message( "steps  ", steps_pp_read );

	if ( period_read != i ) bob.error( "wrong period in file", read_name );
	if ( steps_pp_read != steps_pp ) bob.error( "steps_pp not unique" );

	for( j = 1; j <= traces_read; j++ ) {

	  fread( &position[j], sizeof(float), 1, read_file );  // read position of trace j

	  if (input.Q_fp) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) fp[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_fm) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) fm[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_gp) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) gp[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_gm) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) gm[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_ex) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) ex[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_de) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) de[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_di) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) di[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_jx) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) jx[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_jy) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) jy[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	  if (input.Q_jz) {
	    for( k=0; k<steps_pp; k++ ) vector_read[k]=0;
	    fread( vector_read, sizeof(float), steps_pp, read_file );
	    for( k=0; k<steps_pp; k++ ) jz[j][n_curr+k] += vector_read[k];
	  }
	  else              fseek( read_file, sizeof(float)*steps_pp, 1 );
	}
	fclose( read_file );
	sprintf( read_name, "%s/trace-%d-%d", path, ++domain, i );
      }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////


void trace::transform( parameter &p )
{
  static error_handler bob("trace::transform",errname);
  int i, j;

  bob.message( "steps =", (double)periods*steps_pp, " steps_ft =", (double)ft.steps );

  power_spectrum = dmatrix( 1, traces, 0, ft.steps_half );

  if (input.Q_fp) {
    printf( "transforming fp ... " ); fflush(stdout);
    bob.message("fp");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( fp[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"fp");
    write_traces(p,fp,"fp");
    printf( "done\n" );
  }

  if (input.Q_fm) {
    printf( "transforming fm ... " ); fflush(stdout);
    bob.message("fm");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( fm[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"fm");
    write_traces(p,fm,"fm");
    printf( "done\n" );

    ///////////////////////////////////////////////////////////////
    // class ft can be used to calculate correlations, see the
    // example below
    ///////////////////////////////////////////////////////////////
    //    double mid=10;
    //    double width=2;
    //    FILE   *f;
    //
    //    printf( "correlation fm ... " ); fflush(stdout);
    //    ft.correlation( fm[1], mid, width );
    //
    //    f=fopen( "correlation.dat", "w" );
    //    for(  i=0; i<ft.steps_half; i++ ) {
    //      fprintf( f, "\n %f  %.3e  %.3e", ft.dt*i, ft.local[i], ft.corr[i] );
    //    }
    //    fclose(f);
    //
    //    printf( "done\n" );
    ///////////////////////////////////////////////////////////////
  }

  if (input.Q_gp) {
    printf( "transforming gp ... " ); fflush(stdout);
    bob.message("gp");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( gp[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"gp");
    write_traces(p,gp,"gp");
    printf( "done\n" );
  }

  if (input.Q_gm) {
    printf( "transforming gm ... " ); fflush(stdout);
    bob.message("gm");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( gm[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"gm");
    write_traces(p,gm,"gm");
    printf( "done\n" );
  }

  if (input.Q_ex) {
    printf( "transforming ex ... " ); fflush(stdout);
    bob.message("ex");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( ex[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"ex");
    write_traces(p,ex,"ex");
    printf( "done\n" );
  }

  if (input.Q_ey) {
    printf( "transforming ey ... " ); fflush(stdout);
    bob.message("ey");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( fp[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
      ft.RealFt( fm[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] += ft.power[j];
      for( j=0; j<ft.steps_input; j++ ) ey[i][j] = fp[i][j] + fm[i][j];
    }
    write_transform(p,"ey");
    write_traces(p,ey,"ey");
    printf( "done\n" );
  }

  if (input.Q_bz) {
    printf( "transforming bz ... " ); fflush(stdout);
    bob.message("bz");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( fp[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
      ft.RealFt( fm[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] -= ft.power[j];
      for( j=0; j<ft.steps_input; j++ ) bz[i][j] = fp[i][j] - fm[i][j];
    }
    write_transform(p,"bz");
    write_traces(p,bz,"bz");
    printf( "done\n" );
  }

  if (input.Q_ez) {
    printf( "transforming ez ... " ); fflush(stdout);
    bob.message("ez");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( gp[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
      ft.RealFt( gm[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] += ft.power[j];
      for( j=0; j<ft.steps_input; j++ ) ez[i][j] = gp[i][j] + gm[i][j];
    }
    write_transform(p,"ez");
    write_traces(p,ez,"ez");
    printf( "done\n" );
  }

  if (input.Q_by) {
    printf( "transforming by ... " ); fflush(stdout);
    bob.message("by");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( gp[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
      ft.RealFt( gm[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] -= ft.power[j];
      for( j=0; j<ft.steps_input; j++ ) by[i][j] = gp[i][j] - gm[i][j];
    }
    write_transform(p,"by");
    write_traces(p,by,"by");
    printf( "done\n" );
  }

  if (input.Q_Pi) {
    printf( "transforming Pi ... " ); fflush(stdout);
    bob.message("Pi");
    for( i=1; i<=traces; i++ ) {                       // factor 2 because of E_z and B_y!
      ft.RealFt( fp[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = 2 * ft.power[j];
    }
    write_transform(p,"Pi");
    printf( "done\n" );
  }

  if (input.Q_Pr) {
    printf( "transforming Pr ... " ); fflush(stdout);
    bob.message("Pr");
    for( i=1; i<=traces; i++ ) {                        // factor 2 because of E_y and B_z
      ft.RealFt( fm[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = 2 * ft.power[j];
    }
    write_transform(p,"Pr");
    printf( "done\n" );
  }

  if (input.Q_Sr) {
    printf( "transforming Sr ... " ); fflush(stdout);
    bob.message("Sr");
    for( i=1; i<=traces; i++ ) {                       // factor 2 because of E_z and B_y!
      ft.RealFt( gp[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = 2 * ft.power[j];
    }
    write_transform(p,"Sr");
    printf( "done\n" );
  }

  if (input.Q_Si) {
    printf( "transforming Si ... " ); fflush(stdout);
    bob.message("Si");
    for( i=1; i<=traces; i++ ) {                       // factor 2 because of E_z and B_y!
      ft.RealFt( gm[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = 2 * ft.power[j];
    }
    write_transform(p,"Si");
    printf( "done\n" );
  }

  if (input.Q_de) {
    printf( "transforming de ... " ); fflush(stdout);
    bob.message("de");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( de[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"de");
    write_traces(p,de,"de");
    printf( "done\n" );
  }

  if (input.Q_di) {
    printf( "transforming di ... " ); fflush(stdout);
    bob.message("di");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( di[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"di");
    write_traces(p,di,"di");
    printf( "done\n" );
  }

  if (input.Q_jx) {
    printf( "transforming jx ... " ); fflush(stdout);
    bob.message("jx");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( jx[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"jx");
    write_traces(p,jx,"jx");
    printf( "done\n" );
  }

  if (input.Q_jy) {
    printf( "transforming jy ... " ); fflush(stdout);
    bob.message("jy");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( jy[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"jy");
    write_traces(p,jy,"jy");
    printf( "done\n" );
  }

  if (input.Q_jz) {
    printf( "transforming jz ... " ); fflush(stdout);
    bob.message("jz");
    for( i=1; i<=traces; i++ ) {
      ft.RealFt( jz[i] );
      for( j=0; j<=ft.steps_half; j++ ) power_spectrum[i][j] = ft.power[j];
    }
    write_transform(p,"jz");
    write_traces(p,jz,"jz");
    printf( "done\n" );
  }

  delete_dmatrix(power_spectrum,1,traces,0,ft.steps_half);
}


//////////////////////////////////////////////////////////////////////////////////////////


void trace::write_transform( parameter &p, char* appendix )
{
  static error_handler bob("trace::write_transform",errname);
  int i, j;

  sprintf(power_name, "%s/ft-%s", p.output_path, appendix);
  powerfile.open(power_name);
  if (!powerfile) bob.error("cannot open power_spectrum_file: ", power_name );

  powerfile.precision( 5 );
  powerfile.setf( ios::showpoint | ios::scientific );

  powerfile << "#" << setw(12) << "frequency";
  for( i=1; i<=traces; i++ ) powerfile << setw(13) << position[i];
  powerfile << endl;

  for( j=0; j<=ft.steps_half; j++ ) {
    powerfile << setw(13) << ft.frequency[j];
    for( i=1; i<=traces; i++ ) powerfile << setw(13) << power_spectrum[i][j];
    powerfile << endl;
  }

  powerfile.close();

}


//////////////////////////////////////////////////////////////////////////////////////////


void trace::write_traces( parameter &p, float** input, char* appendix )
{
  static error_handler bob("trace::write_traces",errname);
  double t_start = ft.dt_input * ( period_start - 1 ) * steps_pp;
  int MAX_STEPS  = 5000;
  int sample     = 1 + (int) floor( (double) ft.steps_input / MAX_STEPS );
  int i, j;

//  bob.message("t_start = ", t_start );

  sprintf(trace_name, "%s/ft-%s-trace", p.output_path, appendix);
  tracefile.open(trace_name);
  if (!tracefile) bob.error("cannot open ascii_trace_file: ", trace_name );

  tracefile.precision( 5 );
  tracefile.setf( ios::showpoint | ios::scientific );

  tracefile << "#" << setw(12) << "time";
  for( i=1; i<=traces; i++ ) tracefile << setw(13) << position[i];
  tracefile << endl;

  for( j=0; j<ft.steps_input; j+=sample ) {
    tracefile << setw(13) << t_start + ft.dt_input * j;
    for( i=1; i<=traces; i++ ) tracefile << setw(13) << input[i][j];
    tracefile << endl;
  }

  tracefile.close();
}

//////////////////////////////////////////////////////////////////////////////////////////
//eof





