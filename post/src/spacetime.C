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
#include <spacetime.h>

spacetime::spacetime( parameter &p )
  : input(p),
    ft(   input.periods_x, input.cells_per_wl,     1 ),
    ft2d( input.Q_kw, input.periods_t, input.steps_per_period, 1,
	              input.periods_x, input.cells_per_wl,     1 )
{
  sprintf( errname, "%s/error", p.output_path );
  static error_handler bob("spacetime::Constructor",errname);

  input_path = new char [ filename_size ];
  strcpy(input_path,p.file_path);
  output_path = new char [ filename_size ];
  strcpy(output_path,p.output_path);
}


//////////////////////////////////////////////////////////////////////////////////////////


input_spacetime::input_spacetime( parameter &p )
  : rf()
{
  strcpy( errname, p.errname );
  static error_handler bob("input_spacetime::Constructor",errname);
  char fname[filename_size];

  rf.openinput( p.read_filename );

  t_start   = atoi( rf.setget( "&spacetime", "t_start" ) );
  t_stop    = atoi( rf.setget( "&spacetime", "t_stop" ) );
  x_start   = atof( rf.setget( "&spacetime", "x_start" ) );
  x_stop    = atof( rf.setget( "&spacetime", "x_stop" ) );
  x_offset  = atof( rf.setget( "&spacetime", "x_offset" ) );

  periods_x = (int) floor( x_stop - x_start + 0.5 );
  periods_t = (int) floor( t_stop - t_start + 0.5 );

  average   = atoi( rf.setget( "&spacetime", "smooth" ) );
  size      = atoi( rf.setget( "&spacetime", "imagesize" ) );
  contour_1 = atof( rf.setget( "&spacetime", "contour_1" ) );
  contour_2 = atof( rf.setget( "&spacetime", "contour_2" ) );
  contour_3 = atof( rf.setget( "&spacetime", "contour_3" ) );

  Q_kw = atoi( rf.setget( "&spacetime", "Q_kw" ) );
  Q_kt = atoi( rf.setget( "&spacetime", "Q_kt" ) );

  Q_de = atoi( rf.setget( "&spacetime", "Q_de" ) );
  Q_di = atoi( rf.setget( "&spacetime", "Q_di" ) );
  Q_jx = atoi( rf.setget( "&spacetime", "Q_jx" ) );
  Q_jy = atoi( rf.setget( "&spacetime", "Q_jy" ) );
  Q_jz = atoi( rf.setget( "&spacetime", "Q_jz" ) );
  Q_ex = atoi( rf.setget( "&spacetime", "Q_ex" ) );
  Q_ey = atoi( rf.setget( "&spacetime", "Q_ey" ) );
  Q_ez = atoi( rf.setget( "&spacetime", "Q_ez" ) );
  //  Q_bx = atoi( rf.setget( "&spacetime", "Q_bx" ) );
  Q_bx = 0; // ##
  Q_by = atoi( rf.setget( "&spacetime", "Q_by" ) );
  Q_bz = atoi( rf.setget( "&spacetime", "Q_bz" ) );
  Q_edens = atoi( rf.setget( "&spacetime", "Q_edens" ) );
  Q_de_fi = atoi( rf.setget( "&spacetime", "Q_de_fi" ) );
  Q_de_ii = atoi( rf.setget( "&spacetime", "Q_de_ii" ) );

  C_kw = atof( rf.setget( "&spacetime", "C_kw" ) );
  C_kt = atof( rf.setget( "&spacetime", "C_kt" ) );

  K_cut = atof( rf.setget( "&spacetime", "K_cut" ) );
  W_cut = atof( rf.setget( "&spacetime", "W_cut" ) );

  C_de = atof( rf.setget( "&spacetime", "C_de" ) );
  C_di = atof( rf.setget( "&spacetime", "C_di" ) );
  C_jx = atof( rf.setget( "&spacetime", "C_jx" ) );
  C_jy = atof( rf.setget( "&spacetime", "C_jy" ) );
  C_jz = atof( rf.setget( "&spacetime", "C_jz" ) );
  C_ex = atof( rf.setget( "&spacetime", "C_ex" ) );
  C_ey = atof( rf.setget( "&spacetime", "C_ey" ) );
  C_ez = atof( rf.setget( "&spacetime", "C_ez" ) );
  //  C_bx = atof( rf.setget( "&spacetime", "C_bx" ) );
  C_bx = 0; // ##
  C_by = atof( rf.setget( "&spacetime", "C_by" ) );
  C_bz = atof( rf.setget( "&spacetime", "C_bz" ) );
  C_edens = atof( rf.setget( "&spacetime", "C_edens" ) );
  C_de_fi = atof( rf.setget( "&spacetime", "C_de_fi" ) );
  C_de_ii = atof( rf.setget( "&spacetime", "C_de_ii" ) );

  rf.closeinput();

  sprintf( fname, "%s/lpic.steps", p.file_path );

  rf.openinput( fname );

  cells_per_wl     = atoi( rf.getinput( "spl" ) );
  steps_per_period = atoi( rf.getinput( "spp" ) );

  rf.closeinput();

  bob.message("parameter read");

  save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_spacetime::save( parameter &p )
{
  static error_handler bob("input_spacetime::save",errname);
  ofstream outfile;

  outfile.open(p.save_path_name,ios::app);

  outfile << "spacetime" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "t_start      : " << t_start      << endl;
  outfile << "t_stop       : " << t_stop       << endl;
  outfile << "periods_t    : " << periods_t    << endl;
  outfile << "steps_pp     : " << steps_per_period << endl;
  outfile << "x_start      : " << x_start      << endl;
  outfile << "x_stop       : " << x_stop       << endl;
  outfile << "periods_x    : " << periods_x    << endl;
  outfile << "cells_per_wl : " << cells_per_wl << endl;
  outfile << "x_offset     : " << x_offset     << endl;
  outfile << "average      : " << average      << endl;
  outfile << "size         : " << size         << endl;
  outfile << "contour_1    : " << contour_1    << endl;
  outfile << "contour_2    : " << contour_2    << endl;
  outfile << "contour_3    : " << contour_3    << endl;

  outfile.setf(ios::left);

  outfile << "\n   kw    kt    de    di    jx    jy    jz    ex    ey    ez    ed" <<endl;
  outfile << "Q  ";
  outfile << ":" << setw(5) << Q_kw;
  outfile << ":" << setw(5) << Q_kt;
  outfile << ":" << setw(5) << Q_de;
  outfile << ":" << setw(5) << Q_di;
  outfile << ":" << setw(5) << Q_jx;
  outfile << ":" << setw(5) << Q_jy;
  outfile << ":" << setw(5) << Q_jz;
  outfile << ":" << setw(5) << Q_ex;
  outfile << ":" << setw(5) << Q_ey;
  outfile << ":" << setw(5) << Q_ez;
  //  outfile << ":" << setw(5) << Q_bx;
  outfile << ":" << setw(5) << Q_by;
  outfile << ":" << setw(5) << Q_bz;
  outfile << ":" << setw(5) << Q_edens;
  outfile << ":" << setw(5) << Q_de_fi;
  outfile << ":" << setw(5) << Q_de_ii << endl;

  outfile << "C  ";
  outfile << ":" << setw(5) << C_kw;
  outfile << ":" << setw(5) << C_kt;
  outfile << ":" << setw(5) << C_de;
  outfile << ":" << setw(5) << C_di;
  outfile << ":" << setw(5) << C_jx;
  outfile << ":" << setw(5) << C_jy;
  outfile << ":" << setw(5) << C_jz;
  outfile << ":" << setw(5) << C_ex;
  outfile << ":" << setw(5) << C_ey;
  outfile << ":" << setw(5) << C_ez;
  //  outfile << ":" << setw(5) << C_bx;
  outfile << ":" << setw(5) << C_by;
  outfile << ":" << setw(5) << C_bz;
  outfile << ":" << setw(5) << C_edens;
  outfile << ":" << setw(5) << C_de_fi;
  outfile << ":" << setw(5) << C_de_ii << endl << endl;

  outfile << "K_cut : " << K_cut << endl;
  outfile << "W_cut : " << W_cut << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::select( void )
{
  static error_handler bob("spacetime::select_plot",errname);

  if (input.Q_de) {
    xt_kt_kw( "de", input.C_de, 0, 0 );
    write_idl_header_xt(input.C_de,0, "de", "spacetime-de", "!3n!De!N/n!Dc!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "de", "spacetime-kt-de", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "de", "spacetime-kw-de", "");
  }

  if (input.Q_di) {
    xt_kt_kw( "di", input.C_di, 0, 0 );
    write_idl_header_xt(input.C_di,0, "di", "spacetime-di", "!3n!Di!N/n!Dc!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "di", "spacetime-kt-di", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "di", "spacetime-kw-di", "");
  }

  if (input.Q_jx) {
    xt_kt_kw( "jx", input.C_jx, 1, 0 );
    write_idl_header_xt(input.C_jx,1, "jx", "spacetime-jx", "!3j!Dx!N/(e n!Dc!3 c)");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "jx", "spacetime-kt-jx", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "jx", "spacetime-kw-jx", "");
  }

  if (input.Q_jy) {
    xt_kt_kw( "jy", input.C_jy, 1, 0 );
    write_idl_header_xt(input.C_jy,1, "jy", "spacetime-jy", "!3j!Dy!N/(e n!Dc!3 c)");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "jy", "spacetime-kt-jy", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "jy", "spacetime-kw-jy", "");
  }

  if (input.Q_jz) {
    xt_kt_kw( "jz", input.C_jz, 1, 0 );
    write_idl_header_xt(input.C_jz,1, "jz", "spacetime-jz", "!3j!Dz!N/(e n!Dc!3 c)");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "jz", "spacetime-kt-jz", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "jz", "spacetime-kw-jz", "");
  }

  if (input.Q_ex) {
    xt_kt_kw( "ex", input.C_ex, 1, 0 );
    write_idl_header_xt(input.C_ex,1, "ex", "spacetime-ex", "!3E!Dx!N/E!Dr!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt,1, "ex", "spacetime-kt-ex", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0 , "ex", "spacetime-kw-ex", "");
  }

  if (input.Q_ey) {
    xt_kt_kw( "ey", input.C_ey, 1, 0 );
    write_idl_header_xt(input.C_ey,1, "ey", "spacetime-ey", "!3E!Dy!N/E!Dr!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "ey", "spacetime-kt-ey", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "ey", "spacetime-kw-ey", "");
  }

  if (input.Q_ez) {
    xt_kt_kw( "ez", input.C_ez, 1, 0 );
    write_idl_header_xt(input.C_ez,1, "ez", "spacetime-ez", "!3E!Dz!N/E!Dr!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "ez", "spacetime-kt-ez", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "ez", "spacetime-kw-ez", "");
  }

  if (input.Q_bx) {
    xt_kt_kw( "bx", input.C_bx, 1, 0 );
    write_idl_header_xt(input.C_bx,1, "bx", "spacetime-bx", "!3B!Dx!N/B!Dr!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "bx", "spacetime-kt-bx", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "bx", "spacetime-kw-bx", "");
  }

  if (input.Q_by) {
    xt_kt_kw( "by", input.C_by, 1, 0 );
    write_idl_header_xt(input.C_by,1, "by", "spacetime-by", "!3B!Dy!N/B!Dr!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "by", "spacetime-kt-by", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "by", "spacetime-kw-by", "");
  }

  if (input.Q_bz) {
    xt_kt_kw( "bz", input.C_bz, 1, 0 );
    write_idl_header_xt(input.C_bz,1, "bz", "spacetime-bz", "!3B!Dz!N/B!Dr!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "bz", "spacetime-kt-bz", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "bz", "spacetime-kw-bz", "");
  }

  if (input.Q_edens) {
    xt_kt_kw( "edens", input.C_edens, 0, 0 );
    write_idl_header_xt(input.C_edens,0, "edens", "spacetime-edens", "!3W!N/W!Dr!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "edens", "spacetime-kt-edens", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "edens", "spacetime-kw-edens", "");
  }

  if (input.Q_de_fi) {
    xt_kt_kw( "de_fi", input.C_de_fi, 0, 1 );
    write_idl_header_xt(input.C_de_fi,0, "de_fi", "spacetime-de_fi", "!3n!De!N/n!Dc!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "de_fi", "spacetime-kt-de_fi", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "de_fi", "spacetime-kw-de_fi", "");
  }

  if (input.Q_de_ii) {
    xt_kt_kw( "de_ii", input.C_de_ii, 0, 1 );
    write_idl_header_xt(input.C_de_ii,0, "de_ii", "spacetime-de_ii", "!3n!De!N/n!Dc!3");
    if (input.Q_kt) write_idl_header_kt(input.C_kt, 1, "de_ii", "spacetime-kt-de_ii", "");
    if (input.Q_kw) write_idl_header_kw(input.C_kw, 0, "de_ii", "spacetime-kw-de_ii", "");
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::xt_kt_kw( char *unit, float cut, int sign, int scale_write )
{
  static error_handler bob("spacetime::xt_kt_kw",errname);
  char fname_xt[filename_size];
  char fname_kt[filename_size];
  char fname_kw[filename_size];

  sprintf( fname_xt, "spacetime-%s", unit );
  sprintf( fname_kt, "spacetime-kt-%s", unit );
  sprintf( fname_kw, "spacetime-kw-%s", unit );

  read_input_array_size( fname_xt );       // initializes x_steps_in, t_steps_in
                                           //             spp, spl, fnumber

  printf( "reading %s ...",fname_xt ); fflush(stdout);
  matrix_read  =  fmatrix( 0, t_steps_in-1, 0, x_steps_in-1 ); // array dimensions
  read( fname_xt );
  printf( "done\n" );

  printf( "smoothing ... " ); fflush(stdout);
  smooth( matrix_read );
  printf( "done\n" );

  printf( "scaling to byte ... " ); fflush(stdout);
  matrix_write = ucmatrix( 0, t_steps_in-1, 0, x_steps_in-1 );
  scale( cut, sign, matrix_read, matrix_write, t_steps_in, x_steps_in );
  printf( "done\n" );

  printf( "writing to disk ... " ); fflush(stdout);
  write( fname_xt, scale_write );
  bob.message("check point ktkw");
  delete_ucmatrix( matrix_write, 0, t_steps_in-1, 0, x_steps_in-1 );
  printf( "done\n" );

  if (input.Q_kt) {

    printf( "transforming x,t -> k,t ... " ); fflush(stdout);
    kspace         = fmatrix( 0, t_steps_in-1, 0, ft.steps_half-1 );
    transform_k( matrix_read );
    printf( "done\n" );

    printf( "scaling to byte ... " ); fflush(stdout);
    matrix_write   = ucmatrix( 0, t_steps_in-1, 0, ft.steps_half-1 );
    scale( input.C_kt, 0, kspace, matrix_write, t_steps_in, ft.steps_half );
    delete_fmatrix( kspace, 0, t_steps_in-1, 0, ft.steps_half-1 );
    printf( "done\n" );

    printf( "writing to disk ... " ); fflush(stdout);
    write_transform_k( fname_kt, matrix_write );
    delete_ucmatrix( matrix_write, 0, t_steps_in-1, 0, x_steps_in-1 );
    printf( "done\n" );
  }
  if (input.Q_kw) {

    printf( "transforming x,t -> k,w ... " ); fflush(stdout);
    transform_kw( matrix_read );
    printf( "done\n" );

    printf( "scaling to byte ... " ); fflush(stdout);
    matrix_write   = ucmatrix( 0, ft2d.steps_half_1-1, 0, ft2d.steps_2-1 );
    scale( input.C_kw, 0, ft2d.power, matrix_write, ft2d.steps_half_1, ft2d.steps_2 );
    printf( "done\n" );

    printf( "writing to disk ... " ); fflush(stdout);
    write_transform_kw( fname_kw, matrix_write );
    delete_ucmatrix( matrix_write, 0, ft2d.steps_half_1-1, 0, ft2d.steps_2-1);
    printf( "done\n" );
  }
  delete_fmatrix( matrix_read, 0, t_steps_in-1, 0, x_steps_in-1 );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::read_input_array_size( char *unit )
{
  static error_handler bob("spacetime::read_input_array_size",errname);

  FILE  *file;
  char  fname[ filename_size ];
  float x_start, x_stop;
  int   x_steps;
  int   period;

  fnumber    = 0;
  x_steps_in = 0;
  do                          // read all spacetime file headers in order to determine
    {                         // the dimension of the input array
      sprintf( fname, "%s/%s-%d-%d", input_path, unit, fnumber+1, input.t_start );
      file = fopen( fname, "rb" );
      bob.message( "filename = ",fname);
      if (file) {
	fnumber++;

	fread( &period, sizeof(int), 1, file );
	fread( &spp,  sizeof(int), 1, file );
	fread( &x_start, sizeof(float), 1, file );
	fread( &x_stop,  sizeof(float), 1, file );
	fread( &x_steps, sizeof(int), 1, file );

	fclose( file );

	if (x_steps_in==0 && x_steps>0) x_start_in = 1e-6 * floor( 1e6 * x_start );

	if (x_steps>0 ) {
	  spl = (int) floor( (float) x_steps / (x_stop - x_start) + 0.5 );
	  x_stop_in = 1e-6 * floor( 1e6 * x_stop );
	}

	x_steps_in += x_steps;
      }
    }
  while( file );

  if (fnumber==0) bob.error( "no spacetime files found at time", input.t_start );

  bob.message( "found", fnumber, "domain-spacetime file(s) at time", input.t_start );
  bob.message( "spp     =", spp );
  bob.message( "x_start =", x_start_in );
  bob.message( "x_stop  =", x_stop_in );
  bob.message( "x_steps =", x_steps_in );
  bob.message( "spl     =", spl );

  if ( input.x_start < x_start_in ) input.x_start = x_start_in;
  if ( input.x_stop  > x_stop_in  ) input.x_stop  = x_stop_in;

  t_steps_in = (input.t_stop - input.t_start) * spp;
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::read( char *unit )
{
  static error_handler bob("spacetime::read",errname);

  FILE  *file;
  char  fname[ filename_size ];
  float x_start, x_stop;
  int   x_steps;
  int   ti, fi, *x_steps_previous;
  int   period, period_in, spp_in;

  x_steps_previous = new int [ t_steps_in ];

  for( ti=0; ti<t_steps_in; ti++ ) x_steps_previous[ti] = 0;

  for( period=input.t_start; period<input.t_stop; period++ )
    {
      bob.message("now at time ", period);
      for( fi=1; fi<=fnumber; fi++ )                            // now read the files
	{
	  sprintf( fname, "%s/%s-%d-%d", input_path, unit, fi, period );
	  file = fopen( fname, "rb" );
	  if (!file) bob.error( "Cannot open file", fname );
	  else bob.message( "reading file", fname );

	  fread( &period_in, sizeof(int), 1, file );
	  fread( &spp_in,  sizeof(int), 1, file );

	  if (period!=period_in) bob.error( "wrong period in file", fname );
	  if (spp!=spp_in)       bob.error( "wrong number of time steps in file", fname );

	  for( ti=(period-input.t_start)*spp; ti<(period-input.t_start)*spp+spp; ti++ ) {
	    fread( &x_start, sizeof(float), 1, file );
	    fread( &x_stop,  sizeof(float), 1, file );
	    fread( &x_steps, sizeof(int), 1, file );
	    fread( matrix_read[ti] + x_steps_previous[ti], sizeof(float), x_steps, file );
	    (x_steps_previous[ti]) += x_steps;
	  }

	  fclose( file );
	}
    }

  bob.message( "all spacetime data read!" );

  delete x_steps_previous;
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::scale( float cut, int sign,
		       float **m, unsigned char **mw, int nt, int nx )
  // sign = 1 : signed data
  // sign = 0 : unsigned data ( de, di, edens )
{
  static error_handler bob("spacetime::scale",errname);

  int      INTMAX = 255;
  int INTMAX_HALF = 127;
  int overflow=0;
  float data;
  int ti, xi;

  bob.message( "cut    =", cut );

  for( ti=0; ti<nt; ti++ ) {                                // scale and
    for( xi=0; xi<nx; xi++ )  {                             // convert to unsigned char
      if (sign==0) data = (float) fabs( m[ti][xi]/cut * INTMAX );
      else data = (float) (INTMAX_HALF + INTMAX_HALF * m[ti][xi]/cut);
      if (data>INTMAX) { data=INTMAX; overflow++; }
      if (data<0)      { data=0; }
      mw[ti][xi] = (unsigned char) floor( data );
    }
  }

  bob.message( "overflow =", (float)overflow/(nt*nx) );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::smooth( float **field )
{
  static error_handler bob("spacetime::smooth",errname);

  float *s;
  int   offset = (int) floor( 0.5*( input.average - 1 ) + 0.1 );
  int   ti, xi, xj;

  if ( input.average<3 ) return;

  s = new float [x_steps_in];

  for( ti=0; ti<t_steps_in; ti++ ) {                         // forall times

    for( xi=0; xi<x_steps_in; xi++ ) s[xi] = field[ti][xi];

    for( xi=offset; xi<=x_steps_in-offset; xi++ ) {

      field[ti][xi] = 0;
      for( xj=xi-offset; xj<=xi+offset; xj++ ) field[ti][xi] += s[xj] / input.average;
    }
  }

  delete s;

  bob.message( "smooth =", input.average );
  bob.message( "smoothing done" );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write( char *unit, int scale_write )
{
  static error_handler bob("spacetime::write",errname);

  FILE *file;
  char filename[50];
  float t;
  float x, xmin, xmax;
  int tu, to, xu, xo;
  float data;
  int i, j;
  int k, l, flag, max_flag=0;

  if ( input.x_start < x_start_in ) input.x_start = x_start_in;
  if ( input.x_stop  > x_stop_in  ) input.x_stop  = x_stop_in;

  xmin = floor( (input.x_start - x_start_in) * spl + 0.5 );
  xmax = floor( (input.x_stop - x_start_in) * spl + 0.5 );

  bob.message( "t_start_out =", input.t_start );
  bob.message( "t_stop_out  =", input.t_stop );
  bob.message( "x_start_out =", input.x_start );
  bob.message( "x_stop_out  =", input.x_stop );
  bob.message( "spp         =", spp );
  bob.message( "spl         =", spl );
  bob.message( "t_steps_in  =", t_steps_in );
  bob.message( "xmin        =", xmin );
  bob.message( "xmax        =", xmax );

  sprintf( filename, "%s/%s", output_path, unit );
  file = fopen( filename, "wb" );
  if (!file) bob.error( "cannot open file", filename );

  bob.message("check point 1");
  vector_write = new unsigned char [input.size];

  if (scale_write == 1) // de_fi and de_ii
    {
      bob.message("check point 2");

      for( i=0; i<input.size; i++ )                 // determine max_flag
	{
	  for( j=0; j<input.size; j++ )
	    {
	      t = t_steps_in * i/input.size;           if (t<1) t=1.001;
	      x = xmin + (xmax-xmin) * j/input.size;   if (x<1) x=1.001;
	      to = (int) ceil( t );  tu = (int) floor( t );
	      xo = (int) ceil( x );  xu = (int) floor( x );
	      t = (float) to - t;
	      x = (float) xo - x;

	      flag = 0;
	      for( k =  tu - int(ceil(t_steps_in/input.size));
		   k <= to + int(ceil(t_steps_in/input.size)); k++)
		{
		  for( l =  xu - int(ceil((xmax-xmin)/input.size));
		       l <= xo + int(ceil((xmax-xmin)/input.size)); l++)
		    {
		      if ( k >= 0 && k < t_steps_in && l >= 0 && l < x_steps_in &&
			   matrix_write[k][l]>0) flag ++;
		    }
		}
	      max_flag = flag > max_flag ? flag:max_flag;
	    }
	}
      for( i=0; i<input.size; i++ )             // interpolate to output array dimensions
	{
	  for( j=0; j<input.size; j++ )
	    {
	      t = t_steps_in * i/input.size;           if (t<1) t=1.001;
	      x = xmin + (xmax-xmin) * j/input.size;   if (x<1) x=1.001;
	      to = (int) ceil( t );  tu = (int) floor( t );
	      xo = (int) ceil( x );  xu = (int) floor( x );
	      t = (float) to - t;
	      x = (float) xo - x;

	      flag = 0;
	      for( k =  tu - int(ceil(t_steps_in/input.size));
		   k <= to + int(ceil(t_steps_in/input.size)); k++)
		{
		  for( l =  xu - int(ceil((xmax-xmin)/input.size));
		       l <= xo + int(ceil((xmax-xmin)/input.size)); l++)
		    {
		      if ( k >= 0 && k < t_steps_in && l >= 0 && l < x_steps_in &&
			   matrix_write[k][l]>0) flag ++;
		    }
		}
	      if (flag>0) vector_write[j] = (unsigned char) floor(255*flag/max_flag+0.5);
	      else        vector_write[j] = (unsigned char) 0;
	    }
	  fwrite( vector_write, sizeof(unsigned char), input.size, file );
	}
    }

  if (scale_write == 0)  // all the others
    {
      bob.message("check point 3");

      for( i=0; i<input.size; i++ )          // interpolate to output array dimensions
	{
	  for( j=0; j<input.size; j++ )
	    {
	      t = t_steps_in * i/input.size;           if (t<1) t=1.001;
	      x = xmin + (xmax-xmin) * j/input.size;   if (x<1) x=1.001;
	      to = (int) ceil( t );  tu = (int) floor( t );
	      xo = (int) ceil( x );  xu = (int) floor( x );
	      t = (float) to - t;
	      x = (float) xo - x;
	      data =              t*x*matrix_write[tu][xu] +
		            t*(1.0-x)*matrix_write[tu][xo];
	      data +=       (1.0-t)*x*matrix_write[to][xu] +
		      (1.0-t)*(1.0-x)*matrix_write[to][xo];
	      vector_write[j] = (unsigned char) floor( data + 0.5 );
	    }
	  fwrite( vector_write, sizeof(unsigned char), input.size, file );
	}
    }

  fclose( file );

      bob.message("check point 4");
      delete vector_write;
      bob.message("check point 5");
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_idl_header_xt( float cut, int sign, char *idlname, char *unit,
				     char *axislabel )
{
  static error_handler bob("spacetime::write_idl_header",errname);

  FILE          *file;
  int           i, j;
  unsigned char low;

  char fname[filename_size];
  sprintf( fname, "%s/scale-idl", output_path );
  file = fopen( fname, "wb" );

  for( i=0; i<=255; i++ )                                             // write color table
    {
      low = (unsigned char) i;
      for( j=1; j<=30; j++ ) fwrite( &(low), sizeof(unsigned char), 1, file );
    }
  fclose( file );

  sprintf( fname, "%s/idl_%s.header", output_path, idlname );
  file = fopen( fname, "w" );
  fprintf( file, "pro idl_%s", idlname );
  fprintf( file, "\n\n xoffset = %.4f", input.x_offset );
  fprintf( file, "\n\n xmin    = %.2f - xoffset", input.x_start );
  fprintf( file, "\n xmax    = %.2f - xoffset", input.x_stop );
  fprintf( file, "\n xname   = \"!3x/!7k!3!D0!3\"" );
  fprintf( file, "\n ymin    = %d", input.t_start );
  fprintf( file, "\n ymax    = %d", input.t_stop );
  fprintf( file, "\n yname   = \"!3t/!7s!3\"" );
  fprintf( file, "\n zmax    = %.2e", cut );
  if (sign==0) fprintf( file, "\n zmin    = %.2e", 0.0 );
  else         fprintf( file, "\n zmin    = %.2e", -cut );
  fprintf( file, "\n zname   = \"%s\"", axislabel );
  fprintf( file, "\n\n level1  = %.2e", input.contour_1 );
  fprintf( file, "\n level2  = %.2e", input.contour_2 );
  fprintf( file, "\n level3  = %.2e", input.contour_3 );
  fprintf( file, "\n\n dim     = %d", input.size );
  fprintf( file, "\n\n color   = 4" );
  fprintf( file, "\n\n unit    = \"%s\"\n\n", unit );
  fclose( file );

  sprintf( fname, "%s/idl2ps_%s.header", output_path, idlname );
  file = fopen( fname, "w" );
  fprintf( file, "pro idl2ps_%s", idlname );
  fprintf( file, "\n\n xoffset = %.4f", input.x_offset );
  fprintf( file, "\n\n xmin    = %.2f - xoffset", input.x_start );
  fprintf( file, "\n xmax    = %.2f - xoffset", input.x_stop );
  fprintf( file, "\n xname   = \"!3x/!7k!3!D0!3\"" );
  fprintf( file, "\n ymin    = %d", input.t_start );
  fprintf( file, "\n ymax    = %d", input.t_stop );
  fprintf( file, "\n yname   = \"!3t/!7s!3\"" );
  fprintf( file, "\n zmax    = %.2e", cut );
  if (sign==0) fprintf( file, "\n zmin    = %.2e", 0.0 );
  else         fprintf( file, "\n zmin    = %.2e", -cut );
  fprintf( file, "\n zname   = \"%s\"", axislabel );
  fprintf( file, "\n\n level1  = %.2e", input.contour_1 );
  fprintf( file, "\n level2  = %.2e", input.contour_2 );
  fprintf( file, "\n level3  = %.2e", input.contour_3 );
  fprintf( file, "\n\n dim     = %d", input.size );
  fprintf( file, "\n\n color   = 0" );
  fprintf( file, "\n\n unit    = \"%s\"", unit );
  fprintf( file, "\n\n idl2eps = \"%s.eps\"\n\n", unit );
  fclose( file );

}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::transform_k( float **matrix_read )
{
  static error_handler bob("spacetime::transform",errname);

  int i, j;
  int steps = input.periods_x * input.cells_per_wl;

  bob.message( "steps =", (double)steps, " steps_ft =", (double)ft.steps );

  for( i=0; i<t_steps_in; i++ ) {
    ft.RealFt( matrix_read[i] );
    for( j=0; j<ft.steps_half; j++ ) {
      kspace[i][j] = ft.power[j];
      //      kspace[i][j] = ft.co[j];
      //      kspace[i][j] = ft.si[j];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_transform_k( char *unit, unsigned char **m )
{
  static error_handler bob("spacetime::write_transform_k",errname);

  FILE *file;
  char filename[50];
  float t;
  float x, xmin, xmax, kmax;
  int tu, to, xu, xo;
  float data;
  int i, j;
  int nx=ft.steps_half;

  kmax = ft.df * nx;
  input.K_cut = (int) fabs( input.K_cut );

  xmin = 0;

  if (input.K_cut==0)         { xmax = nx; input.K_cut = kmax; }
  else if (input.K_cut<=kmax) { xmax = (int) floor( 1.0*nx*input.K_cut/kmax ); }
  else                        { xmax = nx; input.K_cut = kmax; }

  sprintf( filename, "%s/%s", output_path, unit );
  file = fopen( filename, "wb" );
  if (!file) bob.error( "cannot open file", filename );

  vector_write = new unsigned char [input.size];

  for( i=0; i<input.size; i++ )                 // interpolate to output array dimensions
    {
      for( j=0; j<input.size; j++ )
	{
	  t = t_steps_in * i/input.size;           if (t<1) t=1.001;
	  x = xmin + (xmax-xmin) * j/input.size;   if (x<1) x=1.001;
	  to = (int) ceil( t );  tu = (int) floor( t );
	  xo = (int) ceil( x );  xu = (int) floor( x );
	  t = (float) to - t;
	  x = (float) xo - x;
	  data =        t*x*m[tu][xu] +       t*(1.0-x)*m[tu][xo];
	  data += (1.0-t)*x*m[to][xu] + (1.0-t)*(1.0-x)*m[to][xo];
	  vector_write[j] = (unsigned char) floor( data + 0.5 );
	}
      fwrite( vector_write, sizeof(unsigned char), input.size, file );
    }

  fclose( file );

  delete vector_write;
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_idl_header_kt( float cut, int sign, char *idlname, char *unit,
				     char *axislabel )
{
  static error_handler bob("spacetime::write_idl_header_kt",errname);

  FILE          *file;
  int           i, j;
  unsigned char low;

  char fname[filename_size];
  sprintf( fname, "%s/scale-idl", output_path );
  file = fopen( fname, "wb" );

  for( i=0; i<=255; i++ )                                             // write color table
    {
      low = (unsigned char) i;
      for( j=1; j<=30; j++ ) fwrite( &(low), sizeof(unsigned char), 1, file );
    }
  fclose( file );

  sprintf( fname, "%s/idl_kt_%s.header", output_path, idlname );
  file = fopen( fname, "w" );
  fprintf( file, "pro idl_kt_%s", idlname );
  fprintf( file, "\n\n xmin    = %.2f", 0.0 );
  fprintf( file, "\n xmax    = %.2f", input.K_cut );
  fprintf( file, "\n xname   = \"!3k/!3k!3!D0!3\"" );
  fprintf( file, "\n ymin    = %d", input.t_start );
  fprintf( file, "\n ymax    = %d", input.t_stop );
  fprintf( file, "\n yname   = \"!3t/!7s!3\"" );
  fprintf( file, "\n zmax    = %.2e", cut );
  if (sign==0) fprintf( file, "\n zmin    = %.2e", 0.0 );
  else         fprintf( file, "\n zmin    = %.2e", -cut );
  fprintf( file, "\n zname   = \"%s\"", axislabel );
  fprintf( file, "\n\n level1  = %.2e", input.contour_1 );
  fprintf( file, "\n level2  = %.2e", input.contour_2 );
  fprintf( file, "\n level3  = %.2e", input.contour_3 );
  fprintf( file, "\n\n dim     = %d", input.size );
  fprintf( file, "\n\n color   = 4" );
  fprintf( file, "\n\n unit    = \"%s\"\n\n", unit );
  fclose( file );

  sprintf( fname, "%s/idl2ps_kt_%s.header", output_path, idlname );
  file = fopen( fname, "w" );
  fprintf( file, "pro idl2ps_kt_%s", idlname );
  fprintf( file, "\n\n xmin    = %.2f", 0.0 );
  fprintf( file, "\n xmax    = %.2f", input.K_cut );
  fprintf( file, "\n xname   = \"!3k/!3k!3!D0!3\"" );
  fprintf( file, "\n ymin    = %d", input.t_start );
  fprintf( file, "\n ymax    = %d", input.t_stop );
  fprintf( file, "\n yname   = \"!3t/!7s!3\"" );
  fprintf( file, "\n zmax    = %.2e", cut );
  if (sign==0) fprintf( file, "\n zmin    = %.2e", 0.0 );
  else         fprintf( file, "\n zmin    = %.2e", -cut );
  fprintf( file, "\n zname   = \"%s\"", axislabel );
  fprintf( file, "\n\n level1  = %.2e", input.contour_1 );
  fprintf( file, "\n level2  = %.2e", input.contour_2 );
  fprintf( file, "\n level3  = %.2e", input.contour_3 );
  fprintf( file, "\n\n dim     = %d", input.size );
  fprintf( file, "\n\n color   = 0" );
  fprintf( file, "\n\n unit    = \"%s\"", unit );
  fprintf( file, "\n\n idl2eps = \"%s.eps\"\n\n", unit );
  fclose( file );

}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::transform_kw( float **matrix_read )
{
  static error_handler bob("spacetime::transform_kw",errname);

  int i, j;
  int steps_x = input.periods_x * input.cells_per_wl;
  int steps_t = input.periods_t * input.steps_per_period;

  bob.message( "steps =", (double)steps_x, " steps_ft =", (double)ft2d.steps_2 );
  bob.message( "steps_t", (double)steps_t, " steps_ft =", (double)ft2d.steps_1 );

  ft2d.RealFt( matrix_read );

  FILE *f;
  f=fopen( "ft2d.dat", "w" );
  for( i=0; i<ft2d.steps_half_1; i+= (int) floor(0.1 * ft2d.steps_half_1) ) {
    fprintf( f, "\n" );
    for( j=0; j<ft2d.steps_half_2; j+=10 ) {
      fprintf( f, "\n %.4e %.3e", ft2d.frequency_2[j],  ft2d.power[i][j] );
    }
  }
  fclose( f );
  // result in ft2d.power[i][j], ft2d.co[i][j], ft2d.si[i][j]
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_transform_kw( char *unit, unsigned char **m )
{
  static error_handler bob("spacetime::write_transform_kw",errname);

  FILE  *file;
  char  filename[50];
  float t, tmin, tmax, wmax;
  float x, xmin, xmax, kmax;
  int   tu, to, xu, xo;
  float data;
  int   i, j;

  wmax = ft2d.df_1 * ft2d.steps_half_1;
  kmax = ft2d.df_2 * ft2d.steps_half_2;

  printf( "wmax = %.3f\n", wmax );
  printf( "kmax = %.3f\n", kmax );

  if (input.K_cut>0 && input.K_cut<=kmax) {
    xmin = ft2d.steps_half_2 - (int) floor(1.0*input.K_cut/ft2d.df_2);
    xmax = ft2d.steps_half_2 + (int) floor(1.0*input.K_cut/ft2d.df_2);
  }
  else {
    xmin = 0;
    xmax = ft2d.steps_2;
  }

  if (input.W_cut>0 && input.W_cut<=wmax) {
    tmin = 0;
    tmax = (int) floor( 1.0*input.W_cut/ft2d.df_1 );
  }
  else {
    tmin = 0;
    tmax = ft2d.steps_half_1;
  }

  sprintf( filename, "%s/%s", output_path, unit );
  file = fopen( filename, "wb" );
  if (!file) bob.error( "cannot open file", filename );

  vector_write = new unsigned char [input.size];

  for( i=0; i<input.size; i++ )                 // interpolate to output array dimensions
    {
      for( j=0; j<input.size; j++ )
	{
	  t = tmin + (tmax-tmin) * i/input.size;   if (t<1) t=1.001;
	  x = xmin + (xmax-xmin) * j/input.size;   if (x<1) x=1.001;
	  to = (int) ceil( t );  tu = (int) floor( t );
	  xo = (int) ceil( x );  xu = (int) floor( x );
	  t = (float) to - t;
	  x = (float) xo - x;
	  data =        t*x*m[tu][xu] +       t*(1.0-x)*m[tu][xo];
	  data += (1.0-t)*x*m[to][xu] + (1.0-t)*(1.0-x)*m[to][xo];
	  vector_write[j] = (unsigned char) floor( data + 0.5 );
	}
      fwrite( vector_write, sizeof(unsigned char), input.size, file );
    }

  fclose( file );

  delete vector_write;
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_idl_header_kw( float cut, int sign, char *idlname, char *unit,
				     char *axislabel )
{
  static error_handler bob("spacetime::write_idl_header_kw",errname);

  FILE          *file;
  int           i, j;
  unsigned char low;

  char fname[filename_size];
  sprintf( fname, "%s/scale-idl", output_path );
  file = fopen( fname, "wb" );

  for( i=0; i<=255; i++ )                                             // write color table
    {
      low = (unsigned char) i;
      for( j=1; j<=30; j++ ) fwrite( &(low), sizeof(unsigned char), 1, file );
    }
  fclose( file );

  sprintf( fname, "%s/idl_kw_%s.header", output_path, idlname );
  file = fopen( fname, "w" );
  fprintf( file, "pro idl_kw_%s", idlname );
  fprintf( file, "\n\n xmin    = %.2f", -input.K_cut );
  fprintf( file, "\n xmax    = %.2f", input.K_cut );
  fprintf( file, "\n xname   = \"!3k/!3k!3!D0!3\"" );
  fprintf( file, "\n ymin    = %.3f", 0.0 );
  fprintf( file, "\n ymax    = %.3f", input.W_cut );
  fprintf( file, "\n yname   = \"!7x/!7x!3!D0\"" );
  fprintf( file, "\n zmax    = %.2e", cut );
  if (sign==0) fprintf( file, "\n zmin    = %.2e", 0.0 );
  else         fprintf( file, "\n zmin    = %.2e", -cut );
  fprintf( file, "\n zname   = \"%s\"", axislabel );
  fprintf( file, "\n\n level1  = %.2e", input.contour_1 );
  fprintf( file, "\n level2  = %.2e", input.contour_2 );
  fprintf( file, "\n level3  = %.2e", input.contour_3 );
  fprintf( file, "\n\n dim     = %d", input.size );
  fprintf( file, "\n\n color   = 4" );
  fprintf( file, "\n\n unit    = \"%s\"\n\n", unit );
  fclose( file );

  sprintf( fname, "%s/idl2ps_kw_%s.header", output_path, idlname );
  file = fopen( fname, "w" );
  fprintf( file, "pro idl2ps_kw_%s", idlname );
  fprintf( file, "\n\n xmin    = %.2f", -input.K_cut );
  fprintf( file, "\n xmax    = %.2f", input.K_cut );
  fprintf( file, "\n xname   = \"!3k/!3k!3!D0!3\"" );
  fprintf( file, "\n ymin    = %.3f", 0.0 );
  fprintf( file, "\n ymax    = %.3f", input.W_cut );
  fprintf( file, "\n yname   = \"!7x/!7x!3!D0\"" );
  fprintf( file, "\n zmax    = %.2e", cut );
  if (sign==0) fprintf( file, "\n zmin    = %.2e", 0.0 );
  else         fprintf( file, "\n zmin    = %.2e", -cut );
  fprintf( file, "\n zname   = \"%s\"", axislabel );
  fprintf( file, "\n\n level1  = %.2e", input.contour_1 );
  fprintf( file, "\n level2  = %.2e", input.contour_2 );
  fprintf( file, "\n level3  = %.2e", input.contour_3 );
  fprintf( file, "\n\n dim     = %d", input.size );
  fprintf( file, "\n\n color   = 0" );
  fprintf( file, "\n\n unit    = \"%s\"", unit );
  fprintf( file, "\n\n idl2eps = \"%s.eps\"\n\n", unit );
  fclose( file );

}


//////////////////////////////////////////////////////////////////////////////////////////
//eof









