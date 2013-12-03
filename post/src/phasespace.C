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
#include <phasespace.h>

phasespace::phasespace( parameter &p )
  : input(p)
{
  sprintf( errname, "%s/error", p.output_path );
  static error_handler bob("phasespace::Constructor",errname);

  dim = 399;  // bins 0...399

  period_start = input.period_start;
  period_stop  = input.period_stop;
  period_step  = input.period_step;
  xmax         = input.xmax;
  xoffset      = input.xoffset;

  Q_el  = input.Q_el;
  Q_ion = input.Q_ion;

  Q_vx = input.Q_vx;
  Q_vy = input.Q_vy;
  Q_vz = input.Q_vz;

  input_path = new char [ filename_size ];
  strcpy(input_path,p.file_path);
  output_path = new char [ filename_size ];
  strcpy(output_path,p.output_path);

  matrix_read  = ucmatrix(0,dim,0,dim);
  matrix_inter = imatrix(0,dim,0,dim);
  matrix_write = ucmatrix(0,dim,0,dim);
}


//////////////////////////////////////////////////////////////////////////////////////////


input_phasespace::input_phasespace( parameter &p )
  : rf()
{
  strcpy( errname, p.errname );
  static error_handler bob("input_phasespace::Constructor",errname);

  rf.openinput( p.read_filename );

  period_start = atof( rf.setget( "&phasespace", "period_start" ) );
  period_stop  = atof( rf.setget( "&phasespace", "period_stop" ) );
  period_step  = atof( rf.setget( "&phasespace", "period_step" ) );
  xmax         = atof( rf.setget( "&phasespace", "xmax" ) );
  xoffset      = atof( rf.setget( "&phasespace", "xoffset" ) );
  Q_el         = atoi( rf.setget( "&phasespace", "Q_el" ) );
  Q_ion        = atoi( rf.setget( "&phasespace", "Q_ion" ) );
  Q_vx         = atoi( rf.setget( "&phasespace", "Q_vx" ) );
  Q_vy         = atoi( rf.setget( "&phasespace", "Q_vy" ) );
  Q_vz         = atoi( rf.setget( "&phasespace", "Q_vz" ) );

  rf.closeinput();

  bob.message("parameter read");

  save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_phasespace::save( parameter &p )
{
  static error_handler bob("input_phasespace::save",errname);
  ofstream outfile;

  outfile.open(p.save_path_name,ios::app);

  outfile << "Phasespace - Plots" << endl;
  outfile << "--------------------------------------------------" << endl;
  outfile << "period_start :" << period_start << endl;
  outfile << "period_stop  :" << period_stop  << endl;
  outfile << "period_step  :" << period_step  << endl;
  outfile << "xmax         :" << xmax         << endl;
  outfile << "xoffset      :" << xoffset      << endl;
  outfile << "Q_el         :" << Q_el << endl;
  outfile << "Q_ion        :" << Q_ion << endl;
  outfile << "Q_vx         :" << Q_vx << endl;
  outfile << "Q_vy         :" << Q_vy << endl;
  outfile << "Q_vz         :" << Q_vz << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void phasespace::concat( void )
{
  static error_handler bob("phasespace::select_plot",errname);  //##
  double time;
  int files;

  if (Q_el){
    for( time=period_start; time<=period_stop; time+=period_step ) {

      if (Q_vx) {
	files=read("phasex", "sp0",time);
	if (files>0) { write("phasex", "sp0",time);     write_idl_header("sp0","x"); }
      }
      if (Q_vy) {
	files=read("phasey", "sp0",time);
	if (files>0) { write("phasey", "sp0",time);     write_idl_header("sp0","y"); }
      }
      if (Q_vz) {
	files=read("phasez", "sp0",time);
	if (files>0) { write("phasez", "sp0",time);     write_idl_header("sp0","z"); }
      }
    }
  }

  if (Q_ion){
    for( time=period_start; time<=period_stop; time+=period_step ) {

      if (Q_vx) {
	files=read("phasex", "sp1",time);
	if (files>0) { write("phasex", "sp1",time);     write_idl_header("sp1","x"); }
      }
      if (Q_vy) {
	files=read("phasey", "sp1",time);
	if (files>0) { write("phasey", "sp1",time);     write_idl_header("sp1","y"); }
      }
      if (Q_vz) {
	files=read("phasez", "sp1",time);
	if (files>0) { write("phasez", "sp1",time);     write_idl_header("sp1","z"); }
      }
    }
  }

  delete_ucmatrix( matrix_read, 0, dim, 0, dim );
  delete_ucmatrix( matrix_write, 0, dim, 0, dim );
  delete_imatrix( matrix_inter, 0, dim, 0, dim );

}


//////////////////////////////////////////////////////////////////////////////////////////


int phasespace::read( char *unit, char *spec, double time )
{
  static error_handler bob("phasespace::read",errname);

  FILE *file;
  char fname[ filename_size ];
  int  fnumber = 0;
  int  INTMAX = 255;
  int  dim1, dim2;
  int  vi, xi;
  int  file_open;
  int  overflow = 0;

  for( vi=0; vi<=dim; vi++ )
    for( xi=0; xi<=dim; xi++ )
      matrix_write[vi][xi]=matrix_inter[vi][xi]=0;

  do
    {
      sprintf( fname, "%s/%s-%d-%s-%.3f", input_path, unit, fnumber+1, spec, time );
      file = fopen( fname, "rb" );
      if (!file) file_open=0;
      else       file_open=1;

      for( vi=0; vi<=dim; vi++ )
	for( xi=0; xi<=dim; xi++ )
	  matrix_read[vi][xi]=0;

      if (file_open) {
	fnumber++;

	fread( &dim1, sizeof(int), 1, file );
	fread( &dim2, sizeof(int), 1, file );

	for( vi=0; vi<=dim; vi++ )
	  fread( matrix_read[vi] + dim1, sizeof(unsigned char), dim2-dim1+1, file );

	for( vi=0; vi<=dim; vi++ ) {
	  for( xi=dim1; xi<=dim2; xi++ ) {

	    matrix_inter[vi][xi] += (int) matrix_read[vi][xi];

	    if (matrix_inter[vi][xi] > INTMAX) {
	      matrix_inter[vi][xi] = INTMAX;
	      overflow++;
	    }
	  }
	}

	fclose( file );
      }
    }
  while( file_open );

  for( vi=0; vi<=dim; vi++ ) {
    for( xi=0; xi<=dim; xi++ ) {
      matrix_write[vi][xi] = (unsigned char) matrix_inter[vi][xi];
    }
  }

  if (fnumber==0) bob.message( "no phasespace file found at time", time );
  else bob.message( "found", fnumber, "phasespace file(s) at time", time );

  bob.message( "overflow =", (float) overflow/(dim*dim) );

  return fnumber;
}


//////////////////////////////////////////////////////////////////////////////////////////


void phasespace::write( char *unit, char *spec, double time )
{
  static error_handler bob("spacetime::write",errname);

  FILE *file;
  char filename[filename_size];
  int vi, xi;
  unsigned char low;

  sprintf( filename, "%s/%s-%s-%.3f", output_path, unit, spec, time );
  file = fopen( filename, "wb" );
  if (!file) bob.error( "cannot open file", filename );

  for( vi=0; vi<=dim; vi++ )
    fwrite( matrix_write[vi], sizeof(unsigned char), dim+1, file );

  fclose( file );

  sprintf( filename, "%s/scale-idl", output_path );
  file = fopen( filename, "wb" );

  for( vi=0; vi<=255; vi++ )                                          // write color table
    {
      low = (unsigned char) vi;
      for( xi=1; xi<=30; xi++ ) fwrite( &(low), sizeof(unsigned char), 1, file );
    }

  fclose( file );
  bob.message("color table written");
}

//////////////////////////////////////////////////////////////////////////////////////////


void phasespace::write_idl_header( char *spec, char *direct )
{
  FILE* file;
  char fname[filename_size];

  sprintf( fname, "%s/idlmovie_%s_%s.header", output_path, direct, spec );
  file = fopen( fname, "w" );
  fprintf( file, "pro idlmovie_%s_%s", direct, spec );
  fprintf( file, "\n\n file0        = \"phase%s-%s-\"", direct, spec );
  fprintf( file, "\n file_begin   = %.3f", period_start );
  fprintf( file, "\n file_end     = %.3f", period_stop );
  fprintf( file, "\n increment    = %.3f", period_step );
  fprintf( file, "\n\n special_file = -100" );
  fprintf( file, "\n\n dimx0   = 400" );
  fprintf( file, "\n dimv0   = 400" );
  fprintf( file, "\n xmax0   = %.3f", xmax );
  fprintf( file, "\n xoffset = %.3f", xoffset );
  fprintf( file, "\n vmax0   = 1.0" );
  fprintf( file, "\n\n cutx   = 0" );
  fprintf( file, "\n dimx   = 400" );
  fprintf( file, "\n cutv   = 0" );
  fprintf( file, "\n dimv   = dimv0 - 2*cutv \n" );
  fclose( file );

}


//////////////////////////////////////////////////////////////////////////////////////////
//eof





