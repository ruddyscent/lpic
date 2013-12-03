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

#include <diagnostic_spacetime.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////


spacetime::spacetime( parameter &p )
  : rf(),
    input(p),
    stepper_de( input.stepper_de, p ),
    stepper_di( input.stepper_di, p ),
    stepper_jx( input.stepper_jx, p ),
    stepper_jy( input.stepper_jy, p ),
    stepper_jz( input.stepper_jz, p ),
    stepper_ex( input.stepper_ex, p ),
    stepper_ey( input.stepper_ey, p ),
    stepper_ez( input.stepper_ez, p ),
    stepper_bx( input.stepper_bx, p ),
    stepper_by( input.stepper_by, p ),
    stepper_bz( input.stepper_bz, p ),
    stepper_edens( input.stepper_edens, p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("spacetime::Constructor",errname);

  name_de = new (char [filename_size]);
  name_di = new (char [filename_size]);
  name_jx = new (char [filename_size]);
  name_jy = new (char [filename_size]);
  name_jz = new (char [filename_size]);
  name_ex = new (char [filename_size]);
  name_ey = new (char [filename_size]);
  name_ez = new (char [filename_size]);
  name_bx = new (char [filename_size]);
  name_by = new (char [filename_size]);
  name_bz = new (char [filename_size]);
  name_edens = new (char [filename_size]);

  stepper_de.t_start += 1;  stepper_de.t_stop += 1;
  stepper_di.t_start += 1;  stepper_di.t_stop += 1;
  stepper_ex.t_start += 1;  stepper_ex.t_stop += 1;
  stepper_ey.t_start += 1;  stepper_ey.t_stop += 1;
  stepper_ez.t_start += 1;  stepper_ez.t_stop += 1;
  stepper_bx.t_start += 1;  stepper_bx.t_stop += 1;
  stepper_by.t_start += 1;  stepper_by.t_stop += 1;
  stepper_bz.t_start += 1;  stepper_bz.t_stop += 1;
  stepper_jx.t_start += 1;  stepper_jx.t_stop += 1;
  stepper_jy.t_start += 1;  stepper_jy.t_stop += 1;
  stepper_jz.t_start += 1;  stepper_jz.t_stop += 1;
  stepper_edens.t_start += 1;  stepper_edens.t_stop += 1;

  if ( input.Q_restart == 0 ){
    output_period_de = (int) floor( input.stepper_de.t_start - 1 + 0.5 );
    output_period_di = (int) floor( input.stepper_di.t_start - 1 + 0.5 );
    output_period_jx = (int) floor( input.stepper_jx.t_start - 1 + 0.5 );
    output_period_jy = (int) floor( input.stepper_jy.t_start - 1 + 0.5 );
    output_period_jz = (int) floor( input.stepper_jz.t_start - 1 + 0.5 );
    output_period_ex = (int) floor( input.stepper_ex.t_start - 1 + 0.5 );
    output_period_ey = (int) floor( input.stepper_ey.t_start - 1 + 0.5 );
    output_period_ez = (int) floor( input.stepper_ez.t_start - 1 + 0.5 );
    output_period_bx = (int) floor( input.stepper_bx.t_start - 1 + 0.5 );
    output_period_by = (int) floor( input.stepper_by.t_start - 1 + 0.5 );
    output_period_bz = (int) floor( input.stepper_bz.t_start - 1 + 0.5 );
    output_period_edens = (int) floor( input.stepper_edens.t_start - 1 + 0.5 );
  }
  else {
    char fname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    stepper_de.t_count = atoi( rf.getinput( "spa.stepper_de.t_count" ) );
    stepper_di.t_count = atoi( rf.getinput( "spa.stepper_di.t_count" ) );
    stepper_jx.t_count = atoi( rf.getinput( "spa.stepper_jx.t_count" ) );
    stepper_jy.t_count = atoi( rf.getinput( "spa.stepper_jy.t_count" ) );
    stepper_jz.t_count = atoi( rf.getinput( "spa.stepper_jz.t_count" ) );
    stepper_ex.t_count = atoi( rf.getinput( "spa.stepper_ex.t_count" ) );
    stepper_ey.t_count = atoi( rf.getinput( "spa.stepper_ey.t_count" ) );
    stepper_ez.t_count = atoi( rf.getinput( "spa.stepper_ez.t_count" ) );
    stepper_bx.t_count = atoi( rf.getinput( "spa.stepper_bx.t_count" ) );
    stepper_by.t_count = atoi( rf.getinput( "spa.stepper_by.t_count" ) );
    stepper_bz.t_count = atoi( rf.getinput( "spa.stepper_bz.t_count" ) );
    stepper_edens.t_count = atoi( rf.getinput( "spa.stepper_edens.t_count" ) );

    output_period_de = atoi( rf.getinput( "spa.output_period_de" ) );
    output_period_di = atoi( rf.getinput( "spa.output_period_di" ) );
    output_period_jx = atoi( rf.getinput( "spa.output_period_jx" ) );
    output_period_jy = atoi( rf.getinput( "spa.output_period_jy" ) );
    output_period_jz = atoi( rf.getinput( "spa.output_period_jz" ) );
    output_period_ex = atoi( rf.getinput( "spa.output_period_ex" ) );
    output_period_ey = atoi( rf.getinput( "spa.output_period_ey" ) );
    output_period_ez = atoi( rf.getinput( "spa.output_period_ez" ) );
    output_period_bx = atoi( rf.getinput( "spa.output_period_bx" ) );
    output_period_by = atoi( rf.getinput( "spa.output_period_by" ) );
    output_period_bz = atoi( rf.getinput( "spa.output_period_bz" ) );
    output_period_edens = atoi( rf.getinput( "spa.output_period_edens" ) );

    rf.closeinput();
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


input_spacetime::input_spacetime( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_spacetime::Constructor",errname);

  rf.openinput( p.input_file_name );

  stepper_de.Q         = atoi( rf.setget( "&de", "Q" ) );
  stepper_de.t_start   = atof( rf.setget( "&de", "t_start" ) );
  stepper_de.t_stop    = atof( rf.setget( "&de", "t_stop" ) );
  stepper_de.t_step    = 0;
  stepper_de.x_start   = atoi( rf.setget( "&de", "x_start" ) );
  stepper_de.x_stop    = atoi( rf.setget( "&de", "x_stop" ) );
  stepper_de.x_step    = 1;       // not used

  stepper_di.Q         = atoi( rf.setget( "&di", "Q" ) );
  stepper_di.t_start   = atof( rf.setget( "&di", "t_start" ) );
  stepper_di.t_stop    = atof( rf.setget( "&di", "t_stop" ) );
  stepper_di.t_step    = 0;
  stepper_di.x_start   = atoi( rf.setget( "&di", "x_start" ) );
  stepper_di.x_stop    = atoi( rf.setget( "&di", "x_stop" ) );
  stepper_di.x_step    = 1;       // not used

  stepper_jx.Q         = atoi( rf.setget( "&jx", "Q" ) );
  stepper_jx.t_start   = atof( rf.setget( "&jx", "t_start" ) );
  stepper_jx.t_stop    = atof( rf.setget( "&jx", "t_stop" ) );
  stepper_jx.t_step    = 0;
  stepper_jx.x_start   = atoi( rf.setget( "&jx", "x_start" ) );
  stepper_jx.x_stop    = atoi( rf.setget( "&jx", "x_stop" ) );
  stepper_jx.x_step    = 1;       // not used

  stepper_jy.Q         = atoi( rf.setget( "&jy", "Q" ) );
  stepper_jy.t_start   = atof( rf.setget( "&jy", "t_start" ) );
  stepper_jy.t_stop    = atof( rf.setget( "&jy", "t_stop" ) );
  stepper_jy.t_step    = 0;
  stepper_jy.x_start   = atoi( rf.setget( "&jy", "x_start" ) );
  stepper_jy.x_stop    = atoi( rf.setget( "&jy", "x_stop" ) );
  stepper_jy.x_step    = 1;       // not used

  stepper_jz.Q         = atoi( rf.setget( "&jz", "Q" ) );
  stepper_jz.t_start   = atof( rf.setget( "&jz", "t_start" ) );
  stepper_jz.t_stop    = atof( rf.setget( "&jz", "t_stop" ) );
  stepper_jz.t_step    = 0;
  stepper_jz.x_start   = atoi( rf.setget( "&jz", "x_start" ) );
  stepper_jz.x_stop    = atoi( rf.setget( "&jz", "x_stop" ) );
  stepper_jz.x_step    = 1;      // not used

  stepper_ex.Q         = atoi( rf.setget( "&ex", "Q" ) );
  stepper_ex.t_start   = atof( rf.setget( "&ex", "t_start" ) );
  stepper_ex.t_stop    = atof( rf.setget( "&ex", "t_stop" ) );
  stepper_ex.t_step    = 0;
  stepper_ex.x_start   = atoi( rf.setget( "&ex", "x_start" ) );
  stepper_ex.x_stop    = atoi( rf.setget( "&ex", "x_stop" ) );
  stepper_ex.x_step    = 1;   // not used

  stepper_ey.Q         = atoi( rf.setget( "&ey", "Q" ) );
  stepper_ey.t_start   = atof( rf.setget( "&ey", "t_start" ) );
  stepper_ey.t_stop    = atof( rf.setget( "&ey", "t_stop" ) );
  stepper_ey.t_step    = 0;
  stepper_ey.x_start   = atoi( rf.setget( "&ey", "x_start" ) );
  stepper_ey.x_stop    = atoi( rf.setget( "&ey", "x_stop" ) );
  stepper_ey.x_step    = 1;    // not used

  stepper_ez.Q         = atoi( rf.setget( "&ez", "Q" ) );
  stepper_ez.t_start   = atof( rf.setget( "&ez", "t_start" ) );
  stepper_ez.t_stop    = atof( rf.setget( "&ez", "t_stop" ) );
  stepper_ez.t_step    = 0;
  stepper_ez.x_start   = atoi( rf.setget( "&ez", "x_start" ) );
  stepper_ez.x_stop    = atoi( rf.setget( "&ez", "x_stop" ) );
  stepper_ez.x_step    = 1;   // not used

  stepper_bx.Q         = atoi( rf.setget( "&bx", "Q" ) );
  stepper_bx.t_start   = atof( rf.setget( "&bx", "t_start" ) );
  stepper_bx.t_stop    = atof( rf.setget( "&bx", "t_stop" ) );
  stepper_bx.t_step    = 0;
  stepper_bx.x_start   = atoi( rf.setget( "&bx", "x_start" ) );
  stepper_bx.x_stop    = atoi( rf.setget( "&bx", "x_stop" ) );
  stepper_bx.x_step    = 1;   // not used

  stepper_by.Q         = atoi( rf.setget( "&by", "Q" ) );
  stepper_by.t_start   = atof( rf.setget( "&by", "t_start" ) );
  stepper_by.t_stop    = atof( rf.setget( "&by", "t_stop" ) );
  stepper_by.t_step    = 0;
  stepper_by.x_start   = atoi( rf.setget( "&by", "x_start" ) );
  stepper_by.x_stop    = atoi( rf.setget( "&by", "x_stop" ) );
  stepper_by.x_step    = 1;    // not used

  stepper_bz.Q         = atoi( rf.setget( "&bz", "Q" ) );
  stepper_bz.t_start   = atof( rf.setget( "&bz", "t_start" ) );
  stepper_bz.t_stop    = atof( rf.setget( "&bz", "t_stop" ) );
  stepper_bz.t_step    = 0;
  stepper_bz.x_start   = atoi( rf.setget( "&bz", "x_start" ) );
  stepper_bz.x_stop    = atoi( rf.setget( "&bz", "x_stop" ) );
  stepper_bz.x_step    = 1;   // not used

  stepper_edens.Q         = atoi( rf.setget( "&edens", "Q" ) );
  stepper_edens.t_start   = atof( rf.setget( "&edens", "t_start" ) );
  stepper_edens.t_stop    = atof( rf.setget( "&edens", "t_stop" ) );
  stepper_edens.t_step    = 0;
  stepper_edens.x_start   = atoi( rf.setget( "&edens", "x_start" ) );
  stepper_edens.x_stop    = atoi( rf.setget( "&edens", "x_stop" ) );
  stepper_edens.x_step    = 1;   // not used

  Q_restart       = atoi( rf.setget( "&restart", "Q"     ) );
  strcpy( restart_file, rf.setget( "&restart", "file"    ) );

  rf.closeinput();

  bob.message("parameter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_spacetime::save( parameter &p )
{
  static error_handler bob("input_spacetime::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "diagnostic spacetime" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "de:" << endl;
  outfile << "   Q             : " << stepper_de.Q       << endl;
  outfile << "   t_start       : " << stepper_de.t_start << endl;
  outfile << "   t_stop        : " << stepper_de.t_stop  << endl;
  outfile << "   x_start       : " << stepper_de.x_start << endl;
  outfile << "   x_stop        : " << stepper_de.x_stop  << endl;
  outfile << "di:" << endl;
  outfile << "   Q             : " << stepper_di.Q       << endl;
  outfile << "   t_start       : " << stepper_di.t_start << endl;
  outfile << "   t_stop        : " << stepper_di.t_stop  << endl;
  outfile << "   x_start       : " << stepper_di.x_start << endl;
  outfile << "   x_stop        : " << stepper_di.x_stop  << endl;
  outfile << "jx:" << endl;
  outfile << "   Q             : " << stepper_jx.Q       << endl;
  outfile << "   t_start       : " << stepper_jx.t_start << endl;
  outfile << "   t_stop        : " << stepper_jx.t_stop  << endl;
  outfile << "   x_start       : " << stepper_jx.x_start << endl;
  outfile << "   x_stop        : " << stepper_jx.x_stop  << endl;
  outfile << "jy:" << endl;
  outfile << "   Q             : " << stepper_jy.Q       << endl;
  outfile << "   t_start       : " << stepper_jy.t_start << endl;
  outfile << "   t_stop        : " << stepper_jy.t_stop  << endl;
  outfile << "   x_start       : " << stepper_jy.x_start << endl;
  outfile << "   x_stop        : " << stepper_jy.x_stop  << endl;
  outfile << "jz:" << endl;
  outfile << "   Q             : " << stepper_jz.Q       << endl;
  outfile << "   t_start       : " << stepper_jz.t_start << endl;
  outfile << "   t_stop        : " << stepper_jz.t_stop  << endl;
  outfile << "   x_start       : " << stepper_jz.x_start << endl;
  outfile << "   x_stop        : " << stepper_jz.x_stop  << endl;
  outfile << "ex:" << endl;
  outfile << "   Q             : " << stepper_ex.Q       << endl;
  outfile << "   t_start       : " << stepper_ex.t_start << endl;
  outfile << "   t_stop        : " << stepper_ex.t_stop  << endl;
  outfile << "   x_start       : " << stepper_ex.x_start << endl;
  outfile << "   x_stop        : " << stepper_ex.x_stop  << endl;
  outfile << "ey:" << endl;
  outfile << "   Q             : " << stepper_ey.Q       << endl;
  outfile << "   t_start       : " << stepper_ey.t_start << endl;
  outfile << "   t_stop        : " << stepper_ey.t_stop  << endl;
  outfile << "   x_start       : " << stepper_ey.x_start << endl;
  outfile << "   x_stop        : " << stepper_ey.x_stop  << endl;
  outfile << "ez:" << endl;
  outfile << "   Q             : " << stepper_ez.Q       << endl;
  outfile << "   t_start       : " << stepper_ez.t_start << endl;
  outfile << "   t_stop        : " << stepper_ez.t_stop  << endl;
  outfile << "   x_start       : " << stepper_ez.x_start << endl;
  outfile << "   x_stop        : " << stepper_ez.x_stop  << endl;
  outfile << "bx:" << endl;
  outfile << "   Q             : " << stepper_bx.Q       << endl;
  outfile << "   t_start       : " << stepper_bx.t_start << endl;
  outfile << "   t_stop        : " << stepper_bx.t_stop  << endl;
  outfile << "   x_start       : " << stepper_bx.x_start << endl;
  outfile << "   x_stop        : " << stepper_bx.x_stop  << endl;
  outfile << "by:" << endl;
  outfile << "   Q             : " << stepper_by.Q       << endl;
  outfile << "   t_start       : " << stepper_by.t_start << endl;
  outfile << "   t_stop        : " << stepper_by.t_stop  << endl;
  outfile << "   x_start       : " << stepper_by.x_start << endl;
  outfile << "   x_stop        : " << stepper_by.x_stop  << endl;
  outfile << "bz:" << endl;
  outfile << "   Q             : " << stepper_bz.Q       << endl;
  outfile << "   t_start       : " << stepper_bz.t_start << endl;
  outfile << "   t_stop        : " << stepper_bz.t_stop  << endl;
  outfile << "   x_start       : " << stepper_bz.x_start << endl;
  outfile << "   x_stop        : " << stepper_bz.x_stop  << endl;
  outfile << "edens:" << endl;
  outfile << "   Q             : " << stepper_edens.Q       << endl;
  outfile << "   t_start       : " << stepper_edens.t_start << endl;
  outfile << "   t_stop        : " << stepper_edens.t_stop  << endl;
  outfile << "   x_start       : " << stepper_edens.x_start << endl;
  outfile << "   x_stop        : " << stepper_edens.x_stop  << endl;
  outfile << "Q_restart        : " << Q_restart       << endl;
  outfile << "restart_file     : " << restart_file    << endl << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::boundaries( float *x_start, float *x_stop, int *x_steps,
			    diagnostic_stepper *stepper, domain *grid )
{
  static error_handler bob("spacetime::boundaries",errname);
  int x1, x2;

  if (stepper->x_start < grid->left->number) x1 = grid->left->number;
  else if (stepper->x_start <= grid->right->number )
                                             x1 = stepper->x_start;
  else                                       x1 = -1;

  if (stepper->x_stop < grid->left->number)  x2 = -1;
  else if (stepper->x_stop <= grid->right->number )
                                             x2 = stepper->x_stop;
  else                                       x2 = grid->right->number;

  if (x2==-1 || x1==-1) {
    *x_steps = 0;
    *x_start = -1;
    *x_stop  = -1;
  }
  else {
    *x_steps = x2 - x1 + 1;
    *x_start = (float) x1 * grid->dx;
    *x_stop  = (float) x2 * grid->dx;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_de( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_de",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_de ++;

    sprintf( name_de, "%s/spacetime-de-%d-%d", p.path, p.domain_number,output_period_de );
    bob.message( "period =", output_period_de, " time_count =", time_out_count );

    file = fopen( name_de, "wb" );
    if (!file) bob.error( "Cannot open file", name_de );

    fwrite( &output_period_de, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_de, "ab" );
  if (!file) bob.error( "Cannot open file", name_de );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_de, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_de.x_start && cell->number <= stepper_de.x_stop ) {
      output = (float) fabs(cell->dens[0]);
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_di( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_di",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_di ++;

    sprintf( name_di, "%s/spacetime-di-%d-%d", p.path, p.domain_number,output_period_di );
    bob.message( "period =", output_period_di, " time_count =", time_out_count );

    file = fopen( name_di, "wb" );
    if (!file) bob.error( "Cannot open file", name_di );

    fwrite( &output_period_di, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_di, "ab" );
  if (!file) bob.error( "Cannot open file", name_di );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_di, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_di.x_start && cell->number <= stepper_di.x_stop ) {
      output = (float) fabs(cell->dens[1]);
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_jx( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_jx",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_jx ++;

    sprintf( name_jx, "%s/spacetime-jx-%d-%d", p.path, p.domain_number,output_period_jx );
    bob.message( "period =", output_period_jx, " time_count =", time_out_count );

    file = fopen( name_jx, "wb" );
    if (!file) bob.error( "Cannot open file", name_jx );

    fwrite( &output_period_jx, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_jx, "ab" );
  if (!file) bob.error( "Cannot open file", name_jx );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_jx, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_jx.x_start && cell->number <= stepper_jx.x_stop ) {
      output = (float) cell->jx;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_jy( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_jy",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_jy ++;

    sprintf( name_jy, "%s/spacetime-jy-%d-%d", p.path, p.domain_number,output_period_jy );
    bob.message( "period =", output_period_jy, " time_count =", time_out_count );

    file = fopen( name_jy, "wb" );
    if (!file) bob.error( "Cannot open file", name_jy );

    fwrite( &output_period_jy, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_jy, "ab" );
  if (!file) bob.error( "Cannot open file", name_jy );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_jy, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_jy.x_start && cell->number <= stepper_jy.x_stop ) {
      output = (float) cell->jy;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_jz( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_jz",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_jz ++;

    sprintf( name_jz, "%s/spacetime-jz-%d-%d", p.path, p.domain_number,output_period_jz );
    bob.message( "period =", output_period_jz, " time_count =", time_out_count );

    file = fopen( name_jz, "wb" );
    if (!file) bob.error( "Cannot open file", name_jz );

    fwrite( &output_period_jz, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_jz, "ab" );
  if (!file) bob.error( "Cannot open file", name_jz );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_jz, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_jz.x_start && cell->number <= stepper_jz.x_stop ) {
      output = (float) cell->jz;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_ex( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_ex",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_ex ++;

    sprintf( name_ex, "%s/spacetime-ex-%d-%d", p.path, p.domain_number,output_period_ex );
    bob.message( "period =", output_period_ex, " time_count =", time_out_count );

    file = fopen( name_ex, "wb" );
    if (!file) bob.error( "Cannot open file", name_ex );

    fwrite( &output_period_ex, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_ex, "ab" );
  if (!file) bob.error( "Cannot open file", name_ex );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_ex, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_ex.x_start && cell->number <= stepper_ex.x_stop ) {
      output = (float) cell->ex;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_ey( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_ey",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_ey ++;

    sprintf( name_ey, "%s/spacetime-ey-%d-%d", p.path, p.domain_number,output_period_ey );
    bob.message( "period =", output_period_ey, " time_count =", time_out_count );

    file = fopen( name_ey, "wb" );
    if (!file) bob.error( "Cannot open file", name_ey );

    fwrite( &output_period_ey, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_ey, "ab" );
  if (!file) bob.error( "Cannot open file", name_ey );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_ey, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_ey.x_start && cell->number <= stepper_ey.x_stop ) {
      output = (float) cell->ey;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_ez( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_ez",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_ez ++;

    sprintf( name_ez, "%s/spacetime-ez-%d-%d", p.path, p.domain_number,output_period_ez );
    bob.message( "period =", output_period_ez, " time_count =", time_out_count );

    file = fopen( name_ez, "wb" );
    if (!file) bob.error( "Cannot open file", name_ez );

    fwrite( &output_period_ez, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_ez, "ab" );
  if (!file) bob.error( "Cannot open file", name_ez );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_ez, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_ez.x_start && cell->number <= stepper_ez.x_stop ) {
      output = (float) cell->ez;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_bx( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_bx",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_bx ++;

    sprintf( name_bx, "%s/spacetime-bx-%d-%d", p.path, p.domain_number,output_period_bx );
    bob.message( "period =", output_period_bx, " time_count =", time_out_count );

    file = fopen( name_bx, "wb" );
    if (!file) bob.error( "Cannot open file", name_bx );

    fwrite( &output_period_bx, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_bx, "ab" );
  if (!file) bob.error( "Cannot open file", name_bx );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_bx, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_bx.x_start && cell->number <= stepper_bx.x_stop ) {
      output = (float) cell->bx;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_by( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_by",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_by ++;

    sprintf( name_by, "%s/spacetime-by-%d-%d", p.path, p.domain_number,output_period_by );
    bob.message( "period =", output_period_by, " time_count =", time_out_count );

    file = fopen( name_by, "wb" );
    if (!file) bob.error( "Cannot open file", name_by );

    fwrite( &output_period_by, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_by, "ab" );
  if (!file) bob.error( "Cannot open file", name_by );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_by, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_by.x_start && cell->number <= stepper_by.x_stop ) {
      output = (float) cell->by;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_bz( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_bz",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_bz ++;

    sprintf( name_bz, "%s/spacetime-bz-%d-%d", p.path, p.domain_number,output_period_bz );
    bob.message( "period =", output_period_bz, " time_count =", time_out_count );

    file = fopen( name_bz, "wb" );
    if (!file) bob.error( "Cannot open file", name_bz );

    fwrite( &output_period_bz, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_bz, "ab" );
  if (!file) bob.error( "Cannot open file", name_bz );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_bz, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_bz.x_start && cell->number <= stepper_bz.x_stop ) {
      output = (float) cell->bz;
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////


void spacetime::write_edens( domain *grid, int time_out_count, parameter &p )
{
  static error_handler bob("spacetime::write_edens",errname);

  struct cell *cell;
  float       output;
  float       x_start, x_stop;
  int         x_steps;
  FILE        *file;

  if ( time_out_count == 1 ) {

    output_period_edens ++;

    sprintf( name_edens, "%s/spacetime-edens-%d-%d",
                          p.path, p.domain_number, output_period_edens );
    bob.message( "period =", output_period_edens, " time_count =", time_out_count );

    file = fopen( name_edens, "wb" );
    if (!file) bob.error( "Cannot open file", name_edens );

    fwrite( &output_period_edens, sizeof(int), 1, file );
    fwrite( &(p.spp), sizeof(int), 1, file );

    fclose( file );
  }

  file = fopen( name_edens, "ab" );
  if (!file) bob.error( "Cannot open file", name_edens );

  boundaries( &x_start, &x_stop, &x_steps, &stepper_edens, grid );

  fwrite( &x_start, sizeof(float), 1, file );
  fwrite( &x_stop, sizeof(float), 1, file );
  fwrite( &x_steps, sizeof(int), 1, file );

  for( cell=grid->left; cell!=grid->rbuf; cell=cell->next ) {
    if (cell->number >= stepper_edens.x_start && cell->number <= stepper_edens.x_stop ) {
      output = (float) ( pow(cell->ex,2) + pow(cell->ey,2) + pow(cell->ez,2) );
      output += (float) ( pow(cell->bx,2) + pow(cell->by,2) + pow(cell->bz,2) );
      fwrite( &output, sizeof(float), 1, file );
    }
  }

  fclose( file );
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof


