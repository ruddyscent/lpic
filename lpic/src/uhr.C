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

#include <uhr.h>

using namespace std;
//////////////////////////////////////////////////////////////////////////////////////////

uhr::uhr( parameter &p, char *name )
  : rf(),
    input(p)
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("uhr::Constructor",errname);

  domain_number = p.domain_number;
  strcpy( uhrname, name );
  strcpy( path, p.path );

  if( input.Q_restart == 0 ) reset();
  else                       restart();
}

//////////////////////////////////////////////////////////////////////////////////////////

input_uhr::input_uhr( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_uhr::Constructor",errname);

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


void input_uhr::save( parameter &p )
{
  static error_handler bob("input_uhr::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "uhr" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Q_restart        : " << Q_restart       << endl;
  outfile << "restart_file     : " << restart_file    << endl;
  outfile << "Q_restart_save   : " << Q_restart_save  << endl;
  outfile << "restart_file_save: " << restart_file_save << endl << endl << endl;;

  outfile.close();

  bob.message("parameter written");
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::start( void )
{
  static error_handler bob("uhr::start",errname);
  start_tics = clock();
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::add( void )
{
  static error_handler bob("uhr::add",errname);
  stop_and_add();
  start_tics = clock();
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::stop_and_add( void )
{
  static error_handler bob("uhr::stop_and_add",errname);
  double h, m, s;

  stop_tics = clock();

  if ( stop_tics - start_tics > 0 )
    {
      tics += ( stop_tics - start_tics );

      sec_cpu += (double)( stop_tics - start_tics ) / CLOCKS_PER_SEC;

      m = 60.0*modf( (double)( stop_tics - start_tics ) / CLOCKS_PER_SEC/3600, &h );
      s = 60.0*modf( m, &m );

      s_cpu += s;
      if ( s_cpu >= 60 ) { s_cpu -= 60; m_cpu += 1; }

      m_cpu += m;
      if ( m_cpu >= 60 ) { m_cpu -= 60; h_cpu += 1; }

      h_cpu += h;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::sys( void )
{
  static error_handler bob("uhr::sys",errname);
  time(&stop_time);

  sec_sys = difftime(stop_time,start_time);

  m_sys = 60*modf( difftime(stop_time,start_time)/3600, &(h_sys) );
  s_sys = 60*modf( m_sys, &(m_sys) );
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::reset( void )
{
  static error_handler bob("uhr::reset",errname);
  start_tics = clock();

  start_time = time( &start_time );

  sec_cpu = sec_sys = 0;

  h_cpu = m_cpu = s_cpu = 0;
  h_sys = m_sys = s_sys = 0;

  tics = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////


void uhr::seconds_cpu( void )
{
  static error_handler bob("uhr::seconds_cpu",errname);
  ofstream f;
  char filename[filename_size];

  double seconds = (double) tics / CLOCKS_PER_SEC;

  cout << " cpu " << setw(7) << seconds << " sec : " << uhrname << " ";
  if ( fabs(seconds-sec_cpu)/sec_cpu > 1e-6 ) cout << sec_cpu;
  cout << endl;

  sprintf( filename, "%s/times-%d", path, domain_number );
  f.open(filename,ios::app);

  f << " cpu " << setw(7) << seconds << " sec : " << uhrname << " ";
  if ( fabs(seconds-sec_cpu)/sec_cpu > 1e-6 ) f << sec_cpu;
  f << endl;

  f.close();
}


//////////////////////////////////////////////////////////////////////////////////////////


void uhr::seconds_sys( void )
{
  static error_handler bob("uhr::seconds_sys",errname);
  ofstream f;
  char filename[filename_size];

  sys();

  cout << " sys " << setw(7) << sec_sys << " sec : " << uhrname << endl;

  sprintf( filename, "%s/times-%d", path, domain_number );
  f.open(filename,ios::app);

  f << " sys " << setw(7) << sec_sys << " sec : " << uhrname << endl;

  f.close();
}


//////////////////////////////////////////////////////////////////////////////////////////

void uhr::restart_save( void )
{
  static error_handler bob("uhr::restart_save",errname);
  ofstream file1;
  char fname[ filename_size ];

  sprintf( fname, "%s/%s-%d-data1", path, input.restart_file_save, domain_number );
  file1.open(fname,ios::app);
  if (!file1) bob.error( "cannot open file", fname );

  file1.precision( 20 );
  file1.setf( ios::showpoint | ios::scientific );

  file1 << uhrname << ".start_tics     = " << start_tics << endl;
  file1 << uhrname << ".stop_tics      = " << stop_tics << endl;
  file1 << uhrname << ".tics           = " << tics << endl;
  file1 << uhrname << ".h_cpu          = " << h_cpu << endl;
  file1 << uhrname << ".m_cpu          = " << m_cpu << endl;
  file1 << uhrname << ".s_cpu          = " << s_cpu << endl;
  file1 << uhrname << ".sec_cpu        = " << sec_cpu << endl;
  file1 << uhrname << ".h_sys          = " << h_sys << endl;
  file1 << uhrname << ".m_sys          = " << m_sys << endl;
  file1 << uhrname << ".s_sys          = " << s_sys << endl;
  file1 << uhrname << ".sec_sys        = " << sec_sys << endl;
  file1 << uhrname << ".start_time     = " << start_time << endl;
  file1 << uhrname << ".stop_time      = " << stop_time << endl << endl;

  file1.close();
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::restart( void )
{
  static error_handler bob("uhr::restart",errname);

  char fname[ filename_size ];
  sprintf( fname, "%s/%s-%d-data1", path, input.restart_file, domain_number );
  rf.openinput(fname);

  char dataname[filename_size];

  sprintf( dataname, "%s.start_tics", uhrname );
  start_tics = atoi( rf.getinput( dataname ) );

  sprintf( dataname, "%s.stop_tics", uhrname );
  stop_tics  = atoi( rf.getinput( dataname ) );

  sprintf( dataname, "%s.tics", uhrname );
  tics       = atoi( rf.getinput( dataname ) );

  sprintf( dataname, "%s.h_cpu", uhrname );
  h_cpu   = atof( rf.getinput( dataname ) );

  sprintf( dataname, "%s.m_cpu", uhrname );
  m_cpu   = atof( rf.getinput( dataname ) );

  sprintf( dataname, "%s.s_cpu", uhrname );
  s_cpu   = atof( rf.getinput( dataname ) );

  sprintf( dataname, "%s.sec_cpu", uhrname );
  sec_cpu = atof( rf.getinput( dataname ) );

  sprintf( dataname, "%s.h_sys", uhrname );
  h_sys   = atof( rf.getinput( dataname ) );

  sprintf( dataname, "%s.m_sys", uhrname );
  m_sys   = atof( rf.getinput( dataname ) );

  sprintf( dataname, "%s.s_sys", uhrname );
  s_sys   = atof( rf.getinput( dataname ) );

  sprintf( dataname, "%s.sec_sys", uhrname );
  sec_sys = atof( rf.getinput( dataname ) );

  sprintf( dataname, "%s.start_time", uhrname );
  start_time = atoi( rf.getinput( dataname ) );

  sprintf( dataname, "%s.stop_time", uhrname );
  stop_time  = atoi( rf.getinput( dataname ) );

  rf.closeinput();
}

//////////////////////////////////////////////////////////////////////////////////////////
//EOF












