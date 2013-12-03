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

//////////////////////////////////////////////////////////////////////////////////////////

uhr::uhr( parameter &p )
{
  sprintf( errname, "%s/error", p.path );
  error_handler bob("uhr::Constructor",errname);

  start_clock = clock();
  start_time  = time( &start_time );

  h_cpu = m_cpu = s_cpu = 0;
  h_sys = m_sys = s_sys = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::proc( void )
{
  double h, m, s;

  stop_clock = clock();

  if ( stop_clock - start_clock > 0 )
    {
      m = 60.0*modf( (double)( stop_clock - start_clock ) / CLOCKS_PER_SEC/3600, &h );
      s = 60.0*modf( m, &m );

      s_cpu += s;
      if ( s_cpu >= 60 ) { s_cpu -= 60; m_cpu += 1; }

      m_cpu += m;
      if ( m_cpu >= 60 ) { m_cpu -= 60; h_cpu += 1; }

      h_cpu += h;
    }

  start_clock = stop_clock;
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::sys( void )
{
  time(&stop_time);

  m_sys = 60*modf( difftime(stop_time,start_time)/3600, &(h_sys) );
  s_sys = 60*modf( m_sys, &(m_sys) );

  if ( m_sys>= 60 ) { m_sys -= 60; h_sys += 1; }
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::init( void )
{
  start_clock = clock();
  start_time  = time( &start_time );

  h_cpu = m_cpu = s_cpu = 0;
  h_sys = m_sys = s_sys = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////

void uhr::exit( void )
{
  proc();
  sys();

  printf( "\n proc " );
  printf( "%02.0f:", h_cpu );
  printf( "%02.0f:", m_cpu );
  printf( "%05.2f ", s_cpu );
  printf( " sys " );
  printf( "%02.0f:", h_sys );
  printf( "%02.0f:", m_sys );
  printf( "%02.0f ", s_sys );
  printf( "\n\n" );
}

//////////////////////////////////////////////////////////////////////////////////////////
//EOF
