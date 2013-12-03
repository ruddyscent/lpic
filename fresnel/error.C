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

#include <error.h>

int error_handler::error_number   = 0;
int error_handler::message_number = 0;
int error_handler::debug_number   = 0;
int error_handler::Q_debug        = 1;
int error_handler::object_number  = 0;

//////////////////////////////////////////////////////////////////////////////////////////

error_handler::error_handler(const char *name, char *error_file_name)
{
  errname = new char [filename_size];
  strcpy(errname,error_file_name);

  errfile.open(errname,ios::app);

  if (!errfile)
    {
      cerr << " Cannot open " << errname << endl;
      exit(1);
    }

  errfile.close();

  my_name = name;
  object_number++;

  debug("");
}

//////////////////////////////////////////////////////////////////////////////////////////

void error_handler::error(char* s1, char* s2, char *s3, char *s4)
{
  error_number++ ;

  errfile.open(errname,ios::app);
  errfile.setf(ios::left);

  errfile << "FAILURE: " << setw(24) << my_name << "       " << s1 << ' ' << s2
    << s3 << s4 << endl;

  errfile.close();

  exit(1);
}

void error_handler::error(char* s1, double d2, char *s3, char *s4)
{
  error_number++ ;

  errfile.open(errname,ios::app);
  errfile.setf(ios::left);

    errfile << "FAILURE: " << setw(24) << my_name << "       " << s1 << ' ' << d2
      << s3 << s4 << endl;

  errfile.close();

  exit(1);
}

//////////////////////////////////////////////////////////////////////////////////////////

void error_handler::message(char *s1, char* s2, char* s3, char* s4)
{
  message_number++ ;

  errfile.open(errname,ios::app);
  errfile.setf(ios::left);

  errfile << setw(33) << my_name << "       "
    << s1   << " " << s2 << " " << s3 << " " << s4 << endl;

  errfile.close();
}

void error_handler::message(char *s1, double d2,
			    char* s3, char* s4)
{
  message_number++ ;

  errfile.open(errname,ios::app);
  errfile.setf(ios::left);

  errfile.precision(12);

  errfile << setw(33) << my_name << "       "
    << s1 << " " << d2 << " " << s3 << " " << s4 <<endl;

  errfile.close();
}

void error_handler::message(char *s1, double d2, char* s3, double d4)
{
  message_number++ ;

  errfile.open(errname,ios::app);
  errfile.setf(ios::left);

  errfile.precision(12);

  errfile << setw(33) << my_name << "       "
    << s1 << " " << d2 << " " << s3 << " " << d4 <<endl;

  errfile.close();
}

//////////////////////////////////////////////////////////////////////////////////////////

void error_handler::debug(char *s1, char* s2, char* s3, char* s4)
{
  if (Q_debug) {
    debug_number++ ;

    errfile.open(errname,ios::app);
    errfile.setf(ios::left);

    errfile << setw(33) << my_name << " DB:" << setw(2) << object_number << " "
      << s1   << " " << s2 << " " << s3 << " " << s4 << endl;

    errfile.close();
  }
}

void error_handler::debug(char *s1, double d2, char* s3, char* s4)
{
  if (Q_debug) {
    debug_number++ ;

    errfile.open(errname,ios::app);
    errfile.setf(ios::left);

    errfile.precision(12);

    errfile << setw(33) << my_name << " DB:" << setw(2) << object_number << " "
      << s1 << " " << d2 << " " << s3 << " " << s4 <<endl;

    errfile.close();
  }
}

void error_handler::debug(char *s1, double d2, char* s3, double d4)
{
  if (Q_debug) {
    debug_number++ ;

    errfile.open(errname,ios::app);
    errfile.setf(ios::left);

    errfile.precision(12);

    errfile << setw(33) << my_name << " DB:" << setw(2) << object_number << " "
      << s1 << " " << d2 << " " << s3 << " " << d4 <<endl;

    errfile.close();
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
//EOF




