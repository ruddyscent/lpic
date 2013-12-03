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

#ifndef ERROR_H
#define ERROR_H

#include <common.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iomanip>
#include <stdlib.h> // <libc.h>

class error_handler {
    static int error_number;
    static int message_number;
    static int Q_debug;
    static int debug_number;
    static int object_number;
    const char *my_name;
    char       *errname;
    std::ofstream   errfile;
public:
    error_handler(const char *, char *error_file_name);
    void error(char* s1,    char*  s2="",
	       char* s3="", char*  s4="");
    void error(char* s1,    double d2,
	       char* s3="", char*  s4="");

    void message(char* m1,
		 char* m2="", char*  m3="", char* m4="");
    void message(char* m1,    double m2,
		 char* m3="", char*  m4="");
    void message(char* m1   , double m2,    char* m3, double m4);

    void debug(char* m1,
	       char* m2="", char*  m3="", char* m4="");
    void debug(char* m1,    double m2,
	       char* m3="", char*  m4="");
    void debug(char* m1   , double m2,    char* m3, double m4);
};

#endif







