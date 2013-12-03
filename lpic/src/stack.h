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

#ifndef STACK_H
#define STACK_H

#include <common.h>
#include <cell.h>
#include <particle.h>
#include <error.h>
#include <parameter.h>

typedef struct stack_member_struct {
  struct stack_member_struct *next;
  struct particle            *part;
  struct cell                *new_cell;
} stack_member;

class stack {

 private:
  char errname[filename_size];

 public:
                   stack( parameter &p );
  void      put_on_stack( struct cell *new_cell, struct particle *part );
  void remove_from_stack( stack_member *member );
  void   insert_particle( struct cell *new_cell, struct particle *part );

  stack_member *zero;
  stack_member *hole;
};

#endif


