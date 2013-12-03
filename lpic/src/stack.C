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

#include <stack.h>

//////////////////////////////////////////////////////////////////////////////////////////

stack::stack( parameter &p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("stack::Constructor", errname );

  zero = new stack_member;
  hole = new stack_member;

  zero->next = hole;
}

//////////////////////////////////////////////////////////////////////////////////////////

void stack::put_on_stack( struct cell *new_cell, struct particle *part )
{
  static error_handler bob("stack::put_on_stack", errname);

  stack_member *new_member;

  new_member = new stack_member;
  if (!new_member) bob.error( "allocation error" );

  new_member->part     = part;
  new_member->new_cell = new_cell;

  new_member->next = zero->next;
  zero->next       = new_member;

#ifdef DEBUG
  if ( (part->x < new_cell->x)  ||  (part->x > new_cell->next->x) )
    bob.error( "particle on stack for wrong cell" );
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////

void stack::remove_from_stack( stack_member *member )
{
  zero->next = member->next;
  delete member;
}

//////////////////////////////////////////////////////////////////////////////////////////

void stack::insert_particle( struct cell *new_cell, struct particle *part )
{
  if (part->prev!=NULL) part->prev->next = part->next; // extract particle from old chain
  else part->cell->first = part->next;
  if (part->next!=NULL) part->next->prev = part->prev;
  else part->cell->last = part->prev;

  if (part->cell->insert == part) part->cell->insert = part->next; // new insert mark!!

  part->prev           = NULL;                    // insert particle into the new chain
  part->next           = new_cell->first;
  new_cell->first      = part;
  if (part->next==NULL) new_cell->last = part;
  else part->next->prev = part;

  part->cell->np[ part->species ] --;
  part->cell->npart               --;
    new_cell->np[ part->species ] ++;
    new_cell->npart               ++;

  part->cell = new_cell;
}




//////////////////////////////////////////////////////////////////////////////////////////
//eof
