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

#include <config.h>

#ifdef LPIC_PARALLEL
#ifdef LPIC_PVM

#include <network.h>

network::network( parameter &p )
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("network::Constructor",errname);

  domain_number = p.domain_number;
  n_domains     = p.n_domains;
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::start_next_task( parameter &p )
{
  static error_handler bob("network::start_next_task",errname);

  if ( p.n_domains > 1 )                // spawn task for following domain
    {                                   // this needs the pvmd running
      tid        = pvm_mytid();
      if (tid<0) bob.error("Maybe you should start the pvm daemon?");

      if (domain_number>1) tid_prev = pvm_parent();
      else          tid_prev = -1;

      bob.message("my tid:   ", tid );
      bob.message("tid_prev: ", tid_prev );

      if (domain_number<p.n_domains)
	{
	  char **arg;
	  arg    = new (char* [3]);
	  arg[0] = new (char [100]);
	  arg[1] = new (char [100]);
	  arg[2] = new (char [100]);
	  sprintf( arg[0], "%d", domain_number+1 );
	  sprintf( arg[1], "%s", p.input_file_name );
	  arg[2] = NULL;

	  bob.message("spawn:   ", p.my_name, arg[0], arg[1] );

	  int numt = pvm_spawn( p.my_name, arg, PvmTaskDefault, "", 1, &tid_next );

	  if (numt!=1) bob.error("cannot start task no.", domain_number+1 );
	  else         bob.message("new task: ", domain_number+1 );
	}
      else tid_next = -1;

      bob.message("tid_next: ", tid_next );
    }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::field( int time_step, domain* grid )
{
  static error_handler bob("network::field",errname);

  if ( domain_number > 1 ) {                        // exchange field copies
    field_send_cpy( grid->left, tid_prev, time_step );
    field_get_cpy( grid->lbuf, tid_prev, time_step );
  }
  if ( domain_number < n_domains ) {
    field_send_cpy( grid->right, tid_next, time_step );
    field_get_cpy( grid->rbuf, tid_next, time_step );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::field_get_cpy( struct cell* cell, int ptid, int time_step )
// recieve from ptid
// store in cell
{
  static error_handler bob("network::field_get",errname);

  int msgtag = time_step;
  double data[9];
                                         // recieve from ptid and store in cell
  pvm_recv( ptid, msgtag );
  pvm_upkdouble( data, 9, 1 );

  cell->fp = data[0];
  cell->gm = data[1];
  cell->fm = data[2];
  cell->gp = data[3];
  cell->ex = data[4];
  cell->ey = data[5];
  cell->ez = data[6];
  cell->by = data[7];
  cell->bz = data[8];
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::field_send_cpy( struct cell* cell, int ptid, int time_step )
// get from cell
// send to ptid
{
  static error_handler bob("network::field_send",errname);

  int msgtag = time_step;
  double data[9];
                                          // send to the next domain
  data[0] = cell->fp;
  data[1] = cell->gm;
  data[2] = cell->fm;
  data[3] = cell->gp;
  data[4] = cell->ex;
  data[5] = cell->ey;
  data[6] = cell->ez;
  data[7] = cell->by;
  data[8] = cell->bz;

  pvm_initsend( PvmDataDefault );
  pvm_pkdouble( data, 9, 1 );
  pvm_send( ptid, msgtag );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::particles( int time_step, domain *grid )
{
  static error_handler bob("network::particles",errname);
  int el_count, ion_count;

  if ( domain_number > 1 ) {                               // exchange particles
    particles_send( grid->lbuf, tid_prev, time_step, &el_count, &ion_count );
    // send particles in lbuf to tid_prev

    grid->n_el   -= el_count;
    grid->n_ion  -= ion_count;
    grid->n_part -= ( el_count + ion_count);

    particles_get( grid->left, tid_prev, time_step, &el_count, &ion_count );
    // get particles from tid_prev into left

    grid->n_el   += el_count;
    grid->n_ion  += ion_count;
    grid->n_part += ( el_count + ion_count);
  }
  if ( domain_number < n_domains ) {
    particles_send( grid->rbuf, tid_next, time_step, &el_count, &ion_count );
    // send particles in rbuf to tid_next

    grid->n_el   -= el_count;
    grid->n_ion  -= ion_count;
    grid->n_part -= ( el_count + ion_count);

    particles_get( grid->right, tid_next, time_step, &el_count, &ion_count );
    // get particles from tid_next into right

    grid->n_el   += el_count;
    grid->n_ion  += ion_count;
    grid->n_part += ( el_count + ion_count);
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::particles_send( struct cell* cell, int ptid, int time_step,
			      int *el_count, int *ion_count )
// cell: take particles from cell
// ptid: send them to tid
{
  static error_handler bob("network::particles_send",errname);

  int msgtag = time_step;
  int npart = cell->npart;
  struct particle *part, *old;

  pvm_initsend( PvmDataDefault );                             // send number of particles
  pvm_pkint( &npart, 1, 1 );
  pvm_send( ptid, msgtag );

  *el_count = *ion_count = 0;

  if ( npart > 0 ) {

    pvm_initsend( PvmDataDefault );   // send particles

    for( part=cell->first; part!=NULL; part=part->next ) {

      pack_particle( part );

      if ( (part->x < cell->x) || (part->x > cell->next->x) )
	bob.error( "particle link to buffer is wrong" );
    }

    part = cell->first;               // delete particles
    do
      { cell->npart --;
	cell->np[part->species] --;
	switch (part->species){       // counters for updating domain's particle
	case 0:                       // numbers grid.n_el, grid.n_ion, grid.n_part
          (*el_count) ++;
          break;
	 case 1:
          (*ion_count) ++;
	  break;
	}
	old  = part;
	part = part->next;
 	delete old;
      }
    while( part!=NULL );
    cell->first=NULL;
    cell->last=NULL;

    pvm_send( ptid, msgtag+1 );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::particles_get( struct cell* cell, int ptid, int time_step,
			     int *el_count, int *ion_count )
// ptid: recieve particles from tid
// cell: put them into this cell
{
  static error_handler bob("network::particles_get",errname);

  int    msgtag = time_step;
  int    i, npart;                         // recieve from next domain
  struct particle *part, *insert_pointer;

  pvm_recv( ptid, msgtag );                // recieve the number of particles to recieve
  pvm_upkint( &npart, 1, 1 );

  *el_count = *ion_count = 0;

  if (npart>0) {

    pvm_recv( ptid, msgtag+1 );            // recieve particles

    insert_pointer = cell->first;

    for( i=0; i<npart; i++ ) {
      part = new( struct particle );       // create a new particle
      if (!part) bob.error("allocation error: part");

      if ( ptid == tid_prev ) {

	unpack_particle( part );

	part->cell    = cell;
	if (insert_pointer!=NULL){         // insert always in front of insert pointer
	    part->next       = insert_pointer;
	    part->prev       = insert_pointer->prev;
	    part->next->prev = part;
	    if (part->prev==NULL) cell->first = part;
	    else part->prev->next = part;
	  }
	else{                              // insert always on bottom
	    part->next    = NULL;
	    part->prev    = cell->last;
	    if (part->prev!=NULL) part->prev->next = part;
	    else cell->first = part;
	    cell->last    = part;
	}
      }
      else if ( ptid == tid_next ) {

	  unpack_particle( part );

	  part->cell    = cell;
	  if (cell->insert!=NULL) part->prev = cell->insert->prev;
	  else part->prev = NULL;
	  part->next    = cell->insert;
	  if (part->prev!=NULL) part->prev->next = part;
	  else cell->first = part;
	  if (part->next!=NULL) part->next->prev = part;
	  else cell->last = part;
      }
      else {
	  bob.error( "ptid neither tid_next nor tid_prev" );
	  exit(-1);
      }

      cell->npart ++;                     // update cell's particle bookkeeping
      cell->np[part->species] ++;
      switch (part->species){             // counters for updating domain's particle
      case 0:                             // numbers grid.n_el, grid.n_ion, grid.n_part
	(*el_count) ++;
	break;
      case 1:
	(*ion_count) ++;
	break;
      }

      if ( (part->x < cell->x) || (part->x > cell->next->x) )
	bob.error( "particle link to new cell in new domain is wrong" );
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current( int time_step, domain* grid )
{
  static error_handler bob("network::current",errname);

  if ( domain_number > 1 ) {
    current_send( grid->Lbuf, tid_prev, time_step );
    // exchange current contributions to cells
    current_get( grid->left, tid_prev, time_step );
    // ""
    current_get_cpy( grid->lbuf, tid_prev, time_step );
    // get a copy of jy and jz into lbuf

    // A copy of currents jy and jz is needed ONLY at the left boundary ( in lbuf )
    // in order to propagate the fields Fplus and Gminus in cell "left".
    // For the propagation of Fminus and Gplus at the right boundary "right",
    // the currents in cell "right" are sufficient!
    // see MPQ-Report 219 p. 22 or propagate::fields in propagate.C
  }
  if ( domain_number < n_domains ) {
    current_send( grid->rbuf, tid_next, time_step );
    // exchange current contributions to cells
    current_get( grid->right->prev, tid_next, time_step );
    // ""
    current_send_cpy( grid->right, tid_next, time_step );
    // send copies of jy and jz to the right __AFTER__ recieving!!
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::density( int time_step, domain* grid )
{
  static error_handler bob("network::density",errname);

  if ( domain_number > 1 ) {
    density_send( grid->lbuf, tid_prev, time_step );
    // exchange density contributions to cells
    density_get( grid->left, tid_prev, time_step );
    // ""
  }
  if ( domain_number < n_domains ) {
    density_send( grid->rbuf, tid_next, time_step );
    // exchange density contributions to cells
    density_get( grid->right, tid_next, time_step );
    // ""
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current_send( struct cell* cell, int ptid, int time_step )
// get currents from cell and cell->next
// send to ptid
{
  static error_handler bob("network::current_send",errname);

  int msgtag = time_step;
  double data[6];
                                          // send to the next domain
  data[0] = cell->jx;
  data[1] = cell->jy;
  data[2] = cell->jz;
  data[3] = cell->next->jx;
  data[4] = cell->next->jy;
  data[5] = cell->next->jz;

  pvm_initsend( PvmDataDefault );
  pvm_pkdouble( data, 6, 1 );
  pvm_send( ptid, msgtag );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current_get( struct cell* cell, int ptid, int time_step )
// recieve from ptid
// add to currents in cell and cell->next
{
  static error_handler bob("network::current_get",errname);

  int msgtag = time_step;
  double data[6];
                                         // recieve from ptid and store in cell
  pvm_recv( ptid, msgtag );
  pvm_upkdouble( data, 6, 1 );
  cell->jx += data[0];
  cell->jy += data[1];
  cell->jz += data[2];
  cell->next->jx += data[3];
  cell->next->jy += data[4];
  cell->next->jz += data[5];
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current_send_cpy( struct cell* cell, int ptid, int time_step )
// send jy, jz from cell to ptid
{
  static error_handler bob("network::current_send",errname);

  int msgtag = time_step;
  double data[2];
                                          // send to the next domain
  data[0] = cell->jy;
  data[1] = cell->jz;

  pvm_initsend( PvmDataDefault );
  pvm_pkdouble( data, 2, 1 );
  pvm_send( ptid, msgtag );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current_get_cpy( struct cell* cell, int ptid, int time_step )
// recieve from ptid copies of jy and jz
// store in cell
{
  static error_handler bob("network::current_get",errname);

  int msgtag = time_step;
  double data[2];
                                         // recieve from ptid and store in cell
  pvm_recv( ptid, msgtag );
  pvm_upkdouble( data, 2, 1 );
  cell->jy = data[0];
  cell->jz = data[1];
}

//////////////////////////////////////////////////////////////////////////////////////////


void network::density_send( struct cell* cell, int ptid, int time_step )
// send densities from cell to ptid
{
  static error_handler bob("network::density_send",errname);

  int msgtag = time_step;
  double data[3];

  data[0] = cell->charge;
  data[1] = cell->dens[0];
  data[2] = cell->dens[1];

  pvm_initsend( PvmDataDefault );
  pvm_pkdouble( data, 3, 1 );
  pvm_send( ptid, msgtag );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::density_get( struct cell* cell, int ptid, int time_step )
// recieve from ptid
// add to density in cell
{
  static error_handler bob("network::density_get",errname);

  int msgtag = time_step;
  double data[3];

  pvm_recv( ptid, msgtag );
  pvm_upkdouble( data, 3, 1 );
  cell->charge += data[0];
  cell->dens[0] += data[1];
  cell->dens[1] += data[2];
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current_1( int time_step, domain* grid )
{
  static error_handler bob("network::current_1",errname);

  if ( domain_number > 1 ) {
    current_get_12( grid->Lbuf, tid_prev, time_step );
    // copy current contributions from previous domain into Lbuf, lbuf, left, left->next
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current_get_12( struct cell* cell, int ptid, int time_step )
// recieve from ptid
// add to currents in cell and cell->next
{
  static error_handler bob("network::current_get_12",errname);

  int msgtag = time_step;
  double data[12];
                                         // recieve from ptid and store in cell
  pvm_recv( ptid, msgtag );
  pvm_upkdouble( data, 12, 1 );
  cell->jx = data[0];
  cell->jy = data[1];
  cell->jz = data[2];
  cell->next->jx = data[3];
  cell->next->jy = data[4];
  cell->next->jz = data[5];
  cell->next->next->jx = data[6];
  cell->next->next->jy = data[7];
  cell->next->next->jz = data[8];
  cell->next->next->next->jx = data[9];
  cell->next->next->next->jy = data[10];
  cell->next->next->next->jz = data[11];
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current_2( int time_step, domain* grid )
{
  static error_handler bob("network::current_2",errname);

  if ( domain_number > 1 ) {
    current_send_12( grid->Lbuf, tid_prev, time_step );
  }
  if ( domain_number < n_domains ) {
    current_send_12( grid->right->prev, tid_next, time_step );
    current_get_12( grid->right->prev, tid_next, time_step );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::current_send_12( struct cell* cell, int ptid, int time_step )
{
  static error_handler bob("network::current_send_12",errname);

  int msgtag = time_step;
  double data[12];
                                          // send to the next domain
  data[0] = cell->jx;
  data[1] = cell->jy;
  data[2] = cell->jz;
  data[3] = cell->next->jx;
  data[4] = cell->next->jy;
  data[5] = cell->next->jz;
  data[6] = cell->next->next->jx;
  data[7] = cell->next->next->jy;
  data[8] = cell->next->next->jz;
  data[9] = cell->next->next->next->jx;
  data[10] = cell->next->next->next->jy;
  data[11] = cell->next->next->next->jz;

  pvm_initsend( PvmDataDefault );
  pvm_pkdouble( data, 12, 1 );
  pvm_send( ptid, msgtag );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::get_part_numbers_from_prev( int *number, int n )
  // in domain #1 : number[i] := 0  (0<=i<n)
  // else         : number[i] = data recieved from previous domain
{
  static error_handler bob("network::get_part_numbers_from_prev",errname);

  int msgtag = domain_number;
  int i;

  if (domain_number > 1) {

      pvm_recv( tid_prev, msgtag );
      pvm_upkint( number, n, 1 );

      bob.message( "recieved  n_el =", number[0], "n_ion =", number[1] );
  }
  else {

      for( i=0; i<n; i++ ) number[i] = 0;

      bob.message( "nothing recieved: domain #", domain_number );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::send_part_numbers_to_next( int *number, int n )
  // in domain #n_doms : -------
  // else              : data sent to next domain = number[i]  (0<=i<n)
{
  static error_handler bob("network::send_part_numbers_to_next",errname);

  int msgtag = domain_number+1;

  if (domain_number < n_domains) {

      pvm_initsend( PvmDataDefault );
      pvm_pkint( number, n, 1 );
      pvm_send( tid_next, msgtag );

      bob.message( "sent  n_el =", number[0], "n_ion =", number[1] );
  }
  else bob.message( "nothing to send: domain #", domain_number );
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::get_total_numbers_from_next( int *number, int n )
  // in domain #n_domains : -------
  // else              : number[i] = data recieved from next domain   (0<=i<n)
{
  static error_handler bob("network::get_total_numbers_from_next",errname);

  int msgtag = domain_number;

  if (domain_number < n_domains) {

      pvm_recv( tid_next, msgtag );
      pvm_upkint( number, n, 1 );

      bob.message( "recieved  n_el =", number[0], "n_ion =", number[1] );
  }
  else bob.message( "nothing to recieve: domain #", domain_number );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::send_total_numbers_to_prev( int *number, int n )
  // in domain #1      : -------
  // else              : data sent to previous domain = number[i]  (0<=i<n)
{
  static error_handler bob("network::send_total_numbers_to_prev",errname);

  int msgtag = domain_number-1;

  if (domain_number > 1) {

      pvm_initsend( PvmDataDefault );
      pvm_pkint( number, n, 1 );
      pvm_send( tid_prev, msgtag );

      bob.message( "sent  n_el =", number[0], "n_ion =", number[1] );
  }
  else bob.message( "nothing to send: domain #", domain_number );
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_get_mesg_from_prev( int* exchange )
  // in domain # 1 : return 0
  // else          : return number of particles to send to ( - )
  //                                         or to recieve from ( + ) previous domain
{
  static error_handler bob("network::reo_get_mesg_from_prev",errname);

  int msgtag = domain_number;

  if (domain_number > 1) {

      pvm_recv( tid_prev, msgtag );
      pvm_upkint( exchange, 1, 1 );

      bob.message( "recieved reo_mesg =", *exchange );
  }
  else {
    bob.message( "no reo_mesg to recieve: domain #", domain_number );
    *exchange = 0;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_send_mesg_to_next( int *exchange )
  // in domain #1      : -------
  // else              : data sent to previous domain = number[i]  (0<=i<n)
{
  static error_handler bob("network::reo_send_mesg_to_next",errname);

  int msgtag = domain_number+1;

  if (domain_number < n_domains) {

      pvm_initsend( PvmDataDefault );
      pvm_pkint( exchange, 1, 1 );
      pvm_send( tid_next, msgtag );

      bob.message( "sent reo_mesg =", *exchange );
  }
  else {
    bob.message( "no reo_mesg to send: domain #", domain_number );
    *exchange = 0;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_from_prev( int *cells_from_prev, int *parts_from_prev )
{
  static error_handler bob("network::reo_from_prev",errname);

  int msgtag = domain_number;

  int data[2];

  if (domain_number > 1) {

      pvm_recv( tid_prev, msgtag );
      pvm_upkint( data, 2, 1 );

      *cells_from_prev = data[0];
      *parts_from_prev = data[1];

      bob.message( "recieve from previous: cells =", data[0], "parts =", data[1] );
  }
  else {
    bob.error( "no previous domain" );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_from_next( int *cells_from_next, int *parts_from_next )
{
  static error_handler bob("network::reo_from_next",errname);

  int msgtag = domain_number;

  int data[2];

  if (domain_number < n_domains) {

      pvm_recv( tid_next, msgtag );
      pvm_upkint( data, 2, 1 );

      *cells_from_next = data[0];
      *parts_from_next = data[1];

      bob.message( "recieve from next: cells =", data[0], "parts =", data[1] );
  }
  else {
    bob.error( "no next domain" );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_recieve_from_prev_and_unpack( int cells_from_prev, int parts_from_prev,
                                  struct cell* firstcell, int *el_count, int *ion_count )
{
  static error_handler bob("network::reo_recieve_from_prev_and_unpack",errname);

  int msgtag = domain_number;
  int i,k;
  int partcount=0;
  struct cell *cell;
  struct particle *part, *part_next;

  *el_count  = 0;
  *ion_count = 0;

  if (cells_from_prev > 0) {

    pvm_recv( tid_prev, msgtag );

    cell = firstcell;
    part = cell->first;

    cell = cell->prev->prev;

    for(i=0;i<cells_from_prev + 2;i++,cell=cell->next) {

      unpack_cell( cell );
      cell->domain  = domain_number;

      for(k=0;k<cell->npart;k++) {

	unpack_particle( part );
	part_next  = part->next;

	if (part->prev!=NULL) part->prev->next  = part->next;
	else                  part->cell->first = part->next;
	if (part->next!=NULL) part->next->prev  = part->prev;
	else                  part->cell->last  = part->prev;
	part->next = NULL;
	part->prev = cell->last;
	if (cell->last!=NULL) cell->last->next = part;
	part->cell       = cell;
	cell->last = part;
	if (part->prev==NULL) cell->first = part;

	switch (part->species){
	case 0:
	  (*el_count) ++;
	  partcount ++;
	  break;
	case 1:
	  (*ion_count) ++;
	  partcount ++;
	  break;
	}
	part = part_next;
      }
    }

    if (partcount!=parts_from_prev) {
      bob.message( "number of particles recieved from prev does" );
      bob.message( "NOT match intended number to receive" );
      bob.error("");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_recieve_from_next_and_unpack( int cells_from_next, int parts_from_next,
                                struct cell* lastcell, int *el_count, int *ion_count )
{
  static error_handler bob("network::reo_recieve_from_next_and_unpack",errname);

  int msgtag = domain_number;
  int i,k;
  int partcount=0;
  struct cell *cell;
  struct particle *part, *part_next;

  *el_count  = 0;
  *ion_count = 0;

  if ( cells_from_next > 0 ) {

    pvm_recv( tid_next, msgtag );

    cell = lastcell;
    part = cell->first;

    for(i=0;i<(cells_from_next - 1);i++, cell=cell->prev);

    for(i=0;i<cells_from_next + 2;i++,cell=cell->next) {

      unpack_cell( cell );
      cell->domain  = domain_number;

      for(k=0;k<cell->npart;k++) {

	unpack_particle( part );
	part_next  = part->next;

	if (part->prev!=NULL) part->prev->next  = part->next;
	else                  part->cell->first = part->next;
	if (part->next!=NULL) part->next->prev  = part->prev;
	else                  part->cell->last  = part->prev;
	part->next = NULL;
	part->prev = cell->last;
	if (cell->last!=NULL) cell->last->next = part;
	part->cell = cell;
	cell->last = part;
	if (part->prev==NULL) cell->first = part;

	switch (part->species){
	case 0:
	  (*el_count) ++;
	  partcount ++;
	  break;
	case 1:
	  (*ion_count) ++;
	  partcount ++;
	  break;
	}

	part = part_next;
      }
    }

    if (partcount!=parts_from_next) {
      bob.message( "number of particles recieved from next does" );
      bob.message( "NOT match intended number to receive" );
      bob.error("");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_to_prev( int cells_to_prev, int parts_to_prev )
  // send  cells_from_prev  and  parts_from_prev  to previous domain
{
  static error_handler bob("network::reo_to_prev",errname);

  int msgtag = domain_number-1;

  int data[2];

  data[0] = cells_to_prev;
  data[1] = parts_to_prev;

  if (domain_number > 1) {

      pvm_initsend( PvmDataDefault );
      pvm_pkint( data, 2, 1 );
      pvm_send( tid_prev, msgtag );

      bob.message( "send to previous: cells =", data[0], "parts =", data[1] );
  }
  else bob.message( "no previous domain" );

}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_to_next( int cells_to_next, int parts_to_next )
  // send  cells_from_next  and  parts_from_next  to next domain
{
  static error_handler bob("network::reo_to_next",errname);

  int msgtag = domain_number+1;

  int data[2];

  data[0] = cells_to_next;
  data[1] = parts_to_next;

  if (domain_number < n_domains) {

      pvm_initsend( PvmDataDefault );
      pvm_pkint( data, 2, 1 );
      pvm_send( tid_next, msgtag );

      bob.message( "send to next: cells =", data[0], "parts =", data[1] );
  }
  else bob.message( "no next domain" );

}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_pack_and_send_to_prev( int cells_to_prev, int parts_to_prev,
                                     struct cell* firstcell )
{
  static error_handler bob("network::reo_pack_and_send_to_prev",errname);

  int msgtag = domain_number-1;
  int i;
  int partcount=0;
  struct cell *cell;
  struct particle *part;

  if ( cells_to_prev > 0 ) {

    pvm_initsend( PvmDataDefault );
    cell = firstcell;

    for(i=0;i<cells_to_prev;i++,cell=cell->next) {
      pack_cell( cell );
      part = cell->first;
      while(part != NULL)
       {
	 pack_particle( part );
	 part = part->next;
	 partcount ++;
       }
      }

    // send two more cells at the right end of the package which
    // will be copied into rbuf and Rbuf in the previous domain
    pack_cell_as_buffer( cell );
    pack_cell_as_buffer( cell->next );

    if (partcount!=parts_to_prev) {
      bob.message( "number of particles sent to prev does" );
      bob.message( "NOT match intended number to send" );
      bob.error("");
    }

    pvm_send( tid_prev, msgtag );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void network::reo_pack_and_send_to_next( int cells_to_next, int parts_to_next,
                                     struct cell* lastcell )
{
  static error_handler bob("network::reo_pack_and_send_to_next",errname);

  int msgtag = domain_number+1;
  int i;
  int partcount=0;
  struct cell *cell;
  struct particle *part;

  if ( cells_to_next > 0 ) {

    pvm_initsend( PvmDataDefault );
    cell = lastcell;

    for(i=0;i<(cells_to_next + 1);i++, cell=cell->prev);

    // also send two more cells at the right end of package which will
    // be copied into Lbuf and lbuf in the next domain
    pack_cell_as_buffer( cell );
    cell = cell->next;
    pack_cell_as_buffer( cell );
    cell = cell->next;

    for(i=0;i<cells_to_next;i++,cell=cell->next) {
      pack_cell( cell );
      part      = cell->first;
      while(part != NULL)
       {
	 pack_particle( part );
	 part      = part->next;
	 partcount ++;
       }
    }

    if (partcount!=parts_to_next) {
      bob.message( "number of particles sent to next does" );
      bob.message( "NOT match intended number to send" );
      bob.error("");
    }

    pvm_send( tid_next, msgtag );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////


void network::pack_particle( struct particle *part )
{
  static error_handler bob("network::pack_particle",errname);

  pvm_pkint( &(part->number), 1, 1 );
  pvm_pkint( &(part->species), 1, 1 );
  pvm_pkint( &(part->fix), 1, 1 );
  pvm_pkdouble( &(part->z), 1, 1 );
  pvm_pkdouble( &(part->m), 1, 1 );
  pvm_pkdouble( &(part->zm), 1, 1 );
  pvm_pkdouble( &(part->x), 1, 1 );
  pvm_pkdouble( &(part->dx), 1, 1 );
  pvm_pkdouble( &(part->igamma), 1, 1 );
  pvm_pkdouble( &(part->ux), 1, 1 );
  pvm_pkdouble( &(part->uy), 1, 1 );
  pvm_pkdouble( &(part->uz), 1, 1 );
  pvm_pkdouble( &(part->n), 1, 1 );
  pvm_pkdouble( &(part->zn), 1, 1 );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::unpack_particle( struct particle *part )
{
  static error_handler bob("network::unpack_particle",errname);

  pvm_upkint( &(part->number), 1, 1 );
  pvm_upkint( &(part->species), 1, 1 );
  pvm_upkint( &(part->fix), 1, 1 );
  pvm_upkdouble( &(part->z), 1, 1 );
  pvm_upkdouble( &(part->m), 1, 1 );
  pvm_upkdouble( &(part->zm), 1, 1 );
  pvm_upkdouble( &(part->x), 1, 1 );
  pvm_upkdouble( &(part->dx), 1, 1 );
  pvm_upkdouble( &(part->igamma), 1, 1 );
  pvm_upkdouble( &(part->ux), 1, 1 );
  pvm_upkdouble( &(part->uy), 1, 1 );
  pvm_upkdouble( &(part->uz), 1, 1 );
  pvm_upkdouble( &(part->n), 1, 1 );
  pvm_upkdouble( &(part->zn), 1, 1 );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::pack_cell( struct cell *cell )
{
  static error_handler bob("network::pack_cell",errname);

  pvm_pkint( &(cell->number), 1, 1 );
  pvm_pkdouble( &(cell->x), 1, 1 );
  pvm_pkdouble( &(cell->charge), 1, 1 );
  pvm_pkdouble( &(cell->jx), 1, 1 );
  pvm_pkdouble( &(cell->jy), 1, 1 );
  pvm_pkdouble( &(cell->jz), 1, 1 );
  pvm_pkdouble( &(cell->ex), 1, 1 );
  pvm_pkdouble( &(cell->ey), 1, 1 );
  pvm_pkdouble( &(cell->ez), 1, 1 );
  pvm_pkdouble( &(cell->bx), 1, 1 );
  pvm_pkdouble( &(cell->by), 1, 1 );
  pvm_pkdouble( &(cell->bz), 1, 1 );
  pvm_pkdouble( &(cell->fp), 1, 1 );
  pvm_pkdouble( &(cell->fm), 1, 1 );
  pvm_pkdouble( &(cell->gp), 1, 1 );
  pvm_pkdouble( &(cell->gm), 1, 1 );
  pvm_pkdouble( (cell->dens), 2, 1 );
  pvm_pkint( (cell->np), 2, 1 );
  pvm_pkint( &(cell->npart), 1, 1 );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::unpack_cell( struct cell *cell )
{
  static error_handler bob("network::unpack_cell",errname);

  pvm_upkint( &(cell->number), 1, 1 );
  pvm_upkdouble( &(cell->x), 1, 1 );
  pvm_upkdouble( &(cell->charge), 1, 1 );
  pvm_upkdouble( &(cell->jx), 1, 1 );
  pvm_upkdouble( &(cell->jy), 1, 1 );
  pvm_upkdouble( &(cell->jz), 1, 1 );
  pvm_upkdouble( &(cell->ex), 1, 1 );
  pvm_upkdouble( &(cell->ey), 1, 1 );
  pvm_upkdouble( &(cell->ez), 1, 1 );
  pvm_upkdouble( &(cell->bx), 1, 1 );
  pvm_upkdouble( &(cell->by), 1, 1 );
  pvm_upkdouble( &(cell->bz), 1, 1 );
  pvm_upkdouble( &(cell->fp), 1, 1 );
  pvm_upkdouble( &(cell->fm), 1, 1 );
  pvm_upkdouble( &(cell->gp), 1, 1 );
  pvm_upkdouble( &(cell->gm), 1, 1 );
  pvm_upkdouble( (cell->dens), 2, 1 );
  pvm_upkint( (cell->np), 2, 1 );
  pvm_upkint( &(cell->npart), 1, 1 );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::pack_cell_as_buffer( struct cell *cell )
  // this cell is packed without particle information since
  // it will be unpacked into a buffer cell
{
  static error_handler bob("network::pack_cell_as_buffer",errname);
  int np[2],npart;
  np[0] = 0;
  np[1] = 0;
  npart = 0;

  pvm_pkint( &(cell->number), 1, 1 );
  pvm_pkdouble( &(cell->x), 1, 1 );
  pvm_pkdouble( &(cell->charge), 1, 1 );
  pvm_pkdouble( &(cell->jx), 1, 1 );
  pvm_pkdouble( &(cell->jy), 1, 1 );
  pvm_pkdouble( &(cell->jz), 1, 1 );
  pvm_pkdouble( &(cell->ex), 1, 1 );
  pvm_pkdouble( &(cell->ey), 1, 1 );
  pvm_pkdouble( &(cell->ez), 1, 1 );
  pvm_pkdouble( &(cell->bx), 1, 1 );
  pvm_pkdouble( &(cell->by), 1, 1 );
  pvm_pkdouble( &(cell->bz), 1, 1 );
  pvm_pkdouble( &(cell->fp), 1, 1 );
  pvm_pkdouble( &(cell->fm), 1, 1 );
  pvm_pkdouble( &(cell->gp), 1, 1 );
  pvm_pkdouble( &(cell->gm), 1, 1 );
  pvm_pkdouble( (cell->dens), 2, 1 );
  pvm_pkint( np, 2, 1 );
  pvm_pkint( &npart, 1, 1 );
}


//////////////////////////////////////////////////////////////////////////////////////////


void network::end_task( void )
{
  printf( "\n end of task in domain #%d\n\n", domain_number );
  if (n_domains>1) pvm_exit();
}


//////////////////////////////////////////////////////////////////////////////////////////
//eof

#endif
#endif
