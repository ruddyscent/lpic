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

#include <propagate.h>

void propagate::particles( domain &grid )
// acceleration according Boris (in Birdsall, Langdon)
{
  static error_handler bob("propagate::particles",errname);

  struct cell *cell;
  struct particle *part;

  // assumes fields of the following domain in cell rbuf

  for( cell=grid.left; cell!=grid.rbuf; cell=cell->next )       // for all cells
    {
      if (cell->npart!=0)
	{
	  part=cell->first;
	  do
	    {
#ifdef DEBUG
	      if (!part) bob.error("segmentation");
	      if (part->x < cell->x || part->x >= cell->x + dx) { // particle in cell ? --
		bob.message( "particle at", part->x );
		bob.message( "    in cell", cell->number, "at", cell->x );
		bob.message( "is linked to wrong cell" );
		bob.error("");
	      }
#endif

	      deposit_charge( cell, part );     // not necessary for the local algorithm
	                                        // charge distribution of the
	                                        // preceeding half time step
	      accelerate_1( cell, part );
	    }
	  while( (part=part->next) );

	  part=cell->first;
	  do part->igamma  = 1.0/sqrt(part->igamma);
	  while( (part=part->next) );

	  part=cell->first;
	  do accelerate_2( cell, part );
	  while( (part=part->next) );

	  part=cell->first;
	  do part->igamma  = 1.0/sqrt(part->igamma);
	  while( (part=part->next) );

	  part=cell->first;
	  do
	    {
	      move( part );                       // move all particles

              has_to_change_cell( cell, part );   // put particles on stack

	      deposit_current( cell, part );      // this step is necessary
	    }
	  while( (part=part->next) );
	}
    }

  do_change_cell( grid ); // particles are removed from stack and linked to their
                          // new cells
                          // cells "Lbuf" and "Rbuf" remain empty

  mask_current( grid ); // makes currents invisible near the box boundaries
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::mask_current( domain &grid )
// this routine sets the currents smoothly equal to zero in
// a region of size 2*MASK at the left and right boundary of the simulation box
// this inhibits artificial radiation at the boundaries

{
  static error_handler bob("propagate::mask_current",errname);

  struct cell *cell;
  int i;

  if (domain_number==1) {
    for( i=1, cell=grid.left; i<=2*MASK; i++, cell=cell->next ) { // the first 2MASK cells
      cell->jx *= mask(i);
      cell->jy *= mask(i);
      cell->jz *= mask(i);
    }
  }

  if (domain_number==n_domains) {
    for( i=1, cell=grid.right; i<=2*MASK; i++, cell=cell->prev ) { // the last 2MASK cells
      cell->jx *= mask(i);
      cell->jy *= mask(i);
      cell->jz *= mask(i);
    }
  }
}


inline double propagate::mask( int i )
{
  if (i<MASK)       return 0;
  else {
    if (i<=2*MASK ) return pow( sin(0.5*PI*(i-MASK)/MASK), 2 );
    else            return 1.0;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::accelerate( struct cell *cell, struct particle *part )
  //
  // this function is currently not used
  //
  // full acceleration according Boris (in Birdsall, Langdon)
{
  static error_handler bob("propagate::accelerate",errname);

  register double zmpidt = part->zm * PI * dt;

  register double w    = weighting(cell,part);
  register double notw = 1.0-w;
  register double ex = w * cell->ex + notw * cell->next->ex;       // interpolate fields
  register double ey = w * cell->ey + notw * cell->next->ey;       // to particle position
  register double ez = w * cell->ez + notw * cell->next->ez;
  register double by = w * cell->by + notw * cell->next->by;
  register double bz = w * cell->bz + notw * cell->next->bz;

  register double ux = part->ux + ex * zmpidt;                     // half acceleration
  register double uy = part->uy + ey * zmpidt;
  register double uz = part->uz + ez * zmpidt;

  register double igamma = (1.0 / sqrt(1.0 + (ux*ux + uy*uy + uz*uz)));


  register double ty = by * zmpidt * igamma;                       // rotation
  register double tz = bz * zmpidt * igamma;
  register double t2 = ty*ty + tz*tz;
  register double sy = 2.0 * ty / ( 1.0 + t2 );
  register double sz = 2.0 * tz / ( 1.0 + t2 );

  register double ux2 =   ux * (1.0-tz*sz-ty*sy) + uy * sz           - uz * sy;
  register double uy2 = - ux * sz                + uy * (1.0-tz*sz)  + uz * ty*sz;
  register double uz2 =   ux * sy                + uy * tz*sy        + uz * (1.0-ty*sy);

  part->ux = ux = ux2 + ex * zmpidt;           // second half acceleration
  part->uy = uy = uy2 + ey * zmpidt;
  part->uz = uz = uz2 + ez * zmpidt;

  part->igamma = (1.0 / sqrt(1.0 + (ux*ux + uy*uy + uz*uz)));

  // we take the square roots for all particles per cell in a separate loop !
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::accelerate_1( struct cell *cell, struct particle *part )
// acceleration according Boris (in Birdsall, Langdon)
{
  static error_handler bob("propagate::accelerate_1",errname);

  register double zmpidt = part->zm * PI * dt;

  register double w     = weighting(cell,part);
  register double notw  = 1.0-w;
  //  register double w0     = weighting_0(cell,part);
  //  register double notw0  = 1.0-w0;
  register double ex = w * cell->ex + notw * cell->next->ex;       // interpolate fields
  register double ey = w * cell->ey + notw * cell->next->ey;       // to particle position
  register double ez = w * cell->ez + notw * cell->next->ez;

  register double ux = part->ux + ex * zmpidt;                     // half acceleration
  register double uy = part->uy + ey * zmpidt;
  register double uz = part->uz + ez * zmpidt;

  part->ux = ux;
  part->uy = uy;
  part->uz = uz;

  part->igamma = 1.0 + ux*ux + uy*uy + uz*uz;

  // we take the inverse square roots for all particles per cell in a separate loop !
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::accelerate_2( struct cell *cell, struct particle *part )
// acceleration according Boris (in Birdsall, Langdon)
{
  static error_handler bob("propagate::accelerate_2",errname);

  register double zmpidt = part->zm * PI * dt;
  register double igamma = part->igamma;

  register double ux = part->ux;
  register double uy = part->uy;
  register double uz = part->uz;

  register double w    = weighting(cell,part);
  register double notw = 1.0-w;
  //  register double w0     = weighting_0(cell,part);
  //  register double notw0  = 1.0-w0;
  register double ex = w * cell->ex + notw * cell->next->ex;      // interpolate fields
  register double ey = w * cell->ey + notw * cell->next->ey;      // to particle position
  register double ez = w * cell->ez + notw * cell->next->ez;
  register double by = w * cell->by + notw * cell->next->by;
  register double bz = w * cell->bz + notw * cell->next->bz;

  register double ty = by * zmpidt * igamma;                          // rotation
  register double tz = bz * zmpidt * igamma;
  register double t2 = ty*ty + tz*tz;
  register double sy = 2.0 * ty / ( 1.0 + t2 );
  register double sz = 2.0 * tz / ( 1.0 + t2 );

  register double ux2 =   ux * (1.0-tz*sz-ty*sy) + uy * sz           - uz * sy;
  register double uy2 = - ux * sz                + uy * (1.0-tz*sz)  + uz * ty*sz;
  register double uz2 =   ux * sy                + uy * tz*sy        + uz * (1.0-ty*sy);

  part->ux = ux = ux2 + ex * zmpidt;                  // second half acceleration
  part->uy = uy = uy2 + ey * zmpidt;
  part->uz = uz = uz2 + ez * zmpidt;

  part->igamma = 1.0 + ux*ux + uy*uy + uz*uz;

  // we take the inverse square roots for all particles per cell in a separate loop !
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::move( struct particle *part )
{
  static error_handler bob("propagate::move",errname);

  if ( part->fix==1 ) part->dx = 0;
  else {

    part->dx = dx * part->ux * part->igamma; // dx = Gamma * dt, see constructor!
    part->x += part->dx;

#ifdef DEBUG
    if ( fabs(part->dx) > dx )
      bob.error( "particle displacement larger than grid spacing!" );
#endif

  }
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::has_to_change_cell( struct cell *cell, struct particle *part )
{
  static error_handler bob("propagate::has_to_change_cell",errname);

  if ( part->x < cell->x )            stk.put_on_stack( cell->prev, part );
  else if ( part->x >= cell->x + dx ) stk.put_on_stack( cell->next, part );
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::do_change_cell( domain &grid )
{
  static error_handler bob("propagate::do_change_cell",errname);

  stack_member    *stack, *stack_delete;
  struct cell     *cell;
  struct cell     *new_cell;
  struct particle *part;

  for( cell=grid.Lbuf; cell!=grid.dummy; cell=cell->next )
    {                                  // keep in mind the pointer to the first particle
      cell->insert = cell->first;      // in cell before inserting additional particles
    }                                  // from stack

  stack = stk.zero->next;
  while( stack != stk.hole )           // now insert particles from stack into new cells
    {                                  // and delete them from stack
      part     = stack->part;
      new_cell = stack->new_cell;

      stk.insert_particle( new_cell, part );

      stack_delete = stack;

      stk.remove_from_stack( stack_delete );

      stack = stk.zero->next;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////


void propagate::reflect_particles( domain &grid )
{
  static error_handler bob("propagate::reflect_particles",errname);

  struct particle *part;

  // stack is empty after particles() has been called

  if ( domain_number == 1 ) {

    for( part=grid.lbuf->first; part!=NULL; part=part->next ) {
      part->ux = - part->ux;
      part->x  -= part->dx;
      bob.message( "re l", part->number, "time", time );

      stk.put_on_stack( grid.left, part );
    }

    do_change_cell( grid );
  }

  if ( domain_number == n_domains ) {

    for( part=grid.rbuf->first; part!=NULL; part=part->next ) {
      part->ux = - part->ux;
      part->x  -= part->dx;
      bob.message( "re r", part->number, "time", time );

      stk.put_on_stack( grid.right, part );
    }

    do_change_cell( grid );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::deposit_charge( struct cell *cell, struct particle *part )
{
  double here, next, prev;                    // contributions to charge of this cell, ...
  register double dist = (part->x - cell->x) * idx;

  if ( dist <= 0.5 ) {
    prev = part->zn * ( 0.5 - dist );
    here = part->zn * ( 0.5 + dist );
    cell->prev->charge              += prev;
    cell->prev->dens[part->species] += prev;
    cell->charge                    += here;
    cell->dens[part->species]       += here;
  }
  else {
    here = part->zn * (  1.5 - dist );
    next = part->zn * ( -0.5 + dist );
    cell->charge                    += here;
    cell->dens[part->species]       += here;
    cell->next->charge              += next;
    cell->next->dens[part->species] += next;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::deposit_current( struct cell *cell, struct particle *part )
// We distinguish six cases:
// first, distinguish former position in the first or second half of the cell
// second, distinguish one one-boundary move and two two-boundary moves
// #-boundary move means contributions to # boundary currents Jx
//
// the currents are calculated from the continuity equation,
// assuming rectangular particle shape and area weighting
// J.Villasenor and O.Buneman, Comp. Phys. Comm. 69 (1992) 306-316
{
  static error_handler bob("propagate::deposit_current",errname);

  register double xm = part->x - part->dx;               // before move
  register double xp = part->x;                          // afterwards
  register double x0 = part->cell->x;                    // former left hand cell boundary

  register double x0p05dx = x0 + 0.5*dx;
  register double x0m05dx = x0 - 0.5*dx;
  register double x0p15dx = x0 + 1.5*dx;

  if ( xm < x0p05dx ) {                       // former position in first half of the cell

    if ( xp < x0m05dx )    left_two_left( cell, part );  // two boundary move to the left

    else {
      if ( xp >= x0p05dx ) left_two_right( cell, part ); // two-boundary move to the right
      else                 left_one ( cell, part );      // one boundary move
    }
  }

  else {                                 // former position in the second half of the cell

    if ( xp > x0p15dx )    right_two_right( cell, part );// two boundary move to the right

    else {
      if ( xp <= x0p05dx ) right_two_left( cell, part ); // two boundary move to the left
      else                 right_one( cell, part );      // one boundary move
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::left_one( struct cell *cell, struct particle *part )
{
  static error_handler bob("propagate::left_one",errname);

  register double xm = part->x - part->dx;               // before move
  register double xp = part->x;                          // afterwards
  register double x0 = part->cell->x;                    // former left hand cell boundary

  register double jx0  = part->zn * (xp-xm)*idx;
  /*
  register double jx = cell->jx;
  register double jy = cell->jy;
  register double jz = cell->jz;
  register double pjy = cell->prev->jy;
  register double pjz = cell->prev->jz;
  */
  register double r_1  = 0.5 * part->zn * ( 1.0 - (xp+xm-2.0*x0)*idx ) * part->igamma;
  register double r0   = 0.5 * part->zn * ( 1.0 + (xp+xm-2.0*x0)*idx ) * part->igamma;

  register double jy_1 = r_1 * part->uy;
  register double jy0  = r0  * part->uy;
  register double jz_1 = r_1 * part->uz;
  register double jz0  = r0  * part->uz;
  /*
  cell->jx = jx + jx0;
  cell->jy = jy + jy0;
  cell->jz = jz + jz0;
  cell->prev->jy = pjy + jy_1;
  cell->prev->jz = pjz + jz_1;
  */
  cell->jx       += jx0;
  cell->jy       += jy0;
  cell->jz       += jz0;
  cell->prev->jy += jy_1;
  cell->prev->jz += jz_1;

#ifdef DEBUG
  r_1 /= part->igamma;
  r0  /= part->igamma;

  if ( fabs(r_1+r0 - part->zn) > TINY || fabs(r_1) > fabs(part->zn)
                                       || fabs(r0) > fabs(part->zn) ) {
    bob.message( "r_1      =", r_1 );
    bob.message( "r0       =", r0 );
    bob.message( "part->zn =", part->zn );
    bob.message( "dif      =", fabs(r_1+r0 - part->zn) );
    bob.message( "part->N  =", part->number );
    bob.error( "density splitting wrong!" );
  }
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::left_two_left(  struct cell *cell, struct particle *part )
{
  static error_handler bob("propagate::left_two_left",errname);

  register double xm = part->x - part->dx;              // before move
  register double xp = part->x;                         // afterwards
  register double x0 = part->cell->x;                   // former left hand cell boundary
  /*
  register double jx = cell->jx;
  register double jy = cell->jy;
  register double jz = cell->jz;
  register double pjx = cell->prev->jx;
  register double pjy = cell->prev->jy;
  register double pjz = cell->prev->jz;
  register double ppjy = cell->prev->prev->jy;
  register double ppjz = cell->prev->prev->jz;
  */
  register double jx_1 = - part->zn * ( 0.5 - (xp-x0+dx)*idx );
  register double jx0  = - part->zn * ( 0.5 + (xm-x0)*idx );

  register double eps  = ( xm - (x0-0.5*dx) ) / ( xm - xp );
  register double r_2  = part->zn * part->igamma * 0.5*(1.0-eps) * ( 0.5 - (xp-x0+dx)*idx );
  register double r0   = part->zn * part->igamma * 0.5*eps * ( 0.5 + (xm-x0)*idx );
  register double r_1  = part->zn * part->igamma - r0 - r_2;

  register double jy_2 = r_2 * part->uy;
  register double jy_1 = r_1 * part->uy;
  register double jy0  = r0  * part->uy;
  register double jz_2 = r_2 * part->uz;
  register double jz_1 = r_1 * part->uz;
  register double jz0  = r0  * part->uz;
  /*
  cell->jx = jx + jx0;
  cell->jy = jy + jy0;
  cell->jz = jz + jz0;
  cell->prev->jx = pjx + jx_1;
  cell->prev->jy = pjy + jy_1;
  cell->prev->jz = pjz + jz_1;
  cell->prev->prev->jy = ppjy + jy_2;
  cell->prev->prev->jz = ppjz + jz_2;
  */
  cell->jx             += jx0;
  cell->jy             += jy0;
  cell->jz             += jz0;
  cell->prev->jx       += jx_1;
  cell->prev->jy       += jy_1;
  cell->prev->jz       += jz_1;
  cell->prev->prev->jy += jy_2;
  cell->prev->prev->jz += jz_2;

#ifdef DEBUG
  r_2 /= part->igamma;
  r_1 /= part->igamma;
  r0  /= part->igamma;

  if ( fabs(r_2+r_1+r0 - part->zn) > TINY || fabs(r_2) > fabs(part->zn) ||
                fabs(r_1) > fabs(part->zn) || fabs(r0) > fabs(part->zn) )   {
    bob.message( "r_2      =", r_2 );
    bob.message( "r_1      =", r_1 );
    bob.message( "r0       =", r0 );
    bob.message( "dif      =", fabs(r_2+r_1+r0 - part->zn) );
    bob.message( "part->zn =", part->zn );
    bob.message( "part->N  =", part->number );
    bob.error( "density splitting wrong!" );
  }
  if ( eps < 0 || eps > 1 )
    bob.error( "eps!" );
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::left_two_right( struct cell *cell, struct particle *part )
{
  static error_handler bob("propagate::left_two_right",errname);

  register double xm = part->x - part->dx;              // before move
  register double xp = part->x;                         // afterwards
  register double x0 = part->cell->x;                   // former left hand cell boundary
  /*
  register double jx = cell->jx;
  register double jy = cell->jy;
  register double jz = cell->jz;
  register double pjy = cell->prev->jy;
  register double pjz = cell->prev->jz;
  register double njx = cell->next->jx;
  register double njy = cell->next->jy;
  register double njz = cell->next->jz;
  */
  register double jx0 = part->zn * ( 0.5 - (xm-x0)*idx );
  register double jx1 = part->zn * ( 0.5 + (xp-x0-dx)*idx );

  register double eps = ( x0 + 0.5*dx - xm ) / ( xp - xm );
  register double r_1 = 0.5*eps * part->zn * part->igamma * ( 0.5 - (xm-x0)*idx );
  register double r1  = 0.5*(1.0-eps) * part->zn * part->igamma * ( 0.5 + (xp-x0-dx)*idx );
  register double r0  = part->zn * part->igamma - r_1 - r1;

  register double jy_1 = r_1 * part->uy;
  register double jy0  = r0  * part->uy;
  register double jy1  = r1  * part->uy;

  register double jz_1 = r_1 * part->uz;
  register double jz0  = r0  * part->uz;
  register double jz1  = r1  * part->uz;
  /*
  cell->prev->jy = pjy + jy_1;
  cell->prev->jz = pjz + jz_1;
  cell->jx       = jx + jx0;
  cell->jy       = jy + jy0;
  cell->jz       = jz + jz0;
  cell->next->jx = njx + jx1;
  cell->next->jy = njy + jy1;
  cell->next->jz = njz + jz1;
  */
  cell->prev->jy += jy_1;
  cell->prev->jz += jz_1;
  cell->jx       += jx0;
  cell->jy       += jy0;
  cell->jz       += jz0;
  cell->next->jx += jx1;
  cell->next->jy += jy1;
  cell->next->jz += jz1;

#ifdef DEBUG
  r_1 /= part->igamma;
  r0  /= part->igamma;
  r1  /= part->igamma;

  if ( fabs(r_1+r0+r1-part->zn)>1e-10 || fabs(r_1) > fabs(part->zn) ||
            fabs(r0) > fabs(part->zn) || fabs(r1) > fabs(part->zn) )   {
    bob.message( "r_1      =", r_1 );
    bob.message( "r0       =", r0 );
    bob.message( "r1       =", r1 );
    bob.message( "part->zn =", part->zn );
    bob.message( "dif      =", fabs(r_1+r0+r1-part->zn) );
    bob.message( "part->N  =", part->number );
    bob.error( "density splitting wrong!" );
  }
  if ( eps < 0 || eps > 1 )
    bob.error( "eps!" );
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::right_one( struct cell *cell, struct particle *part )
{
  static error_handler bob("propagate::right_one",errname);

  register double xm = part->x - part->dx;              // before move
  register double xp = part->x;                         // afterwards
  register double x0 = part->cell->x;                   // former left hand cell boundary
  /*
  register double jy = cell->jy;
  register double jz = cell->jz;
  register double njx = cell->next->jx;
  register double njy = cell->next->jy;
  register double njz = cell->next->jz;
  */
  register double jx1 = part->zn * (xp-xm)*idx;

  register double r0  = 0.5 * part->zn * part->igamma * ( 1.0 - (xp+xm-2.0*(x0+dx))*idx );
  register double r1  = 0.5 * part->zn * part->igamma * ( 1.0 + (xp+xm-2.0*(x0+dx))*idx );

  register double jy0 = r0 * part->uy;
  register double jy1 = r1 * part->uy;

  register double jz0 = r0 * part->uz;
  register double jz1 = r1 * part->uz;
  /*
  cell->jy       = jy + jy0;
  cell->jz       = jz + jz0;
  cell->next->jx = njx + jx1;
  cell->next->jy = njy + jy1;
  cell->next->jz = njz + jz1;
  */

  cell->jy       += jy0;
  cell->jz       += jz0;
  cell->next->jx += jx1;
  cell->next->jy += jy1;
  cell->next->jz += jz1;

#ifdef DEBUG
  r0 /= part->igamma;
  r1 /= part->igamma;

  if ( fabs(r0+r1 - part->zn) > TINY || fabs(r0) > fabs(part->zn) ||
                                         fabs(r1) > fabs(part->zn) )  {
    bob.message( "r0       =", r0 );
    bob.message( "r1       =", r1 );
    bob.message( "part->zn =", part->zn );
    bob.message( "dif      =", fabs(r0+r1 - part->zn) );
    bob.message( "part->N  =", part->number );
    bob.error( "density splitting wrong!" );
  }
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::right_two_right( struct cell *cell, struct particle *part )
{
  static error_handler bob("propagate::right_two_right",errname);

  register double xm = part->x - part->dx;              // before move
  register double xp = part->x;                         // afterwards
  register double x0 = part->cell->x;                   // former left hand cell boundary
  /*
  register double jy = cell->jy;
  register double jz = cell->jz;
  register double njx = cell->next->jx;
  register double njy = cell->next->jy;
  register double njz = cell->next->jz;
  register double nnjx = cell->next->next->jx;
  register double nnjy = cell->next->next->jy;
  register double nnjz = cell->next->next->jz;
  */
  register double jx1 = part->zn * ( 0.5 - (xm-x0-dx)*idx );
  register double jx2 = part->zn * ( 0.5 + (xp-x0-2.0*dx)*idx );

  register double eps = (x0+1.5*dx - xm) / (xp - xm);
  register double r0  = 0.5*eps * part->zn * part->igamma * ( 0.5 - (xm-x0-dx)*idx );
  register double r2  = 0.5*(1.0-eps) * part->zn * part->igamma * ( 0.5 + (xp-x0-2.0*dx)*idx );
  register double r1  = part->zn * part->igamma - r0 - r2;

  register double jy0 = r0 * part->uy;
  register double jy1 = r1 * part->uy;
  register double jy2 = r2 * part->uy;

  register double jz0 = r0 * part->uz;
  register double jz1 = r1 * part->uz;
  register double jz2 = r2 * part->uz;
  /*
  cell->jy             = jy + jy0;
  cell->jz             = jz + jz0;
  cell->next->jx       = njx + jx1;
  cell->next->jy       = njy + jy1;
  cell->next->jz       = njz + jz1;
  cell->next->next->jx = nnjx + jx2;
  cell->next->next->jy = nnjy + jy2;
  cell->next->next->jz = nnjz + jz2;
  */

  cell->jy             += jy0;
  cell->jz             += jz0;
  cell->next->jx       += jx1;
  cell->next->jy       += jy1;
  cell->next->jz       += jz1;
  cell->next->next->jx += jx2;
  cell->next->next->jy += jy2;
  cell->next->next->jz += jz2;

#ifdef DEBUG
  r0 /= part->igamma;
  r1 /= part->igamma;
  r2 /= part->igamma;

  if ( fabs(r0+r1+r2 - part->zn) > TINY || fabs(r0) > fabs(part->zn) ||
               fabs(r1) > fabs(part->zn) || fabs(r2) > fabs(part->zn) )  {
    bob.message( "r0       =", r0 );
    bob.message( "r1       =", r1 );
    bob.message( "r2       =", r2 );
    bob.message( "part->zn =", part->zn );
    bob.message( "dif      =", fabs(r0+r1+r2 - part->zn) );
    bob.message( "part->N  =", part->number );
    bob.error( "density splitting wrong!" );
  }
  if ( eps < 0 || eps > 1 )
    bob.error( "eps!" );
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////


inline void propagate::right_two_left( struct cell *cell, struct particle *part )
{
  static error_handler bob("propagate::right_two_left",errname);

  register double xm = part->x - part->dx;               // before move
  register double xp = part->x;                          // afterwards
  register double x0 = part->cell->x;                    // former left hand cell boundary
  /*
  register double pjy = cell->prev->jy;
  register double pjz = cell->prev->jz;
  register double jx  = cell->jx;
  register double jy  = cell->jy;
  register double jz  = cell->jz;
  register double njx = cell->next->jx;
  register double njy = cell->next->jy;
  register double njz = cell->next->jz;
  */
  register double jx0  = - part->zn * ( 0.5 - (xp-x0)*idx );
  register double jx1  = - part->zn * ( 0.5 + (xm-x0-dx)*idx );

  register double eps  = ( xm - (x0+0.5*dx) ) / ( xm - xp );
  register double r_1  = 0.5*(1.0-eps) * part->zn * part->igamma * ( 0.5 - (xp-x0)*idx );
  register double r1   = 0.5*eps * part->zn * part->igamma * ( 0.5 + (xm-x0-dx)*idx );
  register double r0   = part->zn * part->igamma - r_1 - r1;

  register double jy_1 = r_1 * part->uy;
  register double jy0  = r0  * part->uy;
  register double jy1  = r1  * part->uy;

  register double jz_1 = r_1 * part->uz;
  register double jz0  = r0  * part->uz;
  register double jz1  = r1  * part->uz;
  /*
  cell->prev->jy = pjy + jy_1;
  cell->prev->jz = pjz + jz_1;
  cell->jx       = jx + jx0;
  cell->jy       = jy + jy0;
  cell->jz       = jz + jz0;
  cell->next->jx = njx + jx1;
  cell->next->jy = njy + jy1;
  cell->next->jz = njz + jz1;
  */

  cell->prev->jy += jy_1;
  cell->prev->jz += jz_1;
  cell->jx       += jx0;
  cell->jy       += jy0;
  cell->jz       += jz0;
  cell->next->jx += jx1;
  cell->next->jy += jy1;
  cell->next->jz += jz1;

#ifdef DEBUG
  r_1 /= part->igamma;
  r0  /= part->igamma;
  r1  /= part->igamma;

  if ( fabs(r_1+r0+r1 - part->zn) > TINY || fabs(r_1) > fabs(part->zn) ||
               fabs(r0) > fabs(part->zn)  || fabs(r1) > fabs(part->zn) )   {
    bob.message( " r_1      =", r_1 );
    bob.message( " r0       =", r0 );
    bob.message( " r1       =", r1 );
    bob.message( " part->zn =", part->zn );
    bob.message( " dif      =", fabs(r_1+r0+r1 - part->zn) );
    bob.message( "part->N  =", part->number );
    bob.error( "density splitting wrong!" );
  }

  if ( eps < 0 || eps > 1 )
    bob.error( "eps!" );
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////


inline double propagate::weighting( struct cell *cell, struct particle *part )
//
//  returns contribution to the grid point left of the particle's position
//  (1 - return value) is the contribution to the grid point on its right
//
{
//  zero order weighting
//  return ( 1.0 - floor( (part->x - cell->x) * idx + 0.5 ) );

// linear weighting
  return ( 1.0 - (part->x - cell->x) * idx );
}

inline double propagate::weighting_0( struct cell *cell, struct particle *part )
//
//  returns contribution to the grid point left of the particle's position
//  (1 - return value) is the contribution to the grid point on its right
//
{
  // zero order weighting
  return ( 1.0 - floor( (part->x - cell->x) * idx + 0.5 ) );

  // linear weighting
  //  return ( 1.0 - (part->x - cell->x) * idx );
}


//////////////////////////////////////////////////////////////////////////////////////////
//EOF
