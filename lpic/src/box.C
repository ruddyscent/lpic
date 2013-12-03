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

#include <box.h>

box::box( parameter &p )
#ifdef LPIC_PARALLEL
  : input(p),
    rf(),
    talk(p),      // initialize communication
    grid(p)       // initialize grid, cells, particles
#else
  : input(p),
    grid(p)       // initialize grid, cells, particles
#endif
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("box::Constructor",errname);

  n_domains = input.n_domains;

  n_el   = grid.n_el;    // if there is more than one domain
  n_ion  = grid.n_ion;   // these numbers will be updated
  n_part = grid.n_part;  // in the following

  init_restart(p);                                  // init counter for restart save

#ifdef LPIC_PARALLEL
  talk.start_next_task(p);                          // spawn task for the following domain

  init_reorganize(p);                               // init counter for reorganization

  if(input.Q_restart == 0){
    new_global_particle_numbers(grid,talk);         // introduce GLOBAL particle numbers

    com_total_particle_numbers(grid,talk);          // get total particle numbers

    bob.message( "total particle numbers ", n_el, n_ion, n_part );

    reorganize(grid,talk,0);                        // adjust the size of domain
  }
  else{
    com_total_particle_numbers(grid,talk);          // get total particle numbers

    bob.message( "total particle numbers ", n_el, n_ion, n_part );
  }
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////

input_box::input_box( parameter &p )
  : rf()
{
  sprintf( errname, "%s/error-%d", p.path, p.domain_number );
  static error_handler bob("input_box::Constructor",errname);

  rf.openinput(p.input_file_name);

  Q_restart      = atoi( rf.setget( "&restart", "Q" ) );
  strcpy( restart_file, rf.setget( "&restart", "file" ) );
  Q_restart_save = atoi( rf.setget( "&restart", "Q_save" ) );
  strcpy( restart_file_save, rf.setget( "&restart", "file_save" ) );

  n_domains      = atoi( rf.setget( "&parallel", "N_domains" ) );

  Q_reorganize   = atoi( rf.setget( "&parallel", "Q_reo" ) );
  delta_reo      = atoi( rf.setget( "&parallel", "delta_reo" ) );

  rf.closeinput();

  nsp            = p.nsp;

  bob.message("parameter read");

  if (p.domain_number==1) save(p);
}


//////////////////////////////////////////////////////////////////////////////////////////


void input_box::save( parameter &p )
{
  static error_handler bob("input_box::save",errname);
  ofstream outfile;

  outfile.open(p.outname,ios::app);

  outfile << "box" << endl;
  outfile << "------------------------------------------------------------------" << endl;
  outfile << "Q_restart          : " << Q_restart      << endl;
  outfile << "restart_file       : " << restart_file   << endl;
  outfile << "Q_restart_save     : " << Q_restart_save << endl;
  outfile << "restart_file_save  : " << restart_file_save  << endl;
  outfile << "N_domains          : " << n_domains      << endl;
  outfile << "Q_reorganize       : " << Q_reorganize   << endl;
  outfile << "delta_reo          : " << delta_reo      << endl;
  outfile << "nsp                : " << nsp            << endl << endl << endl;

  outfile.close();

  bob.message("parameter written");
}


//////////////////////////////////////////////////////////////////////////////////////////

#ifdef LPIC_PARALLEL
//////////////////////////////////////////////////////////////////////////////////////////

void box::new_global_particle_numbers( domain &grid, network &talk )
{
  static error_handler bob("box::new_global_particle_numbers",errname);

  int             *number;         // count particles for each species sperately
  struct cell     *cell;
  struct particle *part;

  number = new int [input.nsp];
  if (!number) bob.error("allocation error");

  talk.get_part_numbers_from_prev(number,input.nsp);
  // get accumulated numbers for each species from previous domain,if there is one,
  // otherwise set them to zero

  for( cell=grid.left; cell!=grid.rbuf; cell=cell->next ) // for all cells except buffers
    {
      if (cell->npart!=0)                                 // for occupied cells
	{
	  for( part=cell->first; part!=NULL; part=part->next )
	    {
	      number[part->species]++;
	      part->number = number[part->species];
	    }
	}
    }

  talk.send_part_numbers_to_next(number,input.nsp);
  // send accumulated particle numbers to the following task, if there is one

  delete number;
}


//////////////////////////////////////////////////////////////////////////////////////////


void box::com_total_particle_numbers( domain &grid, network &talk )
{
  static error_handler bob("box::com_total_particle_numbers",errname);

  int             *number;         // count particles for each species sperately
  int             n = input.nsp + 1;

  number = new int [ n ];
  if (!number) bob.error("allocation error");

  grid.count_particles();          // current particle numbers in domain

  talk.get_part_numbers_from_prev(number,n);
  // get accumulated numbers for each species from previous domain, if there is one,
  // otherwise set them to zero

  number[0] += grid.n_el;
  number[1] += grid.n_ion;
  number[2] += grid.n_part;

  talk.send_part_numbers_to_next(number,n);
  // send accumulated particle numbers to the following task, if there is one

  talk.get_total_numbers_from_next(number,n);
  // get the total particle numbers from the following task, if there is one

  n_el   = number[0];
  n_ion  = number[1];
  n_part = number[2];

  talk.send_total_numbers_to_prev(number,n);
  // send the total particle numbers to the previous task, if there is one
}

//////////////////////////////////////////////////////////////////////////////////////////

void box::particle_load( domain &grid )
  // writes to error: deviation of particle number from ideal particle number in percent
  // updates domain's particle numbers
{
  static error_handler bob("box::particle_load",errname);
  double deviation, ideal;

  grid.count_particles();

  ideal = 1.0*n_part/n_domains;
  deviation = 100.0*grid.n_part / ideal - 100.0;

  if ( fabs(deviation) < 1e-10 ) bob.message( "particle balance =  0 %" );
  else {
    if (deviation>0) bob.message( "particle balance =  +", deviation, "%" );
    else             bob.message( "particle balance =  -", fabs(deviation), "%" );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


void box::init_reorganize( parameter &p )
{
  static error_handler bob("box::init_reorganize",errname);
  ofstream file;

  if ( n_domains>1 )  reo.Q_reorganize = input.Q_reorganize;
  else                reo.Q_reorganize = 0;
  reo.delta_reo    = input.delta_reo * p.spp;                     // in time steps

  if ( reo.Q_reorganize ) {
    sprintf( reo.file, "%s/reo-%d", p.path, p.domain_number );
    file.open( reo.file, ios::out );
    file << "# time  left boundary -- right boundary -- number of particles" << endl;
    file.close();
  }

  if( input.Q_restart == 0 ) reo.count_reo = reo.delta_reo;
  else{
    char fname[ filename_size ];
    sprintf( fname, "%s/%s-%d-data1", p.path, input.restart_file, p.domain_number );
    rf.openinput(fname);
    reo.count_reo = atoi( rf.getinput( "reo.count_reo" ) );
    rf.closeinput();
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void box::count_reorganize( void )
{
  static error_handler bob("box::count_reorganize",errname);

  reo.count_reo ++;
}

//////////////////////////////////////////////////////////////////////////////////////////

void box::reorganize( domain &grid, network &talk, double time )
{
  static error_handler bob("box::reorganize",errname);
  ofstream file;

  if ( reo.Q_reorganize ) {

    if ( reo.count_reo == reo.delta_reo ) {

      file.open( reo.file, ios::app );
      file.precision( 3 );
      file.setf( ios::showpoint | ios::scientific );
      file << setw(10) << time
	   << setw(10) << grid.left->number
	   << setw(10) << grid.right->number
	   << setw(10) << grid.n_part << endl;

      particle_load(grid);

      reorganize_f(grid,talk);
      reo.count_reo = 0;

      file << setw(10) << time
	   << setw(10) << grid.left->number
	   << setw(10) << grid.right->number
	   << setw(10) << grid.n_part << endl;
      file.close();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void box::reorganize_f( domain &grid, network &talk )
  // reorganize in forward direction, from the left end of the box to the right end
  // attention: do not use reorganize_f while information is stored in buffer cells
  // (e.g. fields, charge, ....), do use it when the buffer cells are empty, e.g.
  // in the main loop after a whole cycle.
{
  static error_handler bob("box::reorganize_f",errname);

  int request_prev, request_next;
  int cells_from_prev, parts_from_prev;
  int cells_to_prev, parts_to_prev;
  int cells_to_next, parts_to_next;
  int cells_from_next, parts_from_next;
  int el_count, ion_count;

  talk.reo_get_mesg_from_prev( &request_prev );

  if (request_prev > 0) { // recieve cells from previous
    talk.reo_from_prev( &cells_from_prev, &parts_from_prev );
         // get the number of cells and particles which will be to recieve from prev
    grid.reo_alloc_from_prev( cells_from_prev, parts_from_prev );
         // allocate memory for cells and particles and
         // link cells to domain, link all the particles to the first cell (left)
         // and update n_left, n_cells, n_part,
    talk.reo_recieve_from_prev_and_unpack( cells_from_prev, parts_from_prev, grid.left,
                                       &el_count, &ion_count );
         // recieve and unpack the cells and particles from prev,
         // and determine number of electrons/ions recieved
    grid.reo_update_n_el_n_ion( el_count, ion_count );
  }
  if (request_prev < 0) { // send cells to previous
    grid.reo_to_prev( request_prev, &cells_to_prev, &parts_to_prev );
         // determine number of cells and particles actually to be sent to previous domain
    talk.reo_to_prev( cells_to_prev, parts_to_prev );
         // inform previous domain of these numbers
    talk.reo_pack_and_send_to_prev( cells_to_prev, parts_to_prev, grid.left );
         // pack the first cells_to_prev cells of the domain and the particles
         // linked to them and send them to the previous domain
    grid.reo_delete_to_prev( cells_to_prev, parts_to_prev );
         // delete memory which is still allocated by already sent cells and particles
         // and update n_left, n_cells, n_el, n_ion, n_part, lbuf and Lbuf,
         // lbuf->next and the pointer cell->prev of the first occupied cell
  }

  request_next = (int) floor( 1.0*grid.n_part - 1.0*n_part/n_domains );
  // > 0 : domain contains too many particles --> send!
  // < 0 : domain does not contain enough particles --> recieve!

  talk.reo_send_mesg_to_next( &request_next );

  if (request_next > 0) {  // send cells to next
    grid.reo_to_next( request_next, &cells_to_next, &parts_to_next );
         // determine number of cells and particles to be sent to next domain
    talk.reo_to_next( cells_to_next, parts_to_next );
         // inform next domain of these numbers
    talk.reo_pack_and_send_to_next( cells_to_next, parts_to_next, grid.right );
         // pack the last cells_to_next cells of the domain and the particles
         // linked to them and send them to the next domain
    grid.reo_delete_to_next( cells_to_next, parts_to_next );
         // delete memory which is still allocated by already sent cells and particles
         // and update n_right, n_cells, n_el, n_ion, n_part, rbuf and Rbuf,
         // rbuf->prev and the pointer cell->next of the last occupied cell
  }
  if (request_next < 0) {  // recieve cells from next
    talk.reo_from_next( &cells_from_next, &parts_from_next );
         // get the number of cells and particles which will be to recieve from next
    grid.reo_alloc_from_next( cells_from_next, parts_from_next );
         // allocate memory for cells and particles and
         // link cells to domain, link all the particles to the last cell (right)
         // and update n_right, n_cells, n_part,
    talk.reo_recieve_from_next_and_unpack( cells_from_next, parts_from_next, grid.right,
                                       &el_count, &ion_count );
         // recieve and unpack the cells and particles from next,
         // and determine number of electrons/ions recieved
    grid.reo_update_n_el_n_ion( el_count, ion_count );
  }

  if ( request_prev != 0 || request_next != 0 ) {

    bob.message( "after reorganization:" );
    bob.message( "number of cells =", grid.n_cells, " ", grid.n_left, "--", grid.n_right );
    bob.message( "particle numbers ", grid.n_el, grid.n_ion, grid.n_part );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
#endif
//////////////////////////////////////////////////////////////////////////////////////////

void box::init_restart( parameter &p )
{
  static error_handler bob("box::init_restart",errname);

  if( input.Q_restart_save == 0 ) rest.Q_restart_save = 0;
  else                            rest.Q_restart_save = 1;

  rest.delta_rest    = p.spp;
  rest.count_rest    = 0 - 1;              // thus data is saved first in propagate::loop
}                                          // and restart save is performed afterwards,
                                           // restart counter is one timestep behind the
                                           // other counters
//////////////////////////////////////////////////////////////////////////////////////////

void box::count_restart( void )
{
  static error_handler bob("box::count_restart",errname);

  rest.count_rest ++;
}

//////////////////////////////////////////////////////////////////////////////////////////

void box::restart_save( diagnostic &diag, double time, parameter &p,
			uhr &zeit, uhr &zeit_particles, uhr &zeit_fields,
			uhr &zeit_diagnostic)
{
  static error_handler bob("box::restart_save",errname);

  if ( rest.Q_restart_save ) {
    if ( rest.count_rest == rest.delta_rest ) {

      ofstream file1;
      FILE *file;
      char fname[ filename_size ];

      sprintf( fname, "%s/%s-%d-data1", p.path,input.restart_file_save, p.domain_number );
      file1.open(fname,ios::out);
      if (!file1) bob.error( "cannot open file", fname );

      file1.precision( 20 );
      file1.setf( ios::showpoint | ios::scientific );

      file1 << "-- restart file ------"<< endl << endl;
      file1 << "time                = " << time << endl << endl;

#ifdef LPIC_PARALLEL
      file1 << "reo.count_reo       = " << reo.count_reo << endl << endl;
#endif

      file1 << "public_time_steps   = " << diag.public_time_steps << endl << endl;
      file1 << "time_out_count      = " << diag.time_out_count << endl << endl;

      file1 << "ene.flux            = " << diag.ene.flux << endl;
      file1 << "ene.field_0         = " << diag.ene.field_0 << endl;
      file1 << "ene.field_t_0       = " << diag.ene.field_t_0 << endl;
      file1 << "ene.field_l_0       = " << diag.ene.field_l_0 << endl;
      file1 << "ene.kinetic_0       = " << diag.ene.kinetic_0 << endl;
      file1 << "ene.total_0         = " << diag.ene.total_0 << endl;
      file1 << "ene.stepper.t_count = " << diag.ene.stepper.t_count << endl << endl;

      file1 << "flu.stepper.t_count = " << diag.flu.stepper.t_count << endl;
      file1 << "pha_el.stepper.t_count = " << diag.pha_el.stepper.t_count << endl << endl;
      file1 << "pha_ion.stepper.t_count = " << diag.pha_ion.stepper.t_count<<endl << endl;

      file1 << "ref.buf[0]          = " << diag.ref.buf[0] << endl;
      file1 << "ref.buf[1]          = " << diag.ref.buf[1] << endl;
      file1 << "ref.buf[2]          = " << diag.ref.buf[2] << endl;
      file1 << "ref.buf[3]          = " << diag.ref.buf[3] << endl;
      file1 << "ref.stepper.t_count = " << diag.ref.stepper.t_count  << endl << endl;

      file1 << "poi.stepper.t_count = " << diag.poi.stepper.t_count << endl;
      file1 << "sna.stepper.t_count = " << diag.sna.stepper.t_count << endl << endl;

      file1 << "spa.stepper_de.t_count    = " << diag.spa.stepper_de.t_count << endl;
      file1 << "spa.stepper_di.t_count    = " << diag.spa.stepper_di.t_count << endl;
      file1 << "spa.stepper_jx.t_count    = " << diag.spa.stepper_jx.t_count << endl;
      file1 << "spa.stepper_jy.t_count    = " << diag.spa.stepper_jy.t_count << endl;
      file1 << "spa.stepper_jz.t_count    = " << diag.spa.stepper_jz.t_count << endl;
      file1 << "spa.stepper_ex.t_count    = " << diag.spa.stepper_ex.t_count << endl;
      file1 << "spa.stepper_ey.t_count    = " << diag.spa.stepper_ey.t_count << endl;
      file1 << "spa.stepper_ez.t_count    = " << diag.spa.stepper_ez.t_count << endl;
      file1 << "spa.stepper_bx.t_count    = " << diag.spa.stepper_bx.t_count << endl;
      file1 << "spa.stepper_by.t_count    = " << diag.spa.stepper_by.t_count << endl;
      file1 << "spa.stepper_bz.t_count    = " << diag.spa.stepper_bz.t_count << endl;
      file1 << "spa.stepper_edens.t_count = " << diag.spa.stepper_edens.t_count << endl;
      file1 << "spa.output_period_de      = " << diag.spa.output_period_de << endl;
      file1 << "spa.output_period_di      = " << diag.spa.output_period_di << endl;
      file1 << "spa.output_period_jx      = " << diag.spa.output_period_jx << endl;
      file1 << "spa.output_period_jy      = " << diag.spa.output_period_jy << endl;
      file1 << "spa.output_period_jz      = " << diag.spa.output_period_jz << endl;
      file1 << "spa.output_period_ex      = " << diag.spa.output_period_ex << endl;
      file1 << "spa.output_period_ey      = " << diag.spa.output_period_ey << endl;
      file1 << "spa.output_period_ez      = " << diag.spa.output_period_ez << endl;
      file1 << "spa.output_period_bx      = " << diag.spa.output_period_bx << endl;
      file1 << "spa.output_period_by      = " << diag.spa.output_period_by << endl;
      file1 << "spa.output_period_bz      = " << diag.spa.output_period_bz << endl;
      file1 << "spa.output_period_edens   = " << diag.spa.output_period_edens << endl;
      file1 << endl;

      file1 << "tra.stepper.t_count     = " << diag.tra.stepper.t_count << endl;
      file1 << "vel_el.stepper.t_count  = " << diag.vel_el.stepper.t_count << endl<< endl;
      file1 << "vel_ion.stepper.t_count = " << diag.vel_ion.stepper.t_count <<endl<< endl;

      file1.close();

      zeit.restart_save();
      zeit_particles.restart_save();
      zeit_fields.restart_save();
      zeit_diagnostic.restart_save();

      diag.tra.restart_save();

      sprintf( fname, "%s/%s-%d-data2", p.path,input.restart_file_save, p.domain_number );
      file = fopen( fname, "wb" );
      if (!file) bob.error( "cannot open file", fname );

      struct cell *cell;
      struct particle *part;
      int n_cells_check,n_el_check, n_ion_check, n_part_check;

      fwrite( &grid.n_cells, sizeof(int), 1, file );

      n_cells_check = 0;
      n_el_check    = 0;
      n_ion_check   = 0;
      n_part_check  = 0;

      for( cell=grid.Lbuf; cell!=grid.dummy; cell=cell->next )
	{
	  fwrite( &cell->number , sizeof(int), 1, file );
	  fwrite( &cell->x      , sizeof(double), 1, file );
	  fwrite( &cell->charge , sizeof(double), 1, file );
	  fwrite( &cell->jx     , sizeof(double), 1, file );
	  fwrite( &cell->jy     , sizeof(double), 1, file );
	  fwrite( &cell->jz     , sizeof(double), 1, file );
	  fwrite( &cell->ex     , sizeof(double), 1, file );
	  fwrite( &cell->ey     , sizeof(double), 1, file );
	  fwrite( &cell->ez     , sizeof(double), 1, file );
	  fwrite( &cell->bx     , sizeof(double), 1, file );
	  fwrite( &cell->by     , sizeof(double), 1, file );
	  fwrite( &cell->bz     , sizeof(double), 1, file );
	  fwrite( &cell->fp     , sizeof(double), 1, file );
	  fwrite( &cell->fm     , sizeof(double), 1, file );
	  fwrite( &cell->gp     , sizeof(double), 1, file );
	  fwrite( &cell->gm     , sizeof(double), 1, file );
	  fwrite( &(cell->dens[0]), sizeof(double), 1, file );
	  fwrite( &(cell->dens[1]), sizeof(double), 1, file );
	  fwrite( &(cell->np[0])  , sizeof(int), 1, file );
	  fwrite( &(cell->np[1])  , sizeof(int), 1, file );
	  fwrite( &cell->npart    , sizeof(int), 1, file );
	  n_cells_check ++;

	  if (cell->npart!=0){
	    part=cell->first;
	    do{
	      fwrite( &part->number , sizeof(int), 1, file );
	      fwrite( &part->species, sizeof(int), 1, file );
	      fwrite( &part->fix    , sizeof(int), 1, file );
	      fwrite( &part->z      , sizeof(double), 1, file );
	      fwrite( &part->m      , sizeof(double), 1, file );
	      fwrite( &part->zm     , sizeof(double), 1, file );
	      fwrite( &part->x      , sizeof(double), 1, file );
	      fwrite( &part->dx     , sizeof(double), 1, file );
	      fwrite( &part->igamma , sizeof(double), 1, file );
	      fwrite( &part->ux     , sizeof(double), 1, file );
	      fwrite( &part->uy     , sizeof(double), 1, file );
	      fwrite( &part->uz     , sizeof(double), 1, file );
	      fwrite( &part->zn     , sizeof(double), 1, file );

	      switch (part->species){
	      case 0:
		n_el_check   ++;
		n_part_check ++;
		break;
	      case 1:
		n_ion_check  ++;
		n_part_check ++;
		break;
	      }
	    }while( (part=part->next) );
	  }
	}

      fwrite( &grid.n_el  , sizeof(int), 1, file );
      fwrite( &grid.n_ion , sizeof(int), 1, file );
      fwrite( &grid.n_part, sizeof(int), 1, file );
      fclose( file );

      n_cells_check -= 4;
      if( n_cells_check!=grid.n_cells){
	bob.error("n_cells incorrect:");
      }
      if( n_el_check!=grid.n_el){
	bob.error("n_el incorrect");
      }
      if( n_ion_check!=grid.n_ion){
	bob.error("n_ion incorrect");
      }
      if( n_part_check!=grid.n_part){
	bob.error("n_part incorrect");
      }

      rest.count_rest = 0;
      bob.message("restart files stored at time=",time);
    }
  }
}
/////////////////////////////////////////////////////////////////////////////////////////
//eof


