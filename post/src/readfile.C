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

#include <readfile.h>

//////////////////////////////////////////////////////////////////////////////////////////
//
// openinput(file)  has to be called to open/rewind the input file,
// closeinput(file) has to be called to close input file,
// setinput(k)      can be used to set the file pointer to the key 'k'
// getinput(a)      scans the beginning of every line following the current file pointer
//                  position for the desired variable 'a'
// setget(k,a)      resets the file pointer to the key word 'k',
//                  scans the following lines completely for the desired member variable
//                  'a', allowing for variables seperated by commata  ( NAMELIST )
// read_one_line()  reads single lines of the input file skipping blanks and comments
// write_one_line() writes the recently read line to stdout
//
// read_line(narg,arg)   reads columns of a single line into the argument list '**arg'
//                       and the number of arguments into '*narg', skipping comments
//                       return value 1 for successful reading, 0 for end of file
// read_col(d,col,rows)  reads data in column 'col' into the vector 'd'
//                       maximum number of rows to read is specified by 'rows'
//
// copy_file(a,b)        copies file a to file b
// compare_files(a,b)    compare file a with file b
//
//////////////////////////////////////////////////////////////////////////////////////////

readfile::readfile()
{
  already_open = 0;
  buffer = new char [MAX_LINE_LENGTH];
  if(!buffer){ printf( "error in readfile" ); exit(0);}
  result = new char [MAX_LINE_LENGTH];
  if(!result){ printf( "error in readfile" ); exit(0);}
}

//////////////////////////////////////////////////////////////////////////////////////////
/*
   openinput() has to be called to open the file.
   If the file is already open, the file pointer is reset to the begining
*/

void readfile::openinput(char *file)
{
   if(!already_open){
      fd=fopen(file,"r");
      if(fd==NULL){
         printf("readfile::openinput: can't open file %s\n", file);
         exit(1);
      }
      already_open=1;
   }
   else{
      rewind(fd);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////

void readfile::closeinput( void )
{
   if(already_open){
     int value=fclose( fd );
     already_open=0;
     if(value!=0){printf("file not correctly closed");exit(0);}
   }
}

//////////////////////////////////////////////////////////////////////////////////////////
/*
   This funtion moves the file pointer after the line that matches the string
   and returns 1. Blanks and comments ("#" followed by some text) are ignored.
   If not found, 0 is returned. The file has to be opened previouslly by openinput().
*/

int readfile::setinput(char *a)
{
   int m,n;

   n = strlen(a);

   rewind(fd);

   while(read_one_line()){
      m=strlen(buffer);
      if(m==n){
          if(strncmp(buffer,a,n)==0)return(1);
      }
   }

   return(0);
}

//////////////////////////////////////////////////////////////////////////////////////////

char* readfile::setget(char *key, char *a)
{
   int m,i,n,j=0;

   n = strlen(a);
                                                    // reset file pointer to the key word
   if (!setinput(key)) {
     printf( "\n readfile::setget: key word '%s' missing\n", key );
     exit(-1);
     }

   while(read_one_line()){                        // read lines following the key
     m=strlen(buffer);
     if(m>n+1 && !strchr(buffer,38)) {            // length sufficient, '&' not contained
       for(i=0;i<m-n;i++) {                       // scan the line for variable name
	 if(strncmp(buffer+i,a,n)==0){            // if found, write it to result[]
	   if(buffer[n+i]=='='){
	     i++;
	     while(buffer[n+i+j]!=',' && n+i+j<m ) {
	       result[j]=buffer[n+i+j];
	       j++;
	     }
	     result[j]=0;
	     return(result);                        // and return pointer to result
	   }
	 }
       }
     }
   }
   printf(" readfile::setget: can't find name ");   // otherwise: send error message
   for(i=0;i<n;i++)putchar(a[i]);
   printf(" in input file \n\n");
   exit(1);
   return(result);
}

//////////////////////////////////////////////////////////////////////////////////////////

char* readfile::getinput(char *a)
{
   int m,n,i=0,j=0;

   n = strlen(a);

   rewind(fd);

   while(read_one_line()){                        // read lines
     m=strlen(buffer);
     if(m>n+1) {                                  // length sufficient
       for(i=0;i<m-n;i++) {                       // scan the line for variable name
	 if(strncmp(buffer+i,a,n)==0){            // if found, write it to result[]
	   if(buffer[n+i]=='='){
	     i++;
	     while(buffer[n+i+j]!=',' && n+i+j<m ) {
	       result[j]=buffer[n+i+j];
	       j++;
	     }
	     result[j]=0;
	     return(result);                      // and return pointer to result
	   }
	 }
       }
     }
   }
   printf("readfile::getinput: can't find name ");
   for(i=0;i<n;i++)putchar(a[i]);
   printf(" in input file \n");
   exit(1);
   return(result);
}

//////////////////////////////////////////////////////////////////////////////////////////

int readfile::read_one_line( void )
{
   int i=0,c;
   while(i<MAX_LINE_LENGTH){
      c=getc(fd);
      if(c==EOF)return(0);
      else if(c=='\n'){
         buffer[i++]=0;
         return(1);
      }
      else if(c=='#'){
         buffer[i++]=0;
         while(getc(fd)!='\n');
         return(1);
      }
      else if(c!=' '){
         buffer[i++]=c;
      }
   }
   printf("readfile::read_one_line: line too long\n");
   exit(-1);
   return(-1);
}

//////////////////////////////////////////////////////////////////////////////////////////

void readfile::write_one_line( void )
{
  int i=0;

  printf( "\n" );
  while( buffer[i]!=0 ) putchar(buffer[i++]);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


int readfile::read_col( double *data, int col, int rows_max )
{
  int rows_count=0;
  int narg, nmin=MAX_COL, nmax=0;
  char **arg;
  int check;

  if (col==0) { printf( "\n selected colum number %d invalid", col ); exit(-1); }

  arg = cmatrix(MAX_COL,MAX_LINE_LENGTH);

  do
    {
      check = read_line( &narg, arg );
      if (narg>=col) data[rows_count++] = atof( arg[col-1] );
      else
	{
	  data[rows_count++] = 0;
	  printf( "\n warning: %d column(s) in line %d", narg, rows_count );
	}
      if (narg<nmin) nmin=narg;
      if (narg>nmax) nmax=narg;
    }
  while( check>0 && rows_count<rows_max );

  if (nmin==nmax) printf( "\n %d rows, %d cols", rows_count, nmax );
  else            printf( "\n %d rows, %d...%d cols", rows_count, nmin, nmax );

  free_cmatrix(arg);

  return rows_count;
}


//////////////////////////////////////////////////////////////////////////////////////////


int readfile::read_line( int *narg, char **arg )
/*
  read several arguments from one line in a file specified by FILE pointer f
  arguments have to be se seperated by blanks
  a line is terminated with '\n'
  lines of blanks are skipped
  comments are skipped ( anything following a '#' )

  narg: number of arguments that have been read
  arg[]: pointer to strings containing the arguments

  return value: 0, if end of file was reached
		1, if arguments were read successfully
*/
{
  int c, col=0, pos=0, skip=0;

  do
    {
      c=getc(fd);
      arg[col][pos]=0;
      if ( c=='#' )
	{
	  if (pos>0) col++;
	  pos=0;
	  skip++;
	}
      else if ( c==' ' )
	{
	  if (pos>0) col++;
	  pos=0;
	}
      else if ( c=='\n' )
	{
	  if (pos>0) col++;
	  pos=0;
	  *narg=col;
	  if (col>0) return 1;
	  else       skip=0;
	}
      else if ( c==EOF )
	{
	  if (pos>0) col++;
	  *narg=col;
	  return 0;
	}
      else if ( skip==0 ) arg[col][pos++] = c;
    }
  while( col<MAX_COL && pos<MAX_LINE_LENGTH);

  if (col>=MAX_COL) printf( "\n number of columns too large in 'fread_line'\n" );
  else              printf( "\n length of argument too large in 'fread_line'\n" );
  exit(-1);
  return(-1);
}

//////////////////////////////////////////////////////////////////////////////////////////


void readfile::copy_file( char *input, char *output )
{
  FILE *fin, *fout;
  int c, diff;

  fout = fopen( output, "rb" );
  if (fout) {      // file fout exists
    fclose( fout );
    diff=compare_files(input,output);
    if (diff==0) { // files are identical
      printf( "readfile::copy_file: file %s exist \n", output );
      printf( "                     and is equal to %s -> not copied\n", input );
    }
  }
  else {          // file does not exist so far
    fin = fopen( input, "rb" );
    fout = fopen( output, "wb" );
    while( (c=fgetc(fin)) != EOF ) fputc(c,fout);
    fclose( fin );
    fclose( fout );
  }
}


//////////////////////////////////////////////////////////////////////////////////////////


int readfile::compare_files( char *name1, char *name2 )
  // if one or two of these files do not exist: return -1
  // else:  return number of different pairs of characters
{
  FILE *f1, *f2;
  int c1, c2;
  int differences=0;

  f1 = fopen( name1, "rb" );
  f2 = fopen( name2, "rb" );

  if ( (!f1) || (!f2) ) return -1;
  else {
    do {
      c1=fgetc(f1);
      c2=fgetc(f2);
      if (c1!=c2) differences++;
    } while( (c1 != EOF) && (c2 != EOF) );

    fclose( f1 );
    fclose( f2 );

    return differences;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
//
// tools for allocating strings and vectors of strings
//

#define FREE_ARG char*

char** readfile::cmatrix(long nrh, long nch)
/* allocate a char matrix with subscript range m[0..nrh][0..nch] */
{
  long i, nrow=nrh+1, ncol=nch+1;
  char **m;

  /* allocate pointers to rows */
  m=(char **) malloc((size_t)(nrow*sizeof(char*)));
  if (!m) printf("allocation failure 1 in cmatrix()");

  /* allocate rows and set pointers to them */
  m[0]=(char *) malloc((size_t)((nrow*ncol)*sizeof(char)));
  if (!m[0]) printf("allocation failure 2 in cmatrix()");

  for(i=1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void readfile::free_cmatrix( char **m )
/* free a char matrix allocated by cmatrix() */
{
  free((FREE_ARG) (m[0]));
  free((FREE_ARG) (m));
}

char* readfile::cvector(long nch)
/* allocate a char vector with subscript range m[0..nch] */
{
  long ncol=nch+1;
  char *m;

  /* allocate pointers to columns */
  m=(char*) malloc((size_t)(ncol*sizeof(char)));
  if (!m) printf("allocation failure 1 in cvector()");

  /* return pointer to array of pointers to rows */
  return m;
}

void readfile::free_cvector( char *m )
/* free a char vector allocated by cvector() */
{
  free((FREE_ARG) m);
}

//////////////////////////////////////////////////////////////////////////////////////////
//eof







