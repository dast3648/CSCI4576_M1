// Conway's Game of Life
// Main Executable Program
//
// CSCI 4576/5576 High Performance Scientific Computing
// Michael Oberg, modified from code supplied by Dr. Matthew Woitaszek

// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

// Include global variables. Only this file needs the #define
#define __MAIN 
#include "globals.h"
#undef __MAIN

#define verbose 1

// User includes
#include "pprintf.h"
#include "pgm.h"

int main(int argc, char* argv[]) 
{
	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Get the communicator and process information
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	// Print rank and hostname
	MPI_Get_processor_name(my_name, &my_name_len);
	printf("Rank %i is running on %s\n", rank, my_name );

	// Initialize the pretty printer
	init_pprintf( rank );
	pp_set_banner( "main" );
	if( rank==0 )
	pprintf( "Welcome to Conway's Game of Life!\n" );

	// --------------------------------
	// This is where the code goes
	// -------------------------------
	//
	// Determine the partitioning
	//
	// (we fake being a 9x1 block row partitioning)
	//nrows = 9;
	//ncols = 1;
	//if( np != nrows * ncols )
	//{
	  //if( rank==0 )
		//pprintf("Error: %ix%i partitioning requires %i np (%i provided)\n", 
		  //nrows, ncols, nrows * ncols, np );
	  //MPI_Finalize();
	  //return 1;
	//}
	//my_col = 0;
	//my_row = rank;

	// Serial Code
	if ( np == 1 )
	{
		nrows  = 1;
		ncols  = 1;
		my_col = 0;
		my_row = 0;
	}
	// Now, calculate neighbors (N, S, E, W, NW, NE, SW, SE)  
	// ... which means you ...

	/* Read the PGM file. The readpgm() routine reads the PGM file and, based
	 * on the previously set nrows, ncols, my_row, and my_col variables, loads
	 * just the local part of the field onto the current processor. The
	 * variables local_width, local_height, field_width, field_height, as well
	 * as the fields (field_a, field_b) are allocated and filled.
	 */
	if(!readpgm("life.pgm"))
	{
	  if( rank==0 )
		pprintf( "An error occured while reading the pgm file\n" );
	  MPI_Finalize();
	  return 1;
	}

	// primitive bool to discern which field is read from
	// and which field gets written to per iteration.
	// 0: Read from a and write to b.
	// 1: Read from b and write to a.
	int curIter;
	int tolIter = 10;

for (curIter = 0; curIter < tolIter; curIter++)
{
	int i;
	int j;
	int cellNeighbors;
	int	liveCells = 0;

	#if (verbose == 1)
	int deadCells = 0;
	int bornCells = 0;
	int ripCells = 0;
	int survCells = 0;
	#endif
	
	for (i = 1; i < local_height+1; i++)
		for (j = 1; j < local_width+1; j++)
		{
			cellNeighbors = 0;
			
			if ( field_a[ (i-1) * field_width + (j-1) ] ) cellNeighbors++;
			if ( field_a[ (i-1) * field_width + (j+0) ] ) cellNeighbors++;
			if ( field_a[ (i-1) * field_width + (j+1) ] ) cellNeighbors++;
			if ( field_a[ (i+0) * field_width + (j-1) ] ) cellNeighbors++;
		//	if ( field_a[ (i+0) * field_width + (j+0) ] ) cellNeighbors++;
			if ( field_a[ (i+0) * field_width + (j+1) ] ) cellNeighbors++;
			if ( field_a[ (i+1) * field_width + (j-1) ] ) cellNeighbors++;
			if ( field_a[ (i+1) * field_width + (j+0) ] ) cellNeighbors++;
			if ( field_a[ (i+1) * field_width + (j+1) ] ) cellNeighbors++;

			// Cell is dead.
			if (cellNeighbors <= 1 || cellNeighbors >= 4)
			{
				field_b[ i * field_width + j ] = 0;

				// Census bookeeping
				#if (verbose == 1)
				if ( field_a[ i * field_width + j ] ) 
					ripCells++;
				else 	
					deadCells++;
				#endif
			}
			else 
			{
				if (cellNeighbors == 3)
				{
					field_b[ i * field_width + j ] = 1;
				}

				// Census bookeeping
				#if (verbose == 1)
				if ( field_a[ i * field_width + j ] )
					survCells++;
				else 
				{
					if (cellNeighbors == 3)
						bornCells++;
				}
				#endif
			}
		}
	// Count the life forms. Note that we count from [1,1] - [height+1,width+1];
	// we need to ignore the ghost row!
		for( int y=1; y<local_height+1; y++ )
		{
			#if (verbose == 2)
			printf("%d",curIter);
			#endif
	  		for( int x=1; x<local_width+1; x++ )
	  		{
	  			if (field_a[y*field_width+x]) liveCells++;
					field_a[y*field_width+x] = field_b[y*field_width+x];
				#if (verbose == 2)
				printf(",%d", field_a[y*field_width+x]);
				#endif
	  		}
			#if (verbose == 2)
			printf("\n");
			#endif
		}
		#if (verbose == 1)
		if( rank==0 )
		{
			pprintf( "~~ ~  ~   ~ { Iteration: %4d } ~   ~  ~ ~~\n", curIter+1);
		}
		#endif
		//pprintf( "%15i local buggies\n", liveCells);

		int total;
		MPI_Allreduce( &liveCells, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
		if( rank==0 )
		{
			#if (verbose == 1)
	  		pprintf( "%15i total buggies\n", total );
			pprintf( "%15d initial dead cells\n", deadCells);
			pprintf( "%15d surviving cells\n", survCells);
	  		pprintf( "%15d cells born\n", bornCells);
			pprintf( "%15d cells died\n", ripCells);
			pprintf( "~~ ~  ~   ~ { Out of: %7d } ~   ~  ~ ~~\n\n\n", tolIter);
			#endif
		}
	//swap(&field_a, &field_b);
	}

	// Free the fields
	if( field_a != NULL ) free( field_a );
	if( field_b != NULL ) free( field_b );
	// Finalize MPI and terminate
	if( rank==0 )
	  pprintf( "Terminating normally\n" );
	MPI_Finalize();
	return 0;
} 
