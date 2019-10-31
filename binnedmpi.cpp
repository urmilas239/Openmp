#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>
#include "common.h"




//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0,numthreads;
    double davg,dmin, absmin=1.0, absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 

    std::vector<std::vector<int> > bin_map;
    std::vector<std::vector<int> > neighbor_bins;
    std::vector<int> process_bins; 
    std::vector<int> bin_process_map;
    std::vector<int> border_neighbors;

    MPI_Status status; 



    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 5000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );


    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Printout the hosts
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Hello from %s\n", processor_name);
    //

    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;



     int num_of_particles_in_proc;
     int num_of_bins_in_proc;
     int num_of_neighbors_in_proc;


    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    int *partition_sizes= (int*) calloc( n_proc , sizeof(int));
    int *partition_offsets = (int*) calloc( (n_proc) , sizeof(int) );


    MPI_Datatype PARTICLE;
    //Need 6 double and one int. Since we cannot make MPI_Datatype of multiple types using double
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    

    
    

    set_size( n );

    //call the function to set bin variables
    set_bin_count(n);

    // All initialization calls happens in rank 0
     
    bin_map = initialize_bin_vector();
    
    neighbor_bins = initialize_neighbor_bins();


    //init_particles1(n, particles,bin_map);
    init_particles( n, particles);
    bin_particles( n, particles , bin_map);
    process_bins=assign_bins_to_current_process_mpi(n_proc, rank, bin_map, bin_process_map);
    border_neighbors = get_boundary_bins_for_curr_process(process_bins, neighbor_bins);


    std::vector<int> neighbor_list_collect; 
    std::vector<int> particle_list_collect; 
    int number_of_interacting_particles = 0;

    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);
    
     int partition_size_per_rank = 0;
     int npm_size =0;
     int pbm_size = 0;
    //std::cout<<"numthreads:::"<<numthreads<<std::endl;
    for( int step = 0; step < NSTEPS; step++ )
    {
        MPI_Barrier(MPI_COMM_WORLD);
         //printf( ":::::::::::::IN TIME STEP::::::::::::::::::::::::::::::::::::: %d\n" , step);
           
            
        for(i=0;i<process_bins.size();i++)
        {
            //collect neighbor particle for the current bin.
            neighbor_list_collect = neighbor_bins.at(process_bins.at(i));
            neighbor_list_collect.push_back(process_bins.at(i))
            for(int j = 0; j<neighbor_list_collect.size();j++)
            {
                particle_list_collect = bin_map.at(neighbor_list_collect.at(i));
            }

            number_of_interacting_particles = particle_list_collect.size();

            for(int k =0; k<number_of_interacting_particles; k++)
            {
                particles[particle_list_collect.at(k)].ax = particles[particle_list_collect.at(k)].ay = 0;
                    for (int l = 0; l < n; l++ )
                    {
                        apply_force( particles[particle_list_collect.at(k)], particles[particle_list_collect.at(l)],&dmin,&davg,&navg);
                    }
            }


            neighbor_list_collect.clear();
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
           if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }

          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
           if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }

    simulation_time = read_timer( ) - simulation_time;
    
    //printf( "n = %d, simulation time = %g seconds", n, simulation_time);
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //

    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
