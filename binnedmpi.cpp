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
    //int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    //int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    
    int *partition_sizes= (int*) calloc( n_proc , sizeof(int));
    int *partition_offsets = (int*) calloc( (n_proc+1) , sizeof(int) );


    MPI_Datatype PARTICLE;
    //Need 6 double and one int. Since we cannot make MPI_Datatype of multiple types using double
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    

    
    
    int number_of_interacting_particles = 0;
    set_size( n );

    //call the function to set bin variables
    set_bin_count(n);

    // All initialization calls happens in rank 0
     
    bin_map = initialize_bin_vector();
    
    neighbor_bins = initialize_neighbor_bins();


    //init_particles1(n, particles,bin_map);
    init_particles( n, particles);
    bin_particles( n, particles , bin_map);
    process_bins=assign_bins_to_current_process_mpi(n_proc, rank, bin_map, bin_process_map, number_of_interacting_particles);
    border_neighbors = get_boundary_bins_for_curr_process(process_bins, neighbor_bins);
    get_num_of_particles_in_each_process(n_proc, bin_map, &partition_offsets, &partition_sizes);


    //particles on whom apply_force was called in this process
    particle_t *particles_acted_upon = (particle_t*) malloc( number_of_interacting_particles * sizeof(particle_t) );
    int particle_index; 


    std::vector<int> neighbor_list_collect; 
    std::vector<int> particle_list_collect; 
    std::vector<int> particle_list_temp;

    

    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);
    
     int partition_size_per_rank = 0;
     int no_of_particles_in_bin =0;
     int prev_k_index = 0;
     int current_k_index = 0;
    //std::cout<<"numthreads:::"<<numthreads<<std::endl;
    for( int step = 0; step < NSTEPS; step++ )
    {
        MPI_Barrier(MPI_COMM_WORLD);
         //printf( ":::::::::::::IN TIME STEP::::::::::::::::::::::::::::::::::::: %d\n" , step);
           
        particle_index = 0;
        for(int i=0;i<process_bins.size();i++)
        {
            //collect neighbor particle for the current bin.
            neighbor_list_collect = neighbor_bins.at(process_bins.at(i));
            neighbor_list_collect.push_back(process_bins.at(i));
            for(int j =0; j<neighbor_list_collect.size();j++)
            {
               if(neighbor_list_collect.at(j) != -1)
               {
                   particle_list_temp = bin_map.at(neighbor_list_collect.at(j));
                   particle_list_collect.insert(particle_list_collect.end(), particle_list_temp.begin(), particle_list_temp.end());
                   particle_list_temp.clear();
               }

            }

            no_of_particles_in_bin = particle_list_collect.size();

            particle_list_temp = bin_map.at(process_bins.at(i));
            particle_list_collect.insert(particle_list_collect.end(), particle_list_temp.begin(), particle_list_temp.end());

            for(int k =0; k<particle_list_temp.size(); k++)
            {
                current_k_index = particle_list_collect.at(k);
                particles[current_k_index].ax = particles[current_k_index].ay = 0;
                    for (int l = 0; l < no_of_particles_in_bin; l++ )
                    {
                        apply_force( particles[current_k_index], particles[particle_list_collect.at(l)],&dmin,&davg,&navg);


                    }
                    //prev_k_index = current_k_index;
                    particles_acted_upon[particle_index] = particles[current_k_index];
                    particle_index++;
            }


            neighbor_list_collect.clear();
            particle_list_collect.clear();

            //std::cout<<"particle_index:: "<<particle_index<<"  ,number_of_interacting_particles::"<<number_of_interacting_particles<<std::endl;
            if( find_option( argc, argv, "-no" ) == -1 )
            {
              
              MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
              MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
              MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

     
              if (rank == 0){
                //
                // Computing statistical data
                //
                if (rnavg) {
                  absavg +=  rdavg/rnavg;
                  nabsavg++;
                }
                if (rdmin < absmin) absmin = rdmin;
              }
            }


            for( int i = 0; i < number_of_interacting_particles; i++ )
            {
                 move( particles_acted_upon[i] );
            }

           // MPI_Allgatherv( particles_acted_upon, number_of_interacting_particles, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
           }
     }

       simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
   // free( partition_offsets );
    //free( partition_sizes );
    //free( particles );
    //free(particles_acted_upon);
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
    
}
