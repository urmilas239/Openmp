#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "common.h"


std::vector<std::vector<int> > bin_map;
std::vector<std::vector<int> > neighbor_bins;
std::vector<std::vector<int> > process_bins;
std::vector<std::vector<int> > border_neighbors;

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0,numthreads;
    double davg,dmin, absmin=1.0, absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 


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
    
    int n = read_int( argc, argv, "-n", 1000 );

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
    particle_t *particles_to_send;
    particle_t *particles_bin_map;
    particle_t *particles_neighbor_map;
    int *partition_sizes; //= (int*) malloc( n_proc * sizeof(int));
    int *partition_offsets; // = (int*) malloc( (n_proc+1) * sizeof(int) );


    MPI_Datatype PARTICLE;
    //Need 6 double and one int. Since we cannot make MPI_Datatype of multiple types using double
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    



    MPI_Datatype PARTICLE_BIN_MAP;
    MPI_Type_contiguous( 11, MPI_INT, &PARTICLE );
    MPI_Type_commit( &PARTICLE_BIN_MAP );



    MPI_Datatype NEIGHBOR_BIN_MAP;
    MPI_Type_contiguous( 3, MPI_INT, &NEIGHBOR_BIN_MAP );
    MPI_Type_commit( &NEIGHBOR_BIN_MAP );
    
    

    set_size( n );

    //call the function to set bin variables
    set_bin_count(n);
    bin_map = initialize_bin_vector();
    
    neighbor_bins = initialize_neighbor_bins();


    //init_particles1(n, particles,bin_map);
    init_particles( n, particles);
bin_particles( n, particles , bin_map);
            process_bins = assign_bins_to_processes_mpi(n_proc, bin_map);
            border_neighbors = get_boundary_bins(process_bins, neighbor_bins);
    
    


    //Initialize particles and assign bins

        
        //init_particles( n, particles);
        //bin_particles( n, particles , bin_map);
    
    

    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

   particle_t *particles_to_send;
   particle_bin_mapping *pbm;
   neighbor_bin_mapping *n_bins;
	
   

    

    //std::cout<<"numthreads:::"<<numthreads<<std::endl;
    for( int step = 0; step < NSTEPS; step++ )
    {
         //printf( ":::::::::::::IN TIME STEP::::::::::::::::::::::::::::::::::::: %d\n" , step);

         MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0)
        {
            if(bin_map.empty())
            {
                bin_map = initialize_bin_vector();
                bin_particles( n, particles , bin_map);
            }
            

             for(i=0;i<n_proc;i++)
             {
                
                form_particles_array_for_MPI(process_bins.at(i), bin_map, neighbor_bins, particles_to_send, particles_neighbors, pbm, n_bins, particles,partition_sizes,partition_offsets);


                 MPI_Send(sizeof(particles_to_send)/sizeof(particle_t),1,MPI_INT,i,0,MPI_COMM_WORLD);
                MPI_Send(sizeof(pbm)/sizeof(particle_bin_mapping),1,MPI_INT,i,1,MPI_COMM_WORLD);
                MPI_Send(sizeof(n_bins)/sizeof(neighbor_bin_mapping),1,MPI_INT,i,2,MPI_COMM_WORLD);

                MPI_Send(particles_to_send,sizeof(particles_to_send)/sizeof(particle_t),PARTICLE,i,3,MPI_COMM_WORLD);
                MPI_Send(pbm,sizeof(pbm)/sizeof(particle_bin_mapping),PARTICLE_BIN_MAP,i,4,MPI_COMM_WORLD);
                MPI_Send(n_bins,sizeof(n_bins)/sizeof(neighbor_bin_mapping),NEIGHBOR_BIN_MAP,i,5,MPI_COMM_WORLD);

                //send to current_process
             }
         
        }


       
        

        int MPI_Recv(num_of_particles_in_proc, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_Status * status);
        int MPI_Recv(num_of_bins_in_proc, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_Status * status);
        int MPI_Recv(num_of_neighbors_in_proc, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_Status * status);



        int MPI_Recv(particles_in_proc, num_of_particles_in_proc, PARTICLE, 0, 3, MPI_COMM_WORLD, MPI_Status * status);
        int MPI_Recv(bins_in_proc, num_of_bins_in_proc, PARTICLE_BIN_MAP, 0, 4, MPI_COMM_WORLD, MPI_Status * status);
        int MPI_Recv(neighbors_in_proc, num_of_neighbors_in_proc, NEIGHBOR_BIN_MAP, 0, 5, MPI_COMM_WORLD, MPI_Status * status);





        //Add bin ids (in current process and the boundary bins) from to a map - for ease of searching
        map<int, neighbor_bin_mapping> neighbor_map; 


        particle_bin_mapping pb;
        neighbor_bin_mapping nb;
        for(int temp_index = 0;temp_index < num_of_bins_in_proc; temp_index++)
        {
            pb = bins_in_proc[temp_index];
            nb.neighbor_bin_id = pb.bin_id;
            nb.num_particles = pb.num_particles;
            nb.particle_offset = pb.particle_offset;
            neighbor_map.insert ( std::pair<int, neighbor_bin_mapping>(pb.bin_id,nb) );

        }


        for(int temp_index = 0;temp_index < num_of_neighbors_in_proc; temp_index++)
        {
            nb = neighbors_in_proc[temp_index];
            neighbor_map.insert ( std::pair<int, neighbor_bin_mapping>(nb.neighbor_bin_id,nb) );

        }


 
    	navg = 0;
            davg = 0.0;
    	dmin = 1.0;
        //
        //  compute forces
        //
       
          
            int current_bin_index;
            int p_offset;
            int p_count;
            int n_ids[8];
            particle_t *particles_updated;
            //int neighbor_id;
            
            for( int i = 0; i < bins_in_proc; i++ )
            {
                pb = bins_in_proc[i];
                current_bin_index = pb.bin_id;
                p_offset = pb.particle_offset;
                p_count = pb.num_particles;
                n_ids = pb.neighbor_id;

                for(int j = 0; j < p_count ; j++)
                {
                    
                    for(int k =0; k < 8; k++)
                    {

                        if(n_ids[k] != -1)
                        {
                             nb =  neighbor_map.find(n_ids[k]);
                            for(int k =0; k < nb.num_particles; k++)
                            {
                                apply_force( particles_in_proc[p_offset+j], particles_in_proc[nb.particle_offset + k], &dmin, &davg, &navg );
                            }

                        }
                    }
                }
            }
         
       
        

    
           
        
 
        //
        //  move particles
        //
        int temp_offset = 0;
        for( int i = 0; i < bins_in_proc; i++ )
        {
            pb = bins_in_proc[i];
            current_bin_index = pb.bin_id;
            p_offset = pb.particle_offset;
            p_count = pb.num_particles;

            for(int j = 0; j < p_count ; j++)
            {
                move( particles_in_proc[p_offset+j]);
                particles_updated[temp_offset] = particles_in_proc[p_offset+j];
                temp_offset++;
            }
        }
       
       //send results to process in rank 0
       MPI_Send(particles_updated,sizeof(particles_updated)/sizeof(particle_t),PARTICLE,0,0,MPI_COMM_WORLD);

       MPI_Barrier(MPI_COMM_WORLD);


       if(rank == 0)
       {
        //receive particles_in_proc from each process;
        MPI_Allgatherv( particles_updated, sizeof(particles_updated)/sizeof(particle_t), PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );

        //Reform particles array


        //Rebin 
        bin_map.clear();
        
        //bin_particles( n, particles , bin_map); 
        //process_bins.clear();

       }
        
        //if(0)
       // {
         
        
            //
            
            
        //std::cout<<"After rebinning::: "<<omp_get_thread_num()<<std::endl;

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
