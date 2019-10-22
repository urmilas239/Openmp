#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "common.h"
#include "omp.h"


std::vector<std::vector<int> > bin_map;
std::vector<std::vector<int> > neighbor_bins;

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0,numthreads;
    double davg,dmin, absmin=1.0, absavg=0.0;


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
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    


    
    

    set_size( n );

    //call the function to set bin variables
    set_bin_count(n);
    
    neighbor_bins = initialize_neighbor_bins();
    


    //Initialize particles and assign bins

        bin_map = initialize_bin_vector();
        init_particles1(n, particles,bin_map);
        //init_particles( n, particles);
        //bin_particles( n, particles , bin_map);
    
    

    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    omp_set_num_threads(5);

    //int total_bin_count = bin_map.size();
	
    #pragma omp parallel private(dmin)  shared(neighbor_bins, bin_map)
    {


    numthreads = omp_get_num_threads();

    //std::cout<<"numthreads:::"<<numthreads<<std::endl;
    for( int step = 0; step < NSTEPS; step++ )
    {
         //printf( ":::::::::::::IN TIME STEP::::::::::::::::::::::::::::::::::::: %d\n" , step);

	navg = 0;
        davg = 0.0;
	dmin = 1.0;
        //
        //  compute forces
        //

         
           // std::vector<int> particle_ids;
           // std::vector<int> neighbor_bins_list;
            int bin_index;
             //TODO: OpenMP - make particle_ids, neighbor_bins_list private, neighbor_bins;
			 //bin_map , particles- shared
			//std::cout<<"bin_map.size():::"<<bin_map.size()<<std::endl;
            //std::cout<<"neighbor_bins:::"<<neighbor_bins.size()<<std::endl;
            //#pragma omp parallel for firstprivate(neighbor_bins, bin_map)
            #pragma omp for reduction (+:navg) reduction(+:davg)  
            for( int i = 0; i < n; i++ )
            {
                particles[i].ax = particles[i].ay = 0;
                bin_index = compute_bin_index_from_xy(particles[i].x, particles[i].y);

                std::vector<int>  particle_ids ; //= bin_map.at(bin_index);
                std::vector<int> neighbor_bins_list = neighbor_bins.at(bin_index);
                neighbor_bins_list.push_back(bin_index);
				
                for(int k = 0; k<neighbor_bins_list.size();k++)
                {
                    if(neighbor_bins_list.at(k) != -1)
                    {
                        particle_ids = bin_map.at(neighbor_bins_list.at(k));
                         for(int j = 0; j < particle_ids.size(); j++)
                        {
                            apply_force( particles[i], particles[particle_ids.at(j)],&dmin,&davg,&navg);
                        }
                         particle_ids.clear();
                    }
                    
                }
                neighbor_bins_list.clear();
            }
         
        std::cout<<"After FOR::: in thread  "<<omp_get_thread_num()<<std::endl;
        

        #pragma omp barrier 
        {
            std::cout<<"All threads finished gather here  "<<omp_get_thread_num()<<std::endl;
        }
           
        
 
        //
        //  move particles
        //
		#pragma omp for
        for( int i = 0; i < n; i++ ) 
        {
            //move( particles[i]);
            move1( particles[i], bin_map);            
        }
           	


        #pragma omp barrier 
        {
            std::cout<<"All threads finished gather here after move "<<omp_get_thread_num()<<std::endl;
        }
        
        //if(0)
       // {
        /*#pragma omp master
        {
            bin_map = initialize_bin_vector();
            bin_particles( n, particles , bin_map); 
        }*/

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
           #pragma omp master 
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }

          #pragma omp critical
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
} //omp parallel region ends
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
