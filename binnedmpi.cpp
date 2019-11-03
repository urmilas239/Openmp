#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <cmath>
#include <signal.h>
#include <unistd.h>
#include "common.h"

using std::vector;


#define _cutoff 0.01    //Value copied from common.cpp
#define _density 0.0005

double bin_size1, grid_size;
int bin_count;


inline bool cmpf(double A, double B, double epsilon = 0.00000005f)
{
return (fabs(A - B) < epsilon);
}

inline void build_bins(vector<bin_t>& bins, particle_t* particles, int n)
{
    grid_size = sqrt(n * _density);
    bin_size1 = _cutoff;
    bin_count = int(grid_size / bin_size1) + 1; // Should be around sqrt(N/2)

    // printf("Grid Size: %.4lf\n", grid_size);
    // printf("Number of Bins: %d*%d\n", bin_count, bin_count);
      printf("bin_count: %d\n", bin_count);
    // printf("Bin Size: %.2lf\n", bin_size);
    // Increase\Decrease bin_count to be something like 2^k?

    bins.resize(bin_count * bin_count);

    for (int i = 0; i < n; i++)
    {
        int x = int(particles[i].x / bin_size1);
        int y = int(particles[i].y / bin_size1);
        bins.at(x*bin_count + y).push_back(particles[i]);
    }
}



void bin_particle(particle_t& particle, vector<bin_t>& bins)
{
    int x = particle.x / bin_size1;
    int y = particle.y / bin_size1;
    bins.at(x*bin_count + y).push_back(particle);
}


inline void get_neighbors(int i, int j, vector<int>& neighbors)
{
    neighbors.erase(neighbors.begin(), neighbors.end());
   // printf("get_neighbors start");
    int neighbor_count = 0;
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            if (dx == 0 && dy == 0)
                continue;
            if (i + dx >= 0 && i + dx < bin_count && j + dy >= 0 && j + dy < bin_count) {
                int index = (i + dx) * bin_count + j + dy;
                neighbors.push_back((i + dx));
                neighbor_count++;
            }
        }
    }
    //neighbors.resize(8);
    //printf("get_neighbors end");
}




inline int which_process(int bin_id, int x_bins_per_proc, int n_proc )
{
    int process_id = -1;
    for(int i=0; i < n_proc; i++)
    {
        if(bin_id >= i*x_bins_per_proc )//&& bin_id <= (i+1)*x_bins_per_proc)
        {
            process_id = i;
            break;
        }
    }
    int y =  min(bin_id / x_bins_per_proc, n_proc-1);
   // printf("which_process y %d:: process_id::%d\n",y, process_id);
    return y;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    //signal(SIGSEGV, sigsegv);

    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg;

    //
    //  process command line parameters
    //
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

    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


    particle_t *particles = new particle_t[n];

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );

    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);

    vector<int> current_neighbors(8);

    vector<std::vector<particle_t>> bins;
     std::vector<particle_t> receive_updated_particles(2);

    build_bins(bins, particles, n);




    int x_bins_per_proc = bin_count / n_proc;

    // although each worker has all particles, we only access particles within
    // my_bins_start, my_bins_end.

    int my_bins_start = x_bins_per_proc * rank;
    int my_bins_end = x_bins_per_proc * (rank + 1);

    if (rank == n_proc - 1)
        my_bins_end = bin_count;



     int *receive_size = (int*)malloc(n_proc*sizeof (int));
    // printf("worker %d: from %d to %d.\n", rank, my_bins_start, my_bins_end);
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        // if( find_option(argc, argv, "-no" ) == -1 )
        //     if( fsave && (step%SAVEFREQ) == 0 )
        //         save( fsave, n, particles );

        // compute local forces
        for (int i = my_bins_start; i < my_bins_end; ++i) {
            for (int j = 0; j < bin_count; ++j) {
                bin_t& vec = bins[i * bin_count + j];

                for (int k = 0; k < vec.size(); k++)
                    vec[k].ax = vec[k].ay = 0;

                for (int dx = -1; dx <= 1; dx++)   //Search over nearby 8 bins and itself
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        if (i + dx >= 0 && i + dx < bin_count && j + dy >= 0 && j + dy < bin_count)
                        {
                            bin_t& vec2 = bins[(i+dx) * bin_count + j + dy];
                            for (int k = 0; k < vec.size(); k++)
                                for (int l = 0; l < vec2.size(); l++)
                                    apply_force( vec[k], vec2[l], &dmin, &davg, &navg);
                        }
                    }
                }
            }
        }

        if (find_option( argc, argv, "-no" ) == -1) {
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
          if (rank == 0){
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin)
                absmin = rdmin;
          }
        }

        // move, but not rebin
        bin_t local_move;
        std::vector<std::vector<particle_t>> remote_move(n_proc);

        for (int i = my_bins_start; i < my_bins_end; ++i)
        {
            for (int j = 0; j < bin_count; ++j)
            {
                std::vector<particle_t>& bin = bins.at(i * bin_count + j);
                int b_size = bin.size();
                int k = 0;
                for (k=0; k < b_size; ) {
                    particle_t pt = bin[k];
                    move(bin[k]);

                    int x = int(bin[k].x / bin_size1);
                    int y = int(bin[k].y / bin_size1);

                    //If the new bin id is in the current process
                    if (my_bins_start <= x && x < my_bins_end) {

                        //If the new bin id same as current bin id
                        if (x == i && y == j)
                        {
                            //No change dont do anything
                             k++;
                        }

                        else {
                            //moved locally in process
                            local_move.push_back(bin[k]);

                            //Put last bin in at k and resize bin later to remove bin[k]
                            bin[k] = bin[--b_size];
                        }
                    } else {

                        //Which process has the particle moved to?
                        int proc_to_update = which_process(x,x_bins_per_proc, n_proc);
                       // printf("RANK %d, proc_to_update %d , x: %d, x_bins_per_proc:%d, n_proc:%d::\n", rank, proc_to_update, x,  x_bins_per_proc, n_proc);
                        remote_move.at(proc_to_update).push_back(pt);
                        //get neighbors of bin, if neighbors bellong to different bin update them.

                        get_neighbors(x,y, current_neighbors);
                        for(int n_index = 0; n_index < current_neighbors.size(); n_index++)
                        {
                            if(proc_to_update!= which_process(current_neighbors.at(n_index), bin_count,n_proc) )
                            {
                                remote_move.at(which_process(current_neighbors.at(n_index), bin_count,n_proc)).push_back(pt);
                            }
                        }
                        current_neighbors.clear();

                        //remove the particle that has moved.
                        bin[k] = bin[--b_size];
                    }
                }
                bin.resize(k);
            }
        }
            //printf("Moving particles %d::", rank);
        for (int i = 0; i < local_move.size(); ++i) {
            bin_particle(local_move[i], bins);
        }

        if(1)
        {
            int size;
        for (int i = 0; i < remote_move.size(); ++i) {
           if(i != rank)
            {
                 size = 0;
                if(!(remote_move.at(i).empty() ))
                {
                   size = remote_move.at(i).size();

                }
              //  printf("From %d -> to %d, Value sent::: %d\n", rank, i, size);
               // MPI_Isend(&size, 1, MPI_INT, i, rank*10, MPI_COMM_WORLD, &request);
                MPI_Send(&size, 1, MPI_INT, i, rank*10, MPI_COMM_WORLD);

            }
        }

        MPI_Status status;
        //Receive sizes and wait to receive size from particles that only send size != -1
            for(int i =0; i<n_proc; i++)
            {
                if(i != rank)
                {
                   // MPI_Irecv(&receive_size[i], 1, MPI_INT, i, i*10, MPI_COMM_WORLD, &request);
                     MPI_Recv(&receive_size[i], 1, MPI_INT, i, i*10, MPI_COMM_WORLD, &status);
                   // printf("in %d recieve_d value:: %d\n", rank,receive_size[i]);
                   // printf("From %d -> to %d, Value Received::: %d\n\n\n", i, rank,receive_size[i]);

                }

            }
        }

        //Send sizes ends.

        /*if(1)
        {
            //using the last process for gather and update since it could potentially have been allocated less bins than others
            if(rank = n_proc-1)
            {
                //All gather sizes.
                MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

            }
        }*/


        if(1)
        {
            /**Send updated particles*/
                        int size;
                        for(int i =0; i<n_proc; i++)
                        {
                            if(i != rank)
                            {
                                if(!(remote_move.at(i).empty()) && remote_move.at(i).size()>0)
                                {
                                   size = remote_move.at(i).size();
                                    //MPI_Isend(&particles_updates_after_move, particles_updates_after_move.size(), PARTICLE, i, rank*100, MPI_COMM_WORLD, &request);
                                   MPI_Send(remote_move.at(i).data(), size, PARTICLE, i, rank*100, MPI_COMM_WORLD);
                                   remote_move.at(i).erase(remote_move.at(i).begin(), remote_move.at(i).end());
                                   //printf("RANK %d: Sent %d particles to %d: \n", rank, size,i);

                                }
                            }

                        }
                        size = 0;
                        MPI_Status status;


                        /*int prev = receive_size[0] ;
                        int curr;
                        int max = 0;
                        for(int i =1; i<n_proc; i++)
                        {

                            curr = receive_size[i];
                            if(curr > prev)
                            {
                                max = i;
                                prev = curr;
                            }
                        }*/


                        for(int i =0; i<n_proc; i++)
                        {
                            //printf("All sent::: Receiving.....\n");
                            if(i != rank && receive_size[i] > 0)
                            //if(0)
                            {

                                size = receive_size[i];
                               // receive_updated_particles.erase(receive_updated_particles.begin(), receive_updated_particles.end());
                                if(receive_updated_particles.size() < size)
                                {
                                   receive_updated_particles.resize(receive_size[i]+1);
                                }

                                std::vector<particle_t> receive_updated_particles1(size);
                                //MPI_Irecv(&receive_particles, receive_size[i], PARTICLE, i, i*100, MPI_COMM_WORLD, &request);
                                MPI_Recv(receive_updated_particles.data(), size, PARTICLE, i, i*100, MPI_COMM_WORLD, &status);
                                // printf("RANK %d , receive_size[i]::: %d , receive_updated_particles:: %d\n", rank, size,receive_updated_particles.size() );
                                for(int p_index = 0; p_index <receive_updated_particles.size(); p_index++)
                                {
                                    particle_t pt = receive_updated_particles.at(p_index);

                                   // int old_bin_index = compute_bin_index_from_xy(pt.x, pt.y,bin_count);

                                    int x = int(pt.x / bin_size1);
                                    int y = int(pt.y / bin_size1);
                                    int old_bin_index = x*bin_count + y;

                                    //printf("old_bin_index::: %d \n",old_bin_index);

                                      std::vector<particle_t> &bin = bins.at(old_bin_index);
                                      int b_size = bin.size();
                                      //int g ;
                                      for(int g =0; g<b_size;g++)
                                      {
                                        //if((double)bin[g].x == (double)pt.x && (double)bin[g].y == (double)pt.y)
                                        if(cmpf(bin[g].x, pt.x) && cmpf(bin[g].y, pt.y) )
                                        {
                                            //bin[g] = bin[b_size - 1];
                                            //printf(":::MATCH FOUND:::::");
                                            bin.erase(bin.begin()+g);
                                           break;

                                        }

                                      }

                                     //move(pt);
                                    //bin_particle(pt, bins);
                                   // particles[(int) pt.index] = pt;
                                }





                            }

                        }
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
    if( fsave )
        fclose( fsave );

    MPI_Finalize( );

    return 0;
}
