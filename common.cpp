#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <set>
#include <iostream>
#include <algorithm> 
#include "common.h"


double size;
int num_of_bins_x;
int num_of_bins_y;
int total_bin_count;
int bin_density;
double bin_size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005



//
//  timer
//
double read_timer()
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

std::vector<std::vector<int> > initialize_bin_vector()
{
    std::vector<std::vector<int> > bin_map(total_bin_count, std::vector<int>()); 
    return bin_map;
}




//Divides bins(hence particles) across process
std::vector<int > assign_bins_to_current_process_mpi11(int num_of_processes, int current_process, std::vector<std::vector<int> > bin_map, std::vector<int> &bin_process_map, int &num_of_particles_in_curr_process)
{
    //std::cout<<"RANK: "<<current_process<<" assign_bins_to_current_process_mpi START::: "<<std::endl;
    //total bins
    int bin_count = bin_map.size();
    int num_of_particles_in_curr_process_local = 0;
   


    //int num_bins_per_process =  ceil(bin_count/num_of_processes);
    int num_bins_per_process =  ((bin_count + num_of_processes - 1) / num_of_processes);
    std::vector<int> process_bins; 
    int assigned_bin_count = 0;
    int current_process_id = 0;
    int total_particles = 0;
    //std::cout<<"bin_count::: "<<bin_count<<" ,num_bins_per_process ::: "<< num_bins_per_process<< " , num_of_processes::: "<< num_of_processes<< std::endl;

    for(int bin_idex =(num_bins_per_process*current_process); bin_idex < bin_count; bin_idex++)
    {
        //std::cout<<" current_process_id ::: "<< current_process_id<<std::endl;
        bin_process_map.push_back(current_process_id);
        total_particles += bin_map.at(bin_idex).size();
        if(current_process_id == current_process)
        {
            process_bins.push_back(bin_idex);
            num_of_particles_in_curr_process += bin_map.at(bin_idex).size();
           // printf("RANK111 %d assign_bins_to_current_process_mpi :: bin_map.at(bin_idex).size() %d \n",current_process,  bin_map.at(bin_idex).size());

          // printf("RANK111 %d assign_bins_to_current_process_mpi ::total_particles %d \n",current_process, total_particles);
        }
        
        assigned_bin_count++;
        //std::cout<<"assigned_bin_count ::: "<< assigned_bin_count<<std::endl;
        if(assigned_bin_count >= num_bins_per_process)
        {

            //printf("RANK %d assigned %d bins to Process %d\n",current_process, assigned_bin_count, current_process_id);
          // printf("RANK %d assigned %d particles to Process %d\n",current_process, num_of_particles_in_curr_process_local, current_process_id);

            assigned_bin_count = 0;
            current_process_id++;
            num_of_particles_in_curr_process=0;
        }
        
    }


   // printf("RANK %d assigned %d bins to Process %d\n",current_process, assigned_bin_count, current_process_id);
    //printf("RANK %d assigned %d particles to Process %d\n",current_process, num_of_particles_in_curr_process_local, current_process_id);


  //  printf("RANK %d  assign_bins_to_current_process_mpi ::process_bins size %d \n", current_process, process_bins.size());
//
    //printf("RANK %d  assign_bins_to_current_process_mpi ::num_of_particles_in_curr_process %d \n", current_process, num_of_particles_in_curr_process);
   // printf("RANK %d  assign_bins_to_current_process_mpi ::total_particles %d \n", current_process, total_particles);

    // process_bins.push_back(num_of_particles_in_curr_process);
    //std::cout<<"RANK: "<<current_process<<" assign_bins_to_current_process_mpi END::: "<<std::endl;
    return process_bins;
}





//Divides bins(hence particles) across process
std::vector<int > assign_bins_to_current_process_mpi( int num_of_processes, int current_process, std::vector<std::vector<int> > bin_map, std::vector<int> &bin_process_map, int &num_of_particles_in_curr_process)
{
   // std::cout<<"RANK: "<<current_process<<" assign_bins_to_current_process_mpi START::: "<<std::endl;
    //total bins
    int bin_count = bin_map.size();
    int num_of_particles_in_curr_process_local = 0;



    //int num_bins_per_process =  ceil(bin_count/num_of_processes);
    int num_bins_per_process =  ((bin_count + num_of_processes - 1) / num_of_processes);
    std::vector<int> process_bins;
    int assigned_bin_count = 0;
    int current_process_id = 0;
    int total_particles = 0;
    //std::cout<<"bin_count::: "<<bin_count<<" ,num_bins_per_process ::: "<< num_bins_per_process<< " , num_of_processes::: "<< num_of_processes<< std::endl;

    for(int bin_idex =0; bin_idex < bin_count; bin_idex++)
    {

        total_particles += bin_map.at(bin_idex).size();
        num_of_particles_in_curr_process_local += bin_map.at(bin_idex).size();
         bin_process_map.push_back(current_process_id);

        assigned_bin_count++;
        //std::cout<<"assigned_bin_count ::: "<< assigned_bin_count<<std::endl;
        if(current_process_id == current_process)
        {
            process_bins.push_back(bin_idex);

        }


        if(assigned_bin_count >= num_bins_per_process || (bin_count-1) == bin_idex)
        {
            //printf("RANK %d assigned %d bins to Process %d\n",current_process, assigned_bin_count, current_process_id);
            ///printf("RANK %d assigned %d particles to Process %d\n",current_process, num_of_particles_in_curr_process_local, current_process_id);
            //printf("RANK %d assigned %d particles to Process1 %d\n",current_process, num_of_particles_in_curr_process, current_process_id);
            if(current_process_id == current_process)
            {
                 num_of_particles_in_curr_process = num_of_particles_in_curr_process_local;
            }

            assigned_bin_count = 0;
            //num_of_particles_in_curr_process = 0;
            num_of_particles_in_curr_process_local = 0;
            current_process_id++;
            //num_of_particles_in_curr_process=0;
        }

    }
    //printf("RANK %d assigned %d bins to Process %d\n",current_process, assigned_bin_count, current_process_id);
    //printf("RANK %d assigned %d particles to Process %d\n",current_process, num_of_particles_in_curr_process_local, current_process_id);
    //printf("RANK %d assigned %d particles to Process1 %d\n",current_process, num_of_particles_in_curr_process, current_process_id);

   // std::cout<<"RANK: "<<current_process<<" assign_bins_to_current_process_mpi END::: "<<std::endl;
    return process_bins;
}




void get_num_of_particles_in_each_process(int num_of_processes, std::vector<std::vector<int> > bin_map, int** partition_offset, int** partition_sizes)
{

    //std::cout<<"get_num_of_particles_in_each_process START::: "<<std::endl;
    //total bins
    int bin_count = bin_map.size();

    //int num_bins_per_process =  ceil(bin_count/num_of_processes);
    int num_bins_per_process =  ((bin_count + num_of_processes - 1) / num_of_processes);
    std::vector<int > process_bins(num_of_processes); 
    int assigned_bin_count = 0;
    int current_process_id = 0;
    int num_of_particles_in_curr_process = 0;
    //int *partition_offset1 = new int[num_of_processes+1];
    //int *partition_sizes1 = new int[num_of_processes];
    
    for(int bin_idex =0; bin_idex < bin_count; bin_idex++)
    {
        num_of_particles_in_curr_process += bin_map.at(bin_idex).size();
        assigned_bin_count++;
        
        if(assigned_bin_count > num_bins_per_process)
        //if(1)
        {

            *partition_offset[current_process_id+1] = *(partition_offset[current_process_id]) + num_of_particles_in_curr_process + 1;
            *partition_sizes[current_process_id] = num_of_particles_in_curr_process;

            //printf("current_process_id :: %d ,partition_sizes %d, partition_offsets %d \n",current_process_id, *partition_sizes[current_process_id], *partition_offset[current_process_id]);
             assigned_bin_count = 0;
            current_process_id++;
            process_bins.push_back(num_of_particles_in_curr_process);
            num_of_particles_in_curr_process = 0;

        }

    }

    //current_process_id++;
   // printf("num_of_particles_in_curr_process %d num_of_particles_in_curr_process %d \n", num_of_particles_in_curr_process, num_of_particles_in_curr_process);

    process_bins.push_back(num_of_particles_in_curr_process);
    *partition_offset[current_process_id+1] = *(partition_offset[current_process_id]) + num_of_particles_in_curr_process + 1;
    *partition_sizes[current_process_id] = num_of_particles_in_curr_process;

    //printf("\n current_process_id :: %d ,partition_sizes %d, partition_offsets %d \n",current_process_id, *partition_sizes[current_process_id], *partition_offset[current_process_id]);
     assigned_bin_count = 0;

    /*for(int i = 0; i < num_of_processes; i++)
    {
        std::cout<<"i :: "<<i<<" ,partition_sizes "<<*partition_sizes[i]<<", partition_offsets "<<*partition_offset[i]<<std::endl;
    }*/
    //std::cout<<"partition_offsets last   "<<*(partition_offset[num_of_processes])<<std::endl;
    //std::cout<<"get_num_of_particles_in_each_process END::: "<<std::endl;
    //return process_bins;

}



//Divides bins(hence particles) across process
std::vector<std::vector<int> > assign_bins_to_processes_mpi(int num_of_processes, std::vector<std::vector<int> > bin_map)
{
     //std::cout<<"assign_bins_to_processes_mpi START::: "<<std::endl;
    //total bins
    int bin_count = bin_map.size();

    //int num_bins_per_process =  ceil(bin_count/num_of_processes);
    int num_bins_per_process =  ((bin_count + num_of_processes - 1) / num_of_processes);
    std::vector<std::vector<int> > process_bins(num_of_processes, std::vector<int>()); 
    int assigned_bin_count = 0;
    int current_process_id = 0;

    for(int bin_idex =0; bin_idex < bin_count; bin_idex++)
    {
        //std::cout<<" current_process_id ::: "<< current_process_id<<std::endl;
        process_bins.at(current_process_id).push_back(bin_idex);
        assigned_bin_count++;
        //std::cout<<"assigned_bin_count ::: "<< assigned_bin_count<<std::endl;
        if(assigned_bin_count > num_bins_per_process)
        {
            assigned_bin_count = 0;
            current_process_id++;
        }
        
    }
    //std::cout<<"assign_bins_to_processes_mpi END::: "<<std::endl;
    return process_bins;
}


void form_particles_array_for_MPI(std::vector<int>  &bin_ids,
    std::vector<int>  &boundary_bin_ids,
    std::vector<std::vector<int> > &bin_map,
    std::vector<std::vector<int> > &neighbor_bins,
    particle_t *particles_to_send, 
    particle_bin_mapping *pbm,
    neighbor_bin_mapping *n_bins,
    particle_t *particles, 
    int &partition_size_per_rank)
{
    //std::cout<<"form_particles_array_for_MPI START::: "<<std::endl;
    //std::vector<int>  bin_ids = process_bins.at(process_id);
    int num_of_bins_in_curr_proc = bin_ids.size();
    int particles_per_bin;
    std::vector<int>  particle_ids;
    std::vector<int>  neighbor_ids;
    int neighbors_per_bin;
    int bin_id;
    int particle_index = 0;
    int particle_neighbor_index = 0;


    std::vector<int>  particles_in_neighbor;

    int neighbors_id;



    for(int i = 0; i < num_of_bins_in_curr_proc; i++)
    {

        //iterate through bins in bin_ids, get particle id. fetch those from particles_t array and set it into particles_to_send array
        //update pbm to contatin offset and #of particle information for each bin
        bin_id = bin_ids.at(i);
        particle_ids = bin_map.at(bin_id);

        if(particle_ids.empty())
        {
             pbm[i].bin_id = bin_id;
        pbm[i].num_particles = 0;
        pbm[i].particle_offset = particle_index;
            continue;
        }
        particles_per_bin = particle_ids.size();


        pbm[i].bin_id = bin_id;
        pbm[i].num_particles = particles_per_bin;
        pbm[i].particle_offset = particle_index;
        partition_size_per_rank += particles_per_bin;



        for(int j = 0; j<particles_per_bin; j++ )
        {
            particles_to_send[particle_index] = particles[particle_ids.at(j)];
            particle_index++;
        }
        particle_ids.clear();

        
        //get neighbors of the current bin
        neighbor_ids = neighbor_bins.at(bin_id);
        neighbors_per_bin = particle_ids.size(); //number of neighbors
        
        for(int j = 0; j<neighbors_per_bin; j++ )
        {
            neighbors_id = neighbor_ids.at(j);
            pbm[i].neighbor_id[j] = neighbors_id;
            
         }
         

        
    }

    //Assumption num_boundary_bins cannot be zero since all process will have boundary bins
    int num_boundary_bins = boundary_bin_ids.size();
    for(int i =0; i < num_boundary_bins; i ++)
    {

            bin_id = boundary_bin_ids.at(i);
            particles_in_neighbor = bin_map.at(bin_id);
            particles_per_bin = particles_in_neighbor.size();
            n_bins[i].neighbor_bin_id = bin_id;
            n_bins[i].num_particles = particles_per_bin;
            n_bins[i].particle_offset = particle_neighbor_index;
            //neighbor_index++;

            for(int k = 0; k<particles_per_bin; k++ )
            {
                particles_to_send[particle_index] = particles[particles_in_neighbor.at(k)];
                particle_index++;
            }
    }

    //std::cout<<"form_particles_array_for_MPI END::: "<<std::endl;

}



std::vector<int> get_boundary_bins_for_curr_process(std::vector<int > process_bins, std::vector<std::vector<int> > neighbor_bins)
{
     std::vector<int>  neghbors_list;
     std::vector<int> border_neighbors;
    int num_of_bins_in_curr_proc = process_bins.size();
    for(int j =0; j < num_of_bins_in_curr_proc; j++)
    {

        //std::cout<<"i ::: "<< i << " ,j::"<<j<< " ,bin_ids.at(j):::"<< bin_ids.at(j)<< std::endl;
        //neghbors_list.push_back(neighbor_bins.at(bin_ids.at(j)));
        if(process_bins.at(j) != -1)
        {
            std::vector<int>  temp = neighbor_bins.at(process_bins.at(j));
            neghbors_list.insert(std::end(neghbors_list), std::begin(temp), std::end(temp));
        }
        
    }


    for(int k =0; k < neghbors_list.size(); k++)
    {
        if(std::find(process_bins.begin(), process_bins.end(), neghbors_list.at(k)) == process_bins.end())
        {
            //neighbor doesnt exist in the bin list. Add it
            border_neighbors.push_back(neghbors_list.at(k));

        }
    }

    return border_neighbors;

}


std::vector<std::vector<int> > get_boundary_bins(std::vector<std::vector<int> > process_bins, std::vector<std::vector<int> > neighbor_bins)
{
    //std::cout<<"get_boundary_bins START::: "<<std::endl;
    int num_of_processes = process_bins.size();
    std::vector<int>  bin_ids;
    std::vector<int>  neghbors_list;
    std::vector<std::vector<int> > border_neighbors(num_of_processes, std::vector<int>()); 
  
  //std::cout<<"neighbor_bins size ::: "<< neighbor_bins.size()<<std::endl;
  for(int i = 0; i < num_of_processes; i++)
  {
    bin_ids = process_bins.at(i);
    //neghbors_list has the list at all neighbors of bins in bin_ids,
    //Could have repeating values.
    for(int j =0; j < bin_ids.size(); j++)
    {

        //std::cout<<"i ::: "<< i << " ,j::"<<j<< " ,bin_ids.at(j):::"<< bin_ids.at(j)<< std::endl;
        //neghbors_list.push_back(neighbor_bins.at(bin_ids.at(j)));
        if(bin_ids.at(j) != -1)
        {
            std::vector<int>  temp = neighbor_bins.at(bin_ids.at(j));
            neghbors_list.insert(std::end(neghbors_list), std::begin(temp), std::end(temp));
        }
        
    }


    for(int k =0; k < neghbors_list.size(); k++)
    {
        if(std::find(bin_ids.begin(), bin_ids.end(), neghbors_list.at(k)) == bin_ids.end())
        {
            //neighbor doesnt exist in the bin list. Add it
            border_neighbors.at(i).push_back(neghbors_list.at(k));

        }
    }

    bin_ids.clear();
    neghbors_list.clear();

  }
 
   
    
//std::cout<<"get_boundary_bins END::: "<<std::endl;
return border_neighbors;

}



void get_neighbors_in_other_process(int bin_index, int current_rank, std::vector<std::vector<int> > neighbor_bins , std::vector<int> bin_process_map, std::vector<int> &outside_neighbors)
{
    std::vector<int> neighbors = neighbor_bins.at(bin_index);
    for(int i =0; i< neighbors.size(); i++)
    {
        if(bin_process_map[neighbors.at(i)] != current_rank)
        {
            if(std::find(outside_neighbors.begin(), outside_neighbors.end(), bin_process_map[neighbors.at(i)]) == outside_neighbors.end())
            {
                outside_neighbors.push_back(bin_process_map[neighbors.at(i)]);
            }

        }
    }

}




std::vector<std::vector<int> > initialize_neighbor_bins()
{
     //std::cout << ":::IN neighbor_bins::: " << std::endl;
     std::vector<std::vector<int> > neighbor_bins(total_bin_count, std::vector<int>(8, -1));

     //std::cout << ":::neighbor_bins.size()::: " << neighbor_bins.size()<<std::endl;
     int index;
     int last_row = num_of_bins_y*(num_of_bins_y-1);

    

     //#pragma omp parallel
    // {
         for(int i=0; i<total_bin_count;i++)
         { 
            index =0;


           if(i%num_of_bins_y != 0)
           {
            //not first column
              neighbor_bins[i][index] =  i-1;
              index++;

              if(i >= num_of_bins_y)
              {
                //not first row
                neighbor_bins[i][index] =  i-1-num_of_bins_y;
                index++;
              }

              if(i < last_row)
              {
                // Not last row
                neighbor_bins[i][index] =  i-1+num_of_bins_y;
                index++;
              }
              
             
           }


           if((i+1)%num_of_bins_y != 0)
           {
            //not last column
              neighbor_bins[i][index] =  i+1;
              index++;

              if(i >= num_of_bins_y)
              {
                //not first row
                neighbor_bins[i][index] =  i+1-num_of_bins_y;
                index++;
              }

              if(i < last_row)
              {
                neighbor_bins[i][index] =  i+1+num_of_bins_y;
                index++;
              } 
              
           }


           if(i >= num_of_bins_y)
              {
                neighbor_bins[i][index] =  i-num_of_bins_y;
                index++;
              }

              if(i < last_row)
              {
                neighbor_bins[i][index] =  i+num_of_bins_y;
                index++;
              } 
            
         }
     //} //OMP parallel ends

    

      /*
      std::vector<int>  my_neighbors;
        for(int i = 0; i < neighbor_bins.size(); i++)
        {
            std::cout << "Current bin::: " << i << ":: My neighbors:::" << std::endl;
            my_neighbors = neighbor_bins.at(i);
            if(!my_neighbors.empty())
            {
                for(int j = 0; j < my_neighbors.size(); j++)
                {
                    std::cout << "::: " << my_neighbors.at(j) << std::endl;
                }
            }
            

        }
    */

     return neighbor_bins;
}

/*
 * Set number of bins
 */
 void set_bin_count(int n)
 {

    //std::cout<< "n::" <<n<<std::endl;
    //std::cout<< "size::" <<size<<std::endl;

    bin_size = cutoff;
    //std::cout<< "bin_size::" <<bin_size<<std::endl;
    num_of_bins_x = ceil(size/bin_size);
    //std::cout<< "num_of_bins_x::" <<num_of_bins_x<<std::endl;
    num_of_bins_y = num_of_bins_x;
    //std::cout<< "num_of_bins_y::" <<num_of_bins_y<<std::endl;
    total_bin_count = num_of_bins_x*num_of_bins_y;
    //std::cout<< "total_bin_count::" <<total_bin_count<<std::endl;
    bin_density = n/total_bin_count;
    //std::cout<< "bin_density::" <<bin_density<<std::endl;
    
 }


 int compute_bin_index_from_xy(double x, double y)
 {
    
    int bin_index = (floor(x/bin_size)*num_of_bins_y)+floor(y /bin_size);
    return bin_index;

 }


 void remove_particle_from_bin(int old_bin_index, int new_bin_index, int particle_index, std::vector<std::vector<int> > &bin_map )
 {
     std::vector<int> particle_list = bin_map.at(old_bin_index);
     int number_of_particles = particle_list.size();
     for(int i = 0; i< number_of_particles; i++)
     {
        if(particle_list.at(i) == particle_index)
        {
            particle_list.erase(particle_list.begin() + i);
            break;
        }

     }
     //std::remove(particle_list.begin(), particle_list.end(), particle_index); 

     bin_map.erase(bin_map.begin() + old_bin_index);
     bin_map.insert(bin_map.begin()+old_bin_index, particle_list);
     bin_map.at(new_bin_index).push_back(particle_index);

 }


 void mpi_bin_particles(int n, particle_t *p ,  std::vector<std::vector<int> > &bin_map, std::vector<int> &bin_process_map)
 {
     int bin_index;

     for( int i = 0; i < n; i++ ) 
    {
         bin_index = compute_bin_index_from_xy( p[i].x, p[i].y );
         bin_map.at(bin_index).push_back(i);
         
    }
 }


 void bin_particles(int n, particle_t *p ,  std::vector<std::vector<int> > &bin_map)
 {
    

    int bin_index;
    std::vector<int> particle_list;
    //#pragma omp for 
    for( int i = 0; i < n; i++ ) 
    {
        //particle_list.clear();
         bin_index = compute_bin_index_from_xy( p[i].x, p[i].y );
         std::vector<int> &particle_list = bin_map.at(bin_index);
         particle_list.push_back(i);
         //bin_map.at(bin_index).push_back(i);
         
    }

    //std::vector<int> particle_list;
    /*for (int i =0; i < bin_map.size(); i++)
    {
        std::cout<< ":::::::BIN Index:::::::  "<< i << std::endl;
        particle_list = bin_map.at(i);
         for(int j =0; j< particle_list.size(); j++)
         {
            std::cout<< "::::particle Index::  "<< particle_list.at(j) << std::endl;
         }
    }*/

     
 }

//
//  Initialize the particle positions and velocities
//
void init_particles1( int n, particle_t *p ,  std::vector<std::vector<int> > &bin_map)
{
   // std::cout<< "INIT_PARTICLES1 START::"<< std::endl;
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    std::vector<int> particle_list;
    int bin_index = -1;
    for( int i = 0; i < n; i++ ) 
    {
        particle_list.clear();
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);
        //std::cout<< "p[i].x:: " << p[i].x << ",:: p[i].y:: " << p[i].y << std::endl;


        //Assign bin for the particle
         bin_index = compute_bin_index_from_xy( p[i].x, p[i].y );
         bin_map.at(bin_index).push_back(i);
         



        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
        p[i].index = i;
    }
    free( shuffle );
    //std::cout<< "INIT_PARTICLES1 END::"<<  std::endl;
}





void init_particles( int n, particle_t *p)
{
    //std::cout<< "init_particles1 START::"<< std::endl;
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);
        //std::cout<< "p[i].x:: " << p[i].x << ",:: p[i].y::" << p[i].y << std::endl;

        //TODO: Assign bin for the particle
        //bin_map.push_back() vvvv
        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;

        //p[i].index = i;
    }
    free( shuffle );
    //std::cout<< "init_particles1 END::"<< std::endl;
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    if (r2 != 0)
        {
       if (r2/(cutoff*cutoff) < *dmin * (*dmin))
          *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
        
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
    
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //

    //std::cout<< "Move::Before ::::: p.x :: " << p.x << ",:: p.y ::  " << p.y << std::endl;
    int old_bin_index = compute_bin_index_from_xy( p.x, p.y);
    

    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }

    //std::cout<< "Move::After:::::: p.x :: " << p.x << ",:: p.y :: " << p.y << std::endl;
    int new_bin_index = compute_bin_index_from_xy( p.x, p.y);
    //std::cout<< "old_bin_index::  " <<old_bin_index << " ,new_bin_index::  " <<new_bin_index<<std::endl;

}  




void move1( particle_t &p, std::vector<std::vector<int> > &bin_map)
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //

    //Remove the current particle from its old location
    int old_bin_index = compute_bin_index_from_xy( p.x, p.y);
    //std::cout<< "old_bin_index::" <<old_bin_index<<std::endl;
    

    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }

     
     //Insert the particle into its new bin based on its new position cordinates
     int bin_index = compute_bin_index_from_xy( p.x, p.y);
     //std::cout<< "bin_index::" <<bin_index<<std::endl;
     if(old_bin_index != bin_index)
     {
        //std::vector<int> particle_list = bin_map.at(old_bin_index);
        //std::vector<int> particle_list_new  = bin_map.at(bin_index);
        remove_particle_from_bin(old_bin_index, bin_index, p.index, bin_map );
        bin_map.at(bin_index).push_back(p.index);
        
       
     }
      

}


//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
