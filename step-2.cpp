#include <iomanip>
#include <list>
#include <array>
#include <vector>
#include "NBodySimulation.h"
#include <string>

/**
 * You can compile this file with
 *   make step-2-gcc   // Uses the GNU Compiler Collection.
 *   make asisigment-icpc  // Uses the Intel compiler.
 * and run it with
 *   ./step-2-gcc
 *   ./step-2-icpc
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */


/*
Model cells as 3d arrays
cells[x][y][z] = list(p1,p2,p3)--> O(1) space
1. Spawn particles
2. Store them in correct cells based on their positions
3. Calculate forces for each particle
  3.1 Iterate over each cell
  3.2 Get neighboouring cells and their particles
  3.3 Get pairwise forces withing the combined lists
*/

class NBodySimulationMolecularForces : public NBodySimulation {
  public:
  double min_x;
  double max_x;
  double min_y;
  double max_y;
  double min_z;
  double max_z;
  double box_size;
  int NumberOfCellsInX;
  int NumberOfCellsInY;
  int NumberOfCellsInZ;
  const double radius_cutoff = 0.22;
  const double epislon = 0.05;
  const double epislon_dis = 0.000005;

  //Resizeable array
  // std::array<std::array<std::array<std::list<int>>>> cells;
  std::vector<int> ***cells;

  void setUp(int argc, char** argv){
    // Default setup (taken from NBodySimulation class)
    checkInput(argc, argv);
    NumberOfBodies = (argc-4) / 7;

    x = new double[ NumberOfBodies];
    y = new double [NumberOfBodies];
    z = new double [NumberOfBodies];
    vx = new double [NumberOfBodies];
    vy = new double [NumberOfBodies];
    vz = new double [NumberOfBodies];
    mass = new double [NumberOfBodies];
    force0 = new double[NumberOfBodies];
    force1 = new double[NumberOfBodies];
    force2 = new double[NumberOfBodies];

    int readArgument = 1;

    tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
    tFinal       = std::stof(argv[readArgument]); readArgument++;
    timeStepSize = std::stof(argv[readArgument]); readArgument++;

    for (int i=0; i<NumberOfBodies; i++) {
    x[i] = std::stof(argv[readArgument]); readArgument++;
    y[i] = std::stof(argv[readArgument]); readArgument++;
    z[i] = std::stof(argv[readArgument]); readArgument++;

    vx[i] = std::stof(argv[readArgument]); readArgument++;
    vy[i] = std::stof(argv[readArgument]); readArgument++;
    vz[i] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
        std::cerr << "invalid mass for body " << i << std::endl;
        exit(-2);
      }
    } 

    // Create cells
    // Define the variables for max_x and all other minimums to define the limits of the particles
    max_x = std::numeric_limits<double>::lowest();
    min_x = std::numeric_limits<double>::max();
    max_y = std::numeric_limits<double>::lowest();
    min_y = std::numeric_limits<double>::max();
    max_z = std::numeric_limits<double>::lowest();
    min_z = std::numeric_limits<double>::max();


    // this is to find the max and mins from the particles spawned
    for (int i = 0; i < NumberOfBodies; i++ ){
      // std::cout<< "Minimum of distance "<< std::endl;
      max_x = std::max(max_x, x[i]) + epislon;
      min_x = std::min(min_x, x[i]) - epislon;
      max_y = std::max(max_y, y[i]) + epislon;
      min_y = std::min(min_y, y[i]) - epislon;
      max_z = std::max(max_z, z[i]) + epislon; 
      min_z = std::min(min_z, z[i]) - epislon;
    }


    //Defining the size of each cell side
    box_size = 2 * radius_cutoff;

    // Finding the number of cells in the x,y and z
    NumberOfCellsInX = ceil((max_x - min_x) / box_size);
    NumberOfCellsInY = ceil((max_y - min_y) / box_size);
    NumberOfCellsInZ = ceil((max_z - min_z) / box_size);
    cells = new std::vector<int> **[NumberOfCellsInX];

    // Looping over it and defining the datastructure for the cells 
    for (int i = 0; i < NumberOfCellsInX; i++){

      cells[i] = new std::vector<int> *[NumberOfCellsInY];

      for (int j = 0; j < NumberOfCellsInY; j++){

        cells[i][j] = new std::vector<int> [NumberOfCellsInZ];

        for (int k = 0; k < NumberOfCellsInZ; k++){

          cells[i][j][k] = std::vector<int>();
        }
      }
    }
  }

  // This funciton is to populate the cells with the indices of the right particles
  void populate_cells(){
    for (int i = 0; i < NumberOfCellsInX; i++){

      for (int j = 0; j < NumberOfCellsInY; j++){

        for (int k = 0; k < NumberOfCellsInZ; k++){
          
// This here just clears the existing particle indices
          cells[i][j][k].clear();
        }
      }
    }

    // This here is to actually push back or add to the indices of the particles
    for (int i = 0 ; i < NumberOfBodies; i++){

      int x_cell_index = floor((x[i] - min_x) / box_size);

      int y_cell_index = floor((y[i] - min_y) / box_size);

      int z_cell_index = floor((z[i] - min_z) / box_size);

      cells[x_cell_index][y_cell_index][z_cell_index].push_back(i);

    }
  }



  // This is a funciton to collate all particles in the neighbouring cells
  std::vector<int> neighbouring_cells(int x, int y , int z){

  std::vector<int> neighbours;

  // this returns all the particles in the neighbouring cells
  for(int i = std::max(0, x-1); i <= std::min(NumberOfCellsInX -1 ,x + 1);i++){


    for(int j = std::max(0, y-1); j <= std::min(NumberOfCellsInY -1 ,y + 1);j++){


      for(int k = std::max(0, z-1); k <= std::min(NumberOfCellsInZ -1 ,z + 1);k++){

        
        neighbours.insert(neighbours.end(), cells[i][j][k].begin(), cells[i][j][k].end());
        }
      }
    }
    return neighbours;
  }

  // Updating the body function
  void updateBody(){

    // std::cout << 'minDx: ' << std::endl; 
    // std::cout << minDx << std::endl; 

    timeStepCounter++;
    maxV = 0.0;
    minDx  = std::numeric_limits<double>::max();

    // Populating the cells with the right indices
    populate_cells();
    //Calculating forces

    // Defining the datastructure for neighbouring particles
    std::vector<int> neighbours;


    // Then iterating through all cells
    for (int cell_x = 0; cell_x < NumberOfCellsInX; cell_x++){


      for (int cell_y = 0; cell_y < NumberOfCellsInY; cell_y++){


        for (int cell_z = 0; cell_z < NumberOfCellsInZ; cell_z++){

          // Now calculate forces for each particle inside the neighbouring cells
          neighbours = neighbouring_cells(cell_x, cell_y, cell_z );

          for (int neighbour_index : cells[cell_x][cell_y][cell_z]){

              double forcex_temp = 0.0;

              double forcey_temp = 0.0;

              double forcez_temp = 0.0;

            for (int neighbour_index_j = 0; neighbour_index_j < neighbours.size(); neighbour_index_j++){

              if (neighbour_index == neighbour_index_j){

                double dx = x[neighbour_index_j] - x[neighbour_index];

                // Thiss here can happen at the start where particles spawned next to each other
                // This would prevent particles having the same positon
                if (dx == 0){

                  dx = 2 * epislon_dis;

                  x[neighbour_index_j] += epislon_dis;
                  
                  x[neighbour_index] -= epislon_dis;
                }

                // Distance in the y axi
                double dy = y[neighbour_index_j] - y[neighbour_index];

                if (dy == 0){

                  dy = 2 * epislon_dis;

                  y[neighbour_index_j] += epislon_dis;

                  y[neighbour_index] -= epislon_dis;
                }

                // Distance in the z axias
                double dz = z[neighbour_index_j] - z[neighbour_index];

                if (dz == 0){

                  dz = 2 * epislon_dis;

                  z[neighbour_index_j] += epislon_dis;

                  z[neighbour_index] -= epislon_dis;
                }

                // Getting the eucledian distance
                double eucledian_distance = sqrt(dx*dx + dy*dy + dz*dz);
                minDx = std::min(minDx,eucledian_distance);


                // Calculating all the forces now
                double forcex = 0.0;

                double forcey = 0.0;

                double forcez = 0.0;
                

                // if only the distance is less than the cutoff radius then we calclulate the forces
                if (eucledian_distance <= radius_cutoff){
                  forcex = - 10 * (pow((0.1/eucledian_distance),13) - pow((0.1/eucledian_distance),9)) * dx;

                  forcey = - 10 * (pow((0.1/eucledian_distance),13) - pow((0.1/eucledian_distance),9)) * dy;

                  forcez = - 10 * (pow((0.1/eucledian_distance),13) - pow((0.1/eucledian_distance),9)) * dz;

                }
        
                forcex_temp += forcex;

                forcey_temp += forcey;

                forcez_temp += forcez;
              
              }
            }
        force0[neighbour_index] = forcex_temp;

        force1[neighbour_index] = forcey_temp;

        force2[neighbour_index] = forcez_temp;
        }
      }
    }
  }

    //Updating the velocities and positions of the particles
    for (int particle_no = 0; particle_no < NumberOfBodies; particle_no ++){
      maxV = 0.0;
      // This is  to update the velocities and the positions of the particles
      vx[particle_no] = vx[particle_no] + timeStepSize * force0[particle_no] / mass[particle_no];

      vy[particle_no] = vy[particle_no] + timeStepSize * force1[particle_no] / mass[particle_no];

      vz[particle_no] = vz[particle_no] + timeStepSize * force2[particle_no] / mass[particle_no];

      x[particle_no] = x[particle_no] + timeStepSize * vx[particle_no];

      // Finding if the particles are outside the box or not 
      if (x[particle_no] >= max_x){
        
        // spawn them a bit behind the limit
        x[particle_no] = max_x - epislon;

        // Then making the velocity
        vx[particle_no] = vx[particle_no] *  -1;

      }


      if (x[particle_no] <= min_x){

        x[particle_no] = min_x + epislon;

        vx[particle_no] = vx[particle_no] * -1;
      }

      y[particle_no] = y[particle_no] + timeStepSize * vy[particle_no];

      if (y[particle_no] >= max_y){

        y[particle_no] = max_y - epislon;

        vy[particle_no] = vy[particle_no] * -1;

      }

      if (y[particle_no] <= min_y){

        y[particle_no] = min_y  + epislon;

        vy[particle_no] = vy[particle_no]  * -1;

      }


      z[particle_no] = z[particle_no] + timeStepSize * vz[particle_no];

      if (z[particle_no] >= max_z){

        z[particle_no] = max_z - epislon;

        vz[particle_no] = vz[particle_no] * -1;

      }

      if (z[particle_no] <= min_z){

        z[particle_no] = min_z + epislon;

        vz[particle_no] = vz[particle_no] * -1;

      }

      maxV = std::max(maxV, std::sqrt(vx[particle_no] * vx[particle_no] + vy[particle_no] * vy[particle_no] + vz[particle_no] * vz[particle_no]));

    }
    t += timeStepSize;
  }

};
/**
 * Main routine.
 *
 * No major changes are needed in the assignment. You can add initialisation or
 * or remove input checking, if you feel the need to do so. But keep in mind
 * that you may not alter what the program writes to the standard output.
 */
int main (int argc, char** argv) {

  // std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulationMolecularForces nbs;
  nbs.setUp(argc,argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}
