#include <iomanip>
#include "NBodySimulationVectorised.cpp"
#include <chrono>
#include <omp.h>

/**
 * You can compile this file with
 *   make step-4-gcc   // Uses the GNU Compiler Collection.
 *   make asisigment-icpc  // Uses the Intel compiler.
 * and run it with
 *   ./step-4-gcc
 *   ./step-4-icpc
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */



class NBodySimulationParallelised : public NBodySimulationVectorised {
public :

#pragma omp declare simd
void force_calc(){

	// This loop is parallelised and this loop has been done in an non asymetric manner
	#pragma omp parallel for
	for (int i=0; i<NumberOfBodies; i++) {
		// std::cout << omp_get_num_threads() << std::endl;
		double forcex_temp = 0.0;
		double forcey_temp = 0.0;
		double forcez_temp = 0.0;

		#pragma omp simd
		#pragma simd reduction (+:forcex_temp,forcey_temp,forcez_temp)
		// #pragma ivdep
		// #pragma omp simd reduction (+:forcex_temp,forcey_temp,forcez_temp)
		for (int j = 0; j<i; j++) {
			// Distance in the x axiss
			const double dx = x[j]-x[i];
			// Distance in the y axi
			const double dy = y[j]-y[i];
			// Distance in the z axias
			const double dz = z[j]-z[i]; 
			// Getting the eucledian distance
			const double eucledian_distance = sqrt(dx*dx + dy*dy + dz*dz);
			// Min mass of 2 particles for calculating the min distance
			// max_mass_together = std::max(max_mass_together, mass[i] + mass[j]);
			// Cubing the distance
			const double distance3 = eucledian_distance * eucledian_distance * eucledian_distance;
			// This is to get the magnitude of the forces
			const double force_mag = mass[j] * mass[i] / distance3;
			// Gettting the minimum eucledian_distance
  			// minDx2 = std::min(minDx2,eucledian_distance);
			// x,y,z forces acting on particle 0.
			const double forcex = dx * force_mag;
			const double forcey = dy * force_mag;
			const double forcez = dz * force_mag;

			forcex_temp += forcex;
			forcey_temp += forcey;
			forcez_temp += forcez;
		}

		#pragma omp simd
		#pragma simd reduction (+:forcex_temp,forcey_temp,forcez_temp)
		for (int j = i+1; j<NumberOfBodies; j++) {
			// Distance in the x axiss
			const double dx = x[j]-x[i];
			// Distance in the y axi
			const double dy = y[j]-y[i];
			// Distance in the z axias
			const double dz = z[j]-z[i];
			// Getting the eucledian distance
			const double eucledian_distance = sqrt(dx*dx + dy*dy + dz*dz);
			// Min mass of 2 particles for calculating the min distance
			// max_mass_together = std::max(max_mass_together, mass[i] + mass[j]);
			// Cubing the distance
			const double distance3 = eucledian_distance * eucledian_distance * eucledian_distance;
			// This is to get the magnitude of the forces
			const double force_mag = mass[j] * mass[i] / distance3;
			// Gettting the minimum eucledian_distance

			// minDx2 = std::min(minDx2,eucledian_distance);
			// x,y,z forces acting on particle 0.
			const double forcex = dx * force_mag;
			const double forcey = dy * force_mag;
			const double forcez = dz * force_mag;

			forcex_temp += forcex;
			forcey_temp += forcey;
			forcez_temp += forcez;
			}
		
		force0[i] = forcex_temp;
		force1[i] = forcey_temp;
		force2[i] = forcez_temp;
	}
}


void updateBody (){

	// Initialisations 
	timeStepCounter++;
	maxV   = 0.0;
	minDx  = std::numeric_limits<double>::max();

	force_calc();

	const double epislon = 0.001;

	double max_velocity = 0.0;

	#pragma omp simd reduction(max:max_velocity)
	for (int i = 0; i<NumberOfBodies; i++){
		max_velocity = std::max(max_velocity, vx[i]);
		max_velocity = std::max(max_velocity, vy[i]);
		max_velocity = std::max(max_velocity, vz[i]);
	}

	double stable_timestep = (minDx/ 2.0) / (max_velocity + epislon);
	// timeStepSize = stable_timestep;

	#pragma omp simd
	for (int i=0; i<NumberOfBodies; i++){
		x[i] = x[i] + timeStepSize * vx[i];
		y[i] = y[i] + timeStepSize * vy[i];
		z[i] = z[i] + timeStepSize * vz[i];
	}


	#pragma omp simd
	for (int i = 0 ; i < NumberOfBodies; i++){
		vx[i] = vx[i] + timeStepSize * force0[i] / mass[i];
		vy[i] = vy[i] + timeStepSize * force1[i] / mass[i];
		vz[i] = vz[i] + timeStepSize * force2[i] / mass[i];
	}


	#pragma omp simd reduction(max:maxV)
	for (int i = 0; i < NumberOfBodies; i++){
		maxV = std::max(maxV, std::sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]));
	}
	

	// Check for colliding conditions 
	// Performing merging of these masses
	for (int i=0; i<NumberOfBodies; i++){
		for (int j = i + 1; j<NumberOfBodies; j++) {

			double dx = x[j]-x[i];

		// Distance in the y axi
			double dy = y[j]-y[i];

			// Distance in the z axias
			double dz = z[j]-z[i];

			// Getting the eucledian distance
			double eucledian_distance = sqrt(dx*dx + dy*dy + dz*dz);
			minDx = std::min(minDx,eucledian_distance);
			const double mass_added = mass[i] + mass[j];

			if (eucledian_distance <= C * mass_added){
				x[i] = (mass[i] * x[i] + mass[j] * x[j]) / (mass_added);
				y[i] = (mass[i] * y[i] + mass[j] * y[j]) / (mass_added);
				z[i] = (mass[i] * z[i] + mass[j] * z[j]) / (mass_added);
				vx[i] = (mass[i] * vx[i] + mass[j] * vx[j]) / (mass_added);
				vy[i] = (mass[i] * vy[i] + mass[j] * vy[j]) / (mass_added);
				vz[i] = (mass[i] * vz[i] + mass[j] * vz[j]) / (mass_added);
				mass[i] = mass_added;
				NumberOfBodies -= 1;
				x[j] = x[NumberOfBodies];
				y[j] = y[NumberOfBodies];
				z[j] = z[NumberOfBodies];
				vx[j] = vx[NumberOfBodies];
				vy[j] = vy[NumberOfBodies];
				vz[j] = vz[NumberOfBodies];
				mass[j] = mass[NumberOfBodies];
				j -=1;
			} 
		}
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

  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulationParallelised nbs;
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
