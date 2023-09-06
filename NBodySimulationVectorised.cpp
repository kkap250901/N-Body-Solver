#include "NBodySimulation.h"


class NBodySimulationVectorised : public NBodySimulation {
	public :
	void updateBody (){

	// Initialisations 
	timeStepCounter++;
	maxV   = 0.0;
	minDx  = std::numeric_limits<double>::max();

	// Making forces equal to 0 and vectorising this
	#pragma omp simd
	for (int i = 0; i < NumberOfBodies; i++){
		force0[i] = 0.0;
		force1[i] = 0.0;
		force2[i] = 0.0;
	}


	// Again the force calculations here 
	double epislon = 0.001;
	for (int i=0; i<NumberOfBodies; i++) {

		// Vectorising this
		#pragma omp simd
		for (int j = i + 1; j<NumberOfBodies; j++) {
			// Distance in the x axiss
			double dx = x[j]-x[i];
			// Distance in the y axi
			double dy = y[j]-y[i];
			// Distance in the z axias
			double dz = z[j]-z[i];
			// Getting the eucledian distance
			double eucledian_distance = sqrt(dx*dx + dy*dy + dz*dz);
			// Min mass of 2 particles for calculating the min distance
			// Cubing the distance
			const double distance3 = eucledian_distance * eucledian_distance * eucledian_distance;
			// This is to get the magnitude of the forces
			const double force_mag = mass[j] * mass[i] / distance3;
			// Gettting the minimum eucledian_distance
			minDx = std::min(minDx,eucledian_distance);
			// x,y,z forces acting on particle 0.
			double forcex = dx * force_mag;
			double forcey = dy * force_mag;
			double forcez = dz * force_mag;

			force0[i] += forcex;
			force1[i] += forcey;
			force2[i] += forcez;
			force0[j] -= forcex;
			force1[j] -= forcey;
			force2[j] -= forcez;
		}
	}

		//todo edit this with also including the Collision condition
		double max_velocity = 0.0;
		#pragma omp simd reduction(max:max_velocity)
		for (int i = 0; i<NumberOfBodies; i++){
			max_velocity = std::max(max_velocity, vx[i]);
			max_velocity = std::max(max_velocity, vy[i]);
			max_velocity = std::max(max_velocity, vz[i]);
		}

		// Please enable this for the stable timestepping technique
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
			maxV = std::max(maxV, std::sqrt(pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2)));
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