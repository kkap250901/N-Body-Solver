#include "NBodySimulation.h"
#include <iomanip>
#include <fstream>
#include <iostream>

NBodySimulation::NBodySimulation () :
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), y(nullptr), z(nullptr),
  vx(nullptr),vy(nullptr),vz(nullptr), 
  force0(nullptr), force1(nullptr), force2(nullptr),
  mass(nullptr),
  timeStepSize(0), maxV(0), minDx(0), 
  videoFile(nullptr),
  snapshotCounter(0), timeStepCounter(0) {};

NBodySimulation::~NBodySimulation () {
  if (x != nullptr) {
    delete [] x;
  }
  if (y != nullptr) {
    delete [] y;
  }
  if (z != nullptr) {
    delete [] z;
  }
  if (vx != nullptr) {
    delete [] vx;
  }
  if (vy != nullptr) {
    delete [] vy;
  }
  if (vz != nullptr) {
    delete [] vz;
  }
  // Mass memory deleting
  if (mass != nullptr) {
    delete [] mass;
  }
  if (force0 != nullptr){
    delete [] force0;
  }
    if (force1 != nullptr){
    delete [] force1;
  }
    if (force2 != nullptr){
    delete [] force2;
  }
}

void NBodySimulation::checkInput(int argc, char** argv) {
    if (argc==1) {
    std::cerr << "usage: " << std::string(argv[0])
              << " plot-time final-time dt objects" << std::endl
              << " Details:" << std::endl
              << " ----------------------------------" << std::endl
              << "  plot-time:       interval after how many time units to plot."
                 " Use 0 to switch off plotting" << std::endl
              << "  final-time:      simulated time (greater 0)" << std::endl
              << "  dt:              time step size (greater 0)" << std::endl
              << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
              << std::endl
              << "Examples of arguments:" << std::endl
              << "+ One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ One body spiralling around the other" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ Three-body setup from first lecture" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
              << "+ Five-body setup" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
              << std::endl;

    throw -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each body is given by seven entries"
                 " (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments"
                 " (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    throw -2;
  }
}

void NBodySimulation::setUp (int argc, char** argv) {

  checkInput(argc, argv);

  NumberOfBodies = (argc-4) / 7;

  // Changed initialisation for better memory management
  x    = new double [NumberOfBodies];
  y    = new double [NumberOfBodies];
  z    = new double [NumberOfBodies];
  vx    = new double [NumberOfBodies];
  vy    = new double [NumberOfBodies];
  vz    = new double [NumberOfBodies];
  mass = new double [NumberOfBodies];
  force0 = new double [NumberOfBodies];
  force1 = new double [NumberOfBodies];
  force2 = new double [NumberOfBodies];

  C = 10e-2 / NumberOfBodies;

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

  std::cout << "created setup with " << NumberOfBodies << " bodies"
            << std::endl;

  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta
              << " time units" << std::endl;
    tPlot = 0.0;
  }
}


void NBodySimulation::updateBody () {

  // counteer
  timeStepCounter++;
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  // Setting the forces to 0 initially
  for (int i = 0; i < NumberOfBodies; i++){
    force0[i] = 0.0;
    force1[i] = 0.0;
    force2[i] = 0.0;
  }

// Calculating the forces first for each object 
  const double epislon = 0.001;
  for (int i=0; i<NumberOfBodies; i++) {
    // double  minDx2  = std::numeric_limits<double>::max();
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

      // Adding the forces in the arrays and subtracting to decrease the number of calc.
      force0[i] += forcex;
      force1[i] += forcey;
      force2[i] += forcez;
      force0[j] -= forcex;
      force1[j] -= forcey;
      force2[j] -= forcez;
      }
    }

    // This is for calculation of the stable timstep
    double stable_timestep = 0.0;

    // We find the max velocity here
    double max_velocity = 0.0;
    for (int i = 0; i<NumberOfBodies; i++){
      max_velocity = std::max(max_velocity,vx[i]);
      max_velocity = std::max(max_velocity, vy[i]);
      max_velocity = std::max(max_velocity, vz[i]);
    }

    // Stable timstep
    stable_timestep = (minDx/ 2.0) / (max_velocity + epislon);


    // Enable this to enable my stable timestepping scheme please
    // timeStepSize = stable_timestep; 


    // Updating the positions and velocities 
    for (int i=0; i<NumberOfBodies; i++){
      x[i] = x[i] + timeStepSize * vx[i];
      y[i] = y[i] + timeStepSize * vy[i];
      z[i] = z[i] + timeStepSize * vz[i];

      vx[i] = vx[i] + timeStepSize * force0[i] / mass[i];
      vy[i] = vy[i] + timeStepSize * force1[i] / mass[i];
      vz[i] = vz[i] + timeStepSize * force2[i] / mass[i];
      maxV = std::max(maxV, std::sqrt(vx[i] *vx[i]  + vy[i] * vy[i] + vz[i] * vz[i]));
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

/**
 * Check if simulation has been completed.
 */
bool NBodySimulation::hasReachedEnd () {
  return t > tFinal;
}

void NBodySimulation::takeSnapshot () {
  if (t >= tPlot) {
    printParaviewSnapshot();
    printSnapshotSummary();
    tPlot += tPlotDelta;
  }
}


void NBodySimulation::openParaviewVideoFile () {
  videoFile.open("paraview-output/result.pvd");
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
    " version=\"0.1\""
    " byte_order=\"LittleEndian\""
    " compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

void NBodySimulation::closeParaviewVideoFile () {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
  videoFile.close();
}

void NBodySimulation::printParaviewSnapshot () {
  static int counter = -1;
  counter++;
  std::stringstream filename, filename_nofolder;
  filename << "paraview-output/result-" << counter <<  ".vtp";
  filename_nofolder << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\""
    " NumberOfComponents=\"3\""
    " format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i]
        << " "
        << y[i]
        << " "
        << z[i]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  out.close();

  videoFile << "<DataSet timestep=\"" << counter
            << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
            << "\"/>" << std::endl;
}

void NBodySimulation::printSnapshotSummary () {
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t="         << t
            << ",\t dt="        << timeStepSize
            << ",\t v_max="     << maxV
            << ",\t dx_min="    << minDx
            << std::endl;
}

void NBodySimulation::printSummary () {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0] << ", " << y[0] << ", " << z[0] << std::endl;
}

// This is to write the summary 
void NBodySimulation::writeSummary (std::ofstream &myfile) {
  if (NumberOfBodies == 2){
    double vel_1 = sqrt(vx[0] * vx[0] + vy[0] * vy[0] + vz[0] * vz[0]);
    double vel_2 = sqrt(vx[1] * vx[1] + vy[1] * vy[1] + vz[1] * vz[1]);
    myfile << "Velocity of particle 1 : " << vel_1 << " Timestep :" << t << std::endl << "Velocity of particle 2 : " << vel_2 << std::endl;
  }
}
