/*	Author: Dominikus Brian
	"Lennard-Jones Particle Cluster MD simulation"
	2020/July 2/2020
	Version 15.0.2020
*/
/*----------------------Program begin here---------------------- */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include <fstream>
#include <time.h>
using namespace std;

////////////////////////// Global parameters  ////////////////////////////////

const int N = 32;							// N is the number of particle within the system
const int DIM = 3;							// DIM is dimension of the problem
const int tmax = 15000;						// maximum number of timestep or simulation duration
const double rho = 0.8;						// System's volume density
double dt = 0.001;							// Simulation timestep size
double Init_Temp = 1;							// Temperature
const double L = N / rho;					// Box edge size
const double delta_r = 0.2;					// bin width for rdf analysis
const int size_pair = (N * (N - 1)) / 2;	// array size for the total number of pairs 
double Temp = 0;


////////////////////////// Global constants ////////////////////////////////

/*This simulation produce outcome in a MD reduced unit which are defined as follows:
Reduced time, tau = sqrt(epsilon /mass * sigma2) * t
Reduced temperature, T* = kbT/epsilon ; kb is boltzmann constant
Reduced density, rho* = rho * sigma^3
Reduced force , f* = f * sigma /epsilon
Reduced pressure P* = P * sigma^3 /epsilon
*/


const double epsilon = 1;				// Energy scale constant
const double sigma = 1;					// Length scale constant
const double mass = 1;					// Mass of the particle
const double PI = 3.1415926535;			// PI value
const double kB = 1;				// 1.3807 x 10^-23 J /K
ofstream myfile;

/////////////////////////////*Function declaration:*/////////////////////////////////

// Initialize Particle Position
void Config_atoms_lattice(double r[N][DIM]);
void Config_atoms_rand(double r[N][DIM]);

// Initialize Particle Velocity
void Config_velocity(double v[N][DIM]);

// Force Calculation
void Force_Calc(double r[N][DIM], double F[N][DIM], double& E_pot);

// Velocity-Verlet propagator
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& E_pot, double& E_kin);
void Vel_Verlet_first(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& E_pot, double& E_kin);
//Calculate Kinetic Energy
void E_kin_Calc(double v[N][DIM], double& E_kin);

//Set temperature according to velocity
void Config_Temp(double v[N][DIM], double& Temp);

//Periodic Boundary Condition
double PBC(double x1, double x2);
double forcePBC(double x);

//radial distribution function  
void rdf(double r[N][DIM], double& rij, double& count);

/*------------------------Main function begin here------------------------ */

int main() {
	clock_t tstart = clock();
	int i, t;
	double r[N][DIM];				// Particle Position vector in x y z
	double v[N][DIM];				// Particle velocity vector in x y z
	double F[N][DIM];				// Force exerted on atom i, in x y z
	double E_kin = 0;				// Kinetic Energy	
	double E_pot = 0;				// Potential Energy
	double E_tot = 0;				// Total Energy	
	double count = 0;				// Counter for histogram in rdf function
	double rij = 0;					// pair distance
	double Temp = 0;

	cout << "------------------BEGIN SIMULATION----------------------" << endl;

	// Displaying Input parameters and Global constants
	cout << "INPUT PARAMETERS: " << endl;
	cout << "N = " << N << endl;
	cout << "dt = " << dt << endl;
	cout << "tmax = " << tmax << endl;
	cout << "rho = " << rho << endl;
	cout << "dr = " << delta_r << endl;

	cout << "GLOBAL CONSTANTS: " << endl;
	cout << "Epsilon = " << epsilon << endl;
	cout << "Sigma = " << sigma << endl;
	cout << "Mass = " << mass << endl;

	// Export a mdinfo file about the input parameter and global constant employed in the simulation. 
	myfile.open("mdinfo.txt");
	myfile << "INPUT PARAMETERS: " << endl;
	myfile << "N = " << N << endl;
	myfile << "dt = " << dt << endl;
	myfile << "tmax = " << tmax << endl;
	myfile << "rho = " << rho << endl;
	myfile << "dr = " << delta_r << endl;
	myfile << "Temperature = " << Temp << endl;
	myfile << endl;
	myfile << "GLOBAL CONSTANTS: " << endl;
	myfile << "Epsilon = " << epsilon << endl;
	myfile << "Sigma = " << sigma << endl;
	myfile << "Mass = " << mass << endl;
	myfile.close();
	cout << "------------------BEGIN INITIALIZATION----------------------" << endl;

	//Initialize Position
	Config_atoms_rand(r);
	//Initialize Velocity
	Config_velocity(v);
	//Initialize Temperature
	Config_Temp(v, Temp);
	//Initialize Force
	Force_Calc(r, F, E_pot);

	//Export Initial Configurations
	myfile.open("mdinfo.txt");
	myfile << "INPUT PARAMETERS: " << endl;
	myfile << "N = " << N << endl;
	myfile << "dt = " << dt << endl;
	myfile << "tmax = " << tmax << endl;
	myfile << "rho = " << rho << endl;
	myfile << "dr = " << delta_r << endl;
	myfile << "Temperature = " << Temp << endl;
	myfile << endl;
	myfile << "GLOBAL CONSTANTS: " << endl;
	myfile << "Epsilon = " << epsilon << endl;
	myfile << "Sigma = " << sigma << endl;
	myfile << "Mass = " << mass << endl;
	myfile.close();

	//Export initial configuration 
	myfile.open("config.csv");
	myfile << "N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;
	for (i = 0; i < N; i++) {
		myfile << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}
	myfile.close();

	//Set momentum to zero
	for (i = 0; i < N; i++) {
		v[i][0] = 0;
		v[i][1] = 0;
		v[i][2] = 0;
	}
	cout << "------------------INITIALIZATION COMPLETE ----------------------" << endl;
	cout << "------------------BEGIN PROPAGATION ----------------------" << endl;
	//Write Energy value at each time step into properties.csv
	myfile.open("prop.csv");
	myfile << "t" << "," << "Potential" << "," << "Kinetic" << "," << "Total" << "," << "Temperature" << endl;
	

	for (t = 0; t < tmax; t++) {
		if (t == 0) {
			Vel_Verlet_first(r, v, F, E_pot, E_kin);
			E_tot = E_kin + E_pot;
			Temp = ((E_kin / N) * 2 / (3 * kB));
		}
		Vel_Verlet(r, v, F, E_pot, E_kin);
		E_tot = E_kin + E_pot;
		Temp =  ((E_kin / N) * 2 / (3 * kB));

		myfile << t + 1 << "," << E_pot << "," << E_kin << "," << E_tot << "," << Temp << endl; 
		
		if (t == tmax * 0)cout << "The progress is : " << endl;
		if (t == tmax * 0.10)cout << " | 10%";
		if (t == tmax * 0.25)cout << " | 25%";
		if (t == tmax * 0.50)cout << " | 50%";
		if (t == tmax * 0.75)cout << " | 75%";
		if (t == tmax * 0.90)cout << " | 90%";
		if (t == tmax - 1)cout << " | 100% |" << endl;


	}
	myfile.close();
	cout << "------------------PROPAGATION FINISHED ----------------------" << endl;
	cout << "------------------WRITING TRAJECTORY FILES----------------------" << endl;
	myfile.open("rdf.csv");
	for (t = 0; t < tmax; t++) {
		
		if (t % 50 == 0) {
			myfile << "t" << "," << t + 1 << endl;
		
			rdf(r, rij, count);
			
			
		}
		if (t == tmax * 0)cout << "The progress is : " << endl;
		if (t == tmax * 0.10)cout << " | 10%";
		if (t == tmax * 0.25)cout << " | 25%";
		if (t == tmax * 0.50)cout << " | 50%";
		if (t == tmax * 0.75)cout << " | 75%";
		if (t == tmax * 0.90)cout << " | 90%";
		if (t == tmax - 1)cout << " | 100% |" << endl;
	}
	myfile.close();

	//Write Final Configuration to final_config.csv
	myfile.open("finalconfig.csv");
	myfile << "N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;
	for (i = 0; i < N; i++) {
		myfile << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}
	myfile.close();


	cout << "------------------END OF SIMULATION ----------------------" << endl;
	printf("Total simulation time %.2fs\n", (double)(clock() - tstart) / CLOCKS_PER_SEC);


	return 0;
}
/*------------------------Main function end here------------------------ */

/*------------------------Function definition------------------------- */

// Initialize position
void Config_atoms_lattice(double r[N][DIM]) {
	int edge = 2;    // For rho = 0.5 >> edge = 2 gives fill up half of the box, 4 fill up a quarter of the box. 
	int i;
	double ix, iy, iz;

	while ((edge * edge * edge) < N)edge++;

	ix = iy = iz = 0;
	// assigning particle positions
	for (i = 0; i < N; i++) {

		r[i][0] = ix * L / edge;
		r[i][1] = iy * L / edge;
		r[i][2] = iz * L / edge;

		//r[i][0] = ((double)ix + lattice_spacing) * L / edge;
		//r[i][1] = ((double)iy + lattice_spacing) * L / edge;
		//r[i][2] = ((double)iz + lattice_spacing) * L / edge;
		ix++;
		if (ix == edge) {
			ix = 0;
			iy++;
			if (iy == edge) {
				iy = 0;
				iz++;
			}
		}
	}
}

void Config_atoms_rand(double r[N][DIM]) {
	int i, k;

	/// <Random Configuration for many particles>
	/* initialize random seed: */

	srand((unsigned)time(0));

	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {

			/* generate random number between 0.5 to 1.5: */
			r[i][k] = (((double)rand() * (L - 0.1)) / (double)RAND_MAX + 1.0);

		}
	}

}

// Initialize Velocity
void Config_velocity(double v[N][DIM]) {
	
	int i;
	double v_com[DIM] = { 0 };
	double v_mag = sqrt((3 * (N * Temp * kB))/mass);

	// the random number generator generate random number between 1-2
	for (i = 0; i < N; i++) {
	//cout << "v_mag is = " << v_mag << endl;
	v[i][0]= (v_mag * (((double)rand() * (10.0 - 0.1)) / (double)RAND_MAX + 1.0));
	v[i][1] = (v_mag * (((double)rand() * (10.0 - 0.1)) / (double)RAND_MAX + 1.0));
	v[i][2] = (v_mag * (((double)rand() * (10.0 - 0.1)) / (double)RAND_MAX + 1.0));

	//cout << "vx, vy, vz is = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	v_com[0] += v[i][0];
	v_com[1] += v[i][1];
	v_com[2] += v[i][2];
	
	}
	//cout << "v_com x, y, z is = " << v_com[0] << "," << v_com[1] << "," << v_com[2] << endl;
	v_com[0] /= N; //averaged the velocities
	v_com[1] /= N;
	v_com[2] /= N;

	//cout << "v_com x, y, z after averaged is = " << v_com[0] << "," << v_com[1] << "," << v_com[2] << endl;
	for (i = 0; i < N; i++) {
		v[i][0] -= v_com[0];
		v[i][1] -= v_com[1];
		v[i][2] -= v_com[2];

	//cout << "vx, vy, vz after substraction is = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}

}

// Makes total momentum to be zero


// Initialize Temperature
void Config_Temp(double v[N][DIM], double& Temp) {
	int i;
	double v_sq = 0;
	for (i = 0; i < N; i++) {
		v_sq += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
	}
	Temp = (mass * v_sq / (3 * (N * kB)));
}


// Force calculation
void Force_Calc(double r[N][DIM], double F[N][DIM], double& E_pot) {


	// Force,a vector (rij,scalar) = - (1/rij,scalar) * ((d/dr) * U(rij,scalar)) *(rij,vector)
	// ((d/dr) * U(rij,scalar)) = -(48/rij,scalar) *((1/rij12)-0.5*(1/rij6))
	// U(rij,scalar) = 4* (1/rij12 - 1/rij6)
	// rij is given as a condition or taken as an input from the integrator

	int i, j, k, n;
	double rij_vec[DIM] = { 0 };			// vector rij containing the rij_x,rij_y, rij_z--the x,y,z component of rij.
	double rij = 0;							// the magnitude of vector rij (vector r between particle i and j).
	double rij2 = 0;						// squared value of the magnitude of vector rij (vector r between particle i and j).
	double dU_LJ = 0;						// spatial derivative of potential energy dU / dr
	double U_LJ = 0;						// Potential Energy
	const int size_U_LJ = (N * (N - 1)) / 2;			//size for the storage
	double potential[size_U_LJ] = { 0 };				// storage array for the summation of potential energy

	//for simplicity

	const double rcut = (0.5 * L);			// Cut-off Radius
	double sigma2 = sigma * sigma;			// sigma squared
	//double rij1 = 0;						// rij
	double rij6 = 0;						// sigma^6 / rij^6
	double rij12 = 0;						// sigma^12 /rij^12


	E_pot = 0;
	for (n = 0; n < N; n++) {
		for (k = 0; k < DIM; k++) {
			F[n][k] = 0;
			//cout << "are my force all zero ? " << F[n][k] << endl;
		}
	}


	// Looping for pairwise interactions

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {

			//for dimension, where 0 = x,  1 = y = 1, 2 = z

			rij_vec[0] = PBC(r[j][0], r[i][0]);
			rij_vec[1] = PBC(r[j][1], r[i][1]);
			rij_vec[2] = PBC(r[j][2], r[i][2]);


			rij = rij_vec[0] * rij_vec[0] + rij_vec[1] * rij_vec[1] + rij_vec[2] * rij_vec[2];
			//cout << "rij before forcepbc = " << rij << endl;

			forcePBC(rij);
			//cout << "rij after forcepbc = " << rij << endl;

			if (sqrt(rij) < rcut) {
				rij2 = sigma2 / rij;							// sigma^2 / rij^2
				rij6 = rij2 * rij2 * rij2;						// sigma^6 / rij^6
				rij12 = rij6 * rij6;							// sigma^12 / rij^12

				U_LJ = 4 * epsilon * (rij12 - rij6);			// Lennard Jones Potential
				//cout << "U_LJ is " << U_LJ << endl;
				potential[i] = U_LJ;
				E_pot += potential[i];
				//cout <<"the i and j is:"<<i<<","<<j<< " and my U_LJ= " << potential[i][j] << endl;
				//E_pot += potential [i][j];
				//cout << "my E_pot after summation is = " << E_pot << endl;


				//cout << "my E_pot = " << E_pot << endl;
				//computing for spatial derivative of potential energy dU/dr 
				dU_LJ = -48 * epsilon / rij * (rij12 - 0.5 * rij6);
				//cout << "dU_LJ is " << dU_LJ << endl;

				//Force exerted on particle #N 
				//cout << "F before force calc is " << F[i][k] << endl;
				F[i][0] += dU_LJ * rij_vec[0];
				F[i][1] += dU_LJ * rij_vec[1];
				F[i][2] += dU_LJ * rij_vec[2];
				F[j][0] -= dU_LJ * rij_vec[0];
				F[j][1] -= dU_LJ * rij_vec[1];
				F[j][2] -= dU_LJ * rij_vec[2];

			}

			//cout << "my E_pot after summation is = " << E_pot << endl;
		}

	}

}

// Velocity-Verlet Integrator
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& U_LJ, double& E_kin) {

	int i;
	//constants

	double half_dt = 0.5 * dt;					// Half of dt
	double v_half[N][DIM];						// Velocity at half_dt
	//E_kin = 0;
	double kinetic[N];
	E_kin = 0;
	

	//for (i = 0; i < N; i++) {
		//cout << "before the start of verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	//}
	//updating velocity at half time step
	for (i = 0; i < N; i++) {
		//cout << "before verlet v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//cout << "The force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;

		//cout << "Inside the first verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;

		v_half[i][0] = v[i][0] + half_dt * F[i][0] / mass;
		v_half[i][1] = v[i][1] + half_dt * F[i][1] / mass;
		v_half[i][2] = v[i][2] + half_dt * F[i][2] / mass;


		//updating position 

		r[i][0] = r[i][0] + v_half[i][0] * dt;
		r[i][1] = r[i][1] + v_half[i][1] * dt;
		r[i][2] = r[i][2] + v_half[i][2] * dt;

		//cout << "after verlet r values (x,y,z) = " << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;


	//recalculate force
		Force_Calc(r, F, U_LJ);


		//updating velocity at t + dt

			//cout << "Recalculated force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;


		v[i][0] = v_half[i][0] + half_dt * F[i][0] / mass;
		v[i][1] = v_half[i][1] + half_dt * F[i][1] / mass;
		v[i][2] = v_half[i][2] + half_dt * F[i][2] / mass;

		//cout << "after the verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;


		//E_kin = mass * 0.5 * v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
		//cout << "my E_kin = " << E_kin << endl;
		//E_kin += E_kin;
		//cout << "my E_kin_summation = " << E_kin << endl;
		
		//cout << "The E_kin is : " << kinetic[i] << endl;
		

		//velocity scaling for adjusting temperature
		
		v[i][0] = v[i][0] * sqrt(Init_Temp / Temp);
		v[i][1] = v[i][1] * sqrt(Init_Temp / Temp);
		v[i][2] = v[i][2] * sqrt(Init_Temp / Temp);

		//updating kinetic energy 
		kinetic[i] = (mass * 0.5) * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
		//cout << "The E_kin is : " << kinetic[i] << endl;
		E_kin += kinetic[i];
	}
	
	//E_kin_Calc(v,E_kin);
	//cout << "my E_kin_summation = " << E_kin << endl;
}
// Velocity-Verlet Integrator for the first step
void Vel_Verlet_first(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& U_LJ, double& E_kin) {

	int i;
	//constants

	double half_dt = 0.5 * dt;					// Half of dt
	double v_half[N][DIM];						// Velocity at half_dt
	//E_kin = 0;
	double kinetic[N];
	E_kin = 0;

	//for (i = 0; i < N; i++) {
		//cout << "before the start of verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	//}
	//updating velocity at half time step
	for (i = 0; i < N; i++) {
		//cout << "before verlet v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//cout << "The force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;

		//cout << "Inside the first verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;

		v_half[i][0] = v[i][0] + half_dt * F[i][0] / mass;
		v_half[i][1] = v[i][1] + half_dt * F[i][1] / mass;
		v_half[i][2] = v[i][2] + half_dt * F[i][2] / mass;


		//updating position 

		r[i][0] = r[i][0] + v_half[i][0] * dt;
		r[i][1] = r[i][1] + v_half[i][1] * dt;
		r[i][2] = r[i][2] + v_half[i][2] * dt;

		//cout << "after verlet r values (x,y,z) = " << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;


	//recalculate force
		Force_Calc(r, F, U_LJ);


		//updating velocity at t + dt

			//cout << "Recalculated force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;


		v[i][0] = v_half[i][0] + half_dt * F[i][0] / mass;
		v[i][1] = v_half[i][1] + half_dt * F[i][1] / mass;
		v[i][2] = v_half[i][2] + half_dt * F[i][2] / mass;

		//cout << "after the verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;


		//E_kin = mass * 0.5 * v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
		//cout << "my E_kin = " << E_kin << endl;
		//E_kin += E_kin;
		//cout << "my E_kin_summation = " << E_kin << endl;

		//cout << "The E_kin is : " << kinetic[i] << endl;


		//velocity scaling for adjusting temperature
		Config_Temp(v,Temp);
		v[i][0] = v[i][0] * sqrt(Init_Temp / Temp);
		v[i][1] = v[i][1] * sqrt(Init_Temp / Temp);
		v[i][2] = v[i][2] * sqrt(Init_Temp / Temp);

		//updating kinetic energy 
		kinetic[i] = (mass * 0.5) * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
		//cout << "The E_kin is : " << kinetic[i] << endl;
		E_kin += kinetic[i];
	}

	//E_kin_Calc(v,E_kin);
	//cout << "my E_kin_summation = " << E_kin << endl;
}

// Kinetic Energy Calculation
void E_kin_Calc(double v[N][DIM], double& E_kin) {

	int i;
	double kinetic[N];
	E_kin = 0;

	for (i = 0; i < N; i++) {

		//cout << "inside kin loop calc v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		kinetic[i] = mass * 0.5 * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
		//cout << "The E_kin is : " << kinetic[i] << endl;
		E_kin += kinetic[i];
		//cout << "The E_kin after summation is : " << E_kin << endl;

	}


}

// Periodic Boundary Condition for coordinate
double PBC(double x1, double x2) {
	if (x1 - x2 < (-0.5 * L)) {
		return (x1 - x2) + L;
	}

	if (x1 - x2 > (0.5 * L)) {
		return (x1 - x2) - L;
	}
	else {
		return x1 - x2;
	}
	return x1, x2;
}

// Periodic Boundary Condition for force
double forcePBC(double x) {
	if (x < (-0.5 * L)) {
		return  x + L;
	}

	if (x > (0.5 * L)) {
		return x - L;
	}
	else {
		return x;
	}
	return x;
}

void rdf(double r[N][DIM], double& rij, double& count) {
	/*
	g_r  = histogram / rho_bulk *N*(volume of shell = 4*PI *r^2*dr)
	//rho_bulk is Number of particle / box Volume or N/(L*L*L)
	*/
	int i, j, n,k;
	double rij_vec[DIM] = { 0 };
	//double rij = 0;
	const double bin = (L / delta_r);					// array size for the total no of bins.
	const double rcut = (0.5 * L);				// Cut-off Radius
	//double count = 0;
	double histogram[size_pair] = { 0 };
	n = 0;
	double Volume = L * L * L;				//volume of simulation box
	double g_r = 0;
	double Prob = 0;
	double g_r_hist = 0;
	double hn = 0;

	for (n = 2; n < bin + 2; n++) {
		count = 0;
		for (k = 0; k < size_pair; k++) {
			for (i = 0; i < N - 1; i++) {
				for (j = i + 1; j < N; j++) {
					//cout << "my i is " << i << endl;
					//cout << "my j is " << j << endl;
					rij_vec[0] = r[j][0] - r[i][0];
					rij_vec[1] = r[j][1] - r[i][1];
					rij_vec[2] = r[j][2] - r[i][2];
					rij = sqrt(rij_vec[0] * rij_vec[0] + rij_vec[1] * rij_vec[1] + rij_vec[2] * rij_vec[2]);
					//cout << "my rij is " << rij << endl;
					//cout << "my n-1 *delta_r is " << ((n-1.0) * delta_r) << endl;
					//cout << "my n*delta_r is " << (n * delta_r) << endl;
					if (rij > ((n - 1.0) * delta_r) && rij < (n * delta_r)) {

						count++;
						
							//histogram[k] = count;
							//cout <<"my histogram k is "<< histogram[k] << endl;
					}
				}
			}
			//hn += histogram[k];
		}
		
		//hn = hn / N;
		g_r = (count * (Volume)) / (N * 4.0 * PI * (((n - 1.0) * delta_r) * ((n - 1.0) * delta_r)) * delta_r * N);
		Prob = count / size_pair;
		//g_r_hist = (hn  * (Volume)) / (N * 4.0 * PI * (((n - 1.0) * delta_r) * ((n - 1.0) * delta_r)) * delta_r * N);;
		myfile << ((n - 1.0) * delta_r) / sigma << "," << n * delta_r << "," << count << "," << g_r << "," <<g_r_hist<<","<< Prob << endl;
	}

}