/*	Author: Dominikus Brian
	"Lennard-Jones Particle Cluster MD simulation"
	2020/June/25
	Version 12.0.2020
*/
/*----------------------Program begin here---------------------- */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include <fstream>
using namespace std;


////////////////////////// Global parameters  ////////////////////////////////

const int N = 32;						// N is the number of particle within the system
const int DIM = 3;						// DIM is dimension of the problem
const int tmax = 100000;					// maximum number of timestep
const double rho = 1;					// System's volume density
const double L = N / rho;					// Box size

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
const double kb = 1.3807;				// 10^-23 J /K
ofstream myfile;

/////////////////////////////*Function declaration:*/////////////////////////////////

//Configuring position of atoms
void Config_atoms(double r[N][DIM], double v[N][DIM]);
void Config_atoms_2(double r[N][DIM], double v[N][DIM]);
void Config_atoms_rand(double r[N][DIM], double v[N][DIM]);
//void Config_atoms_bcc(double r[N][DIM], double v[N][DIM]);
//void Config_atoms_fcc(double r[N][DIM], double v[N][DIM]);
void Config_atoms_lattice(double r[N][DIM], double v[N][DIM]);

//Calculating force 
void Force_Calc(double r[N][DIM], double F[N][DIM], double& U_LJ);

//Integration with Velocity-Verlet
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& U_LJ);

//Calculate Kinetic Energy
void E_kin_Calc(double v[N][DIM], double& E_kin);

//Periodic Boundary Condition
double PBC(double x1, double x2);
double forcePBC(double x);

//Radial Distribution Function
//void rdf(double r[N][DIM], double& g_r);



/*------------------------Main function begin here------------------------ */

int main() {

	int i, t;

	double r[N][DIM];				// Particle Position vector in x y z
	double v[N][DIM];				// Particle velocity vector in x y z
	double F[N][DIM];				// Force exerted on atom i, in x y z
	double E_kin = 0;				// Kinetic Energy	
	double E_pot = 0;				// Potential Energy
	double E_tot = 0;				// Total Energy	
	double Temp = 0;				// Temperature
	double U_LJ = 0;				// Lennard jones pair potential


	cout << "------------------BEGIN INITIALIZATION ----------------------" << endl;

	//Initialize position and velocity 
	Config_atoms_lattice(r, v);

	//initialize force
	Force_Calc(r, F, E_pot);

	//Write initial configuration to config.csv 
	myfile.open("config.csv");
	myfile << "#N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;
	for (i = 0; i < N; i++) {
		myfile << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}
	myfile.close();

	cout << "------------------BEGIN PROPAGATION ----------------------" << endl;
	//Write Energy value at each time step into properties.csv
	myfile.open("properties.csv");
	myfile << "t" << "," << "Potential" << "," << "Kinetic" << "," << "Total" << "," << "Temperature" << endl;
	for (t = 0; t < tmax; t++) {
		E_pot = 0;
		//cout << "my E_Kin before verlet = " << E_kin << endl;
		//for (i = 0; i < N; i++) {
		//	cout << "velocity before verlet" << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//}
		
		Vel_Verlet(r, v, F, U_LJ);
		E_pot += U_LJ;
		//cout << "my E_Pot = " << E_pot << endl;
		//for (i = 0; i < N; i++) {
		//	cout << "velocity after verlet" << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//}
		//cout << "my E_Kin after verlet= " << E_kin << endl;
		E_kin_Calc(v, E_kin);
		//cout << "my E_Kin  after recalculation = " << E_kin << endl;

		//for (i = 0; i < N; i++) {
		//	cout << "velocity after recalculation" << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//}
		E_tot = E_kin + E_pot;
		Temp = (E_kin / N) * 2 / 3;
		myfile << t + 1 << "," << E_pot << "," << E_kin << "," << E_tot << "," << Temp << endl; \

			if (t == tmax * 0.25)cout << "progress is :| 25%";
		if (t == tmax * 0.50)cout << " | 50%";
		if (t == tmax * 0.75)cout << " | 75%";
		if (t == tmax-1 )cout << " | 100% |" << endl;

	}
	myfile.close();
	cout << "------------------END OF SIMULATION ----------------------" << endl;



	return 0;
}
/*------------------------Main function end here------------------------ */

/*------------------------Function definition------------------------- */

// Functions for specifying initial configurations : 
// 1) Config_Atoms,         // for 3 atoms in a specific arrangement
// 2) Config_atoms_2 ,		// for 2 atoms in a specific arrangement
// 3) Config_atoms_rand		// for N number of atoms with random position and velocity value between 0 and 1
// 'Under construction' 4) Config_atoms_bcc		// for 4*N number of atoms with bcc arrangement
// 'Under construction' 5) Config_atoms_fcc		// for 4*N number of atoms with fcc arrangement
void Config_atoms(double r[N][DIM], double v[N][DIM]) {

	/// <Ordered Configuration for three particles system>
	r[0][0] = 1;
	r[0][1] = 1;
	r[0][2] = 1;
	r[1][0] = 1;
	r[1][1] = 2;
	r[1][2] = 1;
	r[2][0] = 2;
	r[2][1] = 1;
	r[2][2] = 1;

	//assigning velocities
	v[0][0] = 1;
	v[0][1] = 1.5;
	v[0][2] = 2;
	v[1][0] = 1;
	v[1][1] = 1.5;
	v[1][2] = 2;
	v[2][0] = 1;
	v[2][1] = 2;
	v[2][2] = 1.5;
	
	

	
		
			
	

}
void Config_atoms_2(double r[N][DIM], double v[N][DIM]) {
	int i, k;

	/// <Ordered Configuration for two particles system>
	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {
			r[0][k] = 1;
			r[1][k] = 2;
		}
	}

	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {
			v[0][k] = 1;
			v[1][k] = 2;
		}
	}

}
void Config_atoms_rand(double r[N][DIM], double v[N][DIM]) {
	int i, k;

	/// <Random Configuration for many particles>
	/* initialize random seed: */

	srand((unsigned)time(0));

	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {

			/* generate random number between 0.5 to 1.5: */
			r[i][k] = 0.5 + ((double)rand() / RAND_MAX);

		}
	}

	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {

			/* generate random number between 0.5 to 1.5: */
			v[i][k] = 0.5 + ((double)rand() / RAND_MAX);

		}
	}
}
void Config_atoms_lattice(double r[N][DIM], double v[N][DIM]) {
	int edge = 2;    // For rho = 0.5 >> edge = 2 gives fill up half of the box, 4 fill up a quarter of the box. 
	int i, k;
	double ix, iy, iz;

	while ((edge*edge*edge) < N)edge++;

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

	//assigning random particle velocities
	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {
			//v[i][k] = { 0 };
			/* generate random number between 0.5 to 1.5: */
			v[i][k] = ((double)rand() * (10.0 - 1.0)) / (double)RAND_MAX + 1.0;

		}
	}
	;
}

/////////////////////////////THE END OF CONFIG_ATOMS ////////////////////////////////

// Force calculation
void Force_Calc(double r[N][DIM], double F[N][DIM], double &U_LJ) {


	// Force,a vector (rij,scalar) = - (1/rij,scalar) * ((d/dr) * U(rij,scalar)) *(rij,vector)
	// ((d/dr) * U(rij,scalar)) = -(48/rij,scalar) *((1/rij12)-0.5*(1/rij6))
	// U(rij,scalar) = 4* (1/rij12 - 1/rij6)
	// rij is given as a condition or taken as an input from the integrator

	int i, j, k, n;
	double rij_vec[DIM] = { 0 };			// vector rij containing the rij_x,rij_y, rij_z--the x,y,z component of rij.
	double rij = 0;							// the magnitude of vector rij (vector r between particle i and j).
	double rij2 = 0;						// squared value of the magnitude of vector rij (vector r between particle i and j).
	double dU_LJ = 0;						// spatial derivative of potential energy dU / dr
	//double U_LJ = 0;						// Total Potential Energy

	//for simplicity

	const double rcut = 0.5 * L;			// Cut-off Radius
	double sigma2 = sigma * sigma;			// sigma squared
	//double rij1 = 0;						// rij
	double rij6 = 0;						// sigma^6 / rij^6
	double rij12 = 0;						// sigma^12 /rij^12


	//E_pot = 0;
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

			if (rij < rcut * rcut) {
				rij2 = sigma2 / rij;							// sigma^2 / rij^2
				rij6 = rij2 * rij2 * rij2;						// sigma^6 / rij^6
				rij12 = rij6 * rij6;							// sigma^12 / rij^12

				U_LJ = 4 * epsilon * (rij12 - rij6);			// Lennard Jones Potential
				//cout << "U_LJ is " << U_LJ << endl;
				//E_pot += U_LJ;
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


			
		}
	}
	
}


// Velocity-Verlet Integrator
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& U_LJ) {

	int i;
	//constants
	double dt = 0.0001;							// timestep size
	double half_dt = 0.5 * dt;					// Half of dt
	double v_half[N][DIM];						// Velocity at half_dt


	//updating velocity at half time step
	for (i = 0; i < N; i++) {
		//cout << "before verlet v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//cout << "The force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;
		v_half[i][0] = v[i][0] + half_dt * F[i][0] / mass;
		v_half[i][1] = v[i][1] + half_dt * F[i][1] / mass;
		v_half[i][2] = v[i][2] + half_dt * F[i][2] / mass;
		
	}

	//updating position 
	for (i = 0; i < N; i++) {
		
		r[i][0] = r[i][0] + v_half[i][0] * dt;
		r[i][1] = r[i][1] + v_half[i][1] * dt;
		r[i][2] = r[i][2] + v_half[i][2] * dt;
		
		//cout << "after verlet r values (x,y,z) = " << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
	}

	//recalculate force
	Force_Calc(r, F, U_LJ);


	//updating velocity at t + dt
	for (i = 0; i < N; i++) {
		//cout << "Recalculated force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;
		

		v[i][0] = v_half[i][0] + half_dt * F[i][0] / mass;
		v[i][1] = v_half[i][1] + half_dt * F[i][1] / mass;
		v[i][2] = v_half[i][2] + half_dt * F[i][2] / mass;
		
		//cout << "after verlet v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}

}

// Kinetic Energy Calculation
void E_kin_Calc(double v[N][DIM], double& E_kin) {

	int i;

	for (i = 0; i < N; i++) {
	
		//cout << "before kin calc v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		E_kin = 0;
		E_kin = mass * v[i][0] * v[i][0]+ v[i][1] * v[i][1] + v[i][2] * v[i][2];
		E_kin = 0.5 * E_kin;
		
	}

}

// Periodic Boundary Condition 
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

//Radial Distribution Function (RDF)
void rdf(double r[N][DIM], double& g_r) {
	int i, k;
	double delta_r = 0.5;
	double rmax = (0.5 * L) / delta_r;



	for (i = 0; i < rmax; i++) {
		for (k = 0; k < DIM; k++) {
			if (r[i][k] < delta_r * i) {
				//count;//
				;
			}

		}
	}
}