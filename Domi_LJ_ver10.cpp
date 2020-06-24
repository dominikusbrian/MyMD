/*	Author: Dominikus Brian
	"Lennard-Jones Particle Cluster MD simulation"
	2020/June/23
	Version 10.0.2020
*/
/*----------------------Program begin here---------------------- */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include <fstream>
using namespace std;


////////////////////////// Global parameters and constants /////////////////////////

const int N = 3;						// N is the number of particle within the system
const int DIM = 3;						// DIM is dimension of the problem
const int tmax = 10000;					// maximum number of timestep
const double L = 3;						// Box size
const double epsilon = 1;				// Energy scale constant
const double sigma = 1;					// Length scale constant
const double mass = 1;					// Mass of the particle
ofstream myfile;

/////////////////////////////*Function declaration:*/////////////////////////////////

//Configuring position of atoms
void Config_atoms(double r[N][DIM], double v[N][DIM]);
void Config_atoms_2(double r[N][DIM], double v[N][DIM]);
void Config_atoms_rand(double r[N][DIM], double v[N][DIM]);
//void Config_atoms_bcc(double r[N][DIM], double v[N][DIM]);
//void Config_atoms_fcc(double r[N][DIM], double v[N][DIM]);

//Calculating force 
void Force_Calc(double r[N][DIM], double F[N][DIM], double &E_pot);

//Integration with Velocity-Verlet
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double &E_pot);

//Calculate Kinetic Energy
void E_kin_Calc(double v[N][DIM], double &E_kin);

//Periodic Boundary Condition
double PBC(double x1, double x2);

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


	cout << "------------------BEGIN INITIALIZATION ----------------------" << endl;

	//Initialize position and velocity 
	Config_atoms(r, v);

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
		Vel_Verlet(r, v, F, E_pot);
		//cout << "my E_Pot = " << E_pot << endl;
		E_kin_Calc(v, E_kin);
		//cout << "my E_Kin = " << E_kin << endl;
		E_tot = E_kin + E_pot;
		Temp = (E_kin / N) * 2 / 3;
		myfile << t + 1 << "," << E_pot << "," << E_kin << "," << E_tot << "," << Temp << endl;
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
	int i, k;

	/// <Ordered Configuration for three particles system>
	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {
			r[0][k] = 3;
			r[1][k] = 2;
			r[2][k] = 1;
		}
	}

	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {
			v[0][k] = 0.5;
			v[1][k] = 1;
			v[2][k] = 1.5;
		}
	}

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

/////////////////////////////THE END OF CONFIG_ATOMS ////////////////////////////////

// Force calculation
void Force_Calc(double r[N][DIM], double F[N][DIM], double &E_pot) {


	// Force,a vector (rij,scalar) = - (1/rij,scalar) * ((d/dr) * U(rij,scalar)) *(rij,vector)
	// ((d/dr) * U(rij,scalar)) = -(48/rij,scalar) *((1/rij12)-0.5*(1/rij6))
	// U(rij,scalar) = 4* (1/rij12 - 1/rij6)
	// rij is given as a condition or taken as an input from the integrator

	int i, j, k, n;
	double rij_vec[DIM] = { 0 };				// vector rij containing the rij_x,rij_y, rij_z--the x,y,z component of rij.
	double rij = 0;							// the magnitude of vector rij (vector r between particle i and j).
	double rij2 = 0;						// squared value of the magnitude of vector rij (vector r between particle i and j).
	double dU_LJ = 0;						// spatial derivative of potential energy dU / dr
	double U_LJ = 0;						// Total Potential Energy

	//for simplicity

	const double rcut = 10;			// Cut-off Radius
	double rij6;						// rij^6
	double rij12;						// rij^12
	double sigma2 = sigma * sigma;		// sigma squared


	for (n = 0; n < N; n++) {
		for (k = 0; k < DIM; k++) {
			F[n][k] = 0;
		}
	}
	E_pot = 0;

	// Looping for pairwise interactions

	for (i = 0; i < N - 1; i++) {
		for (j = i + 1; j < N; j++) {

			//for dimension loop , where 0 = x,  1 = y = 1, 2 = z
			for (k = 0; k < DIM; k++) {

				rij_vec[k] = PBC(r[j][k], r[i][k]);
				rij = sqrt((rij_vec[0] * rij_vec[0]) + (rij_vec[1] * rij_vec[1]) + (rij_vec[2] * rij_vec[2]));
				if (rij < rcut) {
					rij2 = (sigma2 / (rij * rij));					// sigma^2 / rij^2
					rij6 = rij2 * rij2 * rij2;						// sigma^6 / rij^6
					rij12 = rij6 * rij6;							// sigma^12 / rij^12
					U_LJ = 4 * epsilon * (rij12 - rij6);			// Lennard Jones Potential
					//cout << "U_LJ is " << U_LJ << endl;
					E_pot = U_LJ;

					//computing for spatial derivative of potential energy dU/dr 
					dU_LJ = (-(48 / rij) * (rij12 - (0.5 * rij6)));
					//cout << "dU_LJ is " << dU_LJ << endl;

					//Force exerted on particle #N 
					F[i][k] = F[i][k] + (dU_LJ * rij_vec[k]);
					F[j][k] = F[j][k] - (dU_LJ * rij_vec[k]);

					//cout << "F[i][k] is " << F[i][k] << endl;

				}
				else
				{
					continue;
				}

			}
		}

	}
}




// Velocity-Verlet Integrator
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double &E_pot) {

	int i, k;
	//constants
	double dt = 0.0001;							// timestep size
	double half_dt = 0.5 * dt;					// Half of dt
	double v_half[N][DIM];						// Velocity at half_dt


	//updating velocity at half time step
	for (i = 0; i < N; i++) {
		//cout << "before verlet v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//cout << "The force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;
		for (k = 0; k < DIM; k++) {
			//F[i][k] = { 0.0 };
			v_half[i][k] = v[i][k] + half_dt * (F[i][k] / mass);
		}
	}

	//updating position 
	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {
			r[i][k] = r[i][k] + v_half[i][k] * dt;
		}
		//cout << "after verlet r values (x,y,z) = " << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
	}

	//recalculate force
	Force_Calc(r, F, E_pot);


	//updating velocity at t + dt
	for (i = 0; i < N; i++) {
		//cout << "Recalculated force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;
		for (k = 0; k < DIM; k++) {

			v[i][k] = v_half[i][k] + half_dt * (F[i][k] / mass);
		}
		//cout << "after verlet v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}

}

// Kinetic Energy Calculation
void E_kin_Calc(double v[N][DIM], double &E_kin) {

	int i, k;

	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {
			//cout << "before kin calc v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
			E_kin = 0;
			E_kin += (v[i][k] * v[i][k]);
			E_kin = 0.5 * E_kin;
		}
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

//Radial Distribution Function (RDF)
void rdf(double r[N][DIM], double &g_r) {
	int i,k;
	double delta_r = 0.5;
	double rmax = (0.5 * L)/delta_r;



	for (i = 0; i < rmax; i++) {
		for (k = 0; k < DIM; k++) {
			if ( r[i][k] < delta_r*i) {
				//count;//
				;
			}

		}
	}
}