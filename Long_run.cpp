
/*	Author: Dominikus Brian
	"Lennard-Jones Particle Cluster MD simulation"
	2020/July 10/2020
	NewLJ Version 1.0.2020
*/
/*----------------------Program begin here---------------------- */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include <fstream>
#include <time.h>
#include <random>
using namespace std;
default_random_engine e(23238);

////////////////////////// Global parameters  ////////////////////////////////

const int N = 108;							// N is the number of particle within the system
const int DIM = 3;							// DIM is dimension of the problem
const int tmax_eq = 75000;						// maximum number of timestep or simulation duration
const int tmax_prod = 20000;
const double rho = 0.8;						// System's volume density
double dt = 0.0025;							// Simulation timestep size
double Init_Temp = 1;						// Temperature
const double Volume = N / rho;
double edge = cbrt(Volume);					// Box edge size
const double delta_r = 0.0125;					// bin width for rdf analysis
const int size_pair = (N * (N - 1)) / 2;	// array size for the total number of pairs 
const int size_bin = 200;



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

// file output list.
ofstream myfile;
ofstream otherfile;
ofstream trajfile;
ofstream cvvfile;
ofstream fullprop;

// filestream control
int startcvv = 0;
int vel_scaling_interval = 1000;  // perform velocity scaling every start_vel_scaling steps
int finish_vel_scaling = 50000;
int startrdf = 0;
int write_traj = 500;
int start_write_traj = 10000;


/////////////////////////////*Function declaration:*/////////////////////////////////

// Initialize Particle Position
void Config_atoms_lattice(double r[N][DIM]);
void Config_atoms_rand(double r[N][DIM]);

// Initialize Particle Velocity
void Config_velocity(double v[N][DIM]);

// Scale initial velocity based on Init_Temp
double Config_Temp(double v[N][DIM], double& Temp);

// Force Calculation
void Force_Calc(double r[N][DIM], double F[N][DIM], double& E_pot);

// Velocity-Verlet propagator
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& E_pot, double& E_kin, double& Temp);
void Vel_Verlet_first(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& E_pot, double& E_kin, double& Temp);

//Periodic Boundary Condition
void forcePBC(double& rij);

//radial distribution function  
void rdf(double r[N][DIM], double g_r[size_bin], double g_ave[size_bin], int& divider);

//Velocity scaling
void Vel_scaling(double v[N][DIM], double& E_kin, double& Temp_ave);

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
	int divider = 0;				// Counter for dividing accumulator
	int counter = 0;

	double rij = 0;					// Pair distance
	double Temp = Init_Temp;
	double g_r[size_bin] = { 0 };
	double g_ave[size_bin] = { 0 };	// Accumulated average of rdf for a given frame
	const int N_frame = 300;
	const int storage = 300;
	double vx[storage][N] = { 0 }, vy[storage][N] = { 0 }, vz[storage][N] = { 0 };
	double cvv_dt[N] = { 0 };
	//double cvv_accum = 0;
	double cvv_ave[N_frame] = { 0 };
	double cvv_accum[N_frame] = { 0 };
	double cvv_zero[N_frame] = { 0 };
	int k = 0;
	int j = 0;
	int x = 0;
	double Temp_ave = 0;
	int Delta_t = 1;
	int cvvcount = 1;






	cout << "------------------BEGIN SIMULATION----------------------" << endl;

	// Displaying Input parameters and Global constants
	cout << "INPUT PARAMETERS: " << endl;
	cout << "N = " << N << endl;
	cout << "dt = " << dt << endl;
	cout << "tmax_eq = " << tmax_eq << endl;
	cout << "tmax_prod = " << tmax_prod << endl;
	cout << "rho = " << rho << endl;
	cout << "dr = " << delta_r << endl;
	cout << endl;
	cout << "GLOBAL CONSTANTS: " << endl;
	cout << "Epsilon = " << epsilon << endl;
	cout << "Sigma = " << sigma << endl;
	cout << "Mass = " << mass << endl;
	cout << endl;
	cout << "FILESTREAM CONTROL: " << endl;
	cout << "startcvv = " << startcvv << endl;
	cout << "vel_scaling_interval = " << vel_scaling_interval << endl;
	cout << "finish_vel_scaling = " << finish_vel_scaling << endl;
	cout << "startrdf = " << startrdf << endl;
	cout << "start_write_traj = " << start_write_traj << endl;
	cout << "write_traj = " << write_traj << endl;


	cout << "------------------BEGIN INITIALIZATION----------------------" << endl;

	//Initialize Position
	Config_atoms_lattice(r);
	//Initialize Velocity
	Config_velocity(v);


	//Write INitial Configuration to initconfig.txt
	myfile.open("initconfig.txt");
	for (i = 0; i < N; i++) myfile << r[i][0] << endl;
	for (i = 0; i < N; i++) myfile << r[i][1] << endl;
	for (i = 0; i < N; i++) myfile << r[i][2] << endl;
	for (i = 0; i < N; i++) myfile << v[i][0] << endl;
	for (i = 0; i < N; i++) myfile << v[i][1] << endl;
	for (i = 0; i < N; i++) myfile << v[i][2] << endl;
	myfile.close();


	//Export initial configuration 
	myfile.open("Init_config.csv");
	myfile << "N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;
	for (i = 0; i < N; i++) {
		myfile << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}
	myfile.close();

	//Initialize Force
	Force_Calc(r, F, E_pot);

	// Export a mdinfo file about the input parameter and global constant employed in the simulation. 
	myfile.open("mdinfo.txt");
	myfile << "INPUT PARAMETERS: " << endl;
	myfile << "N = " << N << endl;
	myfile << "dt = " << dt << endl;
	myfile << "tmax_eq = " << tmax_eq << endl;
	myfile << "tmax_prod = " << tmax_prod << endl;
	myfile << "rho = " << rho << endl;
	myfile << "dr = " << delta_r << endl;
	myfile << "Temperature = " << Temp << endl;
	myfile << endl;
	myfile << "GLOBAL CONSTANTS: " << endl;
	myfile << "Epsilon = " << epsilon << endl;
	myfile << "Sigma = " << sigma << endl;
	myfile << "Mass = " << mass << endl;
	myfile << endl;
	myfile << "FILESTREAM CONTROL: " << endl;
	myfile << "startcvv = " << startcvv << endl;
	myfile << "vel_scaling_interval = " << vel_scaling_interval << endl;
	myfile << "finish_vel_scaling = " << finish_vel_scaling << endl;
	myfile << "startrdf = " << startrdf << endl;
	myfile << "start_write_traj = " << start_write_traj << endl;
	myfile << "write_traj = " << write_traj << endl;
	myfile.close();


	//Set momentum to zero
	for (i = 0; i < N; i++) {
		v[i][0] = 0;
		v[i][1] = 0;
		v[i][2] = 0;
	}


	cout << "------------------INITIALIZATION COMPLETE ----------------------" << endl;
	cout << "------------------BEGIN EQUILIBRATION ----------------------" << endl;
	//Write Energy value at each time step into properties.csv

	trajfile.open("traj.csv");
	trajfile << "t" << "N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;


	fullprop.open("fullprop_eq.txt");
	fullprop << "t" << "," << "Potential" << "," << "Kinetic" << "," << "Total" << "," << "Temperature" << endl;

	for (t = 0; t < tmax_eq; t++) {
		//cout << "edge" << edge << endl;
		//cout << "Temp is "<<Temp<<endl;;
		//cout << "the t is : " << t << endl;
		for (i = 0; i < N; i++) {
			//	cout <<"current position " << "at the beginning of t " << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
		}
		if (t == 0) {

			Vel_Verlet_first(r, v, F, E_pot, E_kin, Temp);
			//cout << "first verlet is done " << endl;
		}

		Vel_Verlet(r, v, F, E_pot, E_kin, Temp);
		E_tot = E_kin + E_pot;
		Temp = ((E_kin / N) * 2 / (3 * kB));

		fullprop << t + 1 << "," << E_pot << "," << E_kin << "," << E_tot << "," << Temp << endl;





		// Getting average temperature from previous 100 instants, when t is a multiple of 1000;
		if (t > (1000 * j) - 400 && t <= (1000 * j)) {

			//cout << "Temp is " << Temp << endl;
			Temp_ave += Temp;


			if (t == (1000 * j)) {
				Temp_ave = Temp_ave / 400;
				//cout << "Temp_ave is " << Temp_ave << endl;
				j++;
			}
		}

		if (t > 4000 && t % vel_scaling_interval == 0 && t < finish_vel_scaling) {
			Vel_scaling(v, E_kin, Temp_ave);
		}



		//Positional PBC
		for (i = 0; i < N; i++) {
			while (r[i][0] > (0.5 * edge)) r[i][0] -= edge;
			while (r[i][1] > (0.5 * edge)) r[i][1] -= edge;
			while (r[i][2] > (0.5 * edge)) r[i][2] -= edge;
			while (r[i][0] < (-0.5 * edge)) r[i][0] += edge;
			while (r[i][1] < (-0.5 * edge)) r[i][1] += edge;
			while (r[i][2] < (-0.5 * edge)) r[i][2] += edge;

		}

		// record trajectory every write_traj steps
		if (t > start_write_traj && t % write_traj == 0) {
			for (i = 0; i < N; i++) {
				trajfile << t + 1 << "," << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
			}
		}


		if (t == tmax_eq * 0)cout << "The progress is : " << endl;
		if (t == tmax_eq * 0.10)cout << " | 10%";
		if (t == tmax_eq * 0.25)cout << " | 25%";
		if (t == tmax_eq * 0.50)cout << " | 50%";
		if (t == tmax_eq * 0.75)cout << " | 75%";
		if (t == tmax_eq * 0.90)cout << " | 90%";
		if (t == tmax_eq - 1)cout << " | 100% |" << endl;


	}

	trajfile.close();
	fullprop.close();

	cout << "------------------EQUILIBRATION FINISHED ----------------------" << endl;



	//Write Final Configuration to final_config.csv
	myfile.open("eq_finalconfig.txt");
	for (i = 0; i < N; i++) myfile << r[i][0] << endl;
	for (i = 0; i < N; i++) myfile << r[i][1] << endl;
	for (i = 0; i < N; i++) myfile << r[i][2] << endl;
	for (i = 0; i < N; i++) myfile << v[i][0] << endl;
	for (i = 0; i < N; i++) myfile << v[i][1] << endl;
	for (i = 0; i < N; i++) myfile << v[i][2] << endl;
	myfile.close();

	cout << "------------------BEGIN PRODUCTION ----------------------" << endl;

	//Write Energy value at each time step into properties.csv

	otherfile.open("rdf.csv");
	otherfile << "r" << "," << "r_plus_dr" << "," << "Count" << "," << "g_r" << "," << "gr_ave" << "," << divider << "," << "gr_accum" << endl;

	trajfile.open("traj.csv");
	trajfile << "t" << "N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;

	cvvfile.open("cvv.csv");
	cvvfile << "t" << "," << "cvv_accum" << endl;

	fullprop.open("fullprop_prod.txt");
	fullprop << "t" << "," << "Potential" << "," << "Kinetic" << "," << "Total" << "," << "Temperature" << endl;

	for (t = 0; t < tmax_prod; t++) {


		Vel_Verlet(r, v, F, E_pot, E_kin, Temp);
		E_tot = E_kin + E_pot;
		Temp = ((E_kin / N) * 2 / (3 * kB));


		fullprop << t + 1 << "," << E_pot << "," << E_kin << "," << E_tot << "," << Temp << endl;


		// Radial distribution function (RDF) analysis
		if (t >= startrdf && t % 10 == 0) {
			divider++;
			rdf(r, g_r, g_ave, divider);
		}


		//Positional PBC
		for (i = 0; i < N; i++) {
			while (r[i][0] > (0.5 * edge)) r[i][0] -= edge;
			while (r[i][1] > (0.5 * edge)) r[i][1] -= edge;
			while (r[i][2] > (0.5 * edge)) r[i][2] -= edge;
			while (r[i][0] < (-0.5 * edge)) r[i][0] += edge;
			while (r[i][1] < (-0.5 * edge)) r[i][1] += edge;
			while (r[i][2] < (-0.5 * edge)) r[i][2] += edge;

		}

		// record trajectory every write_traj steps
		if (t > start_write_traj && t % write_traj == 0) {
			for (i = 0; i < N; i++) {
				trajfile << t + 1 << "," << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
			}
		}



		//Velocity auto-correlation function (calculating for Cvv)


		k = t - (N_frame * counter);								// index for storage labeling always 1-N_frame

		//cout << "my k is = " << k << endl;
		//cout << "my counter is " << counter << endl;
		//cout << "my t_corr*counter is " << t_corr * counter << endl;
		for (i = 0; i < N; i++) {

			//cout<< i << " My Velocities are: " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;

	// read the value of v for all atoms
			vx[k][i] = v[i][0];
			vy[k][i] = v[i][1];
			vz[k][i] = v[i][2];
			//cout << "v = " << v[i][0] << endl;
		   // cout << "vx = " << vx[k][i] << endl;
		   // cout << "vy = " << vy[k][i] << endl;
		   // cout << "vz = " << vz[k][i] << endl;
		}


		if (k == N_frame) {

			if (t > startcvv) {

				for (x = 0; x < N_frame; x++) {
					//cvv_accum[0] = 0;
					cvv_ave[x] = 0;
					//cvv_dt[i] = { 0 };
					for (i = 0; i < N; i++) {
						// cout << "vx = " << vx[0][i] << endl;
						//cout << "vy = " << vy[0][i] << endl;
						//cout << "vz = " << vz[0][i] << endl;
						cvv_dt[i] = (vx[0][i] * vx[x][i] + vy[0][i] * vy[x][i] + vz[0][i] * vz[x][i]);
						//cout << "my cvv_dt_i = "<<cvv_dt[i] << endl;
						cvv_ave[x] += cvv_dt[i] / N;
						//cout << "my cvv_ave = " << cvv_ave[x] << endl;
						//cvv_ave[x] = (cvv_ave[x] * (1 / t));
						//cvv_zero[0] += (vx[0][i] * vx[0][i] + vy[0][i] * vy[0][i] + vz[0][i] * vz[0][i]);
						//cout << "my cvv_zero = " << cvv_zero[0]/N << endl;
						//cout << "cvv_dt[i] " << cvv_dt[i] << "," << "cvv_ave[x]" << cvv_ave[x] << "," << "cvv_accum[x] " << cvv_accum[x] << endl;
					}
					cvv_accum[x] += cvv_ave[x];

					//cout << ">>>>>>>>>>>>>>>>>> my cvv_accum[x] = " << cvv_accum[x] << endl;
					cvvfile << t << "," << cvv_accum[x] << "," << cvv_zero[0] / N << "," << k << "," << cvv_accum[x] / cvvcount << endl;

					cvvcount++;

				}


			}
			counter++;
		}



		//cout << "counter " << counter << endl;





		if (t == tmax_prod * 0)cout << "The progress is : " << endl;
		if (t == tmax_prod * 0.10)cout << " | 10%";
		if (t == tmax_prod * 0.25)cout << " | 25%";
		if (t == tmax_prod * 0.50)cout << " | 50%";
		if (t == tmax_prod * 0.75)cout << " | 75%";
		if (t == tmax_prod * 0.90)cout << " | 90%";
		if (t == tmax_prod - 1)cout << " | 100% |" << endl;


	}

	myfile.close();
	otherfile.close();
	trajfile.close();
	cvvfile.close();
	fullprop.close();

	cout << "------------------EQUILIBRATION FINISHED ----------------------" << endl;



	//Write Final Configuration to final_config.csv
	myfile.open("prod_finalconfig.txt");
	for (i = 0; i < N; i++) myfile << r[i][0] << endl;
	for (i = 0; i < N; i++) myfile << r[i][1] << endl;
	for (i = 0; i < N; i++) myfile << r[i][2] << endl;
	for (i = 0; i < N; i++) myfile << v[i][0] << endl;
	for (i = 0; i < N; i++) myfile << v[i][1] << endl;
	for (i = 0; i < N; i++) myfile << v[i][2] << endl;
	myfile.close();
	cout << "------------------END OF SIMULATION ----------------------" << endl;
	printf("Total simulation time %.2fs\n", ((double)clock() - tstart) / CLOCKS_PER_SEC);


	return 0;
}



/*------------------------Main function end here------------------------ */

/*------------------------Function definition------------------------- */

// Initialize position
void Config_atoms_rand(double r[N][DIM]) {
	int i, k;

	/// <Random Configuration for many particles>
	/* initialize random seed: */

	srand((unsigned)time(0));

	for (i = 0; i < N; i++) {
		for (k = 0; k < DIM; k++) {

			/* generate random number between 0.5 to 1.5: */
			r[i][k] = (((double)rand() * (10.0 - 1.0)) / (double)RAND_MAX);

		}
	}
}
void Config_atoms_lattice(double r[N][DIM]) {

	int i;
	double ix, iy, iz;


	ix = iy = iz = 0;
	// assigning particle positions

	for (i = 0; i < N; i++) {

		//r[i][0] = ix *  0.2*edge - (((double)rand() * (0.5 - 0.2)) / (double)RAND_MAX);
		//r[i][1] = iy * 0.2*edge - (((double)rand() * (0.5 - 0.2)) / (double)RAND_MAX);
		//r[i][2] = iz * 0.2* edge - (((double)rand() * (0.5 - 0.2)) / (double)RAND_MAX);
		r[i][0] = ix * 0.2*edge;
		r[i][1] = iy * 0.2* edge;
		r[i][2] = iz * 0.2 * edge;

		ix++;
		if (ix >= 0.8*edge) {
			ix = 0;
			iy++;
			if (iy >= 0.8*edge) {
				iy = 0;
				iz++;
			}
		}
	}

}


// Initialize Velocity
void Config_velocity(double v[N][DIM]) {

	int i;
	double v_com[DIM];  //com is center of mass
	double v_mag = sqrt(3 * Init_Temp * kB / mass);

	for (i = 0; i < DIM; i++) { v_com[i] = 0; }
	// the random number generator generate random number between 1-2
	for (i = 0; i < N; i++) {
		//cout << "v_mag is = " << v_mag << endl;

		//v[i][0] = (v_mag * e()/10000000000);
		//v[i][1] = (v_mag * e()/10000000000);
		//v[i][2] = (v_mag * e()/10000000000);

		v[i][0] = (v_mag * (((double)rand() * (1 - 0.1)) / (double)RAND_MAX));
		v[i][1] = (v_mag * (((double)rand() * (1 - 0.1)) / (double)RAND_MAX));
		v[i][2] = (v_mag * (((double)rand() * (1 - 0.1)) / (double)RAND_MAX));

		//cout << "vx, vy, vz is = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		v_com[0] += (v[i][0] / N);
		v_com[1] += (v[i][1] / N);
		v_com[2] += (v[i][2] / N);

	}

	//cout << "v_com x, y, z after averaged is = " << v_com[0] << "," << v_com[1] << "," << v_com[2] << endl;
	for (i = 0; i < N; i++) {
		v[i][0] -= v_com[0];
		v[i][1] -= v_com[1];
		v[i][2] -= v_com[2];

		//cout << "vx, vy, vz after substraction is = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}

}



double Config_Temp(double v[N][DIM], double& Temp) {
	int i;
	double v_sq[N] = { 0 };
	double v_sq_accum = 0;
	for (i = 0; i < N; i++) {
		v_sq[i] = (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
		v_sq_accum += (v_sq[i] / N);
	}
	Temp = mass * v_sq_accum / (3 * kB);
	//cout << "my v_sq_accum is " << v_sq_accum << endl;
	//cout << "my temp fresh from the oven is " << Temp << endl;
	return Temp;
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
	double potential[N] = { 0 };				// storage array for the summation of potential energy

	//for simplicity

	const double rcut = (0.5 * edge);			// Cut-off Radius
	double sigma2 = sigma * sigma;			// sigma squared
	double rij6 = 0;						// sigma^6 / rij^6
	double rij12 = 0;						// sigma^12 /rij^12
	double no_go_radius = 0.5;

	E_pot = 0;
	for (n = 0; n < N; n++) {
		for (k = 0; k < DIM; k++) {
			F[n][k] = 0;
			//cout << "are my force all zero ? " << F[n][k] << endl;
		}
	}


	// Looping for pairwise interactions
	//for (i = 0; i < N; i++) {
			//cout <<"current position inside force calc " << "," << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
		//}

	// position PBC
	for (i = 0; i < N; i++) {
		//cout << "start position PBC" << endl;
		//cout <<"current position " << "," << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
		
		//while (r[i][0] > (0.5 * edge)) r[i][0] -= edge;
		//while (r[i][1] > (0.5 * edge)) r[i][1] -= edge;
		//while (r[i][2] > (0.5 * edge)) r[i][2] -= edge;
		//while (r[i][0] < (-0.5 * edge)) r[i][0] += edge;
		//while (r[i][1] < (-0.5 * edge)) r[i][1] += edge;
		//while (r[i][2] < (-0.5 * edge)) r[i][2] += edge;


		

		//cout <<"position after PBC " << "," << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
			
	//cout << "finished position PBC" << endl;
		//for (i = 0; i < N; i++) {
			//cout <<"current position " << "," << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
		//}
	
			for (i = 0; i < N - 1; i++) {

				for (j = i + 1; j < N; j++) {

					//for dimension, where 0 = x,  1 = y = 1, 2 = z
					rij_vec[0] = (r[j][0] - r[i][0]);
					rij_vec[1] = (r[j][1] - r[i][1]);
					rij_vec[2] = (r[j][2] - r[i][2]);
					//cout << "rij_Vec[k] before pbc " << rij_vec[0]<<"," <<rij_vec[1]<<","<<rij_vec[2]<< endl;


					//while (rij_vec[0] <= 0.5 * edge && rij_vec[0] > -0.5 * edge) {
					//	rij_vec[0] = rij_vec[0];
					//}
					//while (rij_vec[1] <= 0.5 * edge && rij_vec[1] > -0.5 * edge) {
					//	rij_vec[1] = rij_vec[1];
					//}
					//while (rij_vec[2] <= 0.5 * edge && rij_vec[2] > -0.5 * edge) {
					//	rij_vec[2] = rij_vec[2];
					//}
					
					//cout << "start force PBC" << endl;
					while (rij_vec[0] < (-0.5 * edge)) {
						rij_vec[0] = rij_vec[0] + edge;
					}
					while (rij_vec[1] < (-0.5 * edge)) {
						rij_vec[1] = rij_vec[1] + edge;
					}
					while (rij_vec[2] < (-0.5 * edge)) {
						rij_vec[2] = rij_vec[2] + edge;
					}


					while (rij_vec[0] > (0.5 * edge)) {
						rij_vec[0] = rij_vec[0] - edge;
					}
					while (rij_vec[1] > (0.5 * edge)) {
						rij_vec[1] = rij_vec[1] - edge;
					}
					while (rij_vec[2] > (0.5 * edge)) {
						rij_vec[2] = rij_vec[2] - edge;
					}

					//cout << "finish force PBC 1 " << endl;


					//cout << "rij_vec x,y,z = " << rij_vec[0] << "," << rij_vec[1] << "," << rij_vec[2] << endl;
					rij = sqrt(rij_vec[0] * rij_vec[0] + rij_vec[1] * rij_vec[1] + rij_vec[2] * rij_vec[2]);

					//cout << "rij before forcepbc = " << rij << endl;
					//forcePBC(rij);
					//cout << "rij after forcepbc = " << rij << endl;
					//cout << "finish force PBC 2 " << endl;
					if (rij < rcut) {
						rij2 = sigma2 / (rij * rij);					// sigma^2 / rij^2
						rij6 = rij2 * rij2 * rij2;						// sigma^6 / rij^6
						rij12 = rij6 * rij6;							// sigma^12 / rij^12

						U_LJ = 4 * epsilon * (rij12 - rij6);			// Lennard Jones Potential
						//cout << "U_LJ is " << U_LJ << endl;

						potential[i] = U_LJ;

						E_pot += potential[i];

						//cout <<"the i and j is:"<<i<<","<<j << potential[i]<<"total e_pot "<<E_pot<< endl;
						//E_pot += potential [i][j];
						//cout << "my E_pot after summation is = " << E_pot << endl;


						//cout << "my E_pot = " << E_pot << endl;
						//computing for spatial derivative of potential energy dU/dr 
						dU_LJ = -48 * epsilon / rij * (rij12 - 0.5 * rij6);
						//cout << "dU_LJ is " << dU_LJ << endl;

						//Force exerted on particle #N 
						//cout << "F before force calc is " << F[i][k] << endl;

						//for (k = 0; k < DIM; k++) {
							//F[i][k] += dU_LJ * rij_vec[k];
							//F[j][k] -= dU_LJ * rij_vec[k];

						//}

						F[i][0] += dU_LJ * rij_vec[0];
						F[i][1] += dU_LJ * rij_vec[1];
						F[i][2] += dU_LJ * rij_vec[2];

						F[j][0] -= dU_LJ * rij_vec[0];
						F[j][1] -= dU_LJ * rij_vec[1];
						F[j][2] -= dU_LJ * rij_vec[2];
						//cout << "F after force calc is " << F[i][0]<<","<< F[i][1] << "," << F[i][2] << "," << endl;
					}


					//cout << "i "<<i<<" j "<<j<<" my E_pot after summation is = " << E_pot << endl;
				}

			}
	}

}

// Velocity-Verlet Integrator
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& E_pot, double& E_kin, double& Temp) {

	int i;
	//constants

	double half_dt = 0.5 * dt;					// Half of dt
	double v_half[N][DIM] = { 0 };						// Velocity at half_dt
	//E_kin = 0;
	double kinetic[N];
	E_kin = 0;



	//updating velocity at half time step
	for (i = 0; i < N; i++) {
		//cout << "before verlet v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//cout << "The force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;

		//cout << "Inside the first verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;

		v_half[i][0] = v[i][0] + half_dt * F[i][0] / mass;
		v_half[i][1] = v[i][1] + half_dt * F[i][1] / mass;
		v_half[i][2] = v[i][2] + half_dt * F[i][2] / mass;
	}
	for (i = 0; i < N; i++) {

		//updating position 

		r[i][0] = r[i][0] + v_half[i][0] * dt;
		r[i][1] = r[i][1] + v_half[i][1] * dt;
		r[i][2] = r[i][2] + v_half[i][2] * dt;

		//cout << "after verlet r values (x,y,z) = " << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
	}

	//recalculate force
	Force_Calc(r, F, E_pot);

	for (i = 0; i < N; i++) {
		//updating velocity at t + dt

			//cout << "Recalculated force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;


		v[i][0] = v_half[i][0] + half_dt * F[i][0] / mass;
		v[i][1] = v_half[i][1] + half_dt * F[i][1] / mass;
		v[i][2] = v_half[i][2] + half_dt * F[i][2] / mass;

		//cout << "after the verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}



	//E_kin = 0;
	for (i = 0; i < N; i++) {

		//updating kinetic energy 
		kinetic[i] = ((mass * 0.5) * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]));
		//cout << "The E_kin is : " << kinetic[i] << endl;
		E_kin += kinetic[i];
	}



}

// Velocity-Verlet Integrator for first step (with temperature scaled velocity)
void Vel_Verlet_first(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& E_pot, double& E_kin, double& Temp) {

	int i;
	//constants

	double half_dt = 0.5 * dt;					// Half of dt
	double v_half[N][DIM] = { 0 };						// Velocity at half_dt
	//E_kin = 0;
	double kinetic[N];
	E_kin = 0;

	for (i = 0; i < N; i++) {
		//cout << "The force  before first verlet = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;
	}
	//cout << "starting verlet propagation" << endl;
	//updating velocity at half time step
	for (i = 0; i < N; i++) {
		//cout << "before verlet v values (x,y,z) = " << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	//	cout << "The force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;

		//cout << "Inside the first verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;

		v_half[i][0] = v[i][0] + half_dt * F[i][0] / mass;
		v_half[i][1] = v[i][1] + half_dt * F[i][1] / mass;
		v_half[i][2] = v[i][2] + half_dt * F[i][2] / mass;
	}
	//cout << "v_half is done" << endl;
	for (i = 0; i < N; i++) {

		//updating position 

		r[i][0] = r[i][0] + v_half[i][0] * dt;
		r[i][1] = r[i][1] + v_half[i][1] * dt;
		r[i][2] = r[i][2] + v_half[i][2] * dt;

		//cout << "after verlet r values (x,y,z) = " << r[i][0] << "," << r[i][1] << "," << r[i][2] << endl;
	}
	//cout << "r are updated" << endl;
	//recalculate force
	Force_Calc(r, F, E_pot);

	for (i = 0; i < N; i++) {
		//updating velocity at t + dt

			//cout << "Recalculated force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;


		v[i][0] = v_half[i][0] + half_dt * F[i][0] / mass;
		v[i][1] = v_half[i][1] + half_dt * F[i][1] / mass;
		v[i][2] = v_half[i][2] + half_dt * F[i][2] / mass;

		//cout << "after the verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}
	//cout << "v are updated" << endl;
	Config_Temp(v, Temp);
	for (i = 0; i < N; i++) {
		//kinetic[i] = 0;
		//velocity scaling for adjusting temperature
		//cout << "Temp is "<<Temp<<endl;;
		//cout << "Init_Temp is "<<Init_Temp<<endl;;

		v[i][0] = v[i][0] * sqrt(Init_Temp / Temp);
		v[i][1] = v[i][1] * sqrt(Init_Temp / Temp);
		v[i][2] = v[i][2] * sqrt(Init_Temp / Temp);
		//cout << "after scaling " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;

		//updating kinetic energy 
		kinetic[i] = (mass * 0.5) * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
		//cout << "The E_kin is : " << kinetic[i] << endl;
		E_kin += kinetic[i];

	}
	//cout << "v scaling is done " << endl;
	/*//cout << "Temp is "<<Temp<<endl;;
	if (Temp<(Init_Temp - 0.01) && Temp >(Init_Temp + 0.01)) {
	}*/


}

//void Vel_scaling(double v[N][DIM], double& E_kin, double& Temp) {
void Vel_scaling(double v[N][DIM], double& E_kin, double& Temp_ave) {

	int i;

	//double kinetic[N];
	E_kin = 0;

	for (i = 0; i < N; i++) {
		//kinetic[i] = 0;
		//velocity scaling for adjusting temperature
		//cout << "Temp is "<<Temp<<endl;;
		//cout << "Init_Temp is "<<Init_Temp<<endl;;
		//v[i][0] = v[i][0] * sqrt(Init_Temp / Temp);
		//v[i][1] = v[i][1] * sqrt(Init_Temp / Temp);
		//v[i][2] = v[i][2] * sqrt(Init_Temp / Temp);

		v[i][0] = v[i][0] * sqrt(Init_Temp / Temp_ave);
		v[i][1] = v[i][1] * sqrt(Init_Temp / Temp_ave);
		v[i][2] = v[i][2] * sqrt(Init_Temp / Temp_ave);
		//cout << "after scaling " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;

		//updating kinetic energy 
		//kinetic[i] = (mass * 0.5) * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
		//cout << "The E_kin is : " << kinetic[i] << endl;
		//E_kin += kinetic[i];

	}
}


// Periodic Boundary Condition for force
void forcePBC(double& rij) {
	if (rij > (-0.5) && rij<0)  {
		rij = rij - 0.5;
	}

	if (rij < (0.5) && rij > 0) {
		rij = rij + 0.5;
	}
	

}

void rdf(double r[N][DIM], double g_r[size_bin], double g_ave[size_bin], int& divider) {
	int n, i, j;
	double rij_vec[DIM] = { 0 };
	double rij = 0;
	const double rcut = (0.5 * edge);
	double count = 0;
	double hist[size_bin] = { 0 };


	if (divider == (tmax_prod - startrdf) / 10) {
		otherfile  << 0 << "," << delta_r << "," << 0 << "," << 0 << "," << 0 << "," << divider << "," << 0 << endl;
	}
	for (n = 2; n < size_bin; n++) {
		count = 0;

		for (i = 0; i < N - 1; i++) {
			for (j = i + 1; j < N; j++) {

				rij_vec[0] = r[j][0] - r[i][0];
				rij_vec[1] = r[j][1] - r[i][1];
				rij_vec[2] = r[j][2] - r[i][2];

				while (rij_vec[0] < (-0.5 * edge)) {
					rij_vec[0] = rij_vec[0] + edge;
				}
				while (rij_vec[1] < (-0.5 * edge)) {
					rij_vec[1] = rij_vec[1] + edge;
				}
				while (rij_vec[2] < (-0.5 * edge)) {
					rij_vec[2] = rij_vec[2] + edge;
				}


				while (rij_vec[0] > (0.5 * edge)) {
					rij_vec[0] = rij_vec[0] - edge;
				}
				while (rij_vec[1] > (0.5 * edge)) {
					rij_vec[1] = rij_vec[1] - edge;
				}
				while (rij_vec[2] > (0.5 * edge)) {
					rij_vec[2] = rij_vec[2] - edge;
				}

				rij = sqrt(rij_vec[0] * rij_vec[0] + rij_vec[1] * rij_vec[1] + rij_vec[2] * rij_vec[2]);
				if (rij < rcut) {

					if (rij > ((n - 1.0) * delta_r) && rij < (n * delta_r)) {
						count++;
					}
					else {
						continue;

					}

				}
			}
		}
		hist[n] += count;

		g_r[n] = (hist[n] * (Volume)) / (N * 4.0 * PI * (((n - 1.0) * delta_r) * ((n - 1.0) * delta_r)) * delta_r * N);
		g_ave[n] += g_r[n];
		//g_ave[n] += (g_ave[n] / divider);
		if (divider == (tmax_prod - startrdf) / 10) {
		//if (divider > 3000) {
			otherfile << ((n - 1.0) * delta_r) / sigma << "," << n * delta_r << "," << count << "," << g_r[n] << "," << g_ave[n] << "," << divider << "," << g_ave[n] / divider << endl;
		}
		//}

		//otherfile << "," << ((n - 1.0) * delta_r) / sigma << "," << n * delta_r << "," << count << "," << g_r[n] << endl;
	}
}


