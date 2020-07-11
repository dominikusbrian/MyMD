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
const int tmax = 50000;						// maximum number of timestep or simulation duration
const double rho = 0.8;						// System's volume density
double dt = 0.001;							// Simulation timestep size
double Init_Temp = 1;						// Temperature
const double Volume = N / rho;
double edge = cbrt(Volume);					// Box edge size
const double delta_r = 0.2;					// bin width for rdf analysis
const int size_pair = (N * (N - 1)) / 2;	// array size for the total number of pairs 
const int size_bin = 52;



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
ofstream otherfile;
ofstream trajfile;
ofstream cvvfile;
ofstream cvvfile2;

/////////////////////////////*Function declaration:*/////////////////////////////////

// Initialize Particle Position
void Config_atoms_lattice(double r[N][DIM]);

// Initialize Particle Velocity
void Config_velocity(double v[N][DIM]);

// Scale initial velocity based on Init_Temp
double Config_Temp(double v[N][DIM], double& Temp);

// Force Calculation
void Force_Calc(double r[N][DIM], double F[N][DIM], double& E_pot);

// Velocity-Verlet propagator
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& E_pot, double& E_kin, double& Temp);
void Vel_Verlet_first(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& E_pot, double& E_kin, double& Temp);
//Calculate Kinetic Energy
void E_kin_Calc(double v[N][DIM], double& E_kin);


//Periodic Boundary Condition
double PBC(double x1, double x2);
void forcePBC(double& rij);

//radial distribution function  
void rdf(double r[N][DIM], double g_r[size_bin], double g_ave[size_bin], int& divider);

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
	int write_traj = 250;
	double rij = 0;					// pair distance
	double Temp = Init_Temp;
	double g_r[size_bin] = { 0 };			// Instaneous rdf 
	double g_ave[size_bin] = { 0 };	// Accumulated average of rdf for a given frame
	const double t_corr = 10;
	double Delta_t = 1;
	double duration = t_corr / Delta_t;
	const int storage = 15;
	double vx[storage][N] = { 0 }, vy[storage][N] = { 0 }, vz[storage][N] = { 0 };
	double cvv_dt[N] = { 0 };
	double cvv_accum = 0;
	double cvv_ave[storage] = { 0 };
	int k = 0;
	int startcvv = 0;


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
	Config_atoms_lattice(r);
	//Initialize Velocity
	Config_velocity(v);

	//Initialize Force
	Force_Calc(r, F, E_pot);

	//Export initial configuration 
	myfile.open("config.csv");
	myfile << "N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;
	for (i = 0; i < N; i++) {
		myfile << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}
	myfile.close();

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
	otherfile.open("rdf_t.csv");
	otherfile << "t" << "," << "r" << "," << "r_plus_dr" << "," << "Count" << "," << "g_r" << endl;
	trajfile.open("traj.csv");
	trajfile << "t" << "N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;
	cvvfile.open("cvv.csv");
	cvvfile << "t" << "," << "cvv_accum" << endl;
	cvvfile2.open("cvv2.csv");
	cvvfile2 << "t" << "," << "cvv_ave" << endl;
	
	for (t = 0; t < tmax; t++) {
		//E_tot = 0;
		//cout << "Temp is "<<Temp<<endl;;
		//cout << "the t is : " << t << endl;
		//if (t==0) {
			//Config_Temp(v, Temp);
			//Vel_Verlet_first(r, v, F, E_pot, E_kin, Temp);
			//E_tot = E_kin + E_pot;
			//Temp = ((E_kin / N) * 2 / (3 * kB));

			//myfile << t + 1 << "," << E_pot << "," << E_kin << "," << E_tot << "," << Temp << endl;
		
		//}
		E_kin = 0;
		Vel_Verlet(r, v, F, E_pot, E_kin, Temp);
		
		E_tot = E_kin + E_pot;

		Temp = ((E_kin / N) * 2 / (3 * kB));

		myfile << t+1  << "," << E_pot << "," << E_kin << "," << E_tot << "," << Temp << endl;
		
		

		// Radial distribution function (RDF) analysis
		if (t % 10 == 0) {
			divider++;
			rdf(r, g_r, g_ave, divider);
		}



		//Velocity auto-correlation function (calculating for Cvv)


		//if (t > startcvv) {
		k = t - (t_corr * counter);								// index for storage labeling always 1-10 

			//cout << "my k is = " << k << endl;
			//cout << "my counter is " << counter << endl;
			//cout << "my t_corr*counter is " << t_corr * counter << endl;
		for (i = 0; i < 10; i++) {

			// read the value of v for all atoms
			vx[k][i] = v[i][0];
			vy[k][i] = v[i][1];
			vz[k][i] = v[i][2];
			//cout << "v = " << v[i][0] << endl;
			//cout << "vx = " << vx[k][i] << endl;
			//cout << "vy = " << vy[k][i] << endl;
			//cout << "vz = " << vz[k][i] << endl;

		}


		// a function that call the appropriate velocity for a given time
		//if (t=0) cvv_accum = 

		if (t > (t_corr * counter) + 1) {   //Cvv calculation will be skip for all t=0 and t=1 (the same applies fort =  10,11,20,21, and so on ) 

			if (k > 0 && k <= 10) {
				cvv_ave[k] = 0;
				for (i = 0; i < N; i++) {
					cvv_dt[i] = ((vx[k][i] * vx[k + 1][i] + vy[k][i] * vx[k + 1][i] + vz[k][i] * vx[k + 1][i]) / N);
					cvv_ave[k] += (cvv_dt[i]) / (t - t_corr);
					//cvv_ave[k] = cvv_ave[k] * (1 / (t - t_corr));
					//cout << t << "," << cvv_dt[i] << "," << cvv_ave[k] << "," << cvv_accum << endl;
					//cout << "cvv execution is fine" << endl;
				}
				cvv_accum += (cvv_ave[k]);
				//cout << t << "," << k << "," << cvv_ave[k] << "," << cvv_accum << endl;
				if(t>(startcvv-t_corr) && k>2 && k<10) cvvfile2 << t << "," << cvv_ave[k] << endl;
				//cout << "if k is fine " << endl;
			}
			//cout << t << "," << cvv_accum << endl;

			//if (t>=startcvv+ 20) cvvfile << t << "," << cvv_accum << endl;

		}
		//if (t >= startcvv + 20) cvvfile << t << "," << cvv_accum << endl;
		if (k >= 10) {
			if (t >= startcvv + 20) cvvfile << t << "," << cvv_accum << endl;
			k = 0;
			counter++;

			cvv_accum = 0;
		}

		

		// record trajectory every write_traj steps
		if (t % write_traj == 0) {
			for (i = 0; i < N; i++) {
				trajfile << t + 1 << "," << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
			}
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
	otherfile.close();
	trajfile.close();
	cvvfile.close();

	cout << "------------------PROPAGATION FINISHED ----------------------" << endl;



	//Write Final Configuration to final_config.csv
	myfile.open("finalconfig.csv");
	myfile << "N" << "," << "rx" << "," << "ry" << "," << "rz" << "," << "vx" << "," << "vy" << "," << "vz" << endl;
	for (i = 0; i < N; i++) {
		myfile << i << "," << r[i][0] << "," << r[i][1] << "," << r[i][2] << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}
	myfile.close();

	cout << "------------------END OF SIMULATION ----------------------" << endl;
	printf("Total simulation time %.2fs\n", ((double)clock() - tstart) / CLOCKS_PER_SEC);


	return 0;
}



/*------------------------Main function end here------------------------ */

/*------------------------Function definition------------------------- */

// Initialize position

void Config_atoms_lattice(double r[N][DIM]) {
	 
	int i;
	double ix, iy, iz;


	ix = iy = iz = 0;
	// assigning particle positions

	for (i = 0; i < N; i++) {

		r[i][0] = ix * 0.2* edge;
		r[i][1] = iy * 0.2* edge;
		r[i][2] = iz * 0.2* edge;
		//r[i][0] = ix * edge;
		//r[i][1] = iy * edge;
		//r[i][2] = iz * edge;

		ix++;
		if (ix >= edge) {
			ix = 0;
			iy++;
			if (iy >= edge) {
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
	double v_mag = sqrt(3 *Init_Temp * kB / mass);

	for (i = 0; i < DIM; i++) { v_com[i] = 0; }
	// the random number generator generate random number between 1-2
	for (i = 0; i < N; i++) {
		//cout << "v_mag is = " << v_mag << endl;

		//v[i][0] = (v_mag * e()/10000000000);
		//v[i][1] = (v_mag * e()/10000000000);
		//v[i][2] = (v_mag * e()/10000000000);

		v[i][0] = (v_mag * (((double)rand() * (10 - 0.1)) / (double)RAND_MAX));
		v[i][1] = (v_mag * (((double)rand() * (10 - 0.1)) / (double)RAND_MAX));
		v[i][2] = (v_mag * (((double)rand() * (10 - 0.1)) / (double)RAND_MAX));

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
	n = 0;
	for (i = 0; i < N - 1; i++) {

		for (j = i + 1; j < N; j++) {

			//for dimension, where 0 = x,  1 = y = 1, 2 = z


			rij_vec[0] = (r[j][0] - r[i][0]);
			rij_vec[1] = (r[j][1] - r[i][1]);
			rij_vec[2] = (r[j][2] - r[i][2]);

			//rij_vec[0] = PBC(r[j][0] , r[i][0]);
			//rij_vec[1] = PBC(r[j][1] , r[i][1]);
			//rij_vec[2] = PBC(r[j][2] , r[i][2]);

			//rij_vec[k] = (r[j][k] - r[i][k]);


			//rij_vec[k] = PBC(r[j][k], r[i][k]);

			//cout << "rij_Vec[k] before pbc " << rij_vec[k] << endl;
			//if (rij_vec[k] < (-0.5*edge ))  rij_vec[k] = rij_vec[k] + edge;
			//if (rij_vec[k] > (0.5*edge )) rij_vec[k] = rij_vec[k] - edge ;
			//else {
			//	continue;
			//}
			//cout << "rij_Vec[k] after pbc " << rij_vec[k] << endl;


		//cout << "rij before forcepbc = " << rij << endl;
		//if (rij < (-0.5 * edge))  rij = rij + edge;
		//if (rij > (0.5 * edge)) rij = rij - edge ;
		//else {
		//	continue;
		//}


		//cout << "rij_vec x,y,z = " << rij_vec[0] << "," << rij_vec[1] << "," << rij_vec[2] << endl;
			rij = sqrt(rij_vec[0] * rij_vec[0] + rij_vec[1] * rij_vec[1] + rij_vec[2] * rij_vec[2]);

			//cout << "rij before forcepbc = " << rij << endl;
			forcePBC(rij);
			//cout << "rij after forcepbc = " << rij << endl;

			//if (rij < rcut) {
				rij2 = sigma2 / (rij*rij);					// sigma^2 / rij^2
				rij6 = rij2 * rij2 * rij2;						// sigma^6 / rij^6
				rij12 = rij6 * rij6;							// sigma^12 / rij^12

				U_LJ = 4 * epsilon * (rij12 - rij6);			// Lennard Jones Potential
				//cout << "U_LJ is " << U_LJ << endl;

				potential[i] = U_LJ;

				E_pot += potential[i];
				//n++;
				//cout <<"the i and j is:"<<i<<","<<j<<"my n "<<n<< " and my potential n " << potential[i]<<"total e_pot "<<E_pot<< endl;
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
				//cout << "F after force calc is " << F[i][k] << endl;
			//}


			//cout << "i "<<i<<" j "<<j<<" my E_pot after summation is = " << E_pot << endl;
		}

	}

}

// Velocity-Verlet Integrator
void Vel_Verlet(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& U_LJ, double& E_kin, double& Temp) {

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
	Force_Calc(r, F, U_LJ);

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
		
		//kinetic[i] = 0;
		//velocity scaling for adjusting temperature
		//cout << "Temp is "<<Temp<<endl;;
		//cout << "Init_Temp is "<<Init_Temp<<endl;;
		//v[i][0] = v[i][0] * sqrt(Init_Temp / Temp);
		//v[i][1] = v[i][1] * sqrt(Init_Temp / Temp);
		//v[i][2] = v[i][2] * sqrt(Init_Temp / Temp);
		//cout << "after scaling " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
		//updating kinetic energy 
		kinetic[i] = ((mass * 0.5) * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]));
		//cout << "The E_kin is : " << kinetic[i] << endl;
		E_kin += kinetic[i];
	}
	
	

}

// Velocity-Verlet Integrator for first step (with temperature scaled velocity)
void Vel_Verlet_first(double r[N][DIM], double v[N][DIM], double F[N][DIM], double& U_LJ, double& E_kin, double& Temp) {

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
	Force_Calc(r, F, U_LJ);

	for (i = 0; i < N; i++) {
		//updating velocity at t + dt

			//cout << "Recalculated force is = " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;

		
		v[i][0] = v_half[i][0] + half_dt * F[i][0] / mass;
		v[i][1] = v_half[i][1] + half_dt * F[i][1] / mass;
		v[i][2] = v_half[i][2] + half_dt * F[i][2] / mass;

		//cout << "after the verlet loop " << "," << v[i][0] << "," << v[i][1] << "," << v[i][2] << endl;
	}
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
	/*//cout << "Temp is "<<Temp<<endl;;
	if (Temp<(Init_Temp - 0.01) && Temp >(Init_Temp + 0.01)) {

	}*/


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
	if (x1 - x2 < (-0.5 * edge)) {
		return (x1 - x2) + edge;
	}

	if (x1 - x2 > (0.5 * edge)) {
		return (x1 - x2) - edge;
	}
	else {
		return x1 - x2;
	}
	return x1, x2;
}

// Periodic Boundary Condition for force
void forcePBC(double& rij) {
	if (rij < (-0.5 * edge)) {
		rij = rij + edge;
	}

	if (sqrt(rij) > (0.5 * edge)) {
		rij = rij - edge;
	}
	else {
		rij = rij;
	}
}

void rdf(double r[N][DIM], double g_r[size_bin], double g_ave[size_bin], int& divider) {
	int n, i, j;
	double rij_vec[DIM] = { 0 };
	double rij = 0;
	const double rcut = (0.5 * edge);
	double count = 0;
	//double g_r[52] = { 0 };
	//double g_ave[size_bin] = { 0 };
	double hist[size_bin] = { 0 };
	//int divider = 0;
	//divider++;
	//double Nframe = 0;

	otherfile << "," << 0 << "," << delta_r << "," << 0 << "," << 0 << "," << 0 << "," << divider << "," << 0 << endl;
	for (n = 2; n < size_bin; n++) {
		count = 0;
		
		for (i = 0; i < N - 1; i++) {
			for (j = i + 1; j < N; j++) {

				rij_vec[0] = r[j][0] - r[i][0];
				rij_vec[1] = r[j][1] - r[i][1];
				rij_vec[2] = r[j][2] - r[i][2];
				rij = sqrt(rij_vec[0] * rij_vec[0] + rij_vec[1] * rij_vec[1] + rij_vec[2] * rij_vec[2]);
				//if (rij < rcut) {

					if (rij > ((n - 1.0) * delta_r) && rij < (n * delta_r)) {
						count++;
					}
					else {
						continue;

					}

				///}
			}
		}
		hist[n] += count;
		
		g_r[n] = (hist[n] * (Volume)) / (N * 4.0 * PI * (((n - 1.0) * delta_r) * ((n - 1.0) * delta_r)) * delta_r * N );
		g_ave[n] += g_r[n];
		//g_ave[n] += (g_ave[n] / divider);

		otherfile << "," << ((n - 1.0) * delta_r) / sigma << "," << n * delta_r << "," << count << "," << g_r[n] << "," << g_ave[n] << "," << divider << "," << g_ave[n] / divider << endl;
		//otherfile << "," << ((n - 1.0) * delta_r) / sigma << "," << n * delta_r << "," << count << "," << g_r[n] << endl;
	}
}