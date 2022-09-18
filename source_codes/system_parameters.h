/* fixed system parameters */

#define pi 3.14159265359

#define Nx 128  // number of lattice points
#define Ny 128

#define dx 1.0  // space and time discretization
#define dy 1.0
#define dt 0.01  

#define Nt     1000000 // total number of timesteps
#define nt_int    1000 // timestep interval for checkpoint
#define nt_init      0 // initial timestep

#define A 0.25  // free energy parameters -0.5*A*phi^2 + 0.25*A*phi^4
#define K 1.0  // coefficient of |grad phi|^2
#define D 0.2  // noise strength

double phi_init = 0.6; // initial density
double zeta = 4.0;  // activities
double lambda = 1.0;

double Lx = Nx*dx;
double Ly = Ny*dx;

