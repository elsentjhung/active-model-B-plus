/* fixed system parameters */

#define pi 3.14159265359

#define Nx 128  // number of lattice points
#define Ny 128

#define dx 1.0  // space and time discretization
#define dy 1.0
#define dt 0.01  

#define Ndata  1000  // total number of data sets

#define Nt     17000000 // total number of timesteps
#define nt_int    10000 // timestep interval for checkpoint
#define nt_init      0  // initial timestep

#define A 0.25  // free energy parameters -0.5*A*phi^2 + 0.25*A*phi^4
#define K 1.0  // coefficient of |grad phi|^2
#define D 0.05  // noise strength

double phi0;  // global density
double zeta = -3.0;  // activities
double lambda = -1.5;
double phip = 1.0;  // coexisting densities (depend on zeta and lambda)
double phim = -1.0;

#define REFLECTING_BC 1  // reflecting boundary conditions along y
#define PERIODIC_BC 0  // periodic boundaryt conditions along y

double Lx = Nx*dx;
double Ly = Ny*dx;

#define SAVE_EXP 1  // save data exponentially at nt = 1, 2, 4, 8, 16, ...
#define SAVE_SEQ 0  // save data sequentially at nt = nt_int, 2*nt_int, 3*nt_int, ...
