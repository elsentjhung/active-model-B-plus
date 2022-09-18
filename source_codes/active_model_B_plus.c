/* active model B+ */

/* preprocessor */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "system_parameters.h"
#include "mathematical_functions.h"

/* fields */
double phi[Nx][Ny];  // density field 
double dphi[Nx][Ny];  // changes in phi

double Jx[Nx][Ny], Jy[Nx][Ny];  // total current
double f[Nx][Ny];  // free energy density

/* functions */
void initialize(int ndata);
void calculate_dphi(void);
void update_phi(void);

/* input and output */
FILE *input;
char input_filename[50];
FILE *output1;
char output_filename1[50];
FILE *output2;
char output_filename2[50];
void save_data(int ndata, int nt);

int main(int argc, char** argv) {
    int nt, nt_exp;  // timestep
    if (nt_init == 0) { nt_exp = 1; } 
    else { nt_exp = nt_init; }
    int ndata = atoi(argv[1]);

    srand48((unsigned int) time(NULL) + 100*ndata);
    initialize(ndata);

    for (nt = 0; nt < Nt; nt++) {
#if SAVE_EXP
        if (nt == nt_exp) {
#endif

#if SAVE_SEQ
        if (nt % nt_int == 0) {
#endif
            printf("t = %f \n", nt*dt);
            save_data(ndata, nt);
            nt_exp *= 2;
        }

        calculate_dphi();
        update_phi();                
    }
    return 0;
}

/* algorithm */
void update_phi(void) {
    int i, j;

    for (i = 0; i < Nx; i++) {
#if PERIODIC_BC
        for (j = 0; j < Ny; j++) {
	        phi[i][j] += dt*dphi[i][j];
    	}
#endif

#if REFLECTING_BC
        for (j = 4; j < Ny-4; j++) {
            phi[i][j] += dt*dphi[i][j];
        }

        phi[i][0] = phi[i][4];
        phi[i][1] = phi[i][4];
        phi[i][2] = phi[i][4];
        phi[i][3] = phi[i][4];
        phi[i][Ny-1] = phi[i][Ny-5];
        phi[i][Ny-2] = phi[i][Ny-5];
        phi[i][Ny-3] = phi[i][Ny-5];
        phi[i][Ny-4] = phi[i][Ny-5];
#endif
    }
}
void calculate_dphi(void) {
    int i, j;
    double dphidx[Nx][Ny], dphidy[Nx][Ny];
    double laplacianphi[Nx][Ny];
    double gradphisq[Nx][Ny];
    double mu[Nx][Ny];  // chemical potential
    double Lambdax[Nx][Ny];  // noise current field
    double Lambday[Nx][Ny];  
    double divJ;

    // dphidx, dphidy
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
   	        dphidx[i][j] = diff_x(phi,i,j);
   	        dphidy[i][j] = diff_y(phi,i,j);
	    }
    }

    for (i = 0; i < Nx; i++) {
#if PERIODIC_BC
        for (j = 0; j < Ny; j++) {
#endif

#if REFLECTING_BC
        for (j = 4; j < Ny-4; j++) {
#endif
	        // laplacianphi, gradphisq
	        laplacianphi[i][j] = laplacian(phi,i,j); 
	        gradphisq[i][j] = dphidx[i][j]*dphidx[i][j] + dphidy[i][j]*dphidy[i][j];

	        // chemical potential (last term is AMB chemical potential)
	        mu[i][j] = -A*phi[i][j] + A*phi[i][j]*phi[i][j]*phi[i][j] - K*laplacianphi[i][j] + lambda*gradphisq[i][j];

	        // free energy
	        f[i][j] = -0.5*A*phi[i][j]*phi[i][j] + 0.25*A*phi[i][j]*phi[i][j]*phi[i][j]*phi[i][j] + 0.5*K*gradphisq[i][j];

	        // noise current
	        Lambdax[i][j] = sqrt(2.0*D/(dx*dy*dt))*gaussian_rand();
	        Lambday[i][j] = sqrt(2.0*D/(dx*dy*dt))*gaussian_rand();
	    }
    }

    // current
    for (i = 0; i < Nx; i++) {
#if PERIODIC_BC
        for (j = 0; j < Ny; j++) {
            // total current (last term is AMB+ current)
	        Jx[i][j] = -diff_x1(mu,i,j) + Lambdax[i][j] + zeta*laplacianphi[i][j]*dphidx[i][j];
	        Jy[i][j] = -diff_y1(mu,i,j) + Lambday[i][j] + zeta*laplacianphi[i][j]*dphidy[i][j];
        }
#endif

#if REFLECTING_BC
        for (j = 4; j < Ny-4; j++) {
            // total current (last term is AMB+ current)
	        Jx[i][j] = -diff_x1(mu,i,j) + Lambdax[i][j] + zeta*laplacianphi[i][j]*dphidx[i][j];
	        Jy[i][j] = -diff_y1(mu,i,j) + Lambday[i][j] + zeta*laplacianphi[i][j]*dphidy[i][j];
        }
        Jy[i][0] = 0.0;  // J.n = 0 at j = 4 and j = Ny-5 
        Jy[i][1] = 0.0;
        Jy[i][2] = 0.0;
        Jy[i][3] = 0.0;
        Jy[i][4] = 0.0;
        Jy[i][Ny-1] = 0.0;
        Jy[i][Ny-2] = 0.0;
        Jy[i][Ny-3] = 0.0;
        Jy[i][Ny-4] = 0.0;
        Jy[i][Ny-5] = 0.0;
#endif
    }

    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
	        divJ = diff_x1(Jx,i,j) + diff_y1(Jy,i,j);
  	        dphi[i][j] = -divJ;
	    }
    }
}

/* initialization */
void initialize(int ndata) {
    int i, j;

    init_math();

    // homogenous phi
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            if (j < Ny/2) {
                phi[i][j] = phip;
            } else {
                phi[i][j] = phim;
            }
	    }
    }

    // read data from saved file
    if (nt_init != 0) {
        int i_new, j_new;
        double Jx_new, Jy_new, phi_new;

        sprintf(input_filename, "../data%d/data.txt.%d", ndata, nt_init);
        input = fopen(input_filename, "r");
        for (i = 0; i < Nx; i++) {
            for (j = 0; j < Ny; j++) {
                fscanf(input, "%d %d %lf %lf %lf", &i_new, &j_new, &Jx_new, &Jy_new, &phi_new);
	            phi[i][j] = phi_new;
                Jx[i][j] = Jx_new;
                Jy[i][j] = Jy_new;
	         }
        }
        fclose(input);
    }
}

/* save data at timestep nt to a file */
void save_data(int ndata, int nt) {
    int i, j;
    
    sprintf(output_filename1, "./data%d/data.txt", ndata);
    sprintf(output_filename2, "./data%d/data.txt.%d", ndata, nt);
    output1 = fopen(output_filename1, "w");
    output2 = fopen(output_filename2, "w");
    for (i = 0; i < Nx; i++) {      
        for (j = 0; j < Ny; j++) {
            fprintf(output1, "%d %d %e %e %e \n", 
		        i, j, Jx[i][j], Jy[i][j], phi[i][j]);
	        fprintf(output2, "%d %d %e %e %e \n", 
		        i, j, Jx[i][j], Jy[i][j], phi[i][j]);
	    }
	    fprintf(output1, "\n");
	    fprintf(output2, "\n");
    }
    fclose(output1);
    fclose(output2);    
}

