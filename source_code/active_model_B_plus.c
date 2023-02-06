/* active model B+ */

/* preprocessor */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "system_parameters.h"
#include "mathematical_functions.h"

/* fields */
double phi0;  // global density
double phi[Nx][Ny];  // density field 
double dphi[Nx][Ny];  // changes in phi

double Jx[Nx][Ny], Jy[Nx][Ny];  // total current
double f[Nx][Ny];  // free energy density

/* functions */
void initialize(void);
void calculate_dphi(void);
void update_phi(void);

/* input and output */
FILE *input;
char input_filename[50];
FILE *output1;
char output_filename1[50];
FILE *output2;
char output_filename2[50];
void save_data(int nt);

int main(int argc, char** argv) {
    int nt;  // timestep

    srand48((unsigned int) time(NULL));
    initialize();

    for (nt = 0; nt < Nt; nt++) {
        if (nt % nt_int == 0) {
            printf("t = %f \n", nt*dt);
            save_data(nt);
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
        for (j = 0; j < Ny; j++) {
	        phi[i][j] += dt*dphi[i][j];
    	}
    }
}
void calculate_dphi(void) {
    int i, j;
    double dphidx[Nx][Ny], dphidy[Nx][Ny];
    double laplacianphi[Nx][Ny];
    double gradphisq[Nx][Ny];
    double mueq[Nx][Ny];  // equilibrium chemical potential
    double muact[Nx][Ny];  // active chemical potential
    double Jxeq[Nx][Ny], Jyeq[Nx][Ny];  // equilibrium currents
    double Jxact[Nx][Ny], Jyact[Nx][Ny];  // active currents
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
        for (j = 0; j < Ny; j++) {
	        // laplacianphi, gradphisq
	        laplacianphi[i][j] = laplacian(phi,i,j); 
	        gradphisq[i][j] = dphidx[i][j]*dphidx[i][j] + dphidy[i][j]*dphidy[i][j];

	        // chemical potential
	        mueq[i][j] = -A*phi[i][j] + A*phi[i][j]*phi[i][j]*phi[i][j] - K*laplacianphi[i][j];
            muact[i][j] = lambda*gradphisq[i][j];

	        // free energy
	        f[i][j] = -0.5*A*phi[i][j]*phi[i][j] + 0.25*A*phi[i][j]*phi[i][j]*phi[i][j]*phi[i][j] + 0.5*K*gradphisq[i][j];

	        // noise current
	        Lambdax[i][j] = sqrt(2.0*D/(dx*dy*dt))*gaussian_rand();
	        Lambday[i][j] = sqrt(2.0*D/(dx*dy*dt))*gaussian_rand();
	    }
    }

    // current
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            // equilibrium currents
	        Jxeq[i][j] = -diff_x_8(mueq,i,j); 
	        Jyeq[i][j] = -diff_y_8(mueq,i,j);

            // active currents
            Jxact[i][j] = zeta*laplacianphi[i][j]*dphidx[i][j] - diff_x(muact,i,j);
            Jyact[i][j] = zeta*laplacianphi[i][j]*dphidy[i][j] - diff_y(muact,i,j);

            // total current
            Jx[i][j] = Jxeq[i][j] + Jxact[i][j] + Lambdax[i][j];
            Jy[i][j] = Jyeq[i][j] + Jyact[i][j] + Lambday[i][j];
        }
    }
    
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
	        divJ = diff_x_8(Jxeq,i,j) + diff_y_8(Jyeq,i,j) 
                 + diff_x_8(Lambdax,i,j) + diff_y_8(Lambday,i,j)
                 + diff_x(Jxact,i,j) + diff_y(Jyact,i,j);
  	        dphi[i][j] = -divJ;
	    }
    }
}

/* initialization */
void initialize(void) {
    int i, j;

    init_math();

    // homogenous phi
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            phi[i][j] = phi_init;
	    }
    }

    // read data from saved file
    if (nt_init != 0) {
        int i_new, j_new;
        double Jx_new, Jy_new, phi_new;

        sprintf(input_filename, "../data.txt.%d", nt_init);
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
void save_data(int nt) {
    int i, j;
    
    sprintf(output_filename1, "./data.txt");
    sprintf(output_filename2, "./data.txt.%d", nt);
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

