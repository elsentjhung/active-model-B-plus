/* mathematical functions */

#include <math.h>

int iupa[Nx], idwna[Nx];
int jupa[Ny], jdwna[Ny];

// array iup = [1, 2, 3, 4, ..., Nx-1, 0]
// array iup[iup] = [2, 3, 4, ..., Nx-1, 0, 1]
void init_math(void) {
    int i, j;

    for (i = 0; i < Nx; i++) {
        if (i == Nx - 1) { iupa[i] = 0; } 
        else { iupa[i] = i + 1; }

        if (i == 0) { idwna[i] = Nx - 1; } 
        else { idwna[i] = i - 1; }
    }
    for (j = 0; j < Ny; j++) {
        if (j == Ny - 1) { jupa[j] = 0; } 
        else { jupa[j] = j + 1; }

        if (j == 0) { jdwna[j] = Ny - 1; } 
        else { jdwna[j] = j - 1; }
    }
}

// get random number from normal distribution
double gaussian_rand(void) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;
  
    if (iset == 0) {
        do {
            v1 = 2.0*drand48()-1.0;
            v2 = 2.0*drand48()-1.0;
            rsq = v1*v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0*log(rsq)/rsq);
    
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset = 0;
        return gset;
    }
}

// second order numerical derivative
double diff_x(double phi[Nx][Ny], int i, int j) {
    double dxphi;
    int iup = iupa[i];
    int idwn = idwna[i];
    int jup = jupa[j];
    int jdwn = jdwna[j];

    dxphi = ((phi[iup][jup]  - phi[idwn][jup])*0.1 + 
	         (phi[iup][j]    - phi[idwn][j])*0.3   +
	         (phi[iup][jdwn] - phi[idwn][jdwn])*0.1)/dx;

    return dxphi;
}
double diff_y(double phi[Nx][Ny], int i, int j) {
    double dyphi;
    int iup = iupa[i];
    int idwn = idwna[i];
    int jup = jupa[j];
    int jdwn = jdwna[j];

    dyphi = ((phi[iup][jup]  - phi[iup][jdwn])*0.1 + 
	         (phi[i][jup]    - phi[i][jdwn])*0.3   +
	         (phi[idwn][jup] - phi[idwn][jdwn])*0.1)/dy;

    return dyphi;
}
double laplacian(double phi[Nx][Ny], int i, int j) {
    double laplacianphi;
    int iup = iupa[i];
    int idwn = idwna[i];
    int jup = jupa[j];
    int jdwn = jdwna[j];

    laplacianphi = ((-0.5*phi[idwn][jup]  + 2.0*phi[i][jup]  - 0.5*phi[iup][jup]) + 
	                ( 2.0*phi[idwn][j]    - 6.0*phi[i][j]    + 2.0*phi[iup][j])   +
	                (-0.5*phi[idwn][jdwn] + 2.0*phi[i][jdwn] - 0.5*phi[iup][jdwn]))/(dx*dy);

    return laplacianphi;
}

// eighth order numerical derivative
double diff_x_8(double phi[Nx][Ny], int i, int j) {
    double dxphi;
    int iup = iupa[i];
    int iup2 = iupa[iupa[i]];
    int iup3 = iupa[iupa[iupa[i]]];
    int iup4 = iupa[iupa[iupa[iupa[i]]]];
    int idwn = idwna[i];
    int idwn2 = idwna[idwna[i]];
    int idwn3 = idwna[idwna[idwna[i]]];
    int idwn4 = idwna[idwna[idwna[idwna[i]]]];

    dxphi = (-phi[iup4][j]/280.0  + 4.0*phi[iup3][j]/105.0  - phi[iup2][j]/5.0  + 4.0*phi[iup][j]/5.0  + 
              phi[idwn4][j]/280.0 - 4.0*phi[idwn3][j]/105.0 + phi[idwn2][j]/5.0 - 4.0*phi[idwn][j]/5.0)/dx;

    return dxphi;
}
double diff_y_8(double phi[Nx][Ny], int i, int j) {
    double dyphi;
    int jup = jupa[j];
    int jup2 = jupa[jupa[j]];
    int jup3 = jupa[jupa[jupa[j]]];
    int jup4 = jupa[jupa[jupa[jupa[j]]]];
    int jdwn = jdwna[j];
    int jdwn2 = jdwna[jdwna[j]];
    int jdwn3 = jdwna[jdwna[jdwna[j]]];
    int jdwn4 = jdwna[jdwna[jdwna[jdwna[j]]]];

    dyphi = (-phi[i][jup4]/280.0  + 4.0*phi[i][jup3]/105.0  - phi[i][jup2]/5.0  + 4.0*phi[i][jup]/5.0  + 
              phi[i][jdwn4]/280.0 - 4.0*phi[i][jdwn3]/105.0 + phi[i][jdwn2]/5.0 - 4.0*phi[i][jdwn]/5.0)/dy;

    return dyphi;
}
