/****************************************************************************
    
    Simulation of the hydration layer on an icosahedral virus
    Copyright (C) 2022  Cesar Lopez Pastrana
    
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
****************************************************************************/

#include "icovirus.h"
#include "allvec.h"
#include "rands.h"


int main()
{
    double E0, E;           // Energy before and after spin flip
    double **mean_lattice;     // Mean lattice configuration after equilibration
    double *energymcs;      // Energy per monte carlo step
    long *idum; 
    long foornd; 
    
    FILE *fid_lattice;
    FILE *fid_energymc;
    
    foornd = -time(NULL);
    idum = &foornd;
    
    // PRELIMINARY STEPS   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    
    // Generate the lattice and memory allocation
    struct Lattice lattice;
    struct Virus virus;
    struct Params params;
    
    
    //Assign parameters
    virus.Rc = VIR_R*10.0; // A->nm
    virus.t =  VIR_T*10.0;
    virus.theta0 = THETA0;
    
    params.eps = EPSNN;       
    params.mu = MU;       
    params.bsurf = BSURF;
    params.beta = BETA;      
    
    
    // Construct lattice and draw the virus
    gen_init(&lattice, &virus);
    
    
    // Fill the empty space randomly with water
    for(int i=0; i<lattice.nh; i++){
        for(int j=0; j<lattice.nw; j++){
            if(lattice.lattice[i][j]!=2){
                 lattice.lattice[i][j] = 1;
                if(ran1(idum)>0.5)
                    lattice.lattice[i][j] = 1;
                else
                    lattice.lattice[i][j] = -1;
            }
        }
    }
    
    // Allocate memmory for the mean lattice and fill it with zeros
    mean_lattice = dmatrix(lattice.nh, lattice.nw);
    for(int i=0;i<lattice.nh; i++) {
        for(int j=0; j<lattice.nw; j++)
            mean_lattice[i][j]=0;
    }
    
    
    // Save the initial lattice
    fid_lattice = fopen(FNAME_LATTICE_0, "w");
    for(int i=0; i<lattice.nh; i++){
        for(int j=0; j<lattice.nw; j++){
            fprintf(fid_lattice,"%d\t", lattice.lattice[i][j]);
        }
        fprintf(fid_lattice,"\n");
    }
    fclose(fid_lattice);
    
    
    // Initial energy configuration
    E=0;
    energymcs = malloc(MAX_MCS*sizeof(int *));
    for(int i=0; i<lattice.nh; i++) {
        for(int j=0; j<lattice.nw; j++) {
            if(lattice.lattice[i][j] != 2)
                E += locenergy(&lattice, &params, i, j);
        }
    }
    energymcs[0] = E;
    
    //--------------------------------------------------------------------------------/
    
    // MAIN MONTE CARLO ROUTINE   ++++++++++++++++++++++++++++++++++++++++++++++++++++/
    for(int mcs=1; mcs<MAX_MCS; mcs++)
    {
        for(int i=0; i<lattice.nh; i++){ 
            for(int j=0; j<lattice.nw; j++)
            {
                int x,y;
                x = (int)(lattice.nh)*ran1(idum);
                y = (int)(lattice.nw)*ran1(idum);
                if(lattice.lattice[x][y]!=2) { // Exclude the surfaces
                    
                    E0 = locenergy(&lattice, &params, x, y);
                    lattice.lattice[x][y] *= -1;
                    E = locenergy(&lattice, &params, x, y);

                    if(ran1(idum) > exp(-(E-E0)*params.beta) )
                        lattice.lattice[x][y] *= -1;
                }
            }
        }
        
        E=0;
        #pragma omp parallel for reduction(+:E)
        for(int i=0; i<lattice.nh; i++) {
            for(int j=0; j<lattice.nw; j++) {
                
                // Obtain average after minimisation
                if(mcs>EQ_MCS) {
                    if(lattice.lattice[i][j]==-1)
                        mean_lattice[i][j]++;
                    else if(lattice.lattice[i][j]==2)
                        mean_lattice[i][j]= -2;
                }
                
                // Calculate the total energy of the configuration
                if(lattice.lattice[i][j] != 2)
                    E += locenergy(&lattice, &params, i, j);
            }
        }
        energymcs[mcs] = E;
    }
    
    //--------------------------------------------------------------------------------/
    
    
    // SAVE DATA   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    
    // Save energy evolution (normalised)
    fid_energymc = fopen(FNAME_MCSWEEPS, "w");
    for(int mcs=0; mcs<MAX_MCS; mcs++)
        fprintf(fid_energymc,"%d\t%f\n", mcs, energymcs[mcs]/EPSNN);
    fclose(fid_energymc);

    
    
    // Save one snapshot with the minimal energy configuration configuration
    fid_lattice = fopen(FNAME_LATTICE_MIN, "w");
    for(int i=0; i<lattice.nh; i++){
        for(int j=0; j<lattice.nw; j++){
            if(lattice.lattice[i][j]==1)
                fprintf(fid_lattice,"%d", 0);
            else if(lattice.lattice[i][j]== -1)
                fprintf(fid_lattice,"%d", 1);
            else
                fprintf(fid_lattice,"%d", 2);
            
            if(j<lattice.nw-1)
                fprintf(fid_lattice,"\t");
        }
        fprintf(fid_lattice,"\n");
    }
    fclose(fid_lattice);
    
    
    // Save mean configuration
    fid_lattice = fopen(FNAME_MEAN_LATTICE, "w");
    for(int i=0; i<lattice.nh; i++){
        for(int j=0; j<lattice.nw; j++){
            
            /* All or nothing
            if(mean_lattice[i][j]==-2) {
                fprintf(fid_lattice,"%d", 2);
            } else {
                if( ((double)mean_lattice[i][j])/(MAX_MCS-EQ_MCS-1) > MEAN_LATTICE_THR)
                    fprintf(fid_lattice,"%d", 1);
                else
                    fprintf(fid_lattice,"%d", 0);
            }
            */
            
            if(mean_lattice[i][j]==-2) {
                fprintf(fid_lattice,"%d", 2);
            } else {
                double rho = mean_lattice[i][j]/(MAX_MCS-EQ_MCS-1);
                if(rho >= 0.00 && rho < 0.25 )
                    fprintf(fid_lattice,"%f", 0.25);
                else if(rho >= 0.25 && rho < 0.50 )
                    fprintf(fid_lattice,"%f", 0.50);
                else if(rho >= 0.50 && rho < 0.75 )
                    fprintf(fid_lattice,"%f", 0.75);
                else if(rho >= 0.75)
                    fprintf(fid_lattice,"%f", 1.0);
            }
            
            
            if(j<lattice.nw-1)
                fprintf(fid_lattice,"\t");
        }
        fprintf(fid_lattice,"\n");
    }
   fclose(fid_lattice);
   
   //--------------------------------------------------------------------------------/
    
    free_dmatrix(mean_lattice, lattice.nw);
    free_imatrix(lattice.lattice, lattice.nw);
    free(energymcs);
   
    return 0;
}


/*************************************************************/
/*                         locenergy                         */
/*************************************************************/
double locenergy(struct Lattice *lattice, 
                 struct Params *params, 
                 int y, int x)
{
    
    double E;
    int xl, xr, yt, yb;
    
    // Von Neumann neighborhood with periodic boundary conditions on the x-axis
    yt = y-1;
    yb = y+1;
    xl = x-1;
    xr = x+1;
    
    xl = xl - (int)floor((double)xl/lattice->nw)*(lattice->nw);
    xr = xr - (int)floor((double)xr/lattice->nw)*(lattice->nw);
    yt = yt - (int)floor((double)yt/lattice->nh)*(lattice->nh);
    yb = yb - (int)floor((double)yb/lattice->nh)*(lattice->nh);
    
    // Top neigthbour
    E=0;
    if(lattice->lattice[yt][x] != 2)
        E += -params->eps*lattice->lattice[yt][x]*lattice->lattice[y][x]/4.0;
    else
        E += -params->bsurf*(1.0-lattice->lattice[y][x])/2.0;
        
    // Bottom neigthbour
    if(lattice->lattice[yb][x] != 2)
        E += -params->eps*lattice->lattice[yb][x]*lattice->lattice[y][x]/4.0;
    else
        E += -params->bsurf*(1.0-lattice->lattice[y][x])/2.0;
        
    // Left neigthbour
    if(lattice->lattice[y][xl] != 2)
        E += -params->eps*lattice->lattice[y][xl]*lattice->lattice[y][x]/4.0;    
    else
        E += -params->bsurf*(1.0-lattice->lattice[y][x])/2.0;
        
    // Rigth neightbour
    if(lattice->lattice[y][xr] != 2)
        E += -params->eps*lattice->lattice[y][xr]*lattice->lattice[y][x]/4.0;   
    else
        E += -params->bsurf*(1.0-lattice->lattice[y][x])/2.0;

    // 2.- Chemical potential
    E += (2.0*params->eps + params->mu)*(lattice->lattice[y][x])/2.0;

    return E;
}





/*************************************************************/
/*                          gen_init                         */
/* Draws the virus as icosahedron cross section              */
/*************************************************************/
void gen_init(struct Lattice *lattice, 
              struct Virus *virus)
{
    // The size is 4xRc
    const double FACTR = 6.0;
    lattice->nw = (int)round(FACTR*(virus->Rc + virus->t/2.0)/dw);
    if(lattice->nw % 2 == 0)
        (lattice->nw)++;

    lattice->nh = lattice->nw;
    lattice->lattice = imatrix(lattice->nh, lattice->nw);

    // Draw the virus
    const double dtheta = 2.0*PI/1000.0;
    const int n_theta = (int)round((2.0*PI - virus->theta0)/dtheta);
    double tradius;
    double x,y;
    double theta;
    int xc,yc;
    int nls_thick;
    int center;
    
    nls_thick=(int)round(virus->t/dw);
    center = (int)floor(lattice->nw/2.0) + 1;
    
    for(int h=0; h<nls_thick; h++) { // Iterate over thickness
        tradius = h*dw + virus->Rc - virus->t/2.0;
        for(int i=0; i<=n_theta; i++) { // Iterate over angles
            
            theta = virus->theta0 + i*dtheta;
            /* This can be used to generate multiple cavities*/
            double tcavity = PI/20.0;
            if(theta>PI/3.0-tcavity && theta<PI/3.0 + tcavity)
                continue;
            if(theta>2.0*PI/3.0-tcavity && theta<2.0*PI/3.0 + tcavity)
                continue;
            if(theta>3.0*PI/3.0-tcavity && theta<3.0*PI/3.0 + tcavity)
                continue;
            if(theta>4.0*PI/3.0-tcavity && theta<4.0*PI/3.0 + tcavity)
                continue;
            if(theta>5.0*PI/3.0-tcavity && theta<5.0*PI/3.0 + tcavity)
                continue;
            if(theta>6.0*PI/3.0-tcavity && theta<6.0*PI/3.0 + tcavity)
                continue;
            
            
            icosahedron(tradius, theta, &x, &y);
            
            xc = (int)round(x/dw + (double)center);
            yc = (int)round(y/dw + (double)center);
            lattice->lattice[xc][yc] = 2;

        }
    }
    
    
    // Fill regions that are empty in the middle of the shell
    /*
    for(int i=1; i<lattice->nw-1; i++) {
        for(int j=1; j<lattice->nw-1; j++) {    
            int xl, xr, yu, yd;
            xl = i-1;
            xr = i+1;
            yu = j-1;
            yd = j+1;
            
            int nns=0;
            if( lattice->lattice[yu][i] == 2)
                nns++; 
            if(lattice->lattice[yd][i] == 2)
                nns++;
            if(lattice->lattice[j][xl] == 2)
                nns++;
            if(lattice->lattice[j][xr] == 2)
                nns++;
            if(nns>3)
                lattice->lattice[i][j] = 2;
        }
    }
    */
    
}



