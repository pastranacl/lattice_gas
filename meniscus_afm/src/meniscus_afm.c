/****************************************************************************
    
    Simulation of the hydration layer of an AFM tip
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

#include "meniscus_afm.h"
#include "rands.h"


int main()
{
    double E0, E;
    long *idum; 
    long foornd; 
    int **mean_lattice;
    double *energymcs;
    
    FILE *fid_lattice;
    FILE *fid_energymc;
    
    foornd = -time(NULL);
    idum = &foornd;
    
    // PRELIMINARY STEPS   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    
    // Generate the lattice and memory allocation
    struct Lattice lattice;
    gen_init(&lattice, WIDTH, HEIGHT, SURFACE_THICKNESS, AFM_TIP_RADIUS, AFM_TIP_HEIGTH);
    
    
    // Fill the empty space randomly
    for(int i=0; i<lattice.nh; i++){
        for(int j=0; j<lattice.nw; j++){
            if(lattice.lattice[i][j]!=2){
                if(ran1(idum)>0.5)
                    lattice.lattice[i][j] = 1;
                else
                    lattice.lattice[i][j] = -1;
            }
        }
    }
    
    // Allocate memmory for the mean lattice
    mean_lattice = malloc((lattice.nh)*sizeof(int *));
    for(int i=0;i<(lattice.nh); i++)
        mean_lattice[i] = malloc(lattice.nw*sizeof(int));
    
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
    
    
    // Derived quantities
    struct Params params;
    params.eps = EPSNN;       
    params.mu = MU;       
    params.bsurf = BSURF;
    params.beta = BETA;      
  
    //--------------------------------------------------------------------------------/
    
    
    
    
    // MAIN MONTE CARLO ROUTINE   ++++++++++++++++++++++++++++++++++++++++++++++++++++/
    E=0;
    energymcs = malloc(MAX_MCS*sizeof(int *));

    for(int i=0; i<lattice.nh; i++) {
        for(int j=0; j<lattice.nw; j++) {
            if(lattice.lattice[i][j] != 2)
                E += locenergy(&lattice, &params, i, j);
        }
    }
    energymcs[0] = E;
    
    
    for(int mcs=1; mcs<MAX_MCS; mcs++)
    {
        for(int i=0; i<lattice.nh; i++){ 
            for(int j=0; j<lattice.nw; j++)
            {
                    int x,y;
                    x = (int)(lattice.nh-1)*ran1(idum);
                    y = (int)(lattice.nw-1)*ran1(idum);
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
        for(int i=0; i<lattice.nh; i++) {
            for(int j=0; j<lattice.nw; j++) {
                
                // Obtain average after minimisation
                // TODO: THIS HAS TO BE CHANGED FOR THE RANGE 0 AND 1
                if(mcs>EQ_MCS)
                    mean_lattice[i][j] += lattice.lattice[i][j];
                
                // Calculate the total energy of the configuration
                if(lattice.lattice[i][j] != 2)
                    E += locenergy(&lattice, &params, i, j);
            }
        }
        energymcs[mcs] = E;
        
    }
    
    //--------------------------------------------------------------------------------/
    
    
    // SAVE DATA   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    
    // Save final configuration
    fid_energymc = fopen(FNAME_MCSWEEPS, "w");
    for(int mcs=0; mcs<MAX_MCS; mcs++){
        fprintf(fid_energymc,"%d\t%f\n", mcs, energymcs[mcs]);
    }
    fclose(fid_energymc);

    
    
    // Save minimal energy configuraion configuration
    fid_lattice = fopen(FNAME_LATTICE_MIN, "w");
    for(int i=0; i<lattice.nh; i++){
        for(int j=0; j<lattice.nw; j++){
            fprintf(fid_lattice,"%d\t", lattice.lattice[i][j]);
        }
        fprintf(fid_lattice,"\n");
    }
    fclose(fid_lattice);
   
   //--------------------------------------------------------------------------------/
    
    
    
   //TODO: Liberar memoria!!
    
    
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
    
    
    //TODO: SEEMSTO BE AN ERROR WITH THE BOUNDARY CONDITION PERIODIC ON THE RIGTH SIDE
    
    
    // Top neigthbour
    E=0;
    if(yt>=0) {
       if(lattice->lattice[yt][x] != 2)
           E += -params->eps*lattice->lattice[yt][x]*lattice->lattice[y][x]/4.0;
       else
            E += -params->bsurf*(1.0-lattice->lattice[y][x])/2.0;
    }
        
    // Bottom neigthbour
    if(lattice->lattice[yb][x] != 2 && yb<lattice->nw-1)
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
/* Generates the lattice points with the surface and the AFM */
/* tip. The AFM tip is approximated by a parabola.           */
/*************************************************************/
void gen_init(struct Lattice *lattice, 
              double w,
              double h,
              double surf_thick,
              double tip_radius,
              double dtipsurf)
{

    // 0. Find properties of the lattice, allocate memory and set to zero
    lattice->nw    = (int)round(w/dw);
    lattice->nh    = (int)round(h/dw);
    lattice->nsurf = (int)round(surf_thick/dw);
    
    if(lattice->nw % 2 == 0)
        (lattice->nw)++;
    
    lattice->lattice = malloc((lattice->nh)*sizeof(int *));
    for(int i=0;i<(lattice->nh); i++)
        lattice->lattice[i] = malloc(lattice->nw*sizeof(int));
    
    for(int i=0;i<(lattice->nh); i++){
        for(int j=0;j<(lattice->nw); j++){
             lattice->lattice[i][j] = 1;
        }
    }
    

    // 1. Draw the points of the surface
    int ns = (lattice->nh) - (lattice->nsurf);
    for(int i=ns;i<(lattice->nh); i++){
        for(int j=0;j<(lattice->nw); j++){
             lattice->lattice[i][j] = 2;
        }
    }
    
      
    // 2. Draw the AFM tip
    int x, y, y_parabola;
    int x0 = (int)round(lattice->nw/2.0);
    
    for(int i=0; i<ns; i++){
        for(int j=0;j<(lattice->nw); j++){
            x = (j-x0)*dw;
            y = (ns-1-i)*dw;
            y_parabola = x*x/(2*tip_radius) + dtipsurf;
            if(y>=y_parabola)
                lattice->lattice[i][j] = 2;
        }
    }
}



// TODO: ALLOCATE AND DEALLOCATE FUNCTIONS
