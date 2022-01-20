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



#include "afm_meniscus.h"



int main()
{
    // Generate the lattice
    struct Lattice lattice;
    gen_init(&lattice, 4204.0, 1000.0, dw*20, 100*dw, dw);
    
    
    // Fill the data with water
    for(int i=0; i<lattice.nh; i++){
        for(int j=0; j<lattice.nw; j++){
            if(lattice.lattice[i][j]==0)
                lattice.lattice[i][j]=1;
        }
    }
    
    // Save the initial lattice
    FILE *fid_lattice;
    fid_lattice = fopen("initial_lattice.dat", "w");
    
    for(int i=0; i<lattice.nh; i++){
        for(int j=0; j<lattice.nw; j++){
            fprintf(fid_lattice,"%d\t", lattice.lattice[i][j]);
        }
        fprintf(fid_lattice,"\n");
    }
    fclose(fid_lattice);
    
    
    // MAIN MONTE CARLO ROUTINE ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    double E0, E;
    for(int mcs=0; mcs<MAX_MCS; mcs++){
        
        for(int i=0; i<lattice.nh; i++){ 
            for(int j=0; j<lattice.nw; j++){
                    
                    if(j==2) { // Exclude the surfaces
                        continue;
                    } else {
                        
                        
                        E0 = energy(&lattice, &params, i, j);
                        lattice.lattice[i][j] *= -1;
                        E = energy(&lattice, &params, i, j);
                        
                        
                        //
                        if( exp((E-E0)*beta)< rand())
                            lattice.lattice[i][j] *= -1;
                        
                    }
                    
                    
                    
                    
                    
            }
        }
        
        // Calculate the energy
        double Et=0;
        #pragma omp parallel for atomic
        for(int i=0; i<; i++) {
            for(int j=0; j< j++) {
                Et += locenergy(&lattice, &params, i, j);
            }
        }
        
        
    }
    
    //--------------------------------------------------------------------------------
    
    
    
    // Save final configuration
    
    
    
   // Save energy configuration
   
   
   
   
    return 0;
}


/*************************************************************/
/*                         locenergy                         */
/*************************************************************/
double locenergy(struct Lattice *lattice, 
                 struct Params *params, 
                 int x, int y)
{
    
    double Et=0;
    
    
    // Periodic boundary conditions on the x axis
    xl =
    xr = 
    yt = y-1;
    yb = y+1;
    
    
    // Top neigthbour
    if(yt>=0)
    
        
    // Bottom neigthbour
        
        
    // Left neigthbour
        
        
    // Rigth neightbour
    
        
        
        
    
    // 2. Interaction with surfaces +++++++++
    if(yt>=0) {
        if( lattice[x][yt]==2 || lattice[x][yb]==2 || lattice[xl][y]==2 || lattice[xr][y]==2)
            Et += 0.5*params->b*(lattice[x][y]);
        
    } else {
        if(lattice[x][yb]==2 || lattice[xl][y]==2 || lattice[xr][y]==2)
            Et += 0.5*params->b*(lattice[x][y]);
    }
       
        
        
        
    // 3. Chemical energy due to external source +++++++++
    Et += (2.0*params->eps + params->mu)*(lattice->lattice[i]][j])/2.0:
    
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
             lattice->lattice[i][j] = 0;
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
    int x0 = (int)ceil(lattice->nw/2.0);
    
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
