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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

//git push -u origin main

#define MAX_MCS           3000          // Total number of Monte-Carlo sweeps
#define EQ_MCS            2000          // Number of points to equilibrate and start sampling
#define dw                3.24          // Distance between nodes from LJ minimum of water (A)
#define WIDTH             500.0         // Width of the lattice [A]
#define HEIGHT            300.0         // Heigth of the lattice [A]
#define SURFACE_THICKNESS 20.0          // Thickness of the surface [A]
#define AFM_TIP_RADIUS    70.0          // AFM tip radius [A]
#define AFM_TIP_HEIGTH    15.00          // Heigth to the surface [A]


#define EPSNN   9.0                     // Interaction nearest neighbours [kJ/mol]
#define RH      0.50                    // Relative humidity [%]
#define MU_C    -2.0*EPSNN              // Critical chemical potential [kJ/mol] 
#define BSURF   3.0*EPSNN              // Interaction of water with surface [kJ/mol]
#define R       8.31446261815324e-3     // Ideal gas constant [kJ/mol]
#define T       298.0                   // Temperature [K]
#define BETA    1.0/(R*T)               // (Inverse) Energy of the bath
#define MU      MU_C + R*T*log(RH)      // Chemical potenial energy [kJ/mol] 


#define FNAME_LATTICE_0       "initial_lattice.dat"
#define FNAME_LATTICE_MIN     "minimised_lattice.dat"
#define FNAME_MEAN_LATTICE    "mean_minim_lattice.dat"
#define FNAME_MCSWEEPS        "energy_sweeps.dat"



struct Params
{
  double eps;       // Interaction nearest neighbours [kJ/mol]
  double mu;        // Chemical potenial energy [kJ/mol] 
  double mu_c;      // Critical chemical potential [kJ/mol] 
  double bsurf;     // Interaction of water with surface [kJ/mol]
  double beta;      // 1/kBT parameter []
};


struct Lattice
{
    int nw;         // Number columns
    int nh;         // Number of rows
    int nsurf;      // Points defining the surface
    int **lattice;  // Array for the lattice
};


double locenergy(struct Lattice *lattice, 
                 struct Params *params, 
                 int x, int y);

void gen_init(struct Lattice *lattice, 
              double w,
              double h,
              double surf_thick,
              double tip_radius,
              double dtipsurf);
