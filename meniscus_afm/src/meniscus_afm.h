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

//git push -u origin main

#define MAX_MCS  3000                   // Total number of Monte-Carlo sweeps
#define dw       3.24                   // Distance between nodes from LJ minimum of water (A)
#define WIDTH    2000.0                 // Width of the lattice [A]
#define HEIGHT   500.0                  // Heigth of the lattice [A]
#define SURFACE_THICKNESS 20.0          // Thickness of the surface [A]
#define AFM_TIP_RADIUS    50.0          // AFM tip radius [A]
#define AFM_TIP_HEIGTH    30.0          // Heigth to the surface [A]

#define EPS     9.0                           // Interaction nearest neightbours [kJ/mol]
#define RH      20.0                           // Relative humidity [%]
#define MU_C    -2*EPS                      // Critical chemical potential [kJ/mol] 
#define BSURF   3.0*EPS                     // Interaction of water with surface [kJ/mol]
#define R       8.31446261815324e-3         // Ideal gas constant [kJ/mol]
#define T       298.0                      // Temperature [K]
#define BETA    1.0/(R*T)                  // (Inverse) Energy of the bath
#define MU      MU_C + BETA*log(RH)           // Chemical potenial energy [kJ/mol] 


#define FNAME_LATTICE_0       "initial_lattice.dat"
#define FNAME_LATTICE_MIN     "minimised_lattice.dat"
#define FNAME_MCSWEEPS        "energy_sweeps.dat"



struct Params
{
  double eps;       //
  double mu;        //
  double mu_c;       //
  double bsurf;      //
  double beta;      // 1/kBT
};

struct Lattice
{
    int nw;
    int nh;
    int nsurf;
    int **lattice;
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

 
