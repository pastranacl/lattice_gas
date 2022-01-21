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

#define MAX_MCS  3000                // Total number of Monte-Carlo sweeps
#define R        123                 // Ideal gas constant
#define dw       3.24                // Distance between nodes from LJ minimum of water (A)

#define FNAME_LATTICE_0       "initial_lattice.dat"
#define FNAME_LATTICE_MIN     "minimised_lattice.dat"
#define FNAME_MCSWEEPS        "energy_sweeps.dat"

//git push -u origin main

struct Params
{
  double eps;       //
  double mu;        //
  double b;         //
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

 
