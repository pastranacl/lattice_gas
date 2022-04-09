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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

//git push -u origin main

#define PI 3.14159265358979323846264338327950288419716939937

#define MAX_MCS           4500          // Total number of Monte-Carlo sweeps
#define EQ_MCS            4000          // Number of points to equilibrate before start sampling
#define dw                3.24          // Distance between nodes from LJ minimum of water (A)

#define EPSNN   9.0                     // Interaction nearest neighbours [kJ/mol]
#define RH      0.90                    // Relative humidity [%]
#define MU_C    -2.0*EPSNN              // Critical chemical potential [kJ/mol] 
#define BSURF   3.0*EPSNN               // Interaction of water with surface [kJ/mol]
#define R       8.31446261815324e-3     // Ideal gas constant [kJ/mol]
#define T       298.0                   // Temperature [K]
#define BETA    1.0/(R*T)               // (Inverse) Energy of the bath
#define MU      MU_C + R*T*log(RH)      // Chemical potenial energy [kJ/mol] 

#define VIR_R   25.0                    // Radius of virus [nm]
#define VIR_T   3.0                     // Thickness of the viral shell [nm] 
#define THETA0  0.0                     // Range of theta to plot

#define MEAN_LATTICE_THR   0.50         // Threshold to consider has occupied one spin

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
    int **lattice;  // Array for the lattice
};

struct Virus
{
    double Rc;      // Virus radius capside
    double t;       // Virus shell thickness
    double theta0;  // Cavitiy of the virus
};

void gen_init(struct Lattice *lattice, 
              struct Virus *virus);

double locenergy(struct Lattice *lattice, 
                 struct Params *params, 
                 int x, int y);




/****************************************************************/
/*                          icosahedron                         */
/* Calculates the xy positions of a polar definition to draw a  */
/* a shape that resembles the cross-section of an icosahedron   */ 
/****************************************************************/
inline void icosahedron(double R0, double theta, double *x, double *y)
{
    const double A = 0.03;
    double pf = (1.0 + A*sin(6.0*theta));
    double r;
    
    r = sqrt(R0*R0*pf*pf);
    *x = r*cos(theta);
    *y = r*sin(theta);
}







