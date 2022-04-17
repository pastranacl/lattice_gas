# Lattice gas 

Lattice gas models, interaction water-water and water-surfaces.
(c) 2022 Cesar Lopez Pastrana

## Description
These codes reproduces the lattice-gas model in the presence of surfaces described in [1,2]. 
Two different systems are shown. First, the interaction of water with an AFM tip as in the original works [1,2]. Second, the interaction of water with a viral capsid [3]. 

## Compilation
Simply call `make` on the main folder of each project. The executable files will be generated on the `bin` folder with names `afmeniscus` or `icovirus`.

## Execution
The programs are called simply by their name, with no additional parameters. The exectution produces four files:

1. `energy_sweeps.dat`: Energy as a function of the steps, where each step involves NxM test moves.

2. `initial_lattice.dat`: As the name indicates, simply the original configuration

3. `minimised_lattice.dat`: Snapshot of the configuration at the very last step

4. `mean_minim_lattice.dat`: Average configuration considering multiple steps (configured on the header files). The results produced by this file can be shown by calling the python script in the `bin` folder: `python3 plot_repp.py`



## References

[1] Jang, G. C. Schatz, and M. A. Ratner, The Journal of Chemical Physics 116, 3875 (2002)

[2] J. Jang, G. C. Schatz, and M. A. Ratner, Phys. Rev. Lett. 92, 085504 (2004)

[3] C. Carrasco, M. Douas, R. Miranda, M. Castellanos, P. A. Serena, J. L. Carrascosa, M. G. Mateu, M. I. Marques, and P. J. de Pablo, PNAS 106, 5475 (2009)

