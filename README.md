# Finite-Medium-Splitting-Rates-Using-Non-perturbative-Kernel
Determination of the medium-induced splitting rates using non-perturbatively determined broadening kernel

## Compiling: 
The Makefile has two variable `COLLISION_KERNEL=${LATTICE_EQCD_KERNEL}` and `PROCESS=${GToGG}` which set the broadening kernel used and the process computed. 

Different programs are available to obtain different rates: 
- `make FullRate` creates the executable `FullRate.out`, which computes the full finite medium rate in the output folder "OUTPUT".
- `make Opacity` creates the executable `OpacityRate.out`, which computes the first order Opacity rate in the output folder "Opacity". 
- `make ImprovedOpacity` creates the executable `ImprovedOpacity.out`, which computes the first order Opacity rate in the output folder "OpacityImproved". 
- `make Harmonic` creates the executable `HO.out`, which computes the first order Opacity rate in the output folder "HO". 


## Runing 

To run any executable `exe.out`, run the command `./exe.out -P x -z y` where $P=x$ is the parent particle's energy in units of temperature $[T]$, and $z=y$ is the momentum fraction of the emission with energy $\omega= zP$. 
The rate is written into a file `OUTPUTFolder/Rate-Px-zy.txt` where the first column is the dimensionless time $\tau = \frac{T^2 t}{2Pz(1-z)}T$ and the second column is the rate $\frac{d\Gamma_{a}^{bc}}{dz}(P,z,t)$ in units of $T$. 


## Example with Plot
In addition to the makefile, we provide 3 scripts `GToGG.sh`, `QToGQ.sh` and `GToQQ.sh`. In order to create comparison plots of the non-perturbative broadening kernel using all the different approximation at $P=300T$ and $z=0.25$, follow these steps:

- Run each file successively using `source File.sh`: it computes the radiation rate for all the different approximation and the output is moved to the folders inside `PlotMaking/Process/File`. 
- Then cd to the folder `cd PlotMaking` and run `source MakePlots.sh`
