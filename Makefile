DEBYE_SCREENED_COULOMB_KERNEL=666
LATTICE_EQCD_KERNEL=111
LEADING_ORDER_KERNEL=222
NEXT_TO_LEADING_ORDER_KERNEL=333
MULTIPLE_SOFT_SCATTERING_KERNEL=999


COLLISION_KERNEL=${LATTICE_EQCD_KERNEL}

GToGG=11
QToQG=21
GToQQ=12
PROCESS=${GToGG}

FullRate : 
	g++ -std=c++20 -I/usr/local/include/ -L/usr/local/lib/ -DRATETYPE=${FULLRATE} -DPROCESS=${PROCESS} -DCOLLISION_KERNEL=${COLLISION_KERNEL} -pipe -O3 src/FiniteMedium.cpp -o FullRate.out -lpthread -lm -ldl -lcuba -fopenmp -g -lgsl -lgslcblas 
	mkdir OUTPUT 

OpacityRate : 
	g++ -std=c++20 -DRATETYPE=${OPACITY} -DPROCESS=${PROCESS} -DCOLLISION_KERNEL=${COLLISION_KERNEL} -pipe -O3 src/Opacity.cpp -o OpacityRate.out -lpthread -lm -ldl -lcuba -fopenmp -g -lgsl -lgslcblas 
	mkdir Opacity 

ImprovedOpacity:
	g++ -std=c++20 -DCOLLISION_KERNEL=${COLLISION_KERNEL} -DPROCESS=${PROCESS} -pipe -O3 src/ImprovedOpacity.cpp -o ImprovedOpacity.out -lpthread -lm -ldl -lcuba -fopenmp -g -lgsl -lgslcblas 
	mkdir OpacityImproved

Harmonic:
	g++ -std=c++20 -DCOLLISION_KERNEL=${COLLISION_KERNEL} -DPROCESS=${PROCESS} -pipe -O3 src/HO.cpp -o HO.out -lpthread -lm -ldl -lcuba -fopenmp -g -lgsl -lgslcblas 
	mkdir HO