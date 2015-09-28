FLAGS = -std=c++11 -lm -lstdc++ -I lio -I eigen/ -I raoul/waveblocks -I ben/waveblocks -g -I . -I eigen3-hdf5 -lhdf5 -lhdf5_cpp

harmonic_2D: harmonic_2D.cpp
	clang -o harmonic_2D harmonic_2D.cpp $(FLAGS)

tunneling_1D: tunneling_1D.cpp
	clang -o tunneling_1D tunneling_1D.cpp $(FLAGS)

nonsense_2N_2D: nonsense_2N_2D.cpp
	clang -o nonsense_2N_2D nonsense_2N_2D.cpp $(FLAGS)

nonsense_2N_1D: nonsense_2N_1D.cpp
	clang -o nonsense_2N_1D nonsense_2N_1D.cpp $(FLAGS)
	
nonsense_2N_2D_inhom: nonsense_2N_2D_inhom.cpp
	clang -o nonsense_2N_2D_inhom nonsense_2N_2D_inhom.cpp $(FLAGS)


nonsense_2N_1D_inhom: nonsense_2N_1D_inhom.cpp
	clang -o nonsense_2N_1D_inhom nonsense_2N_1D_inhom.cpp $(FLAGS)
	
all: nonsense_2N_2D nonsense_2N_2D_inhom nonsense_2N_1D_inhom nonsense_2N_1D tunneling_1D harmonic_2D

clean:
	rm -f nonsense_2N_2D nonsense_2N_1D nonsense_2N_1D_inhom nonsense_2N_2D_inhom tunneling_1D harmonic_2D
