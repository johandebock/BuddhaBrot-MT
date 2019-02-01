g++ -Wall -DLINUX -O2 -march=native -fopenmp -IdSFMT -Ipng-linux -o BuddhaBrot-MT-linux BuddhaBrot-MT.cpp dSFMT/*.cpp png-linux/libpng.a png-linux/libz.a -lSDL2
