g++ -Wall -DLINUX -O2 -march=native -fopenmp -IdSFMT -Ipng-linux -o BuddhaBrot-MT-daemon-linux BuddhaBrot-MT-daemon.cpp dSFMT/*.cpp png-linux/libpng.a png-linux/libz.a
