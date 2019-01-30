g++ -Wall -DLINUX -s -static -O3 -march=native -fopenmp -mfpmath=sse -IdSFMT -Ipng-linux -o BuddhaBrot-MT-daemon-linux BuddhaBrot-MT-daemon.cpp dSFMT/*.cpp png-linux/libpng.a png-linux/libz.a
