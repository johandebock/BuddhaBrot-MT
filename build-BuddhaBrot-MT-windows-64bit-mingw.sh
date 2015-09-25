export PATH=".:/usr/local/bin:/mingw/bin:/bin:$PATH"
export PATH="/c/mingw-x86_64-4.9.3-release-win32-sjlj-rt_v4-rev1/bin:$PATH"
g++ -Wall -DWINDOWS -s -static -O3 -mtune=nehalem -fopenmp -mfpmath=sse -IdSFMT -Ipng-mingw-64bit -ISDL-mingw-64bit -o BuddhaBrot-MT-windows-64bit.exe BuddhaBrot-MT.cpp dSFMT/*.cpp png-mingw-64bit/libpng.a png-mingw-64bit/libz.a -lmingw32 SDL-mingw-64bit/libSDL2main.a SDL-mingw-64bit/libSDL2.a -mwindows -lgdi32 -limm32 -lole32 -loleaut32 -lversion -lwinmm
