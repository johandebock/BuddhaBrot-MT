export PATH=".:/usr/local/bin:/mingw/bin:/bin:$PATH"
# http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/4.9.3/threads-win32/sjlj/x86_64-4.9.3-release-win32-sjlj-rt_v4-rev1.7z/download
export PATH="/c/mingw-x86_64-4.9.3-release-win32-sjlj-rt_v4-rev1/bin:$PATH"

g++ -Wall -DWINDOWS -s -static -O3 -march=core2 -fopenmp -mfpmath=sse -IdSFMT -Ipng-windows -ISDL-windows -o BuddhaBrot-MT-windows.exe BuddhaBrot-MT.cpp dSFMT/*.cpp png-windows/libpng.a png-windows/libz.a -lmingw32 SDL-windows/libSDL2main.a SDL-windows/libSDL2.a -mwindows -lgdi32 -limm32 -lole32 -loleaut32 -lversion -lwinmm
