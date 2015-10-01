export PATH=".:/usr/local/bin:/mingw/bin:/bin:$PATH"
# http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/4.9.3/threads-win32/sjlj/i686-4.9.3-release-win32-sjlj-rt_v4-rev1.7z/download
export PATH="/c/mingw-i686-4.9.3-release-win32-sjlj-rt_v4-rev1/bin:$PATH"

# -mconsole for console output
g++ -Wall -DWINDOWS -s -static -O3 -march=core2 -fopenmp -mfpmath=sse -IdSFMT -Ipng-windows-32bit-mingw -ISDL-windows-32bit-mingw -o BuddhaBrot-MT-windows-32bit.exe BuddhaBrot-MT.cpp dSFMT/*.cpp png-windows-32bit-mingw/libpng.a png-windows-32bit-mingw/libz.a -lmingw32 SDL-windows-32bit-mingw/libSDL2main.a SDL-windows-32bit-mingw/libSDL2.a -mwindows -lgdi32 -limm32 -lole32 -loleaut32 -lversion -lwinmm
