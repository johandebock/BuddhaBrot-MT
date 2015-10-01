export PATH=".:/usr/local/bin:/mingw/bin:/bin:$PATH"
# http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/4.9.3/threads-win32/sjlj/i686-4.9.3-release-win32-sjlj-rt_v4-rev1.7z/download
export PATH="/c/mingw-i686-4.9.3-release-win32-sjlj-rt_v4-rev1/bin:$PATH"

mkdir png-windows-32bit-mingw
cd png-windows-32bit-mingw
wget http://prdownloads.sourceforge.net/libpng/zlib-1.2.8.tar.gz?download
wget http://prdownloads.sourceforge.net/libpng/libpng-1.6.18.tar.gz?download
tar -zxf zlib-1.2.8.tar.gz
rm zlib-1.2.8.tar.gz
tar -zxf libpng-1.6.18.tar.gz
rm libpng-1.6.18.tar.gz
mv zlib-1.2.8 zlib
mv libpng-1.6.18 libpng
cd zlib
make -fwin32/Makefile.gcc
cd ..
cd libpng
export LDFLAGS=-L../zlib
export CPPFLAGS=-I../zlib
./configure
make
cd ..

cp zlib/libz.a ./
cp libpng/.libs/libpng16.a libpng.a
cp zlib/zlib.h ./
cp zlib/zconf.h ./
cp libpng/png.h ./
cp libpng/pngconf.h ./
cp libpng/scripts/pnglibconf.h.prebuilt pnglibconf.h
rm -R zlib/
rm -R libpng/
cd ..

mkdir SDL-windows-32bit-mingw
cd SDL-windows-32bit-mingw
wget --no-check-certificate https://www.libsdl.org/release/SDL2-devel-2.0.3-mingw.tar.gz
tar -zxf SDL2-devel-2.0.3-mingw.tar.gz
rm SDL2-devel-2.0.3-mingw.tar.gz
cd SDL2-2.0.3
cp i686-w64-mingw32/lib/libSDL2.a ../
cp i686-w64-mingw32/lib/libSDL2main.a ../
cp i686-w64-mingw32/include/SDL2/*.h ../
cd ..
rm -R SDL2-2.0.3
cd ..
