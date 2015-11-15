mkdir png-linux
cd png-linux
wget http://prdownloads.sourceforge.net/libpng/zlib-1.2.8.tar.gz?download
wget http://prdownloads.sourceforge.net/libpng/libpng-1.6.18.tar.gz?download
mv zlib-1.2.8.tar.gz\?download zlib-1.2.8.tar.gz
mv libpng-1.6.18.tar.gz\?download libpng-1.6.18.tar.gz
tar -zxf zlib-1.2.8.tar.gz
rm zlib-1.2.8.tar.gz
tar -zxf libpng-1.6.18.tar.gz
rm libpng-1.6.18.tar.gz
mv zlib-1.2.8 zlib
mv libpng-1.6.18 libpng
cd zlib
./configure
make
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
