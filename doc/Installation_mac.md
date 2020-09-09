# Installing SQUID on linux machine

### Obtaining Boost Library
1. Download boost library from [http://www.boost.org/](http://www.boost.org/) by running
```
wget https://dl.bintray.com/boostorg/release/1.64.0/source/boost_1_64_0.tar.gz
```
2. Decompress the library into your favorite folder.
```
tar -xzvf boost_1_64_0.tar.gz -C your_favorite_folder
```
3. Open the Makefile under squid in any text editor, modify the path for boost (the first line) to where you decompress boost library:
```
BOOST = <path to your_favorite_folder>/boost_1_64_0
```

### Obtaining BamTools Library
1. Download bamTools library into your favorite folder by running
```
git clone git://github.com/pezmaster31/bamtools.git
```
2. Make sure [CMake](https://cmake.org/) version on our system is >=2.64 by running
```
cmake --version
```
3. Install BamTools
```
cd bamtools
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<bamtools_installation_folder> ..
make
make install
```
4. Open the Makefile under squid again, to change the path for BamTools (the second line). After changing the path for BamTools, the second line in Makefile should look like this
```
BAMTOOLS = <path to bamtools_installation_folder>
```
5. Add BamTools library path to DYLD_LIBRARY_PATH
```
export DYLD_LIBRARY_PATH=<path to your favorite folder>/bamtools/lib:${DYLD_LIBRARY_PATH}
```

### Obtaining GLPK Library
1. Downloading GLPK library from [ftp://ftp.gnu.org/gnu/glpk/](ftp://ftp.gnu.org/gnu/glpk/) by running
```
wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.62.tar.gz
```
2. Decompress GLPK into your favorite folder by running
```
tar -xzvf glpk-4.62.tar.gz -C your_favorite_folder
```
3. Build GLPK library by running
```
cd your_favorate_folder/glpk-4.62
mkdir bin
./configure --prefix=$(pwd)/bin
make
make install 
```
4. Open the Makefile under squid the 3rd time, modify the path for GLPK library on the third line. After modification, the third line should look like
```
GLPK = <path to your_favorite_folder>/glpk-4.62/bin
```
5. Add GLPK to 
```
export DYLD_LIBRARY_PATH=<path to your_favorite_folder>/glpk-4.62/lib:${DYLD_LIBRARY_PATH}
```

### Building SQUID
With the above libraries, you can go into squid folder and build squid simply by running
```
make
```
Then you can find the executable in ./bin under squid folder.
