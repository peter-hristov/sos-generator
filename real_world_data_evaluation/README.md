# Compile

Here's how to compile the code with cmake, using CGAL version 6.0.1 and VTK 9.4.1

Build version
```
cmake -DCGAL_DIR=/home/peter/Projects/libraries/cgal/build -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On ..

```

Release version
```
cmake -DCMAKE_PREFIX_PATH="/home/peter/Projects/libraries/cgal/install;/home/peter/Projects/libraries/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On -DCMAKE_BUILD_TYPE=Release ..
```


make





# Test

make && ./fv99 -f ../data/three-sheet-toy.txt
make && ./fv99 -f ~/Projects/data/reeb-space-test-data/data.vtu
make && ./fv99 -f ~/Projects/data/reeb-space-test-data/torus/torus-factor-50-tets-320.vtu -e 1e-3
make && ./fv99 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-20-300.vtu



# Compiling on Tetralith
cmake -DCMAKE_PREFIX_PATH="/proj/reeb-space-storage/users/x_pethr/librarires/boost_1_88_0/install;/proj/reeb-space-storage/users/x_pethr/librarires/cgal-6.0.1/install;/proj/reeb-space-storage/users/x_pethr/librarires/vtk-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On -DCMAKE_BUILD_TYPE=Release ..


