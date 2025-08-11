# Symbolic Perturbation for Bivariate Jacobi Sets and Reeb spaces
An open-source library for generating evaluation tables for the orientation of points and the order of segments predicate. This is based on the following paper -

```
@article{

}
```

# Required Python Librarires

## Core

```
sympy       1.14.0
numpy       1.16.4
```

## (Optional) Visualisation

```
matplotlib  3.8.4
ipywidgets  8.1.2
```




# Printing Evaluation Tables  

### Point Orientation

```
python3 -m table_generation.printPointOrientationTables
```

### Segment Order
```
python3 -m table_generation.printSegmentOrderTables
```




# Synthetic Data Experiment

In all the following commands replace 100 with the number of iterations that you would like the procedure to run for.

### Point Orientation
```
python3 -m synthetic_data_evaluation.testPointOrientationAll 100
```

### Segment Order
```
python3 -m synthetic_data_evaluation.testSegmentOrderAll.py 100
```




# Real-World Data Experiment




## Installing Dependencies

### CGAL

```
git clone --recursive https://github.com/CGAL/cgal cgal
cd ./cgal
git checkout v6.0.1 

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/librarires/cgal/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4
make install
```


### VTK

```
git clone --recursive https://gitlab.kitware.com/vtk/vtk.git VTK-9.4.1
cd ./VTK-9.4.1
git checkout v9.4.1 

mkdir build install
cd build

cmake -DCMAKE_INSTALL_PREFIX="$projectFolder/librarires/VTK-9.4.1/install" -DCMAKE_BUILD_TYPE="Release" ..
make -j 4
make install
```


### Compiling

## Build
```
cmake -DCMAKE_PREFIX_PATH="$projectFolder/cgal/install;$projectFolder/VTK-9.4.1/install" -DCMAKE_EXPORT_COMPILE_COMMANDS=On -DCMAKE_BUILD_TYPE=Release ..
make
```
