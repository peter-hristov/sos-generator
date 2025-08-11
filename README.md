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

## Point Orientation

```
python3 -m table_generation.printPointOrientationTables
```

## Segment Order
```
python3 -m table_generation.printSegmentOrderTables
```




# Running Synthetic Data Experiment

## Point Orientation
```
python3 -m synthetic_data_evaluation.testPointOrientationAll 100
```

## Segment Order
```
python3 -m synthetic_data_evaluation.testSegmentOrderAll.py 100
```

