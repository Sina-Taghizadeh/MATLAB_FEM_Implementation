# MATLAB_FEM_Implementation
Complete FEM implementation for 2D elasticity problems with truss, triangular and quadrilateral elements
# description
This code provides a complete implementation of the linear finite element method for 2D elasticity problems. The goal of this code is to understand the logic behind well-known finite element software such as Abaqus, ANSYS, etc. The input is an Excel file that must be filled in by the user. The output includes the displacements and forces of each node and the stress components in each triangular and quadrilateral element and the internal force for truss elements; it also shows the overall shape of the structure before and after loading graphically. The program automatically detects the element type based on the filled-in Excel sheet. Automatic triangle mesh generation is also provided for rectangular domains, where only the number of elements in the longitudinal and transverse directions of the rectangle need to be selected in the program. By default, it solves for the plane stress state, but it can be easily changed to plane strain in the code.
# Examples
Below are 3 sample examples solved using this code.
# 1. Truss problem
We have a truss as shown below.

![image](https://github.com/Sina-Taghizadeh/MATLAB_FEM_Implementation/assets/162900845/cf8459f3-973e-47eb-baf4-a105e376ba2c)

To analyze it, we need to enter its properties in the Excel file, which is named "input_sample_truss.xlsx" and exists in the repo. The results include the displacements and forces of each node and the internal force in each element as follows.

![image](https://github.com/Sina-Taghizadeh/MATLAB_FEM_Implementation/assets/162900845/c5d5f95d-0302-401e-9479-8861c3c4e820)

![image](https://github.com/Sina-Taghizadeh/MATLAB_FEM_Implementation/assets/162900845/062ea86b-27b9-43c2-a14c-5e476437eee2)

The overall shape of the truss, before and after loading, is also shown below.

![image](https://github.com/Sina-Taghizadeh/MATLAB_FEM_Implementation/assets/162900845/3ee64c63-fbe2-4c8d-96e0-bf15aa1698a4)

# 2. Beam problem with triangular and quadrilateral elements
Consider the following beam with the following specifications.

![image](https://github.com/Sina-Taghizadeh/MATLAB_FEM_Implementation/assets/162900845/1d1cf5bc-3a3a-46c9-a3d4-8dced8b419ea)
![image](https://github.com/Sina-Taghizadeh/MATLAB_FEM_Implementation/assets/162900845/f8e0c878-8740-4710-9ba4-17a9d08436b2)

The Excel files for solving this problem are named "input_sample_triangle.xlsx" and "input_sample_quadrilateral.xlsx" and exist in the repo. Due to the large number of nodes and elements, we do not show the tables related to them here. The following figures show the beam before and after loading in the triangular and quadrilateral element cases.

![image](https://github.com/Sina-Taghizadeh/MATLAB_FEM_Implementation/assets/162900845/1100cfed-6ec7-4a73-8951-9d17879ea9c7)
![image](https://github.com/Sina-Taghizadeh/MATLAB_FEM_Implementation/assets/162900845/0c84de7a-4c49-47fd-a5fb-f1341d518fa1)

It is worth noting that non-rectangular quadrilateral elements can also be used, but the results for quadrilateral elements seem to be incorrect. I would be happy if someone could correct it and let me know.









