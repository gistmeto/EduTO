# EduTO
A MATLAB code of node-based topology optimization in 3D arbitrary domain for additive manufacturing

https://link.springer.com/article/10.1007/s00158-022-03339-1

Abstract

This paper presents a MATLAB code for node-based topology optimization that can handle a design problem with a three-dimensional (3D) arbitrary-shaped domain. For the meshing of arbitrary geometry, an open-source 3D mesh generator, GMSH, is utilized in this work. Here, a linear four-noded tetrahedral element is utilized due to its advantage in mesh generation. A MATLAB program is composed of three procedures. The pre-processing aims to import mesh and input files into MATLAB workspace. In the main processing, node-based topology optimization is carried out with the well-established three-field projection scheme. The post-processing aims to generate a Computer-Aided Design (CAD) file in an STL format. For this, the zero-level set of filtered density field is utilized to define the boundary of a topology optimization result. From the STL format CAD file, a design result is fabricated using additive manufacturing machines. The effectiveness of the MATLAB code is examined through three design examples including a simply supported beam, bridge, and airplane bearing bracket.
