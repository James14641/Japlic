import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy

# load input data

# load a vtk file as input
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("FEBE_0.vtk")
reader.Update()

# Get the coordinates of nodes in the mesh
nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
