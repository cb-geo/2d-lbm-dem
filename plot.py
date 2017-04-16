import sys
import numpy as np
from vtk import vtkRectilinearGridReader
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt

filename = sys.argv[1]

reader = vtkRectilinearGridReader()
reader.SetFileName(filename)
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

nx, ny, nz = reader.GetOutput().GetDimensions()
data = reader.GetOutput().GetPointData()

grain_pressure = vtk_to_numpy(data.GetArray('grain_pressure')).reshape((ny, nx))
fluid_pressure = vtk_to_numpy(data.GetArray('fluid_pressure')).reshape((ny, nx))

plt.axis('equal')
plt.pcolormesh( np.where(grain_pressure >= 0, grain_pressure, fluid_pressure) )
plt.colorbar()
plt.savefig("plot.png")
