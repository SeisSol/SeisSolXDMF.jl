# SeisSolXDMF

A Julia library for reading SeisSol XDMF outputs (both POSIX and HDF5 are supported).

## Usage

```julia
using SeisSolXDMF

# Parses the file
xdmf = XDMFFile("my-output.xdmf")

# Topology: Integer array of shape (nDims, nSimplices) where nDims is 3 (triangles) or 4 (simplices)
# Geometry: Float array of shape (3, nPoints)
topology, geometry = grid_of(xdmf)

# Returns a float array of nTimesteps elements containing the respective timestamps
timesteps = timesteps_of(xdmf)

# Read the values of variable "u" at timestep 3 into an array. (Timesteps are indexed starting at 1!)
u = data_of(xdmf, 3, "u")
```