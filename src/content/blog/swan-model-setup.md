---
title: "Setting Up SWAN for Coastal Wave Modeling"
description: "Step-by-step guide to installing and configuring the SWAN wave model for coastal applications. Includes practical tips for grid generation and boundary conditions."
pubDate: 2024-10-20
category: "modelling"
tags:
  - "SWAN"
  - "wave modeling"
  - "coastal engineering"
  - "numerical methods"
featured: true
---

## Introduction to SWAN

SWAN (Simulating WAves Nearshore) is a third-generation wave model designed for coastal regions and inland waters. This tutorial will guide you through:

1. Installation and compilation
2. Grid setup
3. Input file configuration
4. Running simulations

## Installation

### Prerequisites

```bash
# Ubuntu/Debian
sudo apt-get install gfortran make

# macOS with Homebrew
brew install gcc make
```

### Compiling SWAN

```bash
# Download SWAN
wget http://swanmodel.sourceforge.net/download/swan4145.tar.gz
tar -xzf swan4145.tar.gz
cd swan4145

# Compile
make config
make install
```

## Creating a Computational Grid

SWAN uses a structured grid. Here's a simple example:

```fortran
CGRID REGULAR -90.0 25.0 0.0 1000 800 0.01 0.01
```

This creates a grid:
- Starting at longitude -90°, latitude 25°
- With 1000 x 800 points
- Grid spacing of 0.01° (~1 km)

## Input File Structure

A basic SWAN input file (`INPUT`) looks like:

```fortran
PROJ 'CoastalWaves' 'v1.0'

SET LEVEL=0.0
MODE NONSTATIONARY

CGRID REGULAR -90.0 25.0 0.0 1000 800 0.01 0.01

INPGRID BOTTOM REGULAR -90.0 25.0 0.0 1000 800 0.01 0.01
READINP BOTTOM 1.0 'bathymetry.dat' 4 0 FREE

BOUND SHAPESPEC JONSWAP PEAK DSPR DEGREES
BOUNDSPEC SIDE WEST CCW CON PAR 2.0 10.0 270.0 20.0

GEN3 KOMEN
FRICTION JONSWAP
BREAKING

OUTPUT OPTIONS TABLE
BLOCK 'COMPGRID' NOHEADER 'output.dat' LAY 3 HSIGN TPS DIR DSPR

COMPUTE 20240101.000000 1 HR 20240131.000000
STOP
```

## Boundary Conditions

### Wind Input

```fortran
INPGRID WIND REGULAR -90.0 25.0 0.0 100 80 0.1 0.1
READINP WIND 1.0 'wind.dat' 4 0 FREE
```

### Wave Spectrum at Boundaries

```fortran
BOUNDSPEC SEGMENT IJ 1 1 1 800 VARIABLE FILE 'boundary_spectra.txt'
```

## Running SWAN

```bash
# Serial execution
./swan.exe INPUT

# Parallel execution with MPI
mpirun -np 4 ./swanrun -input INPUT
```

## Post-Processing

Use Python to analyze SWAN output:

```python
import numpy as np
import matplotlib.pyplot as plt

# Read SWAN table output
data = np.loadtxt('output.dat', skiprows=5)
x, y, hs, tp, dir, dspr = data.T

# Plot significant wave height
plt.scatter(x, y, c=hs, cmap='viridis')
plt.colorbar(label='Hs (m)')
plt.title('Significant Wave Height')
plt.show()
```

## Tips and Best Practices

1. **Grid Resolution**: Balance between accuracy and computational cost
2. **Time Steps**: Use adaptive time stepping for efficiency
3. **Validation**: Always validate against observations
4. **Boundary Conditions**: Ensure proper offshore boundary placement

## Common Issues

### Convergence Problems

If SWAN doesn't converge:
- Check bathymetry for unrealistic values
- Verify boundary conditions
- Adjust numerical parameters (ACCUR, DREL)

### High Computational Time

- Use coarser grids for initial testing
- Enable parallel processing
- Optimize output frequency

## Conclusion

SWAN is a powerful tool for coastal wave modeling. With proper setup and validation, it can provide accurate wave predictions for engineering and research applications.

## Resources

- [SWAN User Manual](http://swanmodel.sourceforge.net/online_doc/swanuse/swanuse.html)
- [SWAN Examples](http://swanmodel.sourceforge.net/examples/examples.htm)
