---
title: "Getting Started with WAVEWATCH III"
description: "A beginner's guide to understanding and using the WAVEWATCH III spectral wave model for ocean wave forecasting and hindcasting."
pubDate: 2024-09-10
category: "modelling"
tags: ["WAVEWATCH III", "wave modeling", "ocean forecasting", "beginner guide"]
featured: false
---

## Introduction

WAVEWATCH III (WW3) is a third-generation spectral wave model developed by NOAA/NCEP. It's widely used for operational wave forecasting and research applications worldwide. This tutorial will help you get started with WW3.

## What is WAVEWATCH III?

WAVEWATCH III solves the spectral action density balance equation for wave propagation. It includes:

- Wind input and dissipation
- Nonlinear wave-wave interactions
- Bottom friction and depth-induced breaking
- Currents and ice coverage effects

## Installation

### Prerequisites

You'll need:
- Linux/Unix environment (or WSL on Windows)
- Fortran compiler (gfortran or ifort)
- NetCDF libraries

### Building WW3

```bash
# Download WW3
wget https://github.com/NOAA-EMC/WW3/archive/refs/tags/7.14.tar.gz
tar -xzf 7.14.tar.gz
cd WW3-7.14

# Set environment
export WWATCH3_DIR=$(pwd)
export WWATCH3_ENV=environment_file

# Compile
./install_ww3_tar

# Select your compiler and options
# Follow the prompts to configure your build
```

## Basic Workflow

### 1. Grid Definition

Create a grid definition file describing your domain:

```fortran
$ WAVEWATCH III Grid preprocessor input file
$
$ Grid name
  'GLOBAL'
$
$ Frequency grid (start, increment, number)
  0.04118 1.1 25
$
$ Directional grid (number of directions)
  24
$
$ Spatial grid
  'RECT' T 'NONE'
  1440 721
  -0.125 -90.0 0.25 0.25
```

This defines:
- Frequency range from 0.041 Hz with 25 bins
- 24 directional bins (15° resolution)
- Global grid at 0.25° resolution

### 2. Input Files

WW3 requires:
- **Wind fields**: Wind speed and direction
- **Ice coverage**: (optional) Sea ice concentration
- **Currents**: (optional) Ocean currents
- **Bathymetry**: Water depths

### 3. Running a Simulation

```bash
# Preprocess grid
ww3_grid

# Generate initial conditions
ww3_strt

# Prepare forcing fields
ww3_prep

# Run the model
ww3_shel

# Extract output
ww3_ounf  # For NetCDF output
ww3_outp  # For point output
```

## Understanding WW3 Output

WW3 produces various output types:

### Gridded Fields

```python
import xarray as xr
import matplotlib.pyplot as plt

# Load WW3 output
ds = xr.open_dataset('ww3_output.nc')

# Plot significant wave height
ds['hs'].isel(time=0).plot()
plt.title('Significant Wave Height')
plt.show()
```

### Key Variables

- `hs`: Significant wave height (m)
- `tp`: Peak period (s)
- `dir`: Mean wave direction (deg)
- `spr`: Directional spread (deg)
- `wnd`: Wind speed (m/s)

## Wave Spectrum

The wave spectrum $S(f, \theta)$ describes energy distribution:

$$
H_s = 4\sqrt{\int_0^\infty \int_0^{2\pi} S(f,\theta) \, d\theta \, df}
$$

## Common Applications

### 1. Operational Forecasting

WW3 is used by:
- NOAA for US wave forecasts
- ECMWF for global wave predictions
- National weather services worldwide

### 2. Climate Studies

- Long-term wave climate analysis
- Extreme wave statistics
- Climate change impacts

### 3. Engineering Design

- Offshore structure design
- Coastal protection planning
- Ship routing optimization

## Tips for Success

1. **Start Simple**: Begin with a small domain and coarse resolution
2. **Validate Output**: Compare with observations (buoys, satellites)
3. **Check Convergence**: Ensure adequate spatial and temporal resolution
4. **Use Appropriate Physics**: Select source terms suitable for your application

## Common Issues

### High CPU Time

- Reduce grid resolution
- Use implicit schemes
- Enable MPI parallelization

### Numerical Instabilities

- Reduce time step
- Check boundary conditions
- Verify input data quality

## Resources

- [WW3 User Manual](https://github.com/NOAA-EMC/WW3/wiki)
- [WW3 GitHub Repository](https://github.com/NOAA-EMC/WW3)
- [WAVEWATCH III Forum](https://groups.google.com/g/ww3)

## Conclusion

WAVEWATCH III is a powerful tool for wave modeling. With proper setup and validation, it can provide accurate wave predictions for various applications. Start with simple cases and gradually increase complexity as you gain experience.

## Next Steps

In the next tutorial, we'll cover:
- Advanced grid nesting techniques
- Coupling WW3 with ocean circulation models
- Post-processing and visualization tools
