---
title: "Setting Up Julia with GPU Support on Ruche HPC"
description: "A practical guide to installing Julia with CUDA support on the Ruche HPC cluster and running GPU-accelerated Oceananigans simulations."
pubDate: 2025-11-06
category: "tutorial"
tags:
  - "Julia"
  - "GPU"
  - "HPC"
  - "CUDA"
  - "Oceananigans"
featured: false
---

## Introduction

To set up Julia with GPU support on the Ruche HPC cluster, I needed to go through several steps to get everything working properly. Here's how I did it, broken down into manageable chunks.

## Requesting GPU Resources

```bash
salloc --nodes=1 --time=00:30:00 -p gpu --gres=gpu:2
```

First, I requested access to GPU nodes through SLURM. This allocates 2 GPUs for 30 minutes, which is plenty of time for installation and testing.

## Loading CUDA and Setting Environment

```bash
module load cuda/12.2.1/gcc-11.2.0

export JULIA_VERSION="1.10.10"
export CUDA_VERSION="12.2"
export HOME=/gpfs/users/dspelaez
export CUDA_PATH=$CUDA_HOME
```

I loaded the CUDA module available on the cluster and set up some essential environment variables. It's important to specify both Julia and CUDA versions to ensure compatibility.

## Configuring Library Paths

```bash
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=$CUDA_HOME/lib64:$LIBRARY_PATH
```

These paths tell Julia where to find the CUDA libraries at runtime. Without this, Julia won't be able to locate the GPU drivers.

## Setting Up Working Directories

```bash
export JULIA_BIN_PATH="$HOME/.local/bin/"
export JULIA_DEPOT_PATH="$HOME/.julia"
mkdir -p "$JULIA_DEPOT_PATH"
```

I created a dedicated directory for Julia packages in my home folder. This keeps everything organised and separate from other software.

Note: Some people prefer moving this to the scratch directory (e.g., `export JULIA_DEPOT_PATH="$SCRATCH/.julia"`) since Julia packages can take up quite a bit of space and home directories often have storage quotas.

## Downloading and Installing Julia

```bash
mkdir -p "$HOME/.local/julia"
cd "$HOME/.local/julia"
wget "https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" -O julia.tar.gz
tar -xzf julia.tar.gz
rm julia.tar.gz
```

I downloaded Julia directly from the official distribution site and extracted it to my local directory. This gives me full control over which version I'm using.

## Creating a Symbolic Link

```bash
mkdir -p "$HOME/.local/bin"
ln -sf "$HOME/.local/julia/julia-${JULIA_VERSION}/bin/julia" "$HOME/.local/bin/julia"
```

Creating a symlink makes it easy to run Julia from anywhere without typing the full path every time.

## Installing Julia Packages

```bash
julia -e 'using Pkg; Pkg.add("Oceananigans")'
julia -e 'using Pkg; Pkg.add("CUDA")'
julia -e 'using CUDA; CUDA.set_runtime_version!(v"12.2")'
julia -e 'using Pkg; Pkg.build("Oceananigans")'
```

This is where the magic happens. I installed Oceananigans first, then CUDA.jl, and explicitly set the CUDA runtime version to match what's loaded on the cluster. Finally, I rebuilt Oceananigans to make sure it's properly configured for GPU use. Why this way? I don't know but it worked 🤷

## Logging into the GPU Node

Now we need to actually log in to the allocated GPU node to test everything:

```bash
# After the salloc command grants you resources, SSH into the node
# Check which node you got assigned (shown in salloc output)
ssh gpu-node-XX  # Replace XX with your assigned node number
```

Once you're on the GPU node, load the CUDA module again and set up your environment paths as before.

## Verifying the Installation

```bash
nvidia-smi
julia -e 'using CUDA; CUDA.versioninfo()'
```

I checked that the GPUs were visible and that Julia could see them correctly. The `nvidia-smi` command shows the GPU status, whilst `CUDA.versioninfo()` confirms Julia's CUDA configuration.

## Testing with a Simulation

```julia
using Oceananigans, CUDA

grid = RectilinearGrid(GPU(), size=(32, 32, 32), extent=(100, 100, 100))
model = NonhydrostaticModel(; grid, tracers=:b, buoyancy=BuoyancyTracer())
simulation = Simulation(model, Δt=1.0, stop_time=10.0)
run!(simulation)

println("Setup complete - simulation ran successfully!")
```

Finally, I ran a quick test simulation on the GPU. Creating a simple 32×32×32 grid with a buoyancy tracer is enough to verify that everything's working properly. If this runs without errors, you're good to go!

## Conclusion

Setting up Julia with GPU support on an HPC cluster requires careful attention to environment variables and library paths. The key is ensuring that the CUDA version loaded on the cluster matches what you tell Julia to use. Once configured, GPU-accelerated simulations with Oceananigans run brilliantly on the cluster.
