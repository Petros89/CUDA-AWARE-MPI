# CUDA-AWARE-MPI

This is an example implementation for the CUDA-AWARE-MPI (non-blocking) for
sending data/buffers between different GPUs and nodes without staggering
through host but instead directly GPU-GPU communications.

The numerical analysis is performed for the 3D conduction on a structured grid
using second-order Finite Differences Scheme and first-order explicit in time.
