# spectral_solver
Solver for nonlinear diffusion equations on disk, using pseudo spectral method with Chebyshev-Fourier series.

It solves N equations formed as:

\partial_t f_i=\sum_{k=1}^N\nabla\cdot(h_{ik}\nabla f_k)+g_i(f_1,f_2,\cdots,f_N)

Now works perfect for linear equations, but still have some instability issues with nonlinear equations. Support for long double and MPI(or openMP) parallel computing.

Comments and documents need to be finished later.
