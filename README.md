# Variational-Method-and-Gradient-Descent
This project explores variational methods and optimization techniques to find the ground state energy of quantum systems, starting with the Quantum Harmonic Oscillator (QHO) and extending to the Lennard-Jones potential for non-attracting diatomic systems.

## Contents
 - $\textbf{C++ Files:}$ Implementation of the variational method and optimization algorithms.
 - $\textbf{Python Jupyter Notebook:}$ Used for plotting and analyzing results.
 - $\textbf{PDF Document:}$ Outlines the theory, reasoning, results, discussion, and conclusions.

 # Project Overview
 ## Quantum Harmonic Oscillator
 ### Variational Method
 - $\textbf{Trial Wavefunction (\psi):}$ A trial wavefunction $\psi$ is used to calculate the energy expectation value of the Hamiltonian.
 - $\textbf{Parameter Optimization:}$ By iterating a parameter of $\psi$ over a set range, we determine the optimal value at the minimum energy. According to the Variational Method, this value is the closest to the true analytic value, which is known for the QHO.

 ### Gradient Descent Optimization
 - $\textbf{Gradient Descent:}$ This method is used to find the same optimized parameter value as obtained through the iterative approach.
 - $\textbf{Efficiency:}$ Gradient Descent converges more quickly and offers an interesting alternative to simple parameter iteration.

 ## Lennard-Jones Potential
 ### New Scenario and Optimization
 - $\textbf{New Trial Wavefunction (\psi):}$ A new trial wavefunction is used for the Lennard-Jones potential.
 - $\textbf{Random Perturbations:}$ $\psi$ is perturbed at random points, and if the energy reduces, the new $\psi$ vector is retained.
 - $\textbf{Convergence:}$ This process is iterated until a convergent minimum energy value is found for the trial $\psi$. By the Variational Principle, this is as close to the analytic value as possible.

 # Dependencies
### C++:
 - Eigen library for linear algebra operations.
### Python:
 - Matplotlib for plotting graphs.
 - Jupyter Notebook for interactive analysis.

# Conclusion
This project demonstrates the application of the Variational Method and Gradient Descent optimization to solve and analyze quantum mechanical systems. The methods are validated against known analytic values for the Quantum Harmonic Oscillator and extended to more complex potentials like the Lennard-Jones potential.
