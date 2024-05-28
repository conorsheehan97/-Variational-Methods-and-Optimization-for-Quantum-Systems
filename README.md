# Variational-Method-and-Gradient-Descent
This project implements the Variational Method, with Gradient Descent in order to optimise hyperparameters for Trial WaveFunctions, estimate the Ground State Energy of two Trial Wavefunctions, and to iteratively improve a Trial Wavefunction.

This project initially tackles the well known potential the Quantum Harmonic Oscillator. We have a Trial Wavefunction Psi which we use to calculate the Energy Expectation Value of the Resultant Hamiltonian. We then iterate a parameter of Psi over a set range, and determine the optimal value at the minimum energy. This value, according to the Variational MEthod is the closest to the true analytic value (which for this potential is also known). 

We then use Gradient Descent optimisation in order to calculate the same optimised parameter value as we had just previously obtained. This converged quicker, while also being an interesting approach to take rather than simply just incrementing over a range. 

Finally, we changed scenario. We investigated a new trial Psi, and a new potential. The Leonnard-Jones potential for non-attracting diatomic systems. We perturbed Psi at random points, checked if the energy reduced, and if so kept the new Psi Vector. This approach was iterated until a convergent minimum energy value was found for our trial Psi. This, by the Variational principle was as close to the analytic value as possible. 

All above code was implemented in C++. Any subsequent graphs were produced using Pythons MAtplotlib library. The PDF attached to this project outlines the theory, reasoning, results, discussion and conclusions from this work.  
