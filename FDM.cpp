#include "FDM.h"

FDM::FDM()
{
	N = 250;
	xmin = -4;
	xmax = 4;
	hbar = m = 1;
} //Default Constructor
FDM::FDM(double m, double hbar, int N, double xmin, double xmax)
{
	this->m = m;
	this->hbar = hbar;
	this->N = N;
	this->xmin = xmin;
	this->xmax = xmax;/*
	bool hbar_input = false;
	std::cout << "\n Please enter the Particle Mass: \n";
	std::cin >> m;
	std::cout << "\n Please enter the minimum x-value: \n";
	std::cin >> xmin;
	std::cout << "\n Please enter the maximum x-value: \n";
	std::cin >> xmax;
	std::cout << "\n Please enter the number of grid points: \n";
	std::cin >> N;
	do {
		char choice;
		std::cout << "\nWould you like the literature hbar value (press l), or default value of 1 (press d)?\n";
		std::cin >> choice;
		switch (choice)
		{
		case 'd':
		case 'D':
			std::cout << "\nGood idea\n";
			hbar = 1;
			break;
		case 'l':
		case 'L':
			std::cout << "\nFeisty one you are.\n";
			hbar = 1.054E-34;
			break;
		default:
			break;
		}
	} while (hbar_input = false);*/
} // Parameterized Constructor
double FDM::wavefunction(double x, double x0)
{
	const double pi = 3.14159;
	double norm_constant = pow(2/(pi*pow(x_0,2)),0.25);
	return norm_constant * exp(-pow(x,2)/pow(x_0,2));
} // Our intial trial Wavefunction
double FDM::trial_wavefunction(double x, double x_0)
{
	double prefix = pow(32 / pi*pow(x_0,6), 0.25);
	return prefix * exp((-1* pow((x - x_0), 2)) / 2);
}; // Our second trial Wavefunction
double FDM::potential(double x, double k, double m, double x_0)
{
	// Our good friend the QHO Potential
	double prefix = 0.5 * k* m ;
	return prefix * pow(x, 2) * pow(wavefunction(x,x_0),2);
}  
double FDM::leonnard_jones_potential(double sigma, double x, double epsilon)
{
	return 4 * epsilon * (pow((sigma/x), 12) - pow((sigma/x), 6));
} //Our Leonnard-Jones potential for dimolecular interaction
void FDM::discretize_x()
{
	deltax = (abs(xmax - xmin)) / (N - 1);
	for (int i = 0; i < N; i++)
	{
		double xn = xmin + i * deltax;
		x_grid.push_back(xn);
	}
} // a method to discretize the x axis
void FDM::set_x0(double x0)
{
	x_0 = x0;
} // a simple setter method 
double FDM::get_x0()
{
	return x_0;
} // anopther simple getter method 
void FDM::evaluate_psi()
{
	for (int i = 0; i < x_grid.size(); ++i)
	{
		wave_function.push_back(wavefunction(x_grid[i], x_0));
	}
} // a method to evaluate the wavefunction and store it in a vector 
void FDM::evaluate_trial_wavefunction()
{
	for (int i = 0; i < x_grid.size(); ++i)
	{
		wave_function.push_back(trial_wavefunction(x_grid[i], x_0));
	}
} // another evaluator method. This and the preceding method could have been implement together, however in the purpose of saving time two methods were used
void FDM::evaluate_potential()
{
	for (int i = 0; i < x_grid.size(); ++i)
	{
		potential_vec.push_back(potential(x_grid[i], 1, 1, x_0));
	}
} // method to evaluate our QHO potential 
void FDM::evaluate_leonnard_jones_potential()
{
	for (int i = 0; i < x_grid.size(); ++i)
	{
		potential_vec.push_back(leonnard_jones_potential(1,x_grid[i],10));
	}
} //a method to evaluate our leonnard jones potential 
void FDM::evaluate_kinetic()
{
	double prefix = (-pow(hbar, 2) / 2 * m);
	for (int i = 0; i < x_grid.size(); ++i)
	{
		if (i == 0)
		{
			kinetic_vec.push_back(0);
		}
		else if (i == N - 1)
		{
			kinetic_vec.push_back(0);
		}
		else 
		{
			kinetic_vec.push_back(wave_function[i]*prefix*(wave_function[i - 1] + wave_function[i + 1] - 2 * wave_function[i]) / (pow(deltax, 2)));
		}
	}
} // Evaluate the kinetic energy using FDM approach
void FDM::evaluate_hamiltonian()
{
	for (int i = 0; i < x_grid.size(); ++i)
	{
		double kinetic_plus_potential = kinetic_vec[i] + potential_vec[i];
		hamiltonian.push_back(kinetic_plus_potential*deltax);
	}
} // Create our Hamiltonian Vector
long double FDM::evaluate_e0()
{
	double sum = 0;
	for (size_t i = 0; i < hamiltonian.size(); i++)
	{
		sum += hamiltonian[i];
	}
	e_0 = sum;
	return sum;
} // And finally, get our energy
long double FDM::gradient_descent(double epsilon)
{
	//This method implements Gradient Descent to find the optimal X_0
	// 
	//Assigning Epsilon
	this->epsilon = epsilon;
	x_0 = x_0 + epsilon;
	//Wavefunction Update
	for (int i = 0; i < x_grid.size(); ++i)
	{
		wave_function_gradient_descent.push_back(wavefunction(x_grid[i], x_0 ));
	}
	//Potential Update
	for (int i = 0; i < x_grid.size(); ++i)
	{
		potential_gradient_descent.push_back(potential(x_grid[i], 1, 1, x_0));
	}
	//Kinetic Update
	double prefix = (-pow(hbar, 2) / 2 * m);
	for (int i = 0; i < x_grid.size(); ++i)
	{
		if (i == 0)
		{
			kinetic_gradient_descent.push_back(0);
		}
		else if (i == N - 1)
		{
			kinetic_gradient_descent.push_back(0);
		}
		else
		{
			kinetic_gradient_descent.push_back(wave_function_gradient_descent[i] * prefix * (wave_function_gradient_descent[i - 1] + wave_function_gradient_descent[i + 1] - 2 * wave_function_gradient_descent[i]) / (pow(deltax, 2)));
		}
	}
	//Hamiltonian Update
	for (int i = 0; i < x_grid.size(); ++i)
	{
		double kinetic_plus_potential = kinetic_gradient_descent[i] + potential_gradient_descent[i];
		hamiltonian_gradient_descent.push_back(kinetic_plus_potential * deltax);
	}
	//Final E0 calculation
	double sum = 0;
	for (size_t i = 0; i < hamiltonian_gradient_descent.size(); i++)
	{
		sum += hamiltonian_gradient_descent[i];
	}
	return sum;
}
void FDM::clearvecs()
{
	// We need to clear the vectors after every energy calculation, as our energy requires a summation of elements. 
	// You see the issue here.  These are the Gradient Descent specific Vectors
	hamiltonian.clear();
	hamiltonian_gradient_descent.clear();
	kinetic_vec.clear();
	kinetic_gradient_descent.clear();
	potential_vec.clear();
	potential_gradient_descent.clear();
	wave_function.clear();
	wave_function_gradient_descent.clear();
} 
void FDM::update_psi()
{
	//this method perturbs our trial wavefunction by an amoun delta below
	double delta = 0.00001;
	// We use the Mersenne Twister for our random graph point to offset
	std::random_device rd;
	std::mt19937 gen(rd());
	int range = this->N; // Any point on Psi can be evaluated
	std::uniform_int_distribution<> dis1(0, range-1);
	int num1 = dis1(gen); // The index of Wavefunction to update
	std::uniform_int_distribution<> dis2(0, 1);
	int plus_minus = (dis2(gen) == 0) ? -1 : 1; //Can either increment up/down a small smidgeon
	long double temp_energy = evaluate_e0(); // Calculating the current energy
	std::vector <double> temp_vec = wave_function; // storing the current vector in case we need to swap back to it
	wave_function[num1] = plus_minus*delta + wave_function[num1]; // perturbing
	normalise_wave_function(); // normalising the new wavefunction
	clear_leonard_jones_vecs(); //clearing our kinetic, potential, and hamiltonian vectors
	evaluate_kinetic();
	evaluate_leonnard_jones_potential();
	evaluate_hamiltonian();
	long double new_energy = evaluate_e0(); //getting the new energy
	//normalise_wave_function();
	//implementing the variational method
	if (new_energy > temp_energy)
	{
		wave_function = temp_vec;
	}
	}
void FDM::clear_leonard_jones_vecs()
{
	//Again, we need to clear our leonnard jones vectors after each iteration
	hamiltonian.clear();
	kinetic_vec.clear();
	potential_vec.clear();
	long double sum = 0;
	for (int i = 0; i < wave_function.size(); ++i)
	{
		sum = sum + wave_function[i];
	}
}
double FDM::get_e0()
{ //simple getter method
	return e_0;
}
void FDM::normalise_wave_function()
{
	// method to normalise wave function. Very important or else you start to get wonky results (as I found out)
	double sum=0;
	for (int i = 0; i < x_grid.size(); ++i)
	{
		sum = sum + wave_function[i];
	}
	for (int i = 0; i < x_grid.size(); ++i)
	{
		wave_function[i] = abs(wave_function[i] / sum);
	}
}
void FDM::wave_func_to_text(int i)
{
	//Method to write wavefunctions to text, just to see how he's looking
	std::string filename = "WaveFunction_at_" + std::to_string(i) + "_iters.txt";

	// Open a file to write the results
	std::ofstream outFile(filename);
	outFile << "X \t" << "Value \n";
	if (!outFile.is_open())
	{
		std::cerr << "Unable to open file: " << filename << std::endl;
		return;
	}

	for (int j = 0; j < wave_function.size(); ++j)
	{
		outFile << x_grid[j] << "\t" << wave_function[j] << "\n";
	}
	outFile.close();
}