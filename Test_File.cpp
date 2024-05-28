#include "FDM.h"

using namespace std;

int main()
{
	// Initial Variational method approach, where we graphed the energy vs parameter in question (x_0)
	/* {
		std::string filename = "Variational_methods.txt";
		double X_0 = 0.5;
		// Open a file to write the results
		std::ofstream outFile(filename);
		if (!outFile.is_open())
		{
			std::cerr << "Unable to open file: " << filename << std::endl;
			return 0;
		}
		outFile << "X_0 Value" << "\t" << "E_0 Value" << "\n";
		for (int i = 0; i <= 1500; ++i)
		{
			double current_X_0 = X_0 + i * 0.001;
			std::unique_ptr<FDM> obj = std::make_unique<FDM>();
			(*obj).discretize_x();
			(*obj).set_x0(current_X_0);
			(*obj).evaluate_psi();
			(*obj).evaluate_kinetic();
			(*obj).evaluate_potential();
			(*obj).evaluate_hamiltonian();
			//cout << "X_0 : " << current_X_0 << " E_0 : " << (*obj).evaluate_e0() << "\n";
			outFile << current_X_0 << "\t" << (*obj).evaluate_e0() << "\n";
		}
		outFile.close();
	}*/

	/* {
		//Second go where we got the same result, but in less time and using the Gradient Descent Approach.
		double epsilon = 0.001;
		double learning_rate = 0.1;
		double tolerance = 0.01;
		double gradient, x_new;
		double current_x0 = 0.5;
		int i = 0;
		std::string filename_new = "Gradient_Descent.txt";
		// Open a file to write the results
		std::ofstream outFile(filename_new);
		if (!outFile.is_open())
		{
			std::cerr << "Unable to open file: " << filename_new << std::endl;
			return 0;
		}
		outFile << "Learning Rate : \t" << learning_rate << "\t Epsilon : \t" << epsilon << "\t Initial Guess : \t" << current_x0 << "\n";
		outFile << "Iter" << "\t" << "Gradient \t" << "X" << "\t" << "E" "\n";
		std::vector <double> energies;
		std::vector <double> x_vals;
		std::unique_ptr<FDM> obj1 = std::make_unique<FDM>();
		(*obj1).discretize_x();
		(*obj1).set_x0(current_x0);
		do
		{
			(*obj1).evaluate_psi();
			(*obj1).evaluate_kinetic();
			(*obj1).evaluate_potential();
			(*obj1).evaluate_hamiltonian();
			long double e_original = (*obj1).evaluate_e0();
			x_vals.push_back((*obj1).get_x0());
			energies.push_back(e_original);
			long double e_new = (*obj1).gradient_descent(epsilon);
			gradient = (e_new - e_original) / epsilon;
			x_new = (*obj1).get_x0() - gradient * learning_rate;
			cout << "Gradient : " << gradient << " X Old : " << (*obj1).get_x0() << " X_new : " << x_new << "\n";
			(*obj1).set_x0(x_new);
			(*obj1).clearvecs();
			
			if (i>0)
			{
				if (x_vals[i] - x_vals[i - 1] < 0.0001)
				{
					outFile << i << "\t" << gradient << "\t" << x_vals[i] << "\t" << energies[i] << "\n";
					cout << "\nFinal E_0 using Gradient Descent and Variational Methods is :" << x_vals[i];
					return EXIT_SUCCESS;
				}
			}
			outFile << i << "\t" << gradient << "\t" <<   x_vals[i] << "\t" << energies[i] << "\n";
			i = i + 1;
			
		} while (abs(x_new - (*obj1).get_x0()) < tolerance);
		outFile.close();
	}*/

	// Our perturbation approach, where we subtly changed Psi, got the energy, and only updated if we got closer to 
	// the true ground state energy, i.e. New_Energy < Old_Energy
	{
		FDM obj(1, 1, 200, 0.7, 5);
		obj.discretize_x();
		obj.set_x0(2.85);
		obj.evaluate_trial_wavefunction();
		obj.normalise_wave_function();
		// Open a file to write the results
		//std::ofstream outFile("Wave_Function_Energy_N_200.txt");
		/*outFile << "Iter \t" << "Energy \n";
		if (!outFile.is_open())
		{
			std::cerr << "Unable to open file: " << "Wave_Function_Energy.txt" << std::endl;
			return 0;
		}*/
		obj.wave_func_to_text(0);
		for (int i = 0; i < 1E6; ++i)
		{
			obj.clear_leonard_jones_vecs();
			obj.evaluate_kinetic();
			obj.evaluate_leonnard_jones_potential();
			obj.evaluate_hamiltonian();
			//long double E_0 = obj.evaluate_e0();
			obj.update_psi();
			if (i % 10000 == 0)
			{
				cout << obj.get_e0() << "\t" << i << "\n";
			//	outFile << i << "\t" << obj.get_e0() << "\n";
			}
			obj.clear_leonard_jones_vecs();
		}
		obj.wave_func_to_text(1000000);
		//outFile.close();

	}
	return 0;
}