#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <random>
#include <string>
#include <cmath>

class FDM
{
private:
	int N;
	double hbar;
	double m;
	double deltax;
	double e_0;
	double xmin, xmax, x_0; 
	double epsilon;
	const double pi = 3.14159;
	std::vector <double> x_grid;
	std::vector <double> wave_function_gradient_descent;
	std::vector <double> kinetic_vec;
	std::vector <double> kinetic_gradient_descent;
	std::vector <double> potential_vec;
	std::vector <double> potential_gradient_descent;
	std::vector <double> hamiltonian;
	std::vector <double> hamiltonian_gradient_descent;
public:
	FDM();
	FDM(double m, double hbar, int N, double xmin, double xmax);
	std::vector <double> wave_function;
	void set_x0(double x0);
	double get_x0();
	double get_e0();
	double wavefunction(double x, double x0);
	double trial_wavefunction(double x, double x0);
	void discretize_x();
	void evaluate_psi();
	void evaluate_trial_wavefunction();
	double potential(double x, double k, double m, double x_0);
	double leonnard_jones_potential(double sigma, double x, double epsilon);
	void evaluate_potential();
	void evaluate_leonnard_jones_potential();
	void evaluate_kinetic();
	void evaluate_hamiltonian(); 
	long double evaluate_e0();
	long double gradient_descent(double epsilon);
	void clearvecs();
	void clear_leonard_jones_vecs();
	void update_psi();
	void normalise_wave_function();
	void wave_func_to_text(int i);
};

