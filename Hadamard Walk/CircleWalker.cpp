#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "Complex.h"
inline int md(int i, int n)
{
	return (i % n + n) % n;
}
int main(int argc, char* argv[])
{
	// Number of step
	unsigned int nstep = 100;
	// Number of location
	unsigned int grid = 200;

	// We can lauch the program with arguments to choose the number of step
	switch (argc)
	{
		case 2:
			nstep = atoi(argv[1]);
			grid = 2 * nstep;
			break;

		default:
			break;
	}

	std::cout << "Computation of a circle Hadamard quantum walk" << std::endl;
	std::cout << "with Tmax = " + std::to_string(nstep) << std::endl;
	std::cout << "and N = " + std::to_string(grid) << std::endl;

	// Definition of the initial state's up composant
	complex initalUp = complex(1.0 / std::sqrt(2.0));

	// Definition of the initial state's down composant
	complex initalDown = complex(0.0,1.0 / std::sqrt(2.0));

	// Creation of the wave function ∈ \mathbb{Z}^{N}\otimes\mathbb{C}^{2}
	complex** state = new complex * [grid];

	// Creation of a temp wave function usefull to apply operator
	complex** tempstate = new complex * [grid];

	for (unsigned int i = 0; i < grid; i++)
	{
		state[i] = new complex[2];
		tempstate[i] = new complex[2];

		for (unsigned int j = 0; j < 2; j++) {

			// Setting the inital state
			state[i][j] = complex(0.0, 0.0);
			tempstate[i][j] = complex(0.0, 0.0);

		}
	}
	// Setting the inital state
	state[grid / 2][0] = initalUp;
	state[grid / 2][1] = initalDown;

	std::cout << "Walking..." << std::endl;
	for (unsigned int n = 0; n < nstep; n++) {
		std::cout << "Step: " + std::to_string(n + 1) << std::endl;
		// Applying the Hardamard coin
		for (unsigned int i = 0; i < grid; i++) {
			// We put in memory the wave funcition's value at the position i
			tempstate[i][0] = state[i][0];
			tempstate[i][1] = state[i][1];

			// We compute the matrix product between the coin and the wave function at the position i
			state[i][0] = 1.0 / std::sqrt(2.0) * (tempstate[i][0] + tempstate[i][1]);
			state[i][1] = 1.0 / std::sqrt(2.0) * (tempstate[i][0] - tempstate[i][1]);

			// We refresh the temp memory to apply the sifht operator correctly
			tempstate[i][0] = state[i][0];
			tempstate[i][1] = state[i][1];
		}

		// Applying the shift operator
		for (unsigned int i = 0; i < grid; i++) {

			// the up composant move form i to i+1
			state[i][0] = tempstate[md(i - 1, grid)][0];

			// the down composant move form i to i-1
			state[i][1] = tempstate[md(i + 1, grid)][1];
		}
	}

	// We release the temp memory
	for (unsigned int i = 0; i < grid; i++) {
		delete[] tempstate[i];
	}
	delete[] tempstate;


	// We write the probability distribution on a file named p.txt
	std::string filename;
	std::ofstream outFile;

	std::cout << "Saving p file's..." << std::endl;
	filename = "p.txt";
	outFile.open(filename.c_str());
	for (unsigned int i = 0; i < grid; i++)
		outFile << state[i][0].abs2() + state[i][1].abs2() << " ";
	outFile.close();

	// We release the wave function memory
	for (unsigned int i = 0; i < grid; i++) {
		delete[] state[i];
	}
	delete[] state;
}