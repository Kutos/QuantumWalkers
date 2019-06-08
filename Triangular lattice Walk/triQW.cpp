#include <iostream> // Console's interactions
#include <string>   // Character chains
#include <math.h>   // Basics mathematical functions
#include <cstdlib>  // System's interactions
#include <stdio.h>  // I/O
#include <fstream>  // I/O
#include <time.h>   // Chronometre

#include "Triangle.h"
#include "Complex.h"

long double const pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513;

// That positive modulo function 
inline int md(int i, int n) { return (i % n + n) % n; }


int main(int argc, char* argv[])
{
	// We just add a chronometre to know how the computation time evolve with the step number
#pragma region Chrono
	float tm;
	clock_t t1, t2;
	t1 = clock();
#pragma endregion

#pragma region Matrix definition

	std::cout << "Matrix creation" << std::endl;
	long double al = 5.0;
	al = std::sqrt(al) / 3.0;
	al = std::acos(al);

	complex S[2][2];
	complex U1[3][2][2];
	complex W1[2][2];

	// S definition
	S[0][0] = exp(complex(0.0, 2.0 * pi / 3.0));
	S[0][1] = 0.0;
	S[1][0] = 0.0;
	S[1][1] = 1.0;

	// U_0 definition
	U1[0][0][0] = std::cos(al / 2.0);
	U1[0][0][1] = -std::sin(al / 2.0);
	U1[0][1][0] = std::sin(al / 2.0);
	U1[0][1][1] = std::cos(al / 2.0);

	// U_1 calculus => matrix product between U_0 and S
	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			U1[1][i][j] = 0.0;

			for (unsigned int k = 0; k < 2; k++)
			{
				U1[1][i][j] += U1[0][i][k] * S[k][j];
			}
		}
	}

	// U_2 calculus => matrix product between U_1 and S
	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			U1[2][i][j] = 0;

			for (unsigned int k = 0; k < 2; k++)
			{
				U1[2][i][j] += U1[1][i][k] * S[k][j];
			}
		}
	}

	// W calculs => matrix product between U_0, S and dagger(U_0)
	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			W1[i][j] = 0;

			for (unsigned int k = 0; k < 2; k++)
			{
				W1[i][j] += U1[1][i][k] * U1[0][j][k].conjugate();
			}
		}
	}

#pragma endregion

	// Number of step
	int nstep = 200;
	// Number of triangle in row and in column
	int grid = 2 * nstep + 3;

	// We can launch the program with arguments to choose the number of step
	switch (argc)
	{
	case 2:
		nstep = atoi(argv[1]);
		grid = 2 * nstep + 3;

	case 3:
		nstep = atoi(argv[1]);
		grid = 2 * nstep + 3;
	}
	std::cout << "Nstep :" + std::to_string(nstep) << std::endl;
	std::cout << "grid :" + std::to_string(grid) << std::endl;

#pragma region Initialisation

	std::cout << "Creation of the initial state" << std::endl;

	// Creation of the wave function ∈ \mathbb{Z}^{N*N}\otimes\mathbb{C}^{2}
	Triangle** T = new Triangle * [grid];

	//Creation of a temp wave function usefull to apply operator
	Triangle** cT = new Triangle * [grid];

	for (unsigned int i = 0; i < grid; i++)
	{
		T[i] = new Triangle[grid];
		cT[i] = new Triangle[grid];
		for (unsigned int j = 0; j < grid; j++)
		{
			T[i][j].setNull();
			cT[i][j].setNull();
		}
	}
	// That choose where is our initial state
	int x0 = grid / 2;
	// Setting the inital state

	T[x0][x0].edge[0][0] = 1 / sqrt(6.0);
	T[x0][x0].edge[0][1] = 1 / sqrt(6.0);

	T[x0][x0].edge[1][0] = 1 / sqrt(6.0);
	T[x0][x0].edge[1][1] = 1 / sqrt(6.0);

	T[x0][x0].edge[2][0] = 1 / sqrt(6.0);
	T[x0][x0].edge[2][1] = 1 / sqrt(6.0);

#pragma endregion

#pragma region Encoding with U
	std::cout << "Encoding the wave function..." << std::endl;
	// Compute the product between U0 and the wave function
	for (int i = 0; i < grid; i++)
	{
		for (int j = 0; j < grid; j++)
		{
			cT[i][j].setNull();
			for (unsigned int k = 0; k < 3; k++)
			{
				for (unsigned int o = 0; o < 2; o++)
				{

					for (unsigned int l = 0; l < 2; l++)
					{
						cT[i][j].edge[k][o] += U1[0][o][l] * T[i][j].edge[k][l];
					}

				}
			}
		}
	}

#pragma endregion

#pragma region The Walk
	int ni, nj;
	for (int n = 0; n < nstep; n++)
	{

		std::cout << "Step " + std::to_string(n) << std::endl;

		// Applying the rotation operator
		for (int i = 0; i < grid; i++)
		{
			for (int j = 0; j < grid; j++)
			{
				T[i][j].setNull();
				for (unsigned int k = 0; k < 3; k++) {
					for (unsigned int cp = 0; cp < 2; cp++) {
						ni = i + (cp==1) * (-((j + 1) % 2) * (k == 0 || k == 2) + (k == 1));
						ni = md(ni, grid);

						nj = j + (cp==1) * ((k == 0) - (k == 2));
						nj = md(nj, grid);

						for (int l = 0; l < 2; l++)
						{
							T[i][j].edge[k][cp] += W1[cp][l] * cT[ni][nj].edge[md(k - 1, 3)][l];
						}
					}
				}
			}
		}
		for (int i = 0; i < grid; i++)
		{
			for (int j = 0; j < grid; j++)
			{
				cT[i][j] = T[i][j];
			}
		}
	}
#pragma endregion

#pragma region Decoding with U
	std::cout << "Decoding the wave function..." << std::endl;
	// Compute the product between dagger(U0) and the wave function
	for (int i = 0; i < grid; i++)
	{
		for (int j = 0; j < grid; j++)
		{
			T[i][j].setNull();
			for (unsigned int k = 0; k < 3; k++) {
				for (unsigned int o = 0; o < 2; o++)
				{

					for (unsigned int l = 0; l < 2; l++)
					{
						T[i][j].edge[k][o] += U1[0][l][o].conjugate() * cT[i][j].edge[k][l];
					}
				}
			}
		}
		delete[] cT[i];
	}
	delete[] cT;
#pragma endregion

#pragma region Out

	std::cout << "Writing outfiles..." << std::endl;

	std::ofstream outFile;
	std::string filename = "./data/x.txt";
	outFile.open(filename.c_str());
	for (int i = 0; i < grid; i++)
	{
		for (int j = 0; j < grid; j++)
		{

			for (unsigned int k = 0; k < 3; k++)
			{
				outFile << 2.0 * i - 0.5 * (k == 0) + 0.5 * (k == 1) + 0.0 * (k == 2) + j % 2 << " ";
			}
		}
	}
	outFile.close();

	filename = "./data/y.txt";
	outFile.open(filename.c_str());
	for (int i = 0; i < grid; i++)
	{
		for (int j = 0; j < grid; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				outFile << std::sqrt(3.0) * (j - 1.0 / 2.0 * (k == 2)) << " ";
			}
		}
	}
	outFile.close();

	double P = 0;
	filename = "./data/p.txt";
	outFile.open(filename.c_str());
	for (int i = 0; i < grid; i++)
	{
		for (int j = 0; j < grid; j++)
		{
			for (unsigned int k = 0; k < 3; k++) {
				outFile << T[i][j].edge[k][0].abs2() + T[i][j].edge[k][1].abs2() << " ";
			}
		}
	}
	outFile.close();

#pragma endregion

#pragma region Chrono
	t2 = clock();
	tm = (float)(t2 - t1) / CLOCKS_PER_SEC;
	printf("temps = %f s\n", tm);
#pragma endregion

	return 0;
}