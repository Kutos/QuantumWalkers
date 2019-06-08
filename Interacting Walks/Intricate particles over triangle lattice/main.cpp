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
inline int md(int i, int n)
{
	return (i % n + n) % n;
}

int main(int argc, char* argv[])
{

#pragma region Chrono
	float tm;
	clock_t t1, t2;
	t1 = clock();
#pragma endregion

#pragma region Matrix definition

	std::cout << "Matrix creation" << std::endl;
	double al = std::acos(std::sqrt(5.0) / 3.0);

	complex U[4][4];
	complex W[4][4];

	complex a, b, c, d;

	a = std::cos(al / 2.0);
	b = -std::sin(al / 2.0);
	c = std::sin(al / 2.0);
	d = std::cos(al / 2.0);

	U[0][0] = a * a;
	U[0][1] = a * b;
	U[1][0] = a * c;
	U[1][1] = a * d;

	U[0][2] = b * a;
	U[0][3] = b * b;
	U[1][2] = b * c;
	U[1][3] = b * d;

	U[2][0] = c * a;
	U[2][1] = c * b;
	U[3][0] = c * c;
	U[3][1] = c * d;

	U[2][2] = d * a;
	U[2][3] = d * b;
	U[3][2] = d * c;
	U[3][3] = d * d;

	a = std::cos(al / 2) * std::cos(al / 2) * exp(complex(0.0, 2.0 * pi / 3.0)) + std::sin(al / 2) * std::sin(al / 2);
	b = std::cos(al / 2) * std::sin(al / 2) * (exp(complex(0.0, 2.0 * pi / 3.0)) - 1.0);
	c = std::cos(al / 2) * std::sin(al / 2) * (exp(complex(0.0, 2.0 * pi / 3.0)) - 1.0);
	d = std::cos(al / 2) * std::cos(al / 2) + std::sin(al / 2) * std::sin(al / 2) * exp(complex(0.0, 2.0 * pi / 3.0));

	std::cout << str_complex(a) << std::endl;
	std::cout << str_complex(b) << std::endl;
	std::cout << str_complex(c) << std::endl;
	std::cout << str_complex(d) << std::endl;

	W[0][0] = a * a;
	W[0][1] = a * b;
	W[1][0] = a * c;
	W[1][1] = a * d;

	W[0][2] = b * a;
	W[0][3] = b * b;
	W[1][2] = b * c;
	W[1][3] = b * d;

	W[2][0] = c * a;
	W[2][1] = c * b;
	W[3][0] = c * c;
	W[3][1] = c * d;

	W[2][2] = d * a;
	W[2][3] = d * b;
	W[3][2] = d * c;
	W[3][3] = d * d;

#pragma endregion

#pragma region Constantes
	int nstep = 20;
	int grid = 2 * nstep + 3;
	double g = pi * 3 / 4;
	switch (argc)
	{
	case 2:
		nstep = atoi(argv[1]);
		grid = 2 * nstep + 3;

	case 3:
		nstep = atoi(argv[1]);
		grid = 2 * nstep + 3;
		if (atof(argv[2]) != 0)
			g = pi / atof(argv[2]);
	}
	std::cout << "Nstep:" + std::to_string(nstep) << std::endl;
	std::cout << "grid:" + std::to_string(grid) << std::endl;
	std::cout << "g: " + std::to_string(g) << std::endl;

	int x_initial = grid / 2;
	complex inter;

#pragma endregion

#pragma region Initialisation

	std::cout << "Creation of the initial state" << std::endl;

	// Allocating memory
	Triangle**** T = new Triangle * **[grid];
	Triangle**** cT = new Triangle * **[grid];

	for (unsigned int i = 0; i < grid; i++)
	{
		T[i] = new Triangle * *[grid];
		cT[i] = new Triangle * *[grid];
		for (unsigned int j = 0; j < grid; j++)
		{
			T[i][j] = new Triangle * [grid];
			cT[i][j] = new Triangle * [grid];
			for (unsigned int u = 0; u < grid; u++)
			{
				T[i][j][u] = new Triangle[grid];
				cT[i][j][u] = new Triangle[grid];
				for (unsigned int v = 0; v < grid; v++)
				{
					T[i][j][u][v] = Triangle();
					cT[i][j][u][v] = Triangle();
				}
			}
		}
	}

	// Setting the initial state
	T[x_initial][x_initial][x_initial][x_initial].edge[0][0][1] = 1 / sqrt(6.0);
	T[x_initial][x_initial][x_initial][x_initial].edge[0][0][2] = 1 / sqrt(6.0);

	T[x_initial][x_initial][x_initial][x_initial].edge[1][1][1] = 1 / sqrt(6.0);
	T[x_initial][x_initial][x_initial][x_initial].edge[1][1][2] = 1 / sqrt(6.0);

	T[x_initial][x_initial][x_initial][x_initial].edge[2][2][1] = 1 / sqrt(6.0);
	T[x_initial][x_initial][x_initial][x_initial].edge[2][2][2] = 1 / sqrt(6.0);

#pragma endregion

#pragma region Encoding with U
	std::cout << "Encoding the wave function..." << std::endl;

	// Encoding the wave function
	for (int i = 0; i < grid; i++)
	{
		for (int j = 0; j < grid; j++)
		{
			for (int u = 0; u < grid; u++)
			{
				for (int v = 0; v < grid; v++)
				{
					cT[i][j][u][v].setNull();
					for (unsigned int k1 = 0; k1 < 3; k1++)
					{
						for (unsigned int k2 = 0; k2 < 3; k2++)
						{
							for (unsigned int o = 0; o < 4; o++)
							{
								for (unsigned int l = 0; l < 4; l++)
								{
									cT[i][j][u][v].edge[k1][k2][o] += U[o][l] * T[i][j][u][v].edge[k1][k2][l];
								}
							}
						}
					}
				}
			}
		}
	}
#pragma endregion

#pragma region The Walk

	int ni, nj, nu, nv;

	for (int n = 0; n < nstep; n++)
	{
		std::cout << "Step " + std::to_string(n) << std::endl;

		// Applying the rotation matrix
		for (int i = 0; i < grid; i++)
		{
			for (int j = 0; j < grid; j++)
			{
				for (int u = 0; u < grid; u++)
				{
					for (int v = 0; v < grid; v++)
					{
						// cT => T
						T[i][j][u][v].setNull();
						for (unsigned int k1 = 0; k1 < 3; k1++)
						{
							for (unsigned int k2 = 0; k2 < 3; k2++)
							{
								inter = (1.0 + (exp(complex(0.0, g)) - 1.0) * (i == u && j == v && k1 == k2));
								for (unsigned int cp = 0; cp < 4; cp++)
								{
									// if k1=0 and (cp=2or3) then ni = i-j%2 and nj = j + 1
									//elif k1=1 and (cp=2or3) then ni = i + 1; nj = j
									//elif k1=2 and (cp=2or3) then ni = i-j%2 and nj = j - 1

									ni = i + (((cp == 2) + (cp == 3)) * (-((j + 1) % 2) * (k1 == 0 || k1 == 2) + (k1 == 1)));
									ni = md(ni, grid);

									nj = j + (((cp == 2) + (cp == 3)) * ((k1 == 0) - (k1 == 2)));
									nj = md(nj, grid);

									nu = u + (((cp == 1) + (cp == 3)) * (-((v + 1) % 2) * (k2 == 0 || k2 == 2) + (k2 == 1)));
									nu = md(nu, grid);

									nv = v + (((cp == 1) + (cp == 3)) * ((k2 == 0) - (k2 == 2)));
									nv = md(nv, grid);

									for (unsigned int l = 0; l < 4; l++)
									{
										T[i][j][u][v].edge[k1][k2][cp] += inter * W[cp][l] * cT[ni][nj][nu][nv].edge[md(k1 - 1, 3)][md(k2 - 1, 3)][l];
									}

								}
							}
						}

					}
				}
			}
		}

		for (int i = 0; i < grid; i++)
		{
			for (int j = 0; j < grid; j++)
			{
				for (int u = 0; u < grid; u++)
				{
					for (int v = 0; v < grid; v++)
					{
						cT[i][j][u][v] = T[i][j][u][v];
					}
				}
			}
		}
	}
#pragma endregion

#pragma region Decoding with U
	std::cout << "Decoding the wave function..." << std::endl;

	// Decoding the wave function
	for (unsigned int i = 0; i < grid; i++)
	{
		for (unsigned int j = 0; j < grid; j++)
		{
			for (unsigned int u = 0; u < grid; u++)
			{
				for (unsigned int v = 0; v < grid; v++)
				{
					T[i][j][u][v].setNull();
					for (unsigned int k1 = 0; k1 < 3; k1++)
					{
						for (unsigned int k2 = 0; k2 < 3; k2++)
						{
							for (unsigned int o = 0; o < 4; o++)
							{
								for (unsigned int l = 0; l < 4; l++)
								{
									T[i][j][u][v].edge[k1][k2][o] += U[l][o].conjugate() * cT[i][j][u][v].edge[k1][k2][l];
								}
							}
						}
					}
				}
			}
		}
	}
#pragma endregion

#pragma region outFiles
	std::cout << "Writing outfiles..." << std::endl;

	std::ofstream outFile;
	std::string filename = "";

	filename = "./data/x.txt";
	outFile.open(filename.c_str());

	for (unsigned int i = 0; i < grid; i++)
	{
		for (unsigned int j = 0; j < grid; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				outFile << (2 * i - 0.5 * (k == 0) + 0.5 * (k == 1) + j % 2) << " ";
			}
		}
	}
	outFile.close();

	filename = "./data/y.txt";
	outFile.open(filename.c_str());

	for (unsigned int i = 0; i < grid; i++)
	{
		for (unsigned int j = 0; j < grid; j++)
		{
			for (unsigned int k = 0; k < 3; k++)
			{
				outFile << std::sqrt(3.0) * (j - 0.5 * (k == 2)) << " ";
			}
		}
	}
	outFile.close();

	filename = "./data/p.txt";
	outFile.open(filename.c_str());
	for (int i = 0; i < grid; i++)
	{
		for (int j = 0; j < grid; j++)
		{
			for (int u = 0; u < grid; u++)
			{
				for (int v = 0; v < grid; v++)
				{
					for (unsigned int k1 = 0; k1 < 3; k1++)
					{
						for (unsigned int k2 = 0; k2 < 3; k2++)
						{
							outFile << T[i][j][u][v].edge[k1][k2][0].abs2() + T[i][j][u][v].edge[k1][k2][1].abs2() + T[i][j][u][v].edge[k1][k2][2].abs2() + T[i][j][u][v].edge[k1][k2][3].abs2() << " ";
						}
					}
				}
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