#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "Complex.h"

int const N = 100;
int const nstep = 50;
double const pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513;

int main(int argc, char* argv[])
{
	double g = pi * 3/4;
	int H[4][4];

#pragma region H def
	H[0][0] = 1;
	H[0][1] = 1;
	H[0][2] = 1;
	H[0][3] = 1;

	H[1][0] = 1;
	H[1][1] = -1;
	H[1][2] = 1;
	H[1][3] = -1;

	H[2][0] = 1;
	H[2][1] = 1;
	H[2][2] = -1;
	H[2][3] = -1;

	H[3][0] = 1;
	H[3][1] = -1;
	H[3][2] = -1;
	H[3][3] = 1;
#pragma endregion

	if (argc > 1) {
		g = pi * atof(argv[1]);
	}
	std::cout << "g: " + std::to_string(g) << std::endl;
	complex*** state = new complex * *[N];
	for (unsigned int i = 0; i < N; i++) {
		state[i] = new complex * [N];
		for (unsigned int j = 0; j < N; j++) {
			state[i][j] = new complex[4];
		}
	}
	complex*** tempState = new complex * *[N];
	for (unsigned int i = 0; i < N; i++) {
		tempState[i] = new complex * [N];
		for (unsigned int j = 0; j < N; j++) {
			tempState[i][j] = new complex[4];
			for (unsigned int k = 0; k < 4; k++) {
				tempState[i][j][k] = 0.0;
			}
		}
	}
	tempState[N / 2][N / 2][0] = complex(1 / sqrt(4.0), 0.0);
	tempState[N / 2][N / 2][1] = complex(1 / sqrt(4.0), 0.0);
	tempState[N / 2][N / 2][2] = complex(1 / sqrt(4.0), 0.0);
	tempState[N / 2][N / 2][3] = complex(1 / sqrt(4.0), 0.0);

	for (unsigned int n = 1; n < nstep; n++)
	{
		// i = N / 2 - n; i <= N / 2 + n; i++
		// j = N / 2 - n; j <= N / 2 + n; j++
		for (unsigned int i = 1; i < N-1; i++) {
			for (unsigned int j = 1; j < N-1; j++) {

				for (unsigned int k = 0; k < 4; k++) state[i][j][k] = 0.0;

				for (unsigned int k = 0; k < 4; k++) {
					state[i][j][0] += (1.0 + double(i == j) * (exp(complex(0.0, g)) - 1.0)) * H[0][k] / 2.0 * tempState[i + 1][j + 1][k];
					state[i][j][1] += (1.0 + double(i == j) * (exp(complex(0.0, g)) - 1.0)) * H[1][k] / 2.0 * tempState[i + 1][j - 1][k];
					state[i][j][2] += (1.0 + double(i == j) * (exp(complex(0.0, g)) - 1.0)) * H[2][k] / 2.0 * tempState[i - 1][j + 1][k];
					state[i][j][3] += (1.0 + double(i == j) * (exp(complex(0.0, g)) - 1.0)) * H[3][k] / 2.0 * tempState[i - 1][j - 1][k];
				}
			}
		}

		for (unsigned int i = 1; i < N - 1; i++) {
			for (unsigned int  j = 1; j < N - 1; j++) {
				for (unsigned int k = 0; k < 4; k++) {
					tempState[i][j][k] = state[i][j][k];
				}
			}
		}

	}
	delete[] state;

	std::string filename = "x.txt";
	std::ofstream outFile;
	outFile.open(filename.c_str());
	for (int i = 1; i < N; i+=2)
	{
		for (int j = 1; j < N; j+=2) {
			outFile << i - N / 2 << " ";
		}
		outFile << std::endl;
	}
	outFile.close();

	filename = "y.txt";
	outFile.open(filename.c_str());
	for (int i = 1; i < N; i+=2)
	{
		for (int j = 1; j < N; j+=2) {
			outFile << j - N / 2 << " ";
		}
		outFile << std::endl;
	}
	outFile.close();

	filename = "p.txt";
	outFile.open(filename.c_str());
	for (unsigned int i = 1; i < N; i+=2)
	{
		for (unsigned int j = 1; j < N; j+=2) {
			outFile << tempState[i][j][0].abs2() + tempState[i][j][1].abs2() + tempState[i][j][2].abs2() + tempState[i][j][3].abs2() << " ";
		}
		outFile << std::endl;
	}
	outFile.close();

	delete[] tempState;
}
