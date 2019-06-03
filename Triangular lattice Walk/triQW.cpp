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
inline int md(int i, int n){return (i % n + n) % n;}


int main(int argc, char *argv[])
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
      S[0][0] = 1.0;
      S[0][1] = 0.0;
      S[1][0] = 0.0;
      S[1][1] = exp(complex(0.0, 2.0 * pi / 3.0));

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
      int nstep = 20;
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
      triangle **T = new triangle *[grid];

      //Creation of a temp wave function usefull to apply operator
      triangle **cT = new triangle *[grid];

      for (unsigned int i = 0; i < grid; i++)
      {
            T[i] = new triangle[grid];
            cT[i] = new triangle[grid];
      }
      // That choose where is our initial state
      int x0 = grid / 2;
      // Setting the inital state
      cT[x0][x0].e0[0] = complex(1.0 / sqrt(6.0), 0.0);
      cT[x0][x0].e1[0] = complex(1.0 / sqrt(6.0), 0.0);
      cT[x0][x0].e2[0] = complex(1.0 / sqrt(6.0), 0.0);
      cT[x0][x0].e0[1] = complex(1.0 / sqrt(6.0), 0.0);
      cT[x0][x0].e1[1] = complex(1.0 / sqrt(6.0), 0.0);
      cT[x0][x0].e2[1] = complex(1.0 / sqrt(6.0), 0.0);

#pragma endregion

#pragma region Encoding with U
      std::cout << "Encoding the wave function..." << std::endl;
      // Compute the product between U0 and the wave function
      for (int i = 0; i < grid; i++)
      {
            for (int j = 0; j < grid; j++)
            {
                  for (unsigned int o = 0; o < 2; o++)
                  {
                        T[i][j].e0[o] = 0.0;
                        T[i][j].e1[o] = 0.0;
                        T[i][j].e2[o] = 0.0;

                        for (unsigned int l = 0; l < 2; l++)
                        {
                              T[i][j].e0[o] += U1[0][o][l] * cT[i][j].e0[l];
                              T[i][j].e1[o] += U1[0][o][l] * cT[i][j].e1[l];
                              T[i][j].e2[o] += U1[0][o][l] * cT[i][j].e2[l];
                        }
                  }
            }
      }

#pragma endregion

#pragma region The Walk
      for (int n = 1; n < nstep; n++)
      {

            std::cout << "Step " + std::to_string(n) << std::endl;

            // Applying the rotation operator
            for (int i = x0 - n; i < x0 + n + 1; i++)
            {
                  for (int j = x0 - n; j < x0 + n + 1; j++)
                  {
                        // Up
                        T[i][j].e1[0] = cT[i][j].e0[0];
                        T[i][j].e2[0] = cT[i][j].e1[0];
                        T[i][j].e0[0] = cT[i][j].e2[0];

                        // Down
                        T[md(i - 1, grid)][j].e1[1] = cT[i][j].e0[1];
                        T[md(i + j % 2, grid)][md(j + 1, grid)].e2[1] = cT[i][j].e1[1];
                        T[md(i + j % 2, grid)][md(j - 1, grid)].e0[1] = cT[i][j].e2[1];
                  }
            }
            // Applying the coin operator
            for (int i = x0 - n - 1; i < x0 + n + 2; i++)
            {
                  for (int j = x0 - n - 1; j < x0 + n + 2; j++)
                  {
                        {
                              cT[i][j] = T[i][j];
                              setNulltriangle(&T[i][j]);
                              for (int o = 0; o < 2; o++)
                              {
                                    for (int l = 0; l < 2; l++)
                                    {
                                          T[i][j].e0[l] += W1[l][o] * cT[i][j].e0[o];
                                          T[i][j].e1[l] += W1[l][o] * cT[i][j].e1[o];
                                          T[i][j].e2[l] += W1[l][o] * cT[i][j].e2[o];
                                    }
                              }

                              cT[i][j] = T[i][j];
                              setNulltriangle(&T[i][j]);
                        }
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
                  for (unsigned int o = 0; o < 2; o++)
                  {
                        T[i][j].e0[o] = 0.0;
                        T[i][j].e1[o] = 0.0;
                        T[i][j].e2[o] = 0.0;

                        for (unsigned int l = 0; l < 2; l++)
                        {
                              T[i][j].e0[o] += U1[0][l][o].conjugate() * cT[i][j].e0[l];
                              T[i][j].e1[o] += U1[0][l][o].conjugate() * cT[i][j].e1[l];
                              T[i][j].e2[o] += U1[0][l][o].conjugate() * cT[i][j].e2[l];
                        }
                  }
            }
            delete[] T[i];
      }
      delete[] T;

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
                  outFile << cT[i][j].e0[0].abs2() + cT[i][j].e0[1].abs2() << " ";
                  outFile << cT[i][j].e1[0].abs2() + cT[i][j].e1[1].abs2() << " ";
                  outFile << cT[i][j].e2[0].abs2() + cT[i][j].e2[1].abs2() << " ";
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