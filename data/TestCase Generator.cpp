// TestCase Generator.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#define MAX 100
using namespace  std;
int main(int argc, char* argv[])
{

	int N;
	int NN = 9;

	for (N = 4; N <= 9; N++) {
		vector<vector<int>>M(N, vector<int>(N));
		srand(time(0));
		for (int it = 1; it < 101; it++) {
			ofstream myfile;
			myfile.open("problem_"+ std::to_string(N)+"_"+std::to_string(it));
			for (int i = 0; i < N; i++) {
				for (int j = i; j < N; j++) {
					if (i == j)
					{
						continue;
					}
					M[i][j] = (rand() % MAX) + 1;
					M[j][i] = (rand() % MAX) + 1;
				}
			}
			myfile << N << endl;
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					myfile << M[i][j];
						if (j != N - 1) {
							myfile << " ";
						}
				}
				if (i != N - 1) {
					myfile << endl;
				}
			}
		}
	}
	
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
