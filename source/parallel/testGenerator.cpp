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