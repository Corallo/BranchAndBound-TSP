#include <iostream>
#include<stdlib.h>
#include <string>
#include <chrono>
#include <thread>

using namespace std;
int main () {
   string file;
   string command;
   int i,j;
   for (int i=8;i<=25;i++){

   	for(j=1;j<=100;j++){

   		//cout<<"Starting problem "+std::to_string(i)+"_"+std::to_string(j)<<endl;
   		file="./TestCase/problem_"+std::to_string(i)+"_"+std::to_string(j);
   		command = "timeout "+std::to_string(i*5)+" mpirun -np 4 ./main "+file;
   		//cout<<command<<endl;
   		system(command.c_str());
   		std::this_thread::sleep_for(std::chrono::milliseconds(1000));
   	}
   }
   

   return(0);
} 