#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <ctime>
#include <string>
#include <chrono>
#include <numeric>
#include <list>
#include <mpi.h>
#include <limits.h>
#include <algorithm>
#define COMUNICATION_SWITCH 1
#define DEBUGLOG 0
#define POLLINGRATE 500
using namespace std;

int checkEndState(int** M, int V);
int checkCost(int** M, int** C, int V);
void saveBestSol(int* bestSol, int** tmp_sol, int V);
int** createNewSol(int** M, int i_new, int j_new, int V, int sol);
list<int> addConstraints(int** M, int V);
int findBestArc(int** M, int** C, int V);
int findSecondBestArc(int** M, int** C, int V);
int ComputeLowerBound(int** M, int V, int** C);
int** createFirstSol(int V);
pair<int, int> findNextChild(int** M, int V);
bool checkFinalSol(int** tmp_sol, int V);
int checkConstrain(int** M, int V);
void TspStep(int** M, int step);
void rollBackChanges(int** M, int V, list<int> changeLog);
void printLog(string message, int** M, int V);
int ParallelTSP(int** M);
void TspExploreToFirstSolution(int** M);

/*******************GLOBAL VARIABLES ****************/
int counter;
int maxDeep;
int** C; 
int V; 
int* Bestsolution;
int bestCost; 
int p; 
int size;
MPI_Request request;
MPI_Status status;
int besttmp;
int bestIsMine;
/*******************GLOBAL VARIABLES ****************/

int main(int argc, char* argv[]) {



    ifstream myfile;
    if (argc != 2) {
        cout << "Error - number of parameter\nUsage: program ./file_distance.txt\nfile_distance must contain the number of vertex and then the distance between them, example:\n4\n0 1 2 3\n1 0 4 5\n2 4 0 6\n3 5 6 0 " << endl;
    }
    /********************************Read file ******************************/
    myfile.open(argv[1]);
    myfile >> V;

    if (V <= 0) {
        cout << "Number of nodes not correct";
        return 1;
    }
    C = createFirstSol(V);
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            myfile >> C[i][j];
        }
    }
    myfile.close();
    /********************************Read file ******************************/
    int** firstSol = createFirstSol(V);
    Bestsolution = (int*)(calloc(V+1, sizeof(int)));
    bestCost = INT_MAX;
    int firstBound = ComputeLowerBound(firstSol, V, C);
    counter=0;
    maxDeep=0;
    bestIsMine=0;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    double duration;

    auto start = chrono::steady_clock::now();

    //cout<<"size:"<<size<<endl;;
    ParallelTSP(firstSol);





    auto end = chrono::steady_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    
    
    cout << endl;
    cout << "Processor: "<<p<<endl;
    cout << "Time elapsed: " << duration << "ms" << endl;
    cout << "State explored counter: "<<counter<<endl;
    cout << "My optimal cost: "<<bestCost<<endl;
    cout << "Deepth explored: "<<maxDeep<<endl;
    
    if(bestIsMine==0){
        bestCost  = INT_MAX;
    }
    MPI_Barrier( MPI_COMM_WORLD );
    int *solList = (int *)malloc(sizeof(int) * size);
    MPI_Allgather(&bestCost, 1, MPI_INT, solList, 1, MPI_INT, MPI_COMM_WORLD);
    int min_idx = min_element(solList,solList+size) - solList;

    if(p==min_idx){
        cout <<endl<<endl<< "Final solution from "<<p<<" :\n" << "Cost:" << bestCost << endl << "Path:" << endl;
        for (int i = 0; i < V + 1; i++) {
            cout << Bestsolution[i] << " ";
        }

    cout<<endl;
    }

    free(solList);
    free(Bestsolution);
    MPI_Finalize();
    return 0;
}

int ParallelTSP(int** M)
{

    TspExploreToFirstSolution(M);
    if(DEBUGLOG){
        if(p==0){
            printLog("Processor "+to_string(p)+" will start from this solution", M, V);
            MPI_Send(&bestCost, 1, MPI_INT, p+1, 0, MPI_COMM_WORLD);
        }else if(p==(size-1)){
            MPI_Recv(&bestCost, 1, MPI_INT, p-1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            printLog("Processor "+to_string(p)+" will start from this solution", M, V);
        }else{
            MPI_Recv(&bestCost, 1, MPI_INT, p-1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            cout<<p<<" "<<size<<endl;
            printLog("Processor "+to_string(p)+" will start from this solution", M, V);
            MPI_Send(&bestCost, 1, MPI_INT, p+1, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier( MPI_COMM_WORLD );
    }
	TspStep(M,0);
    return bestIsMine;
}


void TspExploreToFirstSolution(int** M){

    int choice;
    int p_tmp=p;
    int size_tmp=size;
    for(;size_tmp!=1;p_tmp=p>>1,size_tmp=size_tmp>>1){

        choice=p_tmp%2;
        //printLog("Processor "+to_string(p)+" generated this solution", M, V);
        pair<int, int> idx = findNextChild(M, V);
        int i = idx.first;
        int j = idx.second;
        if(i<0|| i>=V || j<0 || j>=V){
            break;
        }
        M[i][j] = choice==1 ? 1 : -1;
        M[j][i] = choice==1 ? 1 : -1;
        addConstraints(M, V);

    }


}

void TspStep(int** M, int step)
{
    
    counter++; //Global variable, it count the number of "nodes" in the tree visidet for each processor
    //cout<<counter<<endl;
    if(step>maxDeep){
        maxDeep=step; //Global variable, it count the max deepth reached by this processor
    }
    if(COMUNICATION_SWITCH){ //constant ( #define  on top ) to enable/disable comunication between processors
    int flag=0; //flag == 1 if it will recieve a new solution from another processor
        if(step==0){ // first time
            MPI_Irecv(&besttmp, 1, MPI_INT, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &request); //Non blocking recv from all other processor to exange the best sol
            MPI_Test(&request, &flag, &status); //Test if I recived something from someone else 
        }else if(counter%POLLINGRATE==0) { // Polling each 500 state visited.
            MPI_Test(&request, &flag, &status);
        }
        if(flag!=0){ //IF p recived  something
        
        if(DEBUGLOG){cout<<"I am "<<p<<"  today  "<<status.MPI_SOURCE<< "  sent me something: "<<besttmp<<endl;}
        if(DEBUGLOG){cout<<"My current best is "<<bestCost<<endl;}
            if(besttmp<bestCost){ //And the solution recived is better than mine
        if(DEBUGLOG){   cout<<"Updating value!"<<endl;}
                bestCost=besttmp; //Update my solution
                bestIsMine=0; //Flag that says that I am not the processor that has the optimal solution
            }
            MPI_Irecv(&besttmp, 1, MPI_INT, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &request); //p open a new comunication channel, to wait for new sol
        }
    }


    if (checkEndState(M, V) == 1) {          //END STATE
        int tmp_cost = checkCost(M, C, V);

        if (tmp_cost < bestCost && checkFinalSol(M, V)) { //check bestCost and if finalSol is correct 
            bestCost = tmp_cost;
            saveBestSol(Bestsolution, M, V);
            bestIsMine=1;
        if(DEBUGLOG){   cout<<"Processor "<< p << " Found a better solution : "<<bestCost<<endl;}
           if(COMUNICATION_SWITCH){ //Send to all the other p
                for(int i=0;i<size;i++){
                    if(i!=p){
                        MPI_Request new_request;
                        MPI_Isend(&bestCost, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &new_request);
                        //MPI_Send(&bestCost, 1, MPI_INT, i, 0, MPI_COMM_WORLD); 
                    }
                }

            }
            
            //cout<<"I sent "<<bestCost <<" to all other processors"<<endl;
            //MPI_Ibcast(&bestCost, 1, MPI_INT, p, MPI_COMM_WORLD, &found_request);
            //MPI_Wait(&found_request, MPI_STATUS_IGNORE);
        }
        return;
    }

    list<int> changeLogA;

    pair<int, int> idx = findNextChild(M, V);
    int i = idx.first;
    int j = idx.second;

    //child 1
    M[i][j] = 1;
    M[j][i] = 1;
    //newSol = addConstraints(newSol, V);

    const int lowerBoundA = ComputeLowerBound(M, V, C);

    M[i][j] = -1;
    M[j][i] = -1;
    const int lowerBoundB = ComputeLowerBound(M, V, C);

    if (lowerBoundA < lowerBoundB)
    {
        M[i][j] = 1;
        M[j][i] = 1;
        //printLog("Adding arc", M, V);

        changeLogA = addConstraints(M, V);


        //printLog("Adding constraint", M, V);
        if (lowerBoundA < bestCost && checkConstrain(M, V)) {
            TspStep(M,step+1);
        }

        rollBackChanges(M, V, changeLogA);

        //printLog("rolling back constraint", M, V);
        M[i][j] = -1;
        M[j][i] = -1;

        ////printLog("Adding other arc", M, V);
        changeLogA = addConstraints(M, V);

        //printLog("Adding constraints", M, V);
        if (lowerBoundB < bestCost && checkConstrain(M, V)) {
            TspStep(M,step+1);
        }
        rollBackChanges(M, V, changeLogA);

        //printLog("rollingback constraints", M, V);
        M[i][j] = 0;
        M[j][i] = 0;

        //printLog("backtrack arc", M, V);

    }else {

        M[i][j] = -1;
        M[j][i] = -1;

        //printLog("Adding arc", M, V);
        changeLogA = addConstraints(M, V);

        //printLog("Adding constraints", M, V);
        if (lowerBoundB < bestCost && checkConstrain(M, V)) {
            TspStep(M,step+1);
        }
        rollBackChanges(M, V, changeLogA);

        //printLog("rolling back constraints", M, V);
        M[i][j] = 1;
        M[j][i] = 1;

        //printLog("Adding other arc", M, V);
        changeLogA = addConstraints(M, V);

        //printLog("adding constraints", M, V);
        if (lowerBoundA < bestCost && checkConstrain(M, V)) {
            TspStep(M,step+1);
        }
        rollBackChanges(M, V, changeLogA);

        //printLog("rolling back constraints", M, V);
        M[i][j] = 0;
        M[j][i] = 0;

        //printLog("backtrack arc", M, V);
    }


}

void printLog(string message, int** M, int V)
{
    cout << message << endl;
    for (int q = 0; q < V; q++)
    {
        for (int z = 0; z < V; z++)
        {
            printf("%d ", M[q][z]);
        }
        printf("\n");

    }

    printf("\n");
}

pair<int, int> findNextChild(int** M, int V) {
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (M[i][j] == 0) {
                return make_pair(i, j);
            }
        }
    }
    return make_pair(-1, -1);
}

int checkEndState(int** M, int V) {
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (M[i][j] == 0)
                return 0;
        }
    }
    return 1;
}

int checkCost(int** M, int** C, int V) {
    int solCost = 0;
    for (int i = 0; i < V; i++) {
        for (int j = i; j < V; j++) {
            if (M[i][j] == 1) {
                solCost += C[i][j];
            }
        }
    }
    return solCost;
}

bool checkFinalSol(int** tmp_sol, int V) {
    int pos = 0;
    vector<int> visited(V);
    visited[0] = 1;
    for (int k = 0; k < V; k++) {
        for (int i = 0; i < V; i++) {
            if (tmp_sol[pos][i] == 1 && visited[i] == 0) {
                pos = i;
                break;
            }
        }

        visited[pos] = 1;
    }
    for (int i = 0; i < V; i++) {
        if (visited[i] == 0) {
            return false;
        }
    }
    return true;
}

void saveBestSol(int* bestSol, int** tmp_sol, int V) {

    bestSol[0] = 0;
    int pos = 0;
    int prev = -1;
    vector<int> visited(V);
    visited[0] = 1;
    for (int k = 1; k < V; k++) {
        for (int i = 0; i < V; i++) {
            if (tmp_sol[pos][i] == 1 && visited[i] == 0) {
                pos = i;
                break;
            }
        }

        visited[pos] = 1;
        bestSol[k] = pos;
        //prev = best_sol[k-1];

    }
    bestSol[V] = 0;

}

int** createNewSol(int** M, int i_new, int j_new, int V, int sol) {
    int** mat(M);
    mat[i_new][j_new] = sol;
    mat[j_new][i_new] = sol;
    return mat;
}

int** createFirstSol(int V) {
    int** matrix = (int**)malloc(V * sizeof(int*));

    for (int i = 0; i < V; i++) {
        matrix[i] = static_cast<int*>(calloc(V, sizeof(int)));
        matrix[i][i] = -1;
    }
    return matrix;

}

int checkConstrain(int** M, int V) {
    for (int i = 0; i < V; i++) {
        int count_one = 0;
        int count_zero = 0;
        for (int j = 0; j < V; j++) {
            if (M[i][j] == 1) {
                count_one++;
            }
            else if (M[i][j] == 0) {
                count_zero++;
            }
        }

        if (count_one + count_zero < 2) {
            return 0;
        }
        if (count_one > 2)
        {
            return 0;
        }
    }
    return 1;
}

list<int> addConstraints(int** M, int V) {
    int flag = 1;
    list<int> changeLog;
    while (flag) {
        flag = 0;
        for (int i = 0; i < V; i++) {
            int count_one = 0;
            int count_zero = 0;
            for (int j = 0; j < V; j++) {
                if (M[i][j] == 1) {
                    count_one++;
                    if (count_one == 2) { break; }
                }
                else if (M[i][j] == 0) {
                    count_zero++;
                }
            }
            if (count_one == 2) {
                for (int j = 0; j < V; j++) {
                    if (M[i][j] == 0) {
                        M[i][j] = -1;
                        M[j][i] = -1;
                        changeLog.push_back(i * V + j);
                        flag = 1;
                    }
                }
            }
            else if (count_one + count_zero == 2) {
                for (int j = 0; j < V; j++) {
                    if (M[i][j] == 0) {
                        M[i][j] = 1;
                        M[j][i] = 1;
                        changeLog.push_back(i * V + j);
                        flag = 1;
                    }
                }
            }
        }
    }
    return changeLog;
}

void rollBackChanges(int** M, int V, list<int> changeLog)
{

    for (int t : changeLog)
    {
        int i = t / V;
        int j = t % V;
        M[i][j] = 0;
        M[j][i] = 0;
    }
}


int ComputeLowerBound(int** M, int V, int** C) {

    int best = findBestArc(M, C, V);
    int second = findSecondBestArc(M, C, V);

    int lowerBound = best + second;
    if (lowerBound % 2 == 0) {
        lowerBound /= 2;
    }
    else {
        lowerBound = (lowerBound + 1) / 2;
    }
    return lowerBound;
}

int findBestArc(int** M, int** C, int V) {
    int sol = 0;
    for (int i = 0; i < V; i++) {
        int best = INT_MAX;
        for (int j = 0; j < V; j++) {
            if (M[i][j] == 1) {
                best = C[i][j];
                break;
            }
            if (M[i][j] == 0) {
                if (best > C[i][j]) {
                    best = C[i][j];
                }
            }

        }
        sol += best;
    }
    return sol;
}

int findSecondBestArc(int** M, int** C, int V) {
    int sol = 0;


    for (int i = 0; i < V; i++) {
        int best = INT_MAX - 1;
        int secondBest = INT_MAX;
        int flag = 0; //Skip first 1 in matrix

        for (int j = 0; j < V; j++) {
            if (M[i][j] == 1) {
                if (flag == 0) {
                    flag = 1;
                }
                else {
                    secondBest = C[i][j];
                    flag = 0;
                    break;
                }
            }

            if (M[i][j] == 0) {
                if (best > C[i][j]) {
                    secondBest = best;
                    best = C[i][j];
                }
                else if (secondBest > C[i][j]) {
                    secondBest = C[i][j];
                }
            }
        }
        if (flag == 1) {
            sol += best;
        }
        else {
            sol += secondBest;
        }

    }
    return sol;
}

