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
#define COMUNICATION_SWITCH 1
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
void TspStep(int** C, int V, int** M, int* Bestsolution, int& bestCost, int p, int size, MPI_Request request, MPI_Status status,int step,int* besttmp, int* bestIsMine);
void rollBackChanges(int** M, int V, list<int> changeLog);
void printLog(string message, int** M, int V);
int ParallelTSP(int** C, int V, int** M, int* Bestsolution, int& bestCost, int p, int size);
void TspExploreToFirstSolution(int V, int** M, int p,int size);
int counter;
int maxDeep;
int main(int argc, char* argv[]) {



    ifstream myfile;
    int V;
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
    int** C = createFirstSol(V);
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            myfile >> C[i][j];
        }
    }
    myfile.close();
    /********************************Read file ******************************/
    priority_queue<pair<int, int**>> priorityQueue;
    int** firstSol = createFirstSol(V);
    int* finalSol = (int*)(calloc(V, sizeof(int)));
    int finalCost = INT_MAX;
    int firstBound = ComputeLowerBound(firstSol, V, C);
    int size, p;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    double duration;

    auto start = chrono::steady_clock::now();
    counter=0;
    maxDeep=0;
    cout<<"size:"<<size<<endl;;
    int flag=ParallelTSP(C, V, firstSol, finalSol, finalCost,p,size);





    auto end = chrono::steady_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    
    if( flag ){
        cout << "Solution:\n" << "Cost:" << finalCost << endl << "Path:" << endl;
        for (int i = 0; i < V + 1; i++) {
            cout << finalSol[i] << " ";
        }
    }
    cout << endl;
    cout << "Time elapsed: " << duration << "ms" << endl;
    cout << "counter: "<<counter<<endl;
    cout << "Deep: "<<maxDeep<<endl;
    MPI_Finalize();
    return 0;
}

int ParallelTSP(int** C, int V, int** M, int* Bestsolution, int& bestCost,int p, int size)
{

    TspExploreToFirstSolution(V,M,p,size);
    //printLog("Processor "+to_string(p)+" will start from this solution", M, V);
    MPI_Request request;
    MPI_Status status;
    int besttmp;
    int bestIsMine=0;
	TspStep( C,  V, M,  Bestsolution, bestCost, p,size,request,status,0, &besttmp, &bestIsMine);
    return bestIsMine;
}


void TspExploreToFirstSolution(int V, int** M, int p,int size){

    int choice;

    for(;size!=1;p=p>>1,size=size>>1){

        choice=p%2;
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

void TspStep(int** C, int V, int** M, int* Bestsolution, int& bestCost, int p, int size,MPI_Request request, MPI_Status status, int step, int* besttmp, int* bestIsMine)
{
    /*
    cout << "new step:" << endl;
    for(int q=0;q<V;q++)
    {
        for(int z=0;z<V;z++)
        {
            printf("%d ", M[q][z]);
        }
        printf("\n");

    }

    printf("\n");
    printf("\n");

    printf("\n");
    */

    /*
    
    for (int i=0;i<size;i++){
        if(i!=p){
            MPI_Ibcast(besttmp, 1, MPI_INT, i, MPI_COMM_WORLD, &found_request);
        }
    }
    */
    //printLog("Processor "+to_string(p)+" is exploring at step: "+to_string(step), M, V);
    counter++;

    if(step>maxDeep){
        maxDeep=step;
    }
    if(COMUNICATION_SWITCH){
    int flag=0;
        if(step==0){ // first time

            MPI_Irecv(besttmp, 1, MPI_INT, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &request);
            MPI_Test(&request, &flag, &status);
        }else if(counter%500==0) {
            MPI_Test(&request, &flag, &status);
        }
        if(flag!=0){
        
        cout<<"I am "<<p<<"  today  "<<status.MPI_SOURCE<< "  sent me something: "<<*besttmp<<endl;
        cout<<"My current best is "<<bestCost<<endl;
            if(*besttmp<bestCost){
            //cout<<"Updating value!"<<endl;
                bestCost=*besttmp;
                *bestIsMine=0;
            }
            MPI_Irecv(&besttmp, 1, MPI_INT, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &request);
        }
    }


    if (checkEndState(M, V) == 1) {          //END STATE
        int tmp_cost = checkCost(M, C, V);

        if (tmp_cost < bestCost && checkFinalSol(M, V)) {
            bestCost = tmp_cost;
            saveBestSol(Bestsolution, M, V);
            cout<<"Processor "<< p<< " Found a better solution : "<<bestCost<<endl;
           if(COMUNICATION_SWITCH){
                for(int i=0;i<size;i++){
                    if(i!=p){
                        MPI_Request new_request;
                        MPI_Isend(&bestCost, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &new_request);
                        //MPI_Send(&bestCost, 1, MPI_INT, i, 0, MPI_COMM_WORLD); 
                    }
                }
            }
            *bestIsMine=1;
            cout<<"I sent "<<bestCost <<" to all other processors"<<endl;
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
            TspStep(C, V, M, Bestsolution, bestCost,p,size,request,  status, step+1,besttmp,bestIsMine);
        }

        rollBackChanges(M, V, changeLogA);

        //printLog("rolling back constraint", M, V);
        M[i][j] = -1;
        M[j][i] = -1;

        ////printLog("Adding other arc", M, V);
        changeLogA = addConstraints(M, V);

        //printLog("Adding constraints", M, V);
        if (lowerBoundB < bestCost && checkConstrain(M, V)) {
            TspStep(C, V, M, Bestsolution, bestCost,p,size,request,  status, step+1 ,besttmp,bestIsMine);
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
            TspStep(C, V, M, Bestsolution, bestCost,p,size,request,  status, step+1 ,besttmp,bestIsMine);
        }
        rollBackChanges(M, V, changeLogA);

        //printLog("rolling back constraints", M, V);
        M[i][j] = 1;
        M[j][i] = 1;

        //printLog("Adding other arc", M, V);
        changeLogA = addConstraints(M, V);

        //printLog("adding constraints", M, V);
        if (lowerBoundA < bestCost && checkConstrain(M, V)) {
            TspStep(C, V, M, Bestsolution, bestCost,p,size,request,  status, step+1 ,besttmp, bestIsMine);
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

