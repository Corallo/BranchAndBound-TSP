#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <ctime>
#include <string>
#include <chrono>
#include <numeric>
#include <list>

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
void TspStep(int** C, int V, int** M, int* Bestsolution, int& bestCost);
void rollBackChanges(int** M, int V, list<int> changeLog);
void printLog(string message, int** M, int V);

int main(int argc, char* argv[]) {

    int V;
    ifstream myfile;
    ofstream reportFile;

    /*
    if (argc != 2) {
        cout << "Error - number of parameter\nUsage: program ./file_distance.txt\nfile_distance must contain the number of vertex and then the distance between them, example:\n4\n0 1 2 3\n1 0 4 5\n2 4 0 6\n3 5 6 0 " << endl;
    */
    /********************************Read file ******************************/
    reportFile.open("report.txt");
    int a, b;
    for ( a = 21; a <= 25; a++)
    {
        if(a==21)
        {
            b = 19;
        }else{
        
            b = 1;
        }
        vector<double> timeResult(100);
        for ( ; b <= 100; b++)
        {

            string filename = "./TestCase/problem_" + std::to_string(a) + "_" + std::to_string(b);

            myfile.open(filename);
            myfile >> V;

            if (V <= 0) {
                cout << "Number of nodes not correct";
                return 1;
            }
            int** C = createFirstSol(V); //Allocate matrix of cost;

            for (int i = 0; i < V; i++) {
                for (int j = 0; j < V; j++) {
                    myfile >> C[i][j];
                }
            }
            myfile.close();
            /********************************Read file ******************************/
            int** firstSol = createFirstSol(V);
            int* finalSol = (int*)(calloc(V, sizeof(int)));
            int finalCost = INT_MAX;
            int firstBound = ComputeLowerBound(firstSol, V, C);

            double duration;

            auto start = chrono::steady_clock::now();


            TspStep(C, V, firstSol, finalSol, finalCost);



            auto end = chrono::steady_clock::now();
            duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        	
            /*
            cout << "Solution:\n" << "Cost:" << finalCost << endl << "Path:\n" << endl;
            for (int i = 0; i < V + 1; i++) {
                cout <<finalSol[i] << " ";
            }

            cout << endl;
            */
            cout << "Problem_" + std::to_string(a) + "_" + std::to_string(b) << " ";
            cout << "Time elapsed: " << duration << endl;
            
            /*
            ifstream filesol;
                    filesol.open(filename + "_sol");
                    int realCost;
                    filesol >> realCost;
                    if(realCost== finalCost)
                    {
                        reportFile << "Problem " + std::to_string(a) + "_" + std::to_string(b) + " Correct" << endl;
                        cout << "Problem " + std::to_string(a) + "_" + std::to_string(b) + " Correct" << endl;

                    }else
                    {
                        cout << "Problem " + std::to_string(a) + "_" + std::to_string(b) + " UnCorrect ";
                        reportFile << "Problem " + std::to_string(a) + "_" + std::to_string(b) + " Uncorrect ";
                        if (realCost > finalCost)
                        {
                            cout << "Our solution is better!!!!" << endl;
                            reportFile << "Our solution is better!!!!" << endl;
                        }else
                        {
                            cout << "Our solution is worste :(" << endl;
                            reportFile << "Our solution is worste :(" << endl;
                        }
                    }
               */     


        }
        //cout << "Average time for problem of size " + std::to_string(a) + " :";
        double sum_of_elems;
        sum_of_elems = std::accumulate(timeResult.begin(), timeResult.end(), 0);
        //cout << sum_of_elems / 100 << endl;
    }
    return 0;
}
//Alloca best solution (V)
//Aggiungi autoSolver
void TspStep(int** C, int V, int** M, int* Bestsolution, int& bestCost)
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

    if (checkEndState(M, V) == 1) {          //END STATE
        int tmp_cost = checkCost(M, C, V);

        if (tmp_cost < bestCost && checkFinalSol(M, V)) {
            bestCost = tmp_cost;
            saveBestSol(Bestsolution, M, V);
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
            TspStep(C, V, M, Bestsolution, bestCost);
        }

        rollBackChanges(M, V, changeLogA);

        //printLog("rolling back constraint", M, V);
        M[i][j] = -1;
        M[j][i] = -1;

        ////printLog("Adding other arc", M, V);
        changeLogA = addConstraints(M, V);

        //printLog("Adding constraints", M, V);
        if (lowerBoundB < bestCost && checkConstrain(M, V)) {
            TspStep(C, V, M, Bestsolution, bestCost);
        }
        rollBackChanges(M, V, changeLogA);

        //printLog("rollingback constraints", M, V);
        M[i][j] = 0;
        M[j][i] = 0;

        //printLog("backtrack arc", M, V);

    }
    else {

        M[i][j] = -1;
        M[j][i] = -1;

        //printLog("Adding arc", M, V);
        changeLogA = addConstraints(M, V);

        //printLog("Adding constraints", M, V);
        if (lowerBoundB < bestCost && checkConstrain(M, V)) {
            TspStep(C, V, M, Bestsolution, bestCost);
        }
        rollBackChanges(M, V, changeLogA);

        //printLog("rolling back constraints", M, V);
        M[i][j] = 1;
        M[j][i] = 1;

        //printLog("Adding other arc", M, V);
        changeLogA = addConstraints(M, V);

        //printLog("adding constraints", M, V);
        if (lowerBoundA < bestCost && checkConstrain(M, V)) {
            TspStep(C, V, M, Bestsolution, bestCost);
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

