#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <ctime>
#include <string>
#include <chrono>
#include <numeric>

using namespace std;

int checkEndState(vector<vector<int>> M, int V);
int checkCost(vector<vector<int>> M, vector<vector<int>> C, int V);
vector<int> saveBestSol(vector<vector<int>> tmp_sol, int V);
vector<vector<int>> createNewSol(vector<vector<int>> M, int i_new, int j_new, int V, int sol);
vector<vector<int>> addConstraints(vector<vector<int>> M, int V);
int findBestArc(vector<vector<int>> M, vector<vector<int>> C, int V);
int findSecondBestArc(vector<vector<int>> M, vector<vector<int>> C, int V);
int ComputeLowerBound(vector<vector<int>> M, int V, vector<vector<int>> C);
vector<vector<int>> createFirstSol(int V);
pair<int, int> findNextChild(vector<vector<int>> M, int V);
pair<int,vector<int>> Tsp(vector<vector<int>> C, int V, priority_queue<pair<int, vector<vector<int>>>> priorityQueue);
bool checkFinalSol(vector<vector<int>> tmp_sol, int V);
int checkConstrain(vector<vector<int>> M, int V);
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
    for (int a = 4; a <= 9; a++)
    {

        vector<double> timeResult(100);
        for (int b = 1; b <= 100; b++)
        {

            string filename = "../../TestCase/problem_" + std::to_string(a) + "_" + std::to_string(b);

            myfile.open(filename);
            myfile >> V;

            if (V <= 0) {
                cout << "Number of nodes not correct";
                return 1;
            }
            vector<vector<int>> C = createFirstSol(V); //Allocate matrix of cost;

            for (int i = 0; i < V; i++) {
                for (int j = 0; j < V; j++) {
                    myfile >> C[i][j];
                }
            }
            myfile.close();
            /********************************Read file ******************************/
            priority_queue<pair<int, vector<vector<int>>>> priorityQueue;
            vector<vector<int>> firstSol = createFirstSol(V);
            int firstBound = ComputeLowerBound(firstSol, V, C);
            priorityQueue.push(make_pair(-firstBound, firstSol));


            /*
    std::clock_t start;
    double duration;

    start = std::clock();
    */

            auto start = chrono::steady_clock::now();

            pair<int, vector<int>> sol = Tsp(C, V, priorityQueue);
            auto end = chrono::steady_clock::now();
            timeResult[b - 1] = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            /*
            duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

            cout << "Solution:\n" << "Cost:" << sol.first << endl << "Path:\n" << endl;
            for (int i = 0; i < V + 1; i++) {
                cout << sol.second[i] << " ";
            }

            cout << endl;
            cout << "Time elapsed: " << duration << endl;
            */
            /*
            ifstream filesol;
                    filesol.open(filename + "_sol");
                    int realCost;
                    filesol >> realCost;
                    if(realCost== sol.first)
                    {
                        reportFile << "Problem " + std::to_string(a) + "_" + std::to_string(b) + " Correct" << endl;
                        cout << "Problem " + std::to_string(a) + "_" + std::to_string(b) + " Correct" << endl;

                    }else
                    {
                        cout << "Problem " + std::to_string(a) + "_" + std::to_string(b) + " Correct" << endl;

                        reportFile << "Problem " + std::to_string(a) + "_" + std::to_string(b) + " Uncorrect"<<endl;
                    }
                    */


        }
        cout << "Average time for problem of size " + std::to_string(a) +" :";
        double sum_of_elems;
		sum_of_elems = std::accumulate(timeResult.begin(), timeResult.end(), 0);
        cout << sum_of_elems / 100 << endl;
    }
    return 0;
}

pair<int,vector<int>> Tsp(vector<vector<int>> C, int V, priority_queue<pair<int, vector<vector<int>>>> priorityQueue) {

    int best_cost = INT_MAX;
    vector<int> solution(V);
    while (priorityQueue.size() > 0) {
        vector<vector<int>> M = priorityQueue.top().second;
        vector<vector<int>> newSol;
        priorityQueue.pop();

        /***
         * Finché la queue non è vuota:
         *  Prendi una soluzione OK
         *  Controlla se è la soluzione finale OK
         *  se lo è, misura ed eventualmente aggiorna best_sol OK
         *  Trovi i possibili figli Ok
         *  Aggiungi i constraint ai figli Ok
         *  Calcoli i lower bound
         *  aggiungi i figli alla queue
         *
         */

        int tmp_cost;
        int lowerBound;
        if (checkEndState(M, V) == 1) {          //END STATE
            tmp_cost = checkCost(M, C, V);

            if (tmp_cost < best_cost && checkFinalSol(M, V)) {
                best_cost = tmp_cost;
                solution=saveBestSol(M, V);
            }
            continue;
        }
        else { //ELSE search for child

            pair<int, int> idx = findNextChild(M, V);
            int i = idx.first;
            int j = idx.second;
                    if (M[i][j] == 0) { //TAKE first free arc
                        newSol = createNewSol(M, i, j, V, 1);
                        newSol = addConstraints(newSol, V);
                        lowerBound = ComputeLowerBound(newSol, V, C);
                        if (-lowerBound < best_cost && checkConstrain(newSol,V)) {
                            priorityQueue.push(make_pair(-lowerBound, newSol));
                        }

                        newSol = createNewSol(M, i, j, V, -1);
                        newSol = addConstraints(newSol, V);
                        lowerBound = ComputeLowerBound(newSol, V, C);
                        if (-lowerBound < best_cost && checkConstrain(newSol, V)) {
                            priorityQueue.push(make_pair(-lowerBound, newSol));
                        }

                    }
        }
    }
    return make_pair(best_cost,solution);



}


pair<int,int> findNextChild(vector<vector<int>> M, int V) {
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (M[i][j] == 0) {
                return make_pair(i, j);
            }
        }
    }
    return make_pair(-1, -1);
}

int checkEndState(vector<vector<int>> M, int V) {
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (M[i][j] == 0)
                return 0;
        }
    }
    return 1;
}

int checkCost(vector<vector<int>> M, vector<vector<int>> C, int V) {
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

bool checkFinalSol(vector<vector<int>> tmp_sol, int V) {
    int pos = 0;
    vector<int> visited(V);
    visited[0] = 1;
    for(int k=0; k<V;k++) {
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

vector<int> saveBestSol(vector<vector<int>> tmp_sol,  int V) {
    vector<int> best_sol(V + 1);
    best_sol[0]= 0;
    int pos = 0;
    int prev = -1;
    vector<int> visited(V);
    visited[0] = 1;
    for( int k =1; k<V;k++){
            for (int i = 0; i < V; i++) {
                if (tmp_sol[pos][i] == 1 && visited[i]==0) {
                    pos = i;
                    break;
                }
            }

            visited[pos] = 1;
            best_sol[k] = pos;
        //prev = best_sol[k-1];

    } ;
    return best_sol;
}

vector<vector<int>> createNewSol(vector<vector<int>> M, int i_new, int j_new, int V, int sol) {
    vector<vector<int>> mat(M);
    mat[i_new][j_new] = sol;
    mat[j_new][i_new] = sol;
    return mat;
}

vector<vector<int>> createFirstSol(int V) {
    vector<vector<int>> matrix(V, vector<int>(V));
    for (int i = 0; i < V; i++) {
        matrix[i][i] = -1;
    }
    return matrix;

}

int checkConstrain(vector<vector<int>> M, int V){
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
    }
    return 1;
}

vector<vector<int>> addConstraints(vector<vector<int>> M, int V) {
    int flag = 1;
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
                        flag = 1;
                    }
                }
            }
            else if (count_one + count_zero == 2) {
                for (int j = 0; j < V; j++) {
                    if (M[i][j] == 0) {
                        M[i][j] = 1;
                        M[j][i] = 1;
                        flag = 1;
                    }
                }
            }
        }
    }
    return M;
}


int ComputeLowerBound(vector<vector<int>> M, int V, vector<vector<int>> C) {

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

int findBestArc(vector<vector<int>> M, vector<vector<int>> C, int V) {
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

int findSecondBestArc(vector<vector<int>> M, vector<vector<int>> C, int V) {
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

