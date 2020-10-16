#include <iostream>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <bits/stdc++.h>
#include <cmath>
#include <sstream>
#include <string>

#include "streettgraph.h"
#include "parsergraph.h"

using namespace std;

bool compare_vectors1(vector<int>  *v1, vector<int>  *v2)
{
  if(v1->size() != v2->size())
    return false;
  for (int j = 0; j < v1->size(); j++)
    if (v1->at(j) != v1->at(j))
      return false;
  return true;
}

void print1(vector<int>* v)
{
  for (int i = 0; i < v->size(); i++)
    cout<< (*v)[i]<< " ";
  cout<< endl;
}

int main(int argc, char *argv[])
{
    streettgraph G;
    parser P;
    clock_t init1,init2;
    double time1,time2;
    ofstream fout;
    string resultsname = "";

    if (argc != 3) {
        cerr << "Expected two arguments -- mode, filename" << endl;
        return 1;
    }
    istringstream ss(argv[1]);
    int mode;
    if (!(ss >> mode)) {
        cerr << "Invalid number " << argv[1] << endl;
        return 1;
    }
    if (mode < 0 || mode > 6) {
        cerr << "Mode number has to be between 0 and 6" << endl;
        return 1;
    }
    string filename = argv[2];
    int p = 0;
    while (filename[p] != '/') {
        p++;
        if (p == filename.length()-1) {
            p = -1;
            break;
        }
    }
    string shortfilename = filename.substr(p+1);
    if (P.product(filename, G) != 0) { //P.hoa(filename, G)   //P.product(filename, G)
        cerr << "Tried to parse a product file and failed" << endl;
        return 1;
    }

    switch(mode) {
        case 0 :
            resultsname = "results/results0.txt";
            break;
        case 1 :
            resultsname = "results/results1.txt";
            break;
        case 2 :
            resultsname = "results/results2.txt";
            break;
        case 3 :
            resultsname = "results/results3.txt";
            break;
        case 4 :
            resultsname = "results/results4.txt";
            break;
        case 5 :
            resultsname = "results/results5.txt";
            break;
        default :
            resultsname = "results/results.txt";
    }

    fout.open(resultsname, std::ios_base::app);

    fout << setw(40) << "name"        << setw(16) << "states" << setw(16) << "edges" << setw(16) << "pairs" << endl;
    fout << setw(40) << shortfilename << setw(16) << G.n      << setw(16) << G.m     << setw(16) << G.k     << endl;
    fout << setw(40) << "parameter" << setw(16) << "time"  << setw(16) << "timestr"  << setw(16) << "pre" << setw(16) << "post" << setw(16) << "prepost"
         << setw(16) << "mainloopcnt" << setw(16) << "size0cnt" << setw(16) << ">Tcnt" << setw(16) << "<Tcnt" << endl;

    int t = 0;
    switch(mode) {
        case 0 :
            // BASIC
            init1 = clock();
            G.ImprovedSCCFind();
            init2 = clock();
            G.BasicSTROBJ(false);
            time2 = (double)(clock() - init2)/CLOCKS_PER_SEC;
            time1 = (double)(clock() - init1)/CLOCKS_PER_SEC;
            fout << setw(40) << "basic";
            break;
        case 1 :
            // BASIC + BADLOOP
            init1 = clock();
            G.ImprovedSCCFind();
            init2 = clock();
            G.BasicSTROBJ(true);
            time2 = (double)(clock() - init2)/CLOCKS_PER_SEC;
            time1 = (double)(clock() - init1)/CLOCKS_PER_SEC;
            fout << setw(40) << "basicBadLoop";
            break;
        case 2 :
            // IMP 0 (doesn't check)
            t = 0;
            init1 = clock();
            G.ImprovedSCCFind();
            init2 = clock();
            G.ImprovedSTROBJ(t, false);
            time2 = (double)(clock() - init2)/CLOCKS_PER_SEC;
            time1 = (double)(clock() - init1)/CLOCKS_PER_SEC;
            fout << setw(40) << "impr0nocheck";
            break;
        case 3 :
            // IMP 0 (checksedge)
            t = 0;
            init1 = clock();
            G.ImprovedSCCFind();
            init2 = clock();
            G.ImprovedSTROBJ(t, true);
            time2 = (double)(clock() - init2)/CLOCKS_PER_SEC;
            time1 = (double)(clock() - init1)/CLOCKS_PER_SEC;
            fout << setw(40) << "impr0";
            break;
        case 4 :
            // IMP log m (checksedge)
            t = (int) (ceil(   log((double) G.m)   ));
            init1 = clock();
            G.ImprovedSCCFind();
            init2 = clock();
            G.ImprovedSTROBJ(t, true);
            time2 = (double)(clock() - init2)/CLOCKS_PER_SEC;
            time1 = (double)(clock() - init1)/CLOCKS_PER_SEC;
            fout << setw(40) << "logm";
            break;
        case 5 :
            // IMP 4 log m (checksedge)
            t = (int) (ceil(   4 * log((double) G.m)   ));
            init1 = clock();
            G.ImprovedSCCFind();
            init2 = clock();
            G.ImprovedSTROBJ(t, true);
            time2 = (double)(clock() - init2)/CLOCKS_PER_SEC;
            time1 = (double)(clock() - init1)/CLOCKS_PER_SEC;
            fout << setw(40) << "4logm";
            break;
        default :
            fout.close();

            G.ImprovedSCCFind();
            G.BasicSTROBJ(false);
            cout << "BASIC DONE! PREPOST " << G.presteps+G.poststeps << endl << flush;
            vector <int> * ans1 = G.read_answer(G.good_components);
            G.reset_sccs(true);

            G.ImprovedSCCFind();
            G.BasicSTROBJ(true);
            cout << "BASIC + BADLOOP DONE! PREPOST " << G.presteps+G.poststeps << endl << flush;
            vector <int> * ans2 = G.read_answer(G.good_components);
            G.reset_sccs(true);

            t = 0;
            G.ImprovedSCCFind();
            G.ImprovedSTROBJ(t, false);
            cout << "IMPROVED WITHOUT LSS DONE! PREPOST " << G.presteps+G.poststeps << endl << flush;
            vector <int> * ans3 = G.read_answer(G.good_components);
            G.reset_sccs(true);

            t = 0;
            G.ImprovedSCCFind();
            G.ImprovedSTROBJ(t, true);
            cout << "IMPROVED + EDGECHECK WITHOUT LSS DONE! PREPOST " << G.presteps+G.poststeps << endl << flush;
            vector <int> * ans4 = G.read_answer(G.good_components);
            G.reset_sccs(true);

            t = G.n;
            G.ImprovedSCCFind();
            G.ImprovedSTROBJ(t, false);
            cout << "IMPROVED DONE! LSS count: " << G.case3 << endl << flush;
            vector <int> * ans5 = G.read_answer(G.good_components);
            G.reset_sccs(true);

            if (!compare_vectors1(ans1, ans2) || !compare_vectors1(ans1, ans3) ||
                !compare_vectors1(ans1, ans4) || !compare_vectors1(ans1, ans5)) {
                cout << "DIFFERENT ANSWERS" << endl;
                print1(ans1);
                cout << endl;
                print1(ans2);
                cout << endl;
                print1(ans3);
                cout << endl;
                print1(ans4);
                cout << endl;
                print1(ans5);
            } else cout << "SAME ANSWERS" << endl;

            return 0;
    }

    fout << setw(16) << time1 << setw(16) << time2 << setw(16) << G.presteps << setw(16) << G.poststeps
         << setw(16) << G.presteps+G.poststeps << setw(16) << G.mainLoopCounter
         << setw(16) << G.case1 << setw(16) << G.case2 << setw(16) << G.case3 << endl;

    fout.close();
    return 0;
}
