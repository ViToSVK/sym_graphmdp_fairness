#include <iostream>
#include <algorithm>
#include <bits/stdc++.h>
#include <time.h>
#include <fstream>

#include "mec.h"
#include "parser.h"

using namespace std;

bool compare_vectors(vector<vector<int> > *v1, vector<vector<int> > *v2)
{
	if(v1->size() != v2->size())
	{
		return false;
	}
	
	for(int i=0;i<v1->size();i++)
	{
		if (v1->at(i).size() != v2->at(i).size())
			return false;
		for (int j = 0; j < v1->at(i).size(); j++)
		{
			if (v1->at(i)[j] != v1->at(i)[j])
				return false;
		}
	}
	return true;
}

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
    
    if (argc != 4) {
        cerr << "Expected three arguments -- mode, filename, randompercentage" << endl;
        return 1;
    }
    
    istringstream ss(argv[1]);
    int mode;
    if (!(ss >> mode)) {
        cerr << "Invalid number " << argv[1] << endl;
        return 1;
    }
    if (mode < 0 || mode > 5) {
        cerr << "Mode number has to be between 0 and 5" << endl;
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

    istringstream sss(argv[3]);
    int Rpercent;
    if (!(sss >> Rpercent)) {
        cerr << "Invalid number " << argv[3] << endl;
        return 1;
    }
    if (Rpercent < 0 || Rpercent > 50) {
        cerr << "Random percent has to be between 0 and 50" << endl;
        return 1;
    }
    
    mec M;
    parser P;
    if (P.product(filename, M, Rpercent) != 0) { // P.product(filename, M, Rpercent) != 0     P.basic(M) != 0
        cerr << "Tried to parse a product file and failed" << endl;
        return 1;
    }
    
    ofstream fout;
    string resultsname = "results/";
    
    resultsname.append(to_string(Rpercent));
    
    switch(mode) {
        case 0 :
            resultsname.append("results0.txt");
            break;
        case 1 :
            resultsname.append("results1.txt");
            break;
        case 2 :
            resultsname.append("results2.txt");
            break;
        case 3 :
            resultsname.append("results3.txt");
            break;
        case 4 :
            resultsname.append("results4.txt");
            break;
        default :
            resultsname.append("results.txt");
    }
    
    fout.open(resultsname, std::ios_base::app);
    
    fout << setw(40) << "name"        << setw(16) << "states" << setw(16) << "edges" << setw(16) << "pairs" << endl;
    fout << setw(40) << shortfilename << setw(16) << M.n      << setw(16) << M.m     << setw(16) << M.k     << endl;
    fout << setw(40) << "parameter"
         << setw(16) << "MECtime"  << setw(16) << "MECpre"  << setw(16) << "MECpost"
         << setw(16) << "MECprepost" << setw(16) << "MECclassic" << setw(16) << "MEClss"
         << setw(16) << "STRtime" << setw(16) << "STRpre" << setw(16) << "STRpost"
         << setw(16) << "STRprepost" << setw(16) << "STRclassic" << setw(16) << "STRlss" << endl;
    
    switch(mode) {
        case 0 :
            M.streettMDPbasic();
            fout << setw(40) << "basic";
            break;
        case 1 :
            M.threshold = 0;
            M.ImprovedStreettMDP();
            fout << setw(40) << "impr0";
            break;
        case 2:
            M.threshold = (int) (ceil(   log((double) M.m)   ));
            M.ImprovedStreettMDP();
            fout << setw(40) << "logm";
            break;
        case 3:
            M.threshold = (int) (ceil(   sqrt(sqrt((double) M.m ))   ));
            M.ImprovedStreettMDP();
            fout << setw(40) << "sqrt(sqrt(m))";
            break;
        case 4:
            M.threshold = (int) (ceil(   sqrt( ((double) M.m) / (sqrt((double) M.n)) )   )); // (int) (ceil(   sqrt( ((double) M.m) / (log((double) M.m)) )   ));
            M.ImprovedStreettMDP();
            fout << setw(40) << "sqrt(m/(sqrt(n)))"; // "sqrt(m/(log(m)))"
            break;
        default:
            fout.close();
            M.streettMDPbasic();
            cout << "BASIC DONE!" << endl;
            vector<int> *goodB = M.read_answer1();
            sort(goodB->begin(), goodB->end());
            M.reset(2);
            
            M.threshold = 0;
            M.ImprovedStreettMDP();
            cout << "IMPROVED WITHOUT LSS DONE!" << endl;
            vector<int> *goodI = M.read_answer1();
            sort(goodI->begin(), goodI->end());
            M.reset(2);
            
            M.threshold = (int) (ceil(   M.n   ));
            M.ImprovedStreettMDP();
            cout << "IMPROVED DONE! MEC-LSS count: " << M.FMECcountLSS << " STR-LSS count: " << M.STRcountLSS << endl;
            vector<int> *goodIL = M.read_answer1();
            sort(goodIL->begin(), goodIL->end());
            
            if (!compare_vectors1(goodB, goodI) || !compare_vectors1(goodI, goodIL)) {
                cout << "DIFFERENT ANSWERS" << endl;
                print1(goodB);
                cout << endl;
                print1(goodI);
                cout << endl;
                print1(goodIL);
            } else cout << "SAME ANSWERS" << endl;
            return 0;
    }
    
    fout << setw(16) << M.FMECtime  << setw(16) << M.FMECpresteps  << setw(16) << M.FMECpoststeps
         << setw(16) << M.FMECpresteps+M.FMECpoststeps << setw(16) << M.FMECcountClassic << setw(16) << M.FMECcountLSS
         << setw(16) << M.STRtime << setw(16) << M.STRpresteps << setw(16) << M.STRpoststeps
         << setw(16) << M.STRpresteps+M.STRpoststeps << setw(16) << M.STRcountClassic << setw(16) << M.STRcountLSS << endl;
    fout.close();
    
    return 0;
}
