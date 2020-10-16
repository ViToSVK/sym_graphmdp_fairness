#include <iostream>
#include <string>
#include <fstream>

#include "parser.h"

using namespace std;

/* Receives input from the standard input,
 * interprets it as a graph description in the basic format */
int parser::basic(streettgraph &G) {
    int n,k, ec=0;
    cin>>n>> k;
    int p1 = n / 2;
    int p2 = n - p1;
    G.init(p1, p2, k);
    G.n = n;
    int s1,s2;

    while(cin>>s1){
	if(s1==-1){
            break;
	}
	cin>>s2;
	G.add_transition(s1,s2);
	ec++;
    }
    G.m = ec;
    for (int i = 0; i < k; i++)
    {
	int l, u;
	cin>> l>> u;
	for (int j = 0; j < l; j++)
	{
            int s;
            cin>> s;
            G.addLStates(s, i);
	}
	for (int j = 0; j < u; j++)
	{
            int s;
            cin>> s;
            G.addUStates(s, i);
	}
    }
    
    return 0;
}

/* Receives product file provided as an argv,
 * interprets it as a graph description in the basic format */
int parser::product(string filename, streettgraph &G) {
    
    fstream  file;
    file.open(filename.c_str(), ios::in );
    if (!file.is_open()) {
        return 1;
    }
    
    int n=0; // number of states
    int m=0; // number of transitions
    int k=0; // number of pairs
    
    file >> n >> k;
    int p1 = n / 2;
    int p2 = n - p1;
    G.init(p1, p2, k);
    G.n = n;
    int s1,s2;

    while (file >> s1) {
	if (s1==-1) {
            break;
	}
	file >> s2;
	G.add_transition(s1,s2);
	m++;
    }
    G.m = m;
    
    for (int i = 0; i < k; i++) {
	int l, u;
	file >> l >> u;
	for (int j = 0; j < l; j++) {
            int s;
            file >> s;
            G.addLStates(s, i);
	}
	for (int j = 0; j < u; j++) {
            int s;
            file >> s;
            G.addUStates(s, i);
	}
    }
    
    return 0;    
}

/* Receives input from the file provided as an argv,
 * interprets it as a graph description in the HOA format */
int parser::hoa(string filename, streettgraph &G) {
    
    fstream  file;
    file.open(filename.c_str(), ios::in );
    if (!file.is_open()) {
        return 1;
    }
    
    int n=0; // number of states
    int m=0; // number of transitions
    int k=0; // number of pairs
    string line;
    
    bool found = false;
    string prefix("States: ");
    while (!found) {
        getline (file, line);
        if (!line.compare(0,prefix.size(),prefix)) {
            found = true;
            int i = line.size()-1;
            while (line[i] != ' ') i--;
            n = stoi(line.substr(i+1));
            //cout << n << endl;
        }
    }
    
    found = false;
    prefix = "acc-name: Rabin ";
    while (!found) {
        getline (file, line);
        if (!line.compare(0,prefix.size(),prefix)) {
            found = true;
            int i = line.size()-1;
            while (line[i] != ' ') i--;
            k = stoi(line.substr(i+1));
            //cout << k << endl;
        }
    }
    
    int p1 = n / 2;
    int p2 = n - p1;
    G.init(p1, p2, k);
    G.n = n;
    
    found = false;
    prefix = "--BODY--";
    while (!found) {
        getline (file, line);
        if (!line.compare(0,prefix.size(),prefix))
            found = true;
    }
    
    found = false;
    prefix = "--END--";
    int s1, s2, s3;
    while (!found) {
        getline (file, line);
        if (!line.compare(0,prefix.size(),prefix))
            found = true;
        else if (!line.compare(0,7,"State: ")) { // new state
            int i = 0;
            while (line[7+i] != ' ') i++;
            s1 = stoi(line.substr(7,i));
            if (s1 < 0 || s1 >= n) return 1;
            //cout << s1 << endl;
            
            i = line.size() - 1;
            while (line[i] != '{') i--; i++;
            int j = i;
            while (line[j] != '}') j++; j++;
            string numbers = line.substr(i, j-i); // take numbers and '}'
            //cout << numbers << endl;
            
            i = 0; j = 0;
            bool finished = false;
            while (!finished) {
                i = j;
                while (numbers[i] == ' ') i++;
                if (numbers[i] == '}') finished = true;
                else { // a number follows
                    j = i;
                    while (numbers[j] != ' ' && numbers[j] != '}') j++;
                    s3 = stoi(numbers.substr(i, j-i));
                    //cout << "NUMBER: " << s3 << endl;
                
                    if ((int) s3/2 >= k) return 1;
                    if (s3 % 2 == 0) // even -> U_(s3 / 2)
                        G.addUStates(s1, (int) s3/2);
                    else             // odd  -> L_(s3 / 2)
                        G.addLStates(s1, (int) s3/2);
                }
            }

        } else { // new transition
            int i = line.size()-1;
            while (line[i] != ' ') i--; i++;
            s2 = stoi(line.substr(i));
            
            //cout << "Transition from " << s1 << " to " << s2 << endl;
            G.add_transition(s1, s2);
            m++;
        }
    }
    G.m = m;
    
    file.close();
    return 0;
}
