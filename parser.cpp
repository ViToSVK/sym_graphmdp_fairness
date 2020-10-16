#include <iostream>
#include <string>
#include <fstream>

#include "parser.h"

using namespace std;

/* Receives input from the standard input,
 * interprets it as an MDP description in the basic format */
int parser::basic(mec &M) {
	int p1,p2,m,k;// P1states, PRstates, transitions, pairs
        m = 0;
	
	cin >> p1 >> p2 >> k;
	M.init(p1, p2, k);	
        M.n = p1+p2;
        
	int s1,s2;
	while(cin>>s1){
		if(s1==-1){
			break;
		}
		cin>>s2;
		M.add_transition(s1,s2);
		m++;
	}
        M.m = m;
	for (int i = 0; i < k; i++)
	{
		int l, u;
		cin>> l>> u;
		for (int j = 0; j < l; j++)
		{
			int s;
			cin>> s;
			M.addLStates(s, i);
		}
		for (int j = 0; j < u; j++)
		{
			int s;
			cin>> s;
			M.addUStates(s, i);
		}
	}
}

/* Receives product file provided as an argv,
 * interprets it as an MDP description in the basic format */
int parser::product(string filename, mec &M, int Rpercent) {
    
    fstream  file;
    file.open(filename.c_str(), ios::in );
    if (!file.is_open()) {
        return 1;
    }
    
    int n,p1,p2,m,k;// ALLstates, P1states, PRstates, transitions, pairs
    m = 0;
    
    file >> n >> k;
    p2 = (n * Rpercent) / 100;
    p1 = n - p2;
    M.init(p1, p2, k);
    M.n = n;
    
    int s1,s2;
    while (file >> s1) {
	if (s1==-1) {
            break;
	}
	file >> s2;
	M.add_transition(s1,s2);
	m++;
    }
    M.m = m;
    
    for (int i = 0; i < k; i++) {
	int l, u;
	file >> l >> u;
	for (int j = 0; j < l; j++) {
            int s;
            file >> s;
            M.addLStates(s, i);
	}
	for (int j = 0; j < u; j++) {
            int s;
            file >> s;
            M.addUStates(s, i);
	}
    }
    
    return 0;    
}
