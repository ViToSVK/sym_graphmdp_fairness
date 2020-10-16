#ifndef _MEC_H
#define _MEC_H

#include <vector>
#include <queue>
#include <stack>
#include <set>

#include "dd.h"

class mec {
	
	public:
	            
	// Internal Variables		
	DdManager* dd;
	int k;                                  // Number of LU-pairs
	int n, p1, p2;                          // Number of All/P1/PR states
	int m;                                  // Number of transitions
	int no_vars;
	
	int threshold;
        std::set<int> selfloops;
	
	//Storage
	std::vector <bdd_ptr> Lstreett_pair;	// L-pairs
	std::vector <bdd_ptr> Ustreett_pair;	// U-pairs
	bdd_ptr trans;				// Transition Function
	bdd_ptr V1;                             // V1 states
	bdd_ptr VR;                             // VR states
        bdd_ptr universe;			// All states
	std::vector <bdd_ptr> var;		// Variables
	std::vector <bdd_ptr> next_var;
	std::vector <bdd_ptr> state;		// Vector of states
	std::vector <bdd_ptr> next_state;
	bdd_ptr cube;				// set of original state, for existential quantification
	bdd_ptr next_cube;
	
	std::vector <bdd_ptr> allMECs;		// MECs
	std::vector <bdd_ptr> sccs;		// SCCs
	bdd_ptr goodEC;
	bdd_ptr answer;
	
	//Statistical Data
        int presteps, poststeps, countClassic, countLSS;
        int FMECpresteps, FMECpoststeps, FMECcountClassic, FMECcountLSS;
        int STRpresteps, STRpoststeps, STRcountClassic, STRcountLSS;
        double FMECtime, STRtime;
        int fmecaftersccpre, fmecaftersccpost;
	clock_t debugtime, fmecafterscctime;
        bool firstmec;
        
	//Initialization
	void init(int _p1, int _p2, int _k);
	void reset(int level = 0);
	void add_transition(int i, int j);
	void addLStates(int s, int i);
	void addUStates(int s, int i);
	
	// Symbolic Primitive Functions
	bdd_ptr pre(bdd_ptr set_states);
	bdd_ptr post(bdd_ptr set_states);
	bdd_ptr cpre(bdd_ptr set_states, bdd_ptr V);
	bdd_ptr reach(bdd_ptr g, bdd_ptr V);
	bdd_ptr attr(bdd_ptr g, bdd_ptr V);
	bdd_ptr compute_image(DdManager* dd, bdd_ptr set_states, bdd_ptr trans, bdd_ptr cube);
        
	bdd_ptr pick(bdd_ptr B);
	void ImprovedSCCFind(bdd_ptr S, bdd_ptr U, bdd_ptr s);
	void ImprovedSkelFwd(bdd_ptr S, bdd_ptr U, bdd_ptr s, bdd_ptr& FWSet, bdd_ptr& NewSet, bdd_ptr& NewState, bdd_ptr& P);
	void ImprovedSCCFind();
	
	void MECbasic();
	void MECbasic(bdd_ptr input);
	void ImprovedMEC();
	void ImprovedMEC(bdd_ptr input);
	void streettMDPbasic();
	void ImprovedStreettMDP();
	std::vector<std::vector<int> >* read_answer();
	std::vector <int> * read_answer1();
	
	bdd_ptr lockStepSearch(bdd_ptr S, bdd_ptr &H, bdd_ptr &T, std::vector <int> *vh, std::vector <int> *vt, bool &bottomscc);
	std::vector <int>* newsize1(bdd_ptr Src, long threshold, bool& classical);
	std::vector <int>* newsize(bdd_ptr Src, long threshold, bool& classical);
        
};

#endif
