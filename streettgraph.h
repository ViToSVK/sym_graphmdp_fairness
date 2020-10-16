#ifndef _STREETTGRAPH_H
#define _STREETTGRAPH_H

#include <queue>
#include <vector>
#include <set>
#include <stack>

#include "dd.h"

class streettgraph {
    
    public:
    
    // Internal Variables		
    DdManager* dd;
    int k;                              // number of pairs
    int m;				// number of transitions
    int n;                              // number of states
    int no_vars;
    int p1, p2;
    std::set<int> selfloops;
    
    //Storage
    std::vector <bdd_ptr> Lstreett_pair;
    std::vector <bdd_ptr> Ustreett_pair;
    std::vector <bdd_ptr> var;
    std::vector <bdd_ptr> next_var;
    std::vector <bdd_ptr> state;        // each bdd represents a state
    std::vector <bdd_ptr> next_state;
    bdd_ptr universe;			// bdd representing all states
    bdd_ptr cube;
    bdd_ptr next_cube;
    bdd_ptr trans;			// bdd representing the transition function
    bdd_ptr answer;                     // bdd representing the P1 winning region
    bdd_ptr good_components;		// bdd representing the union of good components
    
    std::vector <bdd_ptr> sccs;		// each bdd represents a computed scc
    std::stack <bdd_ptr> candidates;	// stack of good-component candidates
    
    //Statistical Data
    int presteps, poststeps, case1, case2, case3, mainLoopCounter;
    int sym2, sym3, curpre, curpost;
    
    //Initialization
    void init(int _p1, int _p2, int _k);
    void reset_sccs(bool comp = false);
    void add_transition(int i, int j);
    void addLStates(int s, int i);
    void addUStates(int s, int i);
    
    bdd_ptr pre(bdd_ptr set_states);
    bdd_ptr post(bdd_ptr set_states);
    bdd_ptr compute_image(bdd_ptr set_states, bdd_ptr trans, bdd_ptr cube);
    bdd_ptr reach(bdd_ptr g, bdd_ptr V);
    
    bdd_ptr pick(bdd_ptr B);
    void ImprovedSCCFind();
    void ImprovedSCCFind(bdd_ptr S, bdd_ptr U, bdd_ptr s);
    void ImprovedSkelFwd(bdd_ptr S, bdd_ptr U, bdd_ptr s, bdd_ptr& FWSet, bdd_ptr& NewSet, bdd_ptr& NewState, bdd_ptr& P);
    
    void temp_print(std::vector<int>* v);
    std::vector<int>* read_answer(bdd_ptr& setToCount);
    std::vector<int>* newsize(bdd_ptr Src, int threshold, bool& classical);
    std::vector<int>* newsize1(bdd_ptr Src, int threshold, bool& classical);
    
    bdd_ptr lockStepSearch(bdd_ptr S, bdd_ptr& H, bdd_ptr& T, std::vector<int>* vh, std::vector<int>* vt, bool& bottomscc);
    
    void BasicSTROBJ(bool badloop);
    void ImprovedSTROBJ(int threshold, bool checkedge);
    
};

#endif
