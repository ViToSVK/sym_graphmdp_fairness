#include <vector>
#include <stack>
#include <algorithm>
#include <bits/stdc++.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <queue>

#include "mec.h"

using namespace std;

/* The non-symbolic way of computing size:
   This is not used in the code
   It is just to provide a better understanding of what the symbolic procedure below does */
std::vector<int>* mec::newsize1(bdd_ptr Src, long threshold, bool& classical) {
  std::vector<int>* ans = new std::vector<int>();
  bdd_ptr temp;
  classical = false;
  for(int i=0;i<state.size(); i++) {
    temp = bdd_and(dd,Src,state[i]);
    if(temp != Cudd_ReadLogicZero(dd)) {
      ans->push_back(i);
      if(ans->size() > threshold) { // was ans->size() * ans->size() > threshold
        classical = true;
        break;
      }
    }
    Cudd_RecursiveDeref(dd,temp);
  }
  return ans;
}

/* Symbolic Size Computation, checks whether the size of Src is more than threshold
   if it is at most threshold, returns the states enumerated and sets classical to false, else sets classical to true */
std::vector<int>* mec::newsize(bdd_ptr Src, long threshold, bool& classical) {
  DdGen *gen;
  int *cube;
  CUDD_VALUE_TYPE value;
  int global_cnt = 0;
  std::vector<int> position2;
  std::vector<int>* ans = new std::vector<int>();
  int index;
  bool is_player_random = false;
  classical = false;

  Cudd_ForeachCube(dd, Src, gen, cube, value) {
    index = 0;
    position2.clear();
    for(int i=0;i<no_vars;i++) {
      if(cube[i] == 2) {
        position2.push_back(i);
      }
      else if(cube[i] == 1) {
        if(i == 0) {
          is_player_random = true;
        }
        else
          index += (1 << (i-1));
      }
    }

    global_cnt += (1 << position2.size());
    if(global_cnt > threshold) { // was global_cnt * global_cnt > threshold
      classical = true;
      break;
    }

    for(int i=0; i < (1 << position2.size()); i++) {
      int currindex = index;
      int j=i;
      int bit_no = 0;
      while(j) {
        int c = j%2;
        j>>=1;
        if(c)
          currindex += (1 << (position2[bit_no]-1));
        bit_no++;
      }
      if(is_player_random)
        currindex += p1;
      if(currindex >= p1+p2) {    // ERROR
        std::cerr << "Cube: ";
        for(int x=0;x<no_vars; x++)
          std::cerr << cube[x] << " ";
        std::cerr << "Index: " << index << std::endl;
        std::cerr << "Position2: ";
        for(int x=0;x<position2.size(); x++)
          std::cerr << position2[x] << " ";
        std::cerr << std::endl;
        std::cerr << "Assignment(i): " << i << std::endl;
        std::cerr << "Currindex: " << currindex << std::endl;
        std::cerr << "p1,p2: "<< p1 << "," << p2 << std::endl;
      }
      assert(currindex < p1+p2);
      ans->push_back(currindex);
    }
  }
  return ans;
}

/* Symbolic Helper Function for PRE and POST: Computing Image */
bdd_ptr mec::compute_image(DdManager* dd, bdd_ptr set_states, bdd_ptr trans, bdd_ptr cube) {
  bdd_ptr result,temp;
  result = bdd_and(dd,set_states,trans);
  temp = bdd_forsome(dd,result,cube);
  Cudd_RecursiveDeref(dd,result);
  result = temp;
  return result;
}

/* Symbolic Primitive: PRE */
bdd_ptr mec::pre(bdd_ptr set_states) {
  presteps++;

  bdd_ptr result = set_states;
  Cudd_Ref(result);
  bdd_ptr temp;
  for(int i=0;i<var.size();i++) {
      /* compose replaces var with next_var */
      temp = bdd_compose(dd,result,next_var[i],var[i]->index);
      Cudd_RecursiveDeref(dd,result);
      result = temp;
  }

  temp = compute_image(dd,result,trans,next_cube);
  Cudd_RecursiveDeref(dd,result);
  result = temp;
  return result;
}

/* Symbolic Primitive: POST */
bdd_ptr mec::post(bdd_ptr set_states) {
  poststeps++;

  bdd_ptr result = set_states;
  Cudd_Ref(result);
  bdd_ptr temp;
  temp = compute_image(dd,result,trans,cube);
  Cudd_RecursiveDeref(dd,result);
  result = temp;

  for(int i=0;i<next_var.size();i++) {
      /* compose replaces var with next_var */
      temp = bdd_compose(dd,result,var[i],next_var[i]->index);
      Cudd_RecursiveDeref(dd,result);
      result = temp;
  }
  return result;
}

/* CPRE */
bdd_ptr mec::cpre(bdd_ptr set_states, bdd_ptr V) {

    // ViTo 25.9.
    bdd_ptr temp1, temp2, temp3;

    temp2 = bdd_not(dd,set_states);
    temp3 = bdd_and(dd,V,temp2);
    Cudd_RecursiveDeref(dd,temp2);
    // temp3 = V \ Z

    temp2 = pre(temp3);
    Cudd_RecursiveDeref(dd,temp3);
    // temp2 = Pre(V \ Z)

    temp3 = bdd_and(dd,V1,temp2);
    Cudd_RecursiveDeref(dd,temp2);
    // temp3 = V1 intersect (Pre(V \ Z))

    temp2 = bdd_not(dd,temp3);
    Cudd_RecursiveDeref(dd,temp3);
    // temp2 = Not ( V1 intersect (Pre(V \ Z)) )

    temp1 = pre(set_states);
    // temp1 = Pre(Z)

    bdd_ptr result = bdd_and(dd,temp1,temp2);
    Cudd_RecursiveDeref(dd,temp1);
    Cudd_RecursiveDeref(dd,temp2);
    // result = Pre(Z) \ ( V1 intersect (Pre(V \ Z)) )

    return result;
}

/* Reverse Reachability */
bdd_ptr mec::reach(bdd_ptr g, bdd_ptr V){
  bdd_ptr Rnew,Rold,Re,preRold;
  Rnew = bdd_and(dd,g,V);
  Rold = bdd_one(dd);

  do{
    Cudd_RecursiveDeref(dd,Rold);
    Rold = Rnew;
    preRold = pre(Rold);
    Re = bdd_and(dd,preRold,V);
    Cudd_RecursiveDeref(dd,preRold);
    Rnew = bdd_or(dd,Rold,Re);
    Cudd_RecursiveDeref(dd,Re);
  }while(Rold != Rnew);

  Cudd_RecursiveDeref(dd,Rold);
  return Rnew;
}

/* Attractor */
bdd_ptr mec::attr(bdd_ptr g, bdd_ptr V){
  bdd_ptr Rold,Rnew,Re,cpreRold;
  Rnew = g;
  Cudd_Ref(Rnew);
  Rold = bdd_one(dd);

  do{
    Cudd_RecursiveDeref(dd,Rold);
    Rold = Rnew;
    cpreRold = cpre(Rold,V);
    Re = bdd_and(dd,cpreRold,V);
    Cudd_RecursiveDeref(dd,cpreRold);
    Rnew = bdd_or(dd,Rold,Re);
    Cudd_RecursiveDeref(dd,Re);
  }while(Rold != Rnew);

  Cudd_RecursiveDeref(dd,Rold);
  return Rnew;
}

/* Pick subroutine - This is inexpensive and hence not a symbolic step */
bdd_ptr mec::pick(bdd_ptr B) {
  DdGen *gen;
  int *cube;
  CUDD_VALUE_TYPE value;

  gen = Cudd_FirstCube(dd, B, &cube, &value);
  bdd_ptr ans = bdd_one(dd);
  bdd_ptr temp,temp2;

  for(int i=0;i<no_vars;i++) {
    if(cube[i] == 0 || cube[i] == 2) {
      temp = bdd_new_var_with_index(dd,i);
      Cudd_Ref(temp);
      temp2 = bdd_not(dd,temp);
      Cudd_RecursiveDeref(dd,temp);
      temp = bdd_and(dd,ans,temp2);
      Cudd_RecursiveDeref(dd,ans);
      Cudd_RecursiveDeref(dd,temp2);
      ans = temp;       // ans = ans & ~var(i)
    }
    else if(cube[i] == 1) {
      temp = bdd_new_var_with_index(dd,i);
      Cudd_Ref(temp);
      temp2 = bdd_and(dd,ans,temp);
      Cudd_RecursiveDeref(dd,ans);
      Cudd_RecursiveDeref(dd,temp);
      ans = temp2;        // ans = ans & var(i)
    }
    else
      assert(0 && "Cube has integer other than 0,1,2");
  }
  return ans;
}

/* Initialization */
void mec::init(int _p1, int _p2, int _k) {

  p1 = _p1;
  p2 = _p2;
  k = _k;
  presteps = 0; poststeps = 0; countClassic = 0; countLSS = 0;
        FMECpresteps = 0; FMECpoststeps = 0; FMECcountClassic = 0; FMECcountLSS = 0;
        STRpresteps = 0; STRpoststeps = 0; STRcountClassic = 0; STRcountLSS = 0;
        selfloops.clear();

  bdd_ptr temp,temp2;

  dd = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0);
  if(dd == NULL) {
    std::cerr<<"DdManager failed to initialize"<<std::endl;
  }

  no_vars = 1 + (int)std::ceil(std::log((double)std::max(p1,p2))  /   std::log(2)  );
  int cnt;
  for(cnt=0;cnt<no_vars;cnt++) {
    var.push_back(bdd_new_var_with_index(dd,cnt));
  }
  for(int j=0;j<no_vars;j++,cnt++) {
    next_var.push_back(bdd_new_var_with_index(dd,cnt));
  }

  V1 = bdd_zero(dd);
  VR = bdd_zero(dd);

  for(int i=0;i<p1;i++) {
    bdd_ptr s = bdd_not(dd,var[0]);
    int k = i;
    int j;
    for(j=1;k;j++) {
      int c = k%2;
      k = k>>1;
      assert(j < var.size());

      if(c == 0)
        temp = bdd_not(dd,var[j]);
      else {
        temp = var[j];
        Cudd_Ref(temp);
      }

      temp2 = bdd_and(dd,s,temp);
      Cudd_RecursiveDeref(dd,s);
      Cudd_RecursiveDeref(dd,temp);
      s = temp2;
    }

    while(j < var.size()) {
      temp = bdd_not(dd,var[j]);
      temp2 = bdd_and(dd,s,temp);
      Cudd_RecursiveDeref(dd,s);
      Cudd_RecursiveDeref(dd,temp);
      s = temp2;
      j++;
    }
    temp = bdd_or(dd, V1, s);
    Cudd_RecursiveDeref(dd, V1);
    V1 = temp;
    state.push_back(s);
  }

  for(int i=0;i<p2;i++) {
    bdd_ptr s = var[0];
    Cudd_Ref(s);
    int k = i;
    int j;
    for(j=1;k;j++) {
      int c = k%2;
      k = k>>1;

      if(c == 0)
        temp = bdd_not(dd,var[j]);
      else {
        temp = var[j];
        Cudd_Ref(temp);
      }

      temp2 = bdd_and(dd,s,temp);
      Cudd_RecursiveDeref(dd,s);
      Cudd_RecursiveDeref(dd,temp);
      s = temp2;
    }

    while(j < var.size()) {
      temp = bdd_not(dd,var[j]);
      temp2 = bdd_and(dd,s,temp);
      Cudd_RecursiveDeref(dd,s);
      Cudd_RecursiveDeref(dd,temp);
      s = temp2;
      j++;
    }
    temp = bdd_or(dd, VR, s);
    Cudd_RecursiveDeref(dd, VR);
    VR = temp;
    state.push_back(s);
  }
  universe = bdd_zero(dd);
  for(int i=0;i<state.size(); i++) {
    temp = bdd_or(dd,universe,state[i]);
    Cudd_RecursiveDeref(dd,universe);
    universe = temp;
  }

  for(int i=0;i<p1;i++) {
    bdd_ptr s = bdd_not(dd,next_var[0]);
    int k = i;
    int j;
    for(j=1;k;j++) {
      int c = k%2;
      k = k>>1;

      if(c == 0)
        temp = bdd_not(dd,next_var[j]);
      else {
        temp = next_var[j];
        Cudd_Ref(temp);
      }

      temp2 = bdd_and(dd,s,temp);
      Cudd_RecursiveDeref(dd,s);
      Cudd_RecursiveDeref(dd,temp);
      s = temp2;
    }

    while(j < next_var.size()) {
      temp = bdd_not(dd,next_var[j]);
      temp2 = bdd_and(dd,s,temp);
      Cudd_RecursiveDeref(dd,s);
      Cudd_RecursiveDeref(dd,temp);
      s = temp2;
      j++;
    }
    next_state.push_back(s);
  }

  for(int i=0;i<p2;i++) {
    bdd_ptr s = next_var[0];
    Cudd_Ref(s);
    int k = i;
    int j;
    for(j=1;k;j++) {
      int c = k%2;
      k = k>>1;

      if(c == 0)
        temp = bdd_not(dd,next_var[j]);
      else {
        temp = next_var[j];
        Cudd_Ref(temp);
      }

      temp2 = bdd_and(dd,s,temp);
      Cudd_RecursiveDeref(dd,s);
      Cudd_RecursiveDeref(dd,temp);
      s = temp2;
    }

    while(j < next_var.size()) {
      temp = bdd_not(dd,next_var[j]);
      temp2 = bdd_and(dd,s,temp);
      Cudd_RecursiveDeref(dd,s);
      Cudd_RecursiveDeref(dd,temp);
      s = temp2;
      j++;
    }
    next_state.push_back(s);
  }

  cube = bdd_one(dd);
  for(int i=0;i<var.size();i++) {
    temp = bdd_and(dd,cube,var[i]);
    Cudd_RecursiveDeref(dd,cube);
    cube = temp;
  }

  next_cube = bdd_one(dd);
  for(int i=0;i<next_var.size();i++) {
    temp = bdd_and(dd,next_cube,next_var[i]);
    Cudd_RecursiveDeref(dd,next_cube);
    next_cube = temp;
  }
  for (int i = 0; i < k; i++)
  {
    temp = bdd_zero(dd);
    Lstreett_pair.push_back(temp);
    temp = bdd_zero(dd);
    Ustreett_pair.push_back(temp);
  }
  trans = bdd_zero(dd);
  goodEC = bdd_zero(dd);
  return;
}

/* Resetting SCC Information to run different algorithm */
void mec::reset(int level) {
  if (level > 1)
  {
    Cudd_RecursiveDeref(dd, goodEC);
    //Cudd_RecursiveDeref(dd, answer);
    goodEC = bdd_zero(dd);
  }
  if (level > 0)
  {
    for (int i = 0; i < allMECs.size(); i++)
      Cudd_RecursiveDeref(dd, allMECs[i]);
    allMECs.clear();
  }
  for(int i=0;i<sccs.size();i++)
    Cudd_RecursiveDeref(dd,sccs[i]);
  sccs.clear();
  return;
}

std::vector<std::vector<int> >* mec::read_answer()
{
  std::vector<std::vector<int> >* ans = new std::vector<std::vector<int> >();
  bdd_ptr temp;
  for (int j = 0; j < allMECs.size(); j++)
  {
    std::vector<int>* single_components = new std::vector<int>();
    for(int i=0;i<state.size();i++)
    {
      temp = bdd_and(dd, state[i], allMECs[j]);
      if (temp != Cudd_ReadLogicZero(dd))
        single_components->push_back(i);
      Cudd_RecursiveDeref(dd, temp);
    }
    sort(single_components->begin(), single_components->end());
    ans->push_back(*single_components);
  }
  sort(ans->begin(), ans->end());
  return ans;
}

std::vector <int> * mec::read_answer1()
{
  std::vector <int> * ans = new std::vector<int> ();
  bdd_ptr temp;
  for(int i=0;i<state.size();i++)
  {
    temp = bdd_and(dd, state[i], goodEC);
    if (temp != Cudd_ReadLogicZero(dd))
      ans->push_back(i);
    Cudd_RecursiveDeref(dd, temp);
  }
  return ans;
}


/* Adding Transitions */
void mec::add_transition(int i, int j) {
  bdd_ptr res,temp;
  res = bdd_and(dd,state[i],next_state[j]);
  temp = bdd_or(dd,trans,res);
  Cudd_RecursiveDeref(dd,trans);
  Cudd_RecursiveDeref(dd,res);
  trans = temp;
        if (i == j) selfloops.insert(i);
  return;
}


/* ImprovedSCCFind Procedure - Our New Improved Linear Algorithm */
void mec::ImprovedSCCFind() {
  Cudd_Ref(universe);
  ImprovedSCCFind(universe,bdd_zero(dd),bdd_zero(dd));
  return;
}

/* ImprovedSCCFind Algorithm */
void mec::ImprovedSCCFind(bdd_ptr S, bdd_ptr U, bdd_ptr s) {
  bdd_ptr FWSet, NewSet, NewState, P;
  bdd_ptr SCC;
  bdd_ptr temp,temp2,temp3;
  bdd_ptr Spr, Upr, spr;

  if(S == Cudd_ReadLogicZero(dd)) {return;}

  if(U == Cudd_ReadLogicZero(dd)) {s = pick(S);}

  ImprovedSkelFwd(S, U, s, FWSet, NewSet, NewState, P);
  if(P != Cudd_ReadLogicZero(dd)) {SCC = P; Cudd_Ref(SCC);}
  else {SCC = s; Cudd_Ref(SCC);}

  while(1) {
    temp = pre(SCC);
    temp2 = bdd_and(dd,temp,FWSet);
    Cudd_RecursiveDeref(dd,temp);
    temp = bdd_not(dd,SCC);
    temp3 = bdd_and(dd,temp2,temp);
    Cudd_RecursiveDeref(dd,temp2);
    Cudd_RecursiveDeref(dd,temp);

    if(temp3 == Cudd_ReadLogicZero(dd)) {
      Cudd_RecursiveDeref(dd,temp3);
      break;
    }
    temp = bdd_or(dd,SCC,temp3);
    Cudd_RecursiveDeref(dd,temp3);
    Cudd_RecursiveDeref(dd,SCC);
    SCC = temp;
  }
  sccs.push_back(SCC);

  temp = bdd_not(dd,FWSet);
  Spr = bdd_and(dd,S,temp);       // Spr = S \ FWSet
  Cudd_RecursiveDeref(dd,temp);
  temp = bdd_not(dd,SCC);         // temp = ~SCC
  Upr = bdd_and(dd,U,temp);       // Upr = U \ SCC

  temp2 = bdd_and(dd,SCC,U);
  temp3 = pre(temp2);
  Cudd_RecursiveDeref(dd,temp2);
  temp2 = bdd_and(dd,U,temp);
  spr = bdd_and(dd,temp2,temp3);  // spr = Pre(SCC inter U) inter (U \ SCC)

  ImprovedSCCFind(Spr, Upr, spr);

  Spr = bdd_and(dd,FWSet,temp);
  Upr = bdd_and(dd,NewSet,temp);
  spr = bdd_and(dd,NewState,temp);

  ImprovedSCCFind(Spr, Upr, spr);

  Cudd_RecursiveDeref(dd,temp);
  Cudd_RecursiveDeref(dd,FWSet);
  Cudd_RecursiveDeref(dd,NewSet);
  Cudd_RecursiveDeref(dd,NewState);
  Cudd_RecursiveDeref(dd,S);
  Cudd_RecursiveDeref(dd,U);
  Cudd_RecursiveDeref(dd,s);
  Cudd_RecursiveDeref(dd,P);
  return;
}

/* Helper Function ImprovedSkelFwd - For ImprovedSCCFind */
void mec::ImprovedSkelFwd(bdd_ptr S, bdd_ptr U, bdd_ptr s, bdd_ptr& FWSet, bdd_ptr& NewSet, bdd_ptr& NewState, bdd_ptr& P) {
  std::stack<bdd_ptr> stack;
  bdd_ptr L;
  bdd_ptr temp, temp2;
  L = s; Cudd_Ref(L);
  FWSet = bdd_zero(dd);     // FWSet = phi

  while(L != Cudd_ReadLogicZero(dd)) {
    stack.push(L);
    temp = bdd_or(dd,FWSet,L);
    Cudd_RecursiveDeref(dd,FWSet);
    FWSet = temp;                   // FWSet = FWSet union L

    temp = post(L);
    temp2 = bdd_and(dd,temp,S);     // post(L)
    Cudd_RecursiveDeref(dd,temp);
    temp = bdd_not(dd,FWSet);
    L = bdd_and(dd,temp2,temp); // post(L) \ FWSet
                Cudd_RecursiveDeref(dd,temp);   // ViTo 29.9.
    Cudd_RecursiveDeref(dd,temp2);
  }
  Cudd_RecursiveDeref(dd,L);

  P = bdd_and(dd,FWSet,U);
  L = stack.top();
        Cudd_Ref(L); // ViTo 29.9.
  stack.pop();

  NewState = pick(L); NewSet = NewState; Cudd_Ref(NewSet);
  Cudd_RecursiveDeref(dd,L);

  while(!stack.empty()) {
    L = stack.top();
    stack.pop();
    temp = bdd_and(dd,L,P);
    if(temp != Cudd_ReadLogicZero(dd)) {  // There must be an intersection
                    Cudd_RecursiveDeref(dd,temp); // ViTo 29.9
                    break;
    }
    Cudd_RecursiveDeref(dd,temp);

    temp = pre(NewSet);
    temp2 = bdd_and(dd,temp,L);
    Cudd_RecursiveDeref(dd,L);
    Cudd_RecursiveDeref(dd,temp);

    temp = pick(temp2);
    Cudd_RecursiveDeref(dd,temp2);
    temp2 = bdd_or(dd,NewSet,temp);
    Cudd_RecursiveDeref(dd,NewSet);
    Cudd_RecursiveDeref(dd,temp);
    NewSet = temp2;
  }
  return;

}

void mec::addLStates(int s, int i)
{
  bdd_ptr temp;
  temp = bdd_or(dd, state[s], Lstreett_pair[i]);
  Cudd_RecursiveDeref(dd, Lstreett_pair[i]);
  Lstreett_pair[i] = temp;

  return;
}

void mec::addUStates(int s, int i)
{
  bdd_ptr temp;
  temp = bdd_or(dd, state[s], Ustreett_pair[i]);
  Cudd_RecursiveDeref(dd, Ustreett_pair[i]);
  Ustreett_pair[i] = temp;

  return;
}

bdd_ptr mec::lockStepSearch(bdd_ptr S, bdd_ptr &H, bdd_ptr &T, std::vector <int> *vh, std::vector<int> *vt, bool &bottomscc)
{
  std::queue <bdd_ptr> Q, Q2; // Q for Heads, Q2 for Tails
  bdd_ptr temp, temp2, temp3, C;
  for (int i = 0; i < vh->size(); i++)
  {
    temp = state[vh->at(i)];
    Cudd_Ref(temp);
    Q.push(temp);
  }
  for (int i = 0; i < vt->size(); i++)
  {
    temp = state[vt->at(i)];
    Cudd_Ref(temp);
    Q2.push(temp);
  }
  while (true)
  {
    int i = 0;
    while (i < vh->size()) {
      C = Q.front();  Q.pop();
      Cudd_Ref(C); // ViTo 29.9.
      temp2 = pre(C);
      temp3 = bdd_or(dd, temp2, C);
      temp = bdd_and(dd, temp3, S);  // temp = C'
      Cudd_RecursiveDeref(dd, temp2);
      Cudd_RecursiveDeref(dd, temp3);

      std::vector <int> *aux;
      temp2 = bdd_and(dd, temp, H);
      bool classical;
      aux = newsize(temp2, 1, classical);
      Cudd_RecursiveDeref(dd, temp2);
      delete aux;

      if (classical) {
        temp3 = bdd_not(dd, state[vh->at(i)]);
        temp2 = bdd_and(dd, temp3, H);
        Cudd_RecursiveDeref(dd, temp3);
        Cudd_RecursiveDeref(dd, H);
        H = temp2;
        vh->erase(vh->begin() + i);
        // Q.push(C); Do not push anything
        // i++; Do not increase i, since there is a new element at pos i now
      } else {
        temp2 = bdd_not(dd, C);
        temp3 = bdd_and(dd, temp, temp2); // C' \ C
        if (temp3 == Cudd_ReadLogicZero(dd)) {
          Cudd_RecursiveDeref(dd, temp2);
          Cudd_RecursiveDeref(dd, temp3);
          temp2 = bdd_not(dd, temp);
          temp3 = bdd_and(dd, temp2, C); // C \ C'
          if (temp3 == Cudd_ReadLogicZero(dd)) {
            Cudd_RecursiveDeref(dd, temp2);
            Cudd_RecursiveDeref(dd, temp3);
            bottomscc = false;
            return C;
          }
        }
        Cudd_RecursiveDeref(dd, temp2);
        Cudd_RecursiveDeref(dd, temp3);
        Q.push(temp);
        i++;
      }
    }
    i = 0;
    while (i < vt->size()) {
      C = Q2.front(); Q2.pop();
      Cudd_Ref(C); // ViTo 29.9.
      temp2 = post(C);
      temp3 = bdd_or(dd, temp2, C);
      temp = bdd_and(dd, temp3, S); // temp = C'
      Cudd_RecursiveDeref(dd, temp2);
      Cudd_RecursiveDeref(dd, temp3);

      std::vector <int> *aux;
      temp2 = bdd_and(dd, temp, T);
      bool classical;
      aux = newsize(temp2, 1, classical);
      Cudd_RecursiveDeref(dd, temp2);
                  delete aux;

      if (classical) {
        temp3 = bdd_not(dd, state[vt->at(i)]);
        temp2 = bdd_and(dd, temp3, T);
        Cudd_RecursiveDeref(dd, temp3);
        Cudd_RecursiveDeref(dd, T);
        T = temp2; // was H = temp2
        vt->erase(vt->begin() + i);
        // Q2.push(C);
        // i++;
      } else {
        temp2 = bdd_not(dd, C);
        temp3 = bdd_and(dd, temp, temp2); // C' \ C
        if (temp3 == Cudd_ReadLogicZero(dd)) {
          Cudd_RecursiveDeref(dd, temp2);
          Cudd_RecursiveDeref(dd, temp3);
          temp2 = bdd_not(dd, temp);
          temp3 = bdd_and(dd, temp2, C); // C \ C'
          if (temp3 == Cudd_ReadLogicZero(dd)) {
            Cudd_RecursiveDeref(dd, temp2);
            Cudd_RecursiveDeref(dd, temp3);
            bottomscc = true;
            return C;
          }
        }
        Cudd_RecursiveDeref(dd, temp2);
        Cudd_RecursiveDeref(dd, temp3);
        Q2.push(temp);
        i++;
      }
    }
  }
}

void mec::MECbasic()
{
  Cudd_Ref(universe);
  MECbasic(universe);
  return;
}

void mec::MECbasic(bdd_ptr input)
{
  ImprovedSCCFind(input, bdd_zero(dd), bdd_zero(dd));
  if (firstmec) {
    fmecafterscctime = clock();
    fmecaftersccpre = presteps;
    fmecaftersccpost = poststeps;
  }

  std::queue <bdd_ptr> Q;
  bdd_ptr temp, temp2, temp3, C, D;
  for (int i = 0; i < sccs.size(); i++)
  {
    temp = sccs[i];
    Cudd_Ref(temp);
    Q.push(temp);
  }
  while (!Q.empty())
  {
    //cout << "start " << flush;
    //debugtime = clock();
    C = Q.front();  Q.pop();
    //Cudd_Ref(C); // ViTo 29.9. commented 2.10.

    if (firstmec) {
      bool classical;
      std::vector<int> *aux = newsize(C, 1, classical);
      if (!classical && aux->size()==1) {
        if (selfloops.count((*aux)[0]) == 0) {
          Cudd_RecursiveDeref(dd, C);
          delete aux;
          continue;
        }
      }
      delete aux;
    }

    D = bdd_and(dd, C, VR);
    temp2 = bdd_not(dd, C);
    temp3 = bdd_and(dd, universe, temp2);
    Cudd_RecursiveDeref(dd, temp2);
    temp2 = pre(temp3);
    temp = bdd_and(dd, D, temp2);
    Cudd_RecursiveDeref(dd, temp3);
    Cudd_RecursiveDeref(dd, temp2);
    Cudd_RecursiveDeref(dd, D);
    D = temp;
    if (D != Cudd_ReadLogicZero(dd))
    {
      //cout << "routnonzero(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << flush;
      //debugtime = clock();
      countClassic++;
      temp2 = attr(D, universe);
      temp3 = bdd_not(dd, temp2);
      temp = bdd_and(dd, temp3, C);
      Cudd_RecursiveDeref(dd, C);
      Cudd_RecursiveDeref(dd, temp2);
      Cudd_RecursiveDeref(dd, temp3);
      C = temp;
      reset();
      //cout << "removed_attr(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << flush;
      //debugtime = clock();
      ImprovedSCCFind(C,bdd_zero(dd),bdd_zero(dd));
      //cout << "comp_sccs(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << endl << flush;
      //debugtime = clock();
      for (int i = 0; i < sccs.size(); i++)
      {
        temp = sccs[i];
        Cudd_Ref(temp);
        Q.push(temp);
      }
    }
    else
    {
      //cout << "routzero(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << flush;
      //debugtime = clock();
      temp2 = post(C);
      temp3 = bdd_and(dd, temp2, C);
      if (temp3 != Cudd_ReadLogicZero(dd))
        allMECs.push_back(C);
      else
      //Cudd_RecursiveDeref(dd, C); commented ViTo 29.9.
      Cudd_RecursiveDeref(dd, temp2);
      Cudd_RecursiveDeref(dd, temp3);
      //cout << "checked_edge(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << endl << flush;
      //debugtime = clock();
    }
    Cudd_RecursiveDeref(dd, D);
  }

  return;
}

void mec::ImprovedMEC()
{
  Cudd_Ref(universe);
  ImprovedMEC(universe);
  return;
}

void mec::ImprovedMEC(bdd_ptr input)
{
  ImprovedSCCFind(input, bdd_zero(dd), bdd_zero(dd));
  if (firstmec) {
    fmecafterscctime = clock();
    fmecaftersccpre = presteps;
    fmecaftersccpost = poststeps;
  }

  std::queue <bdd_ptr> Q;
  bdd_ptr temp, temp2, temp3, S, A, D, T, C, H;
  for (int i = 0; i < sccs.size(); i++)
  {
    temp = sccs[i];
    Cudd_Ref(temp);
    Q.push(temp);
    temp2 = bdd_zero(dd);
    Q.push(temp2);
  }
  while (!Q.empty())
  {
    S = Q.front();  Q.pop();
    //Cudd_Ref(S); // ViTo 29.9. commented 2.10.
    T = Q.front();  Q.pop();
    //Cudd_Ref(T); // ViTo 29.9. commented 2.10.

    if (firstmec) {
      bool classical;
      std::vector<int> *aux = newsize(S, 1, classical);
      if (!classical && aux->size()==1) {
        if (selfloops.count((*aux)[0]) == 0) {
          Cudd_RecursiveDeref(dd, S);
          Cudd_RecursiveDeref(dd, T);
          delete aux;
          continue;
        }
      }
      delete aux;
    }

    D = bdd_and(dd, S, VR);
    temp2 = bdd_not(dd, S);
    temp3 = bdd_and(dd, universe, temp2);
    Cudd_RecursiveDeref(dd, temp2);
    temp2 = pre(temp3);
    temp = bdd_and(dd, temp2, D);
    Cudd_RecursiveDeref(dd, temp2);
    Cudd_RecursiveDeref(dd, temp3);
    Cudd_RecursiveDeref(dd, D);
    D = temp;
    A = attr(D, universe);
    temp2 = bdd_not(dd, A);
    temp = bdd_and(dd, temp2, S);
    Cudd_RecursiveDeref(dd, temp2);
    Cudd_RecursiveDeref(dd, S);
    S = temp;

    if (firstmec) {
      bool classical;
      std::vector<int> *aux = newsize(S, 1, classical);
      if (!classical) {
        if ((aux->size()==0) || (aux->size()==1 && selfloops.count((*aux)[0]) == 0)) {
          Cudd_RecursiveDeref(dd, S);
          Cudd_RecursiveDeref(dd, T);
          delete aux;
          continue;
        }
      }
      delete aux;
    }

    temp2 = pre(A);
    temp3 = bdd_or(dd, temp2, T);
    Cudd_RecursiveDeref(dd, temp2);
    temp = bdd_and(dd, temp3, S);
    Cudd_RecursiveDeref(dd, temp3);
    Cudd_RecursiveDeref(dd, T);
    T = temp;

    Cudd_RecursiveDeref(dd, D);
    Cudd_RecursiveDeref(dd, A);

    if (T == Cudd_ReadLogicZero(dd)) {
      Cudd_RecursiveDeref(dd, T);
      allMECs.push_back(S);
      continue;
    }

    temp2 = post(S);
    temp3 = bdd_and(dd, temp2, S);
    bool hasedge = (temp3 != Cudd_ReadLogicZero(dd));
    Cudd_RecursiveDeref(dd, temp2);
    Cudd_RecursiveDeref(dd, temp3);

    if (hasedge) {
      bool classical;
      std::vector<int>* pos = new std::vector<int>(), *pos2;
      pos2 = newsize(T, threshold, classical);
      if (classical) // more than threshold, do classical
      {
        countClassic++;
        reset();
        ImprovedSCCFind(S,bdd_zero(dd),bdd_zero(dd));
        for (int i = 0; i < sccs.size(); i++)
        {
          temp = sccs[i];
          Cudd_Ref(temp);
          Q.push(temp);
          temp2 = bdd_zero(dd);
          Q.push(temp2);
        }
        Cudd_RecursiveDeref(dd, T);
      }
      else // at most threshold, do LSS
      {
        countLSS++;
        bool bottomscc; // ignore in mec subroutine
        H = bdd_zero(dd);
        C = lockStepSearch(S, H, T, pos, pos2, bottomscc);
        delete pos; delete pos2;
        Cudd_RecursiveDeref(dd, H);

        bool trivialbscc = false;
        bool classical;
        std::vector<int> *aux = newsize(S, 1, classical);
        if (!classical && aux->size()==1) {
          if (selfloops.count((*aux)[0]) == 0) {
            trivialbscc = true;
          }
        }
        delete aux;

        if (!trivialbscc) {
          allMECs.push_back(C);
        }

        /*
        temp = post(C);
        temp2 = bdd_and(dd, C, temp);
        Cudd_RecursiveDeref(dd, temp);
        if (temp2 != Cudd_ReadLogicZero(dd)) {
          allMECs.push_back(C);
        }
        Cudd_RecursiveDeref(dd, temp2);
        */

        temp2 = bdd_not(dd, C);
        temp = bdd_and(dd, temp2, S);
        Cudd_RecursiveDeref(dd, temp2);
        Cudd_RecursiveDeref(dd, S);
        S = temp;

        if (firstmec) {
          bool classical;
          std::vector<int> *aux = newsize(S, 1, classical);
          if (!classical) {
            if ((aux->size()==0) || (aux->size()==1 && selfloops.count((*aux)[0]) == 0)) {
              Cudd_RecursiveDeref(dd, S);
              Cudd_RecursiveDeref(dd, T);
              delete aux;
              continue;
            }
          }
          delete aux;
        }

        temp2 = pre(C);
        temp3 = bdd_or(dd, temp2, T);
        temp = bdd_and(dd, temp3, S);
        Cudd_RecursiveDeref(dd, temp2);
        Cudd_RecursiveDeref(dd, temp3);
        Cudd_RecursiveDeref(dd, T);
        T = temp;
        Q.push(S);
        Q.push(T);
      }
    }
  }

  return;
}

void mec::streettMDPbasic()
{
  presteps = 0;
  poststeps = 0;
  countClassic = 0;
  countLSS = 0;

  firstmec = true;
  MECbasic();
  firstmec = false;

  FMECtime = (double)(clock() - fmecafterscctime)/CLOCKS_PER_SEC;
  FMECpresteps = presteps - fmecaftersccpre;
  FMECpoststeps = poststeps - fmecaftersccpost;
  FMECcountClassic = countClassic;
  FMECcountLSS = countLSS;

  presteps = 0;
  poststeps = 0;
  countClassic = 0;
  countLSS = 0;

  //cout << "MECDONE" << endl << flush;
  //debugtime = clock();

  clock_t timeinit = clock();

  std::queue <bdd_ptr> Q;
  bdd_ptr temp, temp2, temp3, X, B;
  for (int i = 0; i < allMECs.size(); i++)
  {
    temp = allMECs[i];
    Cudd_Ref(temp);
    Q.push(temp);
  }
  while (!Q.empty())
  {
    //cout << "Sstart " << flush;
    //debugtime = clock();
    X = Q.front();  Q.pop();
    //Cudd_Ref(X); // ViTo 29.9. commented 2.10.
    // Here the check for size 1 no selfloop is useless,
    // since candidates are computed by ALLMEC and those
    // never satisfy the above condition
    B = bdd_zero(dd);
    for (int i = 0; i < k; i++)
    {
      temp2 = bdd_and(dd, X, Ustreett_pair[i]);
      if (temp2 == Cudd_ReadLogicZero(dd))
      {
        temp3 = bdd_and(dd, X, Lstreett_pair[i]);
        temp = bdd_or(dd, temp3, B);
        Cudd_RecursiveDeref(dd, B);
        B = temp;
        Cudd_RecursiveDeref(dd, temp3);
      }
      Cudd_RecursiveDeref(dd, temp2);
    }
    //cout << "Bad(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << flush;
    //debugtime = clock();
    if (B != Cudd_ReadLogicZero(dd))
    {
      countClassic++;
      temp2 = attr(B, X);
      temp3 = bdd_not(dd, temp2);
      temp = bdd_and(dd, temp3, X);
      Cudd_RecursiveDeref(dd, temp2);
      Cudd_RecursiveDeref(dd, temp3);
      Cudd_RecursiveDeref(dd, X);
      X = temp;
      reset(1);
      //cout << "nonzeroWILLDOMEC(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << flush;
      //debugtime = clock();
      MECbasic(X);
      for (int i = 0; i < allMECs.size(); i++) {
        temp = allMECs.at(i);
        Cudd_Ref(temp); // ViTo 2.10.
        Q.push(temp);
      }
      //cout << "INSIDEMECDONE(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << endl << flush;
      //debugtime = clock();
    }
    else
    {
      temp = bdd_or(dd, goodEC, X);
      Cudd_RecursiveDeref(dd, goodEC);
      Cudd_RecursiveDeref(dd, X); // ViTo 2.10.
      goodEC = temp;
      //cout << "zero(" << ((double)(clock() - debugtime)/CLOCKS_PER_SEC) << ") " << endl << flush;
      //debugtime = clock();
    }
    Cudd_RecursiveDeref(dd, B);
  }
  //answer = reach(goodEC, universe);

  STRtime = (double)(clock() - timeinit)/CLOCKS_PER_SEC;
  STRpresteps = presteps;
  STRpoststeps = poststeps;
  STRcountClassic = countClassic;
  STRcountLSS = countLSS;

  return;
}

void mec::ImprovedStreettMDP()
{
  presteps = 0;
  poststeps = 0;
  countClassic = 0;
  countLSS = 0;

  firstmec = true;
  ImprovedMEC();
  firstmec = false;

  FMECtime = (double)(clock() - fmecafterscctime)/CLOCKS_PER_SEC;
  FMECpresteps = presteps - fmecaftersccpre;
  FMECpoststeps = poststeps - fmecaftersccpost;
  FMECcountClassic = countClassic;
  FMECcountLSS = countLSS;

  presteps = 0;
  poststeps = 0;
  countClassic = 0;
  countLSS = 0;

  clock_t timeinit = clock();

  std::queue<bdd_ptr> Q;
  bdd_ptr temp, temp2, temp3, B, S, T, H, A, D, As, Ac, C;
  for (int i = 0; i < allMECs.size(); i++)
  {
    temp = allMECs[i];
    Cudd_Ref(temp);
    Q.push(temp);
    temp2 = bdd_zero(dd);
    Q.push(temp2);
    temp3 = bdd_zero(dd);
    Q.push(temp3);
  }
  while (!Q.empty())
  {
    S = Q.front();  Q.pop();
    H = Q.front();  Q.pop();
    T = Q.front();  Q.pop();
    bool classical;
    std::vector<int> *aux = newsize(S, 1, classical);
    if (!classical && aux->size()==1) {
      if (selfloops.count((*aux)[0]) == 0) {
        Cudd_RecursiveDeref(dd, S);
        Cudd_RecursiveDeref(dd, H);
        Cudd_RecursiveDeref(dd, T);
        delete aux;
        continue;
      }
    }
    delete aux;
    B = bdd_zero(dd);
    for (int i = 0; i < k; i++)
    {
      temp2 = bdd_and(dd, S, Ustreett_pair[i]);
      if (temp2 == Cudd_ReadLogicZero(dd))
      {
        temp3 = bdd_and(dd, S, Lstreett_pair[i]);
        temp = bdd_or(dd, temp3, B);
        Cudd_RecursiveDeref(dd, B);
        Cudd_RecursiveDeref(dd, temp3);
        B = temp;
      }
      Cudd_RecursiveDeref(dd, temp2);
    }
    bool done = false;
    while (B != Cudd_ReadLogicZero(dd))
    {
      A = attr(B, S);
      temp2 = bdd_not(dd, A);
      temp = bdd_and(dd, temp2, S);
      Cudd_RecursiveDeref(dd, S);
      Cudd_RecursiveDeref(dd, temp2);
      S = temp;

      bool classical;
      std::vector<int> *aux = newsize(S, 1, classical);
      if (!classical) {
        if (aux->size() == 0 || (aux->size() == 1 && selfloops.count((*aux)[0]) == 0)) {
          Cudd_RecursiveDeref(dd, S);
          Cudd_RecursiveDeref(dd, H);
          Cudd_RecursiveDeref(dd, T);
          Cudd_RecursiveDeref(dd, A);
          delete aux;
          done = true;
          break;
        }
      }
      delete aux;

      temp2 = post(A);
      temp3 = bdd_or(dd, temp2, H);
      temp = bdd_and(dd, temp3, S);
      Cudd_RecursiveDeref(dd, H);
      Cudd_RecursiveDeref(dd, temp2);
      Cudd_RecursiveDeref(dd, temp3);
      H = temp;

      temp2 = pre(A);
      temp3 = bdd_or(dd, temp2, T);
      temp = bdd_and(dd, temp3, S);
      Cudd_RecursiveDeref(dd, T); // was Cudd_RecursiveDeref(dd, H)
      Cudd_RecursiveDeref(dd, temp2);
      Cudd_RecursiveDeref(dd, temp3);
      T = temp; // was H = temp

      Cudd_RecursiveDeref(dd, A);
      Cudd_RecursiveDeref(dd, B);
      B = bdd_zero(dd);
      for (int i = 0; i < k; i++)
      {
        temp2 = bdd_and(dd, S, Ustreett_pair[i]);
        if (temp2 == Cudd_ReadLogicZero(dd))
        {
          temp3 = bdd_and(dd, S, Lstreett_pair[i]);
          temp = bdd_or(dd, temp3, B);
          Cudd_RecursiveDeref(dd, B);
          Cudd_RecursiveDeref(dd, temp3);
          B = temp;
        }
        Cudd_RecursiveDeref(dd, temp2);
      }
    }
    Cudd_RecursiveDeref(dd, B);
    if (done) continue;

    temp2 = post(S);
    temp3 = bdd_and(dd, temp2, S);
    bool containsedge = temp3 != Cudd_ReadLogicZero(dd);
    Cudd_RecursiveDeref(dd, temp2);
    Cudd_RecursiveDeref(dd, temp3);
    if (containsedge)
    {
      bool classical, classical2;
      std::vector<int>* pos, *pos2;
      pos = newsize(H, threshold, classical);
      pos2 = newsize(T, threshold, classical2);
      if (!classical && !classical2 && pos->size() + pos2->size() == 0)
      {
        temp = bdd_or(dd, goodEC, S);
        Cudd_RecursiveDeref(dd, goodEC);
        Cudd_RecursiveDeref(dd, H);
        Cudd_RecursiveDeref(dd, T);
        goodEC = temp;
      }
      else if (classical || classical2 || ((pos->size() + pos2->size()) > threshold))
      {  // more than threshold, do classical
        countClassic++;
        reset();
        temp = S;
        Cudd_Ref(temp);
        ImprovedSCCFind(temp,bdd_zero(dd),bdd_zero(dd));
        if (sccs.size() == 1)
        {
          temp = bdd_or(dd, goodEC, S);
          Cudd_RecursiveDeref(dd, goodEC);
          Cudd_RecursiveDeref(dd, H);
          Cudd_RecursiveDeref(dd, T);
          goodEC = temp;
        }
        else
        {
          for (int i = 0; i < sccs.size(); i++)
          {
            C = sccs[i];
            Cudd_Ref(C);

            bool classical;
            std::vector<int> *aux = newsize(C, 1, classical);
            if (!classical && aux->size()==1) {
              if (selfloops.count((*aux)[0]) == 0) {
                Cudd_RecursiveDeref(dd, C);
                delete aux;
                continue;
              }
            }
            delete aux;

            temp = bdd_not(dd, C);
            temp2 = bdd_and(dd, temp, S);
            temp3 = pre(temp2);
            Cudd_RecursiveDeref(dd, temp);
            Cudd_RecursiveDeref(dd, temp2);
            temp2 = bdd_and(dd, temp3, VR);
            temp = bdd_and(dd, temp2, C);
            //Cudd_RecursiveDeref(dd, D);  // why here? Caused segfault
            Cudd_RecursiveDeref(dd, temp2);
            Cudd_RecursiveDeref(dd, temp3);
            D = temp;
            A = attr(D, C);

            temp2 = bdd_not(dd, A);
            temp = bdd_and(dd, temp2, C);
            Cudd_RecursiveDeref(dd, C);
            Cudd_RecursiveDeref(dd, temp2);
            C = temp;

            aux = newsize(C, 1, classical);
            if (!classical) {
              if (aux->size() == 0 || (aux->size() == 1 && selfloops.count((*aux)[0]) == 0)) {
                Cudd_RecursiveDeref(dd, C);
                Cudd_RecursiveDeref(dd, A);
                Cudd_RecursiveDeref(dd, D);
                delete aux;
                continue;
              }
            }
            delete aux;

            temp2 = post(A);
            temp = bdd_and(dd, temp2, C);
            Cudd_RecursiveDeref(dd, H);
            Cudd_RecursiveDeref(dd, temp2);
            H = temp;
            Cudd_Ref(H); // ViTo 4.10.
            temp2 = pre(A);
            temp = bdd_and(dd, temp2, C);
            Cudd_RecursiveDeref(dd, T);
            Cudd_RecursiveDeref(dd, temp2);
            T = temp;
            Cudd_Ref(T); // ViTo 4.10.
            Q.push(C);
            Q.push(H);
            Q.push(T);
            Cudd_RecursiveDeref(dd, A);
            Cudd_RecursiveDeref(dd, D);
          }
        }
      }
      else
      {  // at most threshold, do LSS
        countLSS++;
        bool bottomscc;
        C = lockStepSearch(S, H, T, pos, pos2, bottomscc);
        delete pos; delete pos2;
        temp = bdd_not(dd, C);
        temp3 = bdd_and(dd, temp, S); // S \ C
        Cudd_RecursiveDeref(dd, temp);
        temp = bdd_not(dd, S);
        temp2 = bdd_and(dd, temp, C); // C \ S
        if (temp2 == Cudd_ReadLogicZero(dd) && temp3 == Cudd_ReadLogicZero(dd))
        {
          temp = bdd_or(dd, goodEC, S);
          Cudd_RecursiveDeref(dd, goodEC);
          goodEC = temp;
          Cudd_RecursiveDeref(dd, C);
          Cudd_RecursiveDeref(dd, temp2);
          Cudd_RecursiveDeref(dd, temp3);
        }
        else
        {
          // ViTo 4.10. added optimizations
          Cudd_RecursiveDeref(dd, temp2);
          Cudd_RecursiveDeref(dd, temp3);

          bool disregardC = false;
          bool disregardS = false;
          bool derefAc = false;
          bool derefAs = false;

          std::vector<int> *aux = newsize(C, 1, classical);
          if (!classical)
            if (aux->size()==1 && selfloops.count((*aux)[0]) == 0)
              disregardC = true;
          delete aux;

          if (!bottomscc && !disregardC) {
            // note: we can skip this is disregardC
            // since Ac is only used to modify C anyway
            // lines 33 and 34
            temp = bdd_not(dd, C);
            temp2 = bdd_and(dd, temp, S);
            temp3 = pre(temp2);
            Cudd_RecursiveDeref(dd, temp);
            Cudd_RecursiveDeref(dd, temp2);
            temp2 = bdd_and(dd, temp3, VR);
                                          Cudd_RecursiveDeref(dd, temp3);
            temp = bdd_and(dd, temp2, C);
            Cudd_RecursiveDeref(dd, temp2);
            D = temp;
            Ac = attr(D, C); // Ac works, others deref.
            Cudd_RecursiveDeref(dd, D);
            derefAc = true;
          }

          if (bottomscc) {
            // line 35
            As = attr(C, S);
            derefAs = true;
          } else As = bdd_or(dd, C, bdd_zero(dd)); // copy of C

          if (!bottomscc && !disregardC) {
            // line 36
            temp2 = bdd_not(dd, Ac);
            temp = bdd_and(dd, temp2, C);
            Cudd_RecursiveDeref(dd, C);
            Cudd_RecursiveDeref(dd, temp2);
            C = temp;

            // update disregardC
            aux = newsize(C, 1, classical);
            if (!classical)
              if (aux->size()==0 || (aux->size()==1 && selfloops.count((*aux)[0]) == 0))
                disregardC = true;
            delete aux;
          }

          //always do this
          // line 37
          temp2 = bdd_not(dd, As);
          temp = bdd_and(dd, temp2, S);
          Cudd_RecursiveDeref(dd, S);
          Cudd_RecursiveDeref(dd, temp2);
          S = temp;

          // update disregardS
          aux = newsize(S, 1, classical);
          if (!classical)
            if (aux->size()==0 || (aux->size()==1 && selfloops.count((*aux)[0]) == 0))
              disregardS = true;
          delete aux;

          if (!disregardC) {
            // push C
            Q.push(C);

            // line 38
            if (derefAc)
              temp2 = post(Ac);
            else
              temp2 = bdd_zero(dd);
            temp3 = bdd_or(dd, temp2, H);
            Cudd_RecursiveDeref(dd, temp2);
            temp = bdd_and(dd, temp3, C);
            Cudd_RecursiveDeref(dd, temp3);
            Q.push(temp);

            // line 39
            if (derefAc)
              temp2 = pre(Ac);
            else
              temp2 = bdd_zero(dd);
            temp3 = bdd_or(dd, temp2, T);
            Cudd_RecursiveDeref(dd, temp2);
            temp = bdd_and(dd, temp3, C);
            Cudd_RecursiveDeref(dd, temp3);
            Q.push(temp);
          }

          if (!disregardS) {
            // push S
            Q.push(S);

            // line 40
            temp2 = post(As);
            temp3 = bdd_or(dd, temp2, H);
            Cudd_RecursiveDeref(dd, temp2);
            temp = bdd_and(dd, temp3, S);
            Cudd_RecursiveDeref(dd, temp3);
            Q.push(temp);

            // line 41
            if (bottomscc)
              temp2 = pre(As);
            else
              temp2 = bdd_zero(dd); // if C is top, then As = C, and pre(C) = empty
            temp3 = bdd_or(dd, temp2, T);
            Cudd_RecursiveDeref(dd, temp2);
            temp = bdd_and(dd, temp3, S);
            Cudd_RecursiveDeref(dd, temp3);
            Q.push(temp);
          }

          if (derefAc) Cudd_RecursiveDeref(dd, Ac);
          Cudd_RecursiveDeref(dd, As);
        }
      }
    }
  }
  //answer = reach(goodEC, universe);

  STRtime = (double)(clock() - timeinit)/CLOCKS_PER_SEC;
  STRpresteps = presteps;
  STRpoststeps = poststeps;
  STRcountClassic = countClassic;
  STRcountLSS = countLSS;

  return;
}
