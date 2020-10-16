#include <vector>
#include <stack>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cassert>
#include <queue>
#include <set>

#include "streettgraph.h"

using namespace std;

/* The non-symbolic way of computing size:
   This is not used in the code
   It is just to provide a better understanding of what the symbolic procedure below does */
std::vector<int>* streettgraph::newsize1(bdd_ptr Src, int threshold, bool& classical) {
  std::vector<int>* ans = new std::vector<int>();
  bdd_ptr temp;
  classical = false;
  for(int i=0;i<state.size(); i++) {
    temp = bdd_and(dd,Src,state[i]);
    if(temp != Cudd_ReadLogicZero(dd)) {
      ans->push_back(i);
      if(1ll * ans->size() > threshold) {
        classical = true;
        break;
      }
    }
    Cudd_RecursiveDeref(dd,temp);
  }
  return ans;
}

/* Symbolic Size Computation, checks whether the size of Src is more than threshold
   if it is at most the threshold, returns the states enumerated and sets classical to false, else sets classical to true */
std::vector<int>* streettgraph::newsize(bdd_ptr Src, int threshold, bool& classical) {
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
    if(global_cnt > threshold) { // 1ll * global_cnt
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

/* Pick subroutine - This is inexpensive and hence not a symbolic step */
bdd_ptr streettgraph::pick(bdd_ptr B) {
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
      ans = temp;                       // ans = ans & ~var(i)
    }
    else if(cube[i] == 1) {
      temp = bdd_new_var_with_index(dd,i);
      Cudd_Ref(temp);
      temp2 = bdd_and(dd,ans,temp);
      Cudd_RecursiveDeref(dd,ans);
      Cudd_RecursiveDeref(dd,temp);
      ans = temp2;                      // ans = ans & var(i)
    }
    else
      assert(0 && "Cube has integer other than 0,1,2");
  }
  return ans;
}

/* Symbolic Pre Operation */
bdd_ptr streettgraph::pre(bdd_ptr set_states) {
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

  temp = compute_image(result,trans,next_cube);
  Cudd_RecursiveDeref(dd,result);
  result = temp;
  return result;
}

/* Symbolic Post Operation */
bdd_ptr streettgraph::post(bdd_ptr set_states) {
  poststeps++;

  bdd_ptr result = set_states;
  Cudd_Ref(result);
  bdd_ptr temp;
  temp = compute_image(result,trans,cube);
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

/* Image Computation */
bdd_ptr streettgraph::compute_image(bdd_ptr set_states, bdd_ptr trans, bdd_ptr cube) {
  bdd_ptr result,temp;
  result = bdd_and(dd,set_states,trans);
  temp = bdd_forsome(dd,result,cube);
  Cudd_RecursiveDeref(dd,result);
  result = temp;
  return result;
}

/* Initialization */
void streettgraph::init(int _p1, int _p2, int _k) {
  presteps = poststeps = 0;
  case1 = case2 = case3 = 0;
  mainLoopCounter = 0;
  p1 = _p1;
  p2 = _p2;
  k = _k;
        selfloops.clear();

  bdd_ptr temp,temp2,temp3;
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

  trans = bdd_zero(dd);
  for (int i = 0; i < k; i++)
  {
    temp = bdd_zero(dd);
    Lstreett_pair.push_back(temp);
    temp = bdd_zero(dd);
    Ustreett_pair.push_back(temp);
  }
  good_components = bdd_zero(dd);
  return;
}

/* Resetting SCC Information to run different algorithm */
void streettgraph::reset_sccs(bool comp) {
  if (comp)
        {
    //Cudd_RecursiveDeref(dd, answer);
    Cudd_RecursiveDeref(dd, good_components);
    good_components = bdd_zero(dd);
    presteps = poststeps = 0;
  }
  for(int i=0;i<sccs.size();i++)
    Cudd_RecursiveDeref(dd,sccs[i]);
  sccs.clear();
  return;
}

/* Adding Transitions */
void streettgraph::add_transition(int i, int j) {
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
void streettgraph::ImprovedSCCFind() {
  Cudd_Ref(universe);
  ImprovedSCCFind(universe,bdd_zero(dd),bdd_zero(dd));
  return;
}

/* ImprovedSCCFind Algorithm */
void streettgraph::ImprovedSCCFind(bdd_ptr S, bdd_ptr U, bdd_ptr s) {
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
void streettgraph::ImprovedSkelFwd(bdd_ptr S, bdd_ptr U, bdd_ptr s, bdd_ptr& FWSet, bdd_ptr& NewSet, bdd_ptr& NewState, bdd_ptr& P) {
  std::stack<bdd_ptr> stack;
  bdd_ptr L;
  bdd_ptr temp, temp2, temp3;
  L = s; Cudd_Ref(L);
  FWSet = bdd_zero(dd);     // FWSet = phi

  while(L != Cudd_ReadLogicZero(dd)) {
    stack.push(L);
    temp = bdd_or(dd,FWSet,L);
    Cudd_RecursiveDeref(dd,FWSet);
    FWSet = temp;     // FWSet = FWSet union L

    temp = post(L);
    temp2 = bdd_and(dd,temp,S); // post(L) // QUESTION: Is this required? Where else is this required?
    Cudd_RecursiveDeref(dd,temp);
    temp = bdd_not(dd,FWSet);
    L = bdd_and(dd,temp2,temp); // post(L) \ FWSet
    Cudd_RecursiveDeref(dd,temp2);
  }
  Cudd_RecursiveDeref(dd,L);

  P = bdd_and(dd,FWSet,U);
  L = stack.top();
  stack.pop();

  NewState = pick(L); NewSet = NewState; Cudd_Ref(NewSet);
  Cudd_RecursiveDeref(dd,L);

  while(!stack.empty()) {
    L = stack.top();
    stack.pop();
    temp = bdd_and(dd,L,P);
    if(temp != Cudd_ReadLogicZero(dd)) {  // There must be an intersection
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

/* Reverse Reachability */
bdd_ptr streettgraph::reach(bdd_ptr g, bdd_ptr V){
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
  } while(Rold != Rnew);

  Cudd_RecursiveDeref(dd,Rold);
  return Rnew;
}

std::vector <int> * streettgraph::read_answer(bdd_ptr &setToCount)
{
  std::vector <int> * ans = new std::vector<int> ();
  bdd_ptr temp;
  for(int i=0;i<state.size();i++)
  {
    temp = bdd_and(dd, state[i], setToCount);
    if (temp != Cudd_ReadLogicZero(dd))
      ans->push_back(i);
    Cudd_RecursiveDeref(dd, temp);
  }
  return ans;
}

void streettgraph::temp_print(std::vector<int>* v)
{
    if (true) { // change to false once output undesired
  //std::cout<<"here are the vertices of good component:\n";
  for (int i = 0; i < v->size(); i++)
    std::cout<< (*v)[i]<< " ";
  std::cout<< std::endl;
    }
}

void streettgraph::addLStates(int s, int i)
{
  bdd_ptr temp;
  temp = bdd_or(dd, state[s], Lstreett_pair[i]);
  Cudd_RecursiveDeref(dd, Lstreett_pair[i]);
  Lstreett_pair[i] = temp;
  return;
}

void streettgraph::addUStates(int s, int i)
{
  bdd_ptr temp;
  temp = bdd_or(dd, state[s], Ustreett_pair[i]);
  Cudd_RecursiveDeref(dd, Ustreett_pair[i]);
  Ustreett_pair[i] = temp;
  return;
}

bdd_ptr streettgraph::lockStepSearch(bdd_ptr S, bdd_ptr& H, bdd_ptr& T, std::vector <int>* vh, std::vector<int>* vt, bool& bottomscc)
{
    std::queue <bdd_ptr> Q, Q2; // Q for Heads, Q2 for Tails
    bdd_ptr temp, temp2, temp3, C;
    for (int i = 0; i < vh->size(); i++) {
      temp = state[vh->at(i)];
      Cudd_Ref(temp);
      Q.push(temp);
    }
    for (int i = 0; i < vt->size(); i++) {
      temp = state[vt->at(i)];
      Cudd_Ref(temp);
      Q2.push(temp);
    }

    while (true) {

        int i = 0;
        while (i < vh->size()) {
            C = Q.front();  Q.pop();
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

            if (classical) {
                temp3 = bdd_not(dd, state[vh->at(i)]);
                temp2 = bdd_and(dd, temp3, H);
                Cudd_RecursiveDeref(dd, temp3);
                Cudd_RecursiveDeref(dd, H);
                H = temp2;
                vh->erase(vh->begin() + i);
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

            if (classical) {
                temp3 = bdd_not(dd, state[vt->at(i)]);
                temp2 = bdd_and(dd, temp3, T);
                Cudd_RecursiveDeref(dd, temp3);
                Cudd_RecursiveDeref(dd, T);
                T = temp2;
                vt->erase(vt->begin() + i);
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

void streettgraph::BasicSTROBJ(bool badloop)
{
    presteps = poststeps = 0;
    case1 = case2 = case3 = 0;
    mainLoopCounter = 0;

    bdd_ptr temp, temp2, temp3, C, B;

    std::queue <bdd_ptr> Q;
    for (int i = 0; i < sccs.size(); i++) {
        temp = sccs[i];
        Cudd_Ref(temp);
        candidates.push(temp);
    }

    while (!candidates.empty()) {
      C = candidates.top(); candidates.pop();
      bool classical;
      std::vector<int>* aux = newsize(C, 1, classical);
      if (!classical && aux->size()==1) {
            if (selfloops.count((*aux)[0]) == 0) {
                Cudd_RecursiveDeref(dd, C);
                delete aux;
                continue;
            }
      }
      delete aux;

      bool removedsome = false;
      bool removedlastiter = false;
      bool done = false;

      do {
          removedlastiter = false;
          B = bdd_zero(dd);
          for (int i = 0; i < k; i++) {
              temp2 = bdd_and(dd, Ustreett_pair[i], C);
              if (temp2 == Cudd_ReadLogicZero(dd)) {
                  temp3 = bdd_and(dd, Lstreett_pair[i], C);
                  temp = bdd_or(dd, temp3, B);
                  Cudd_RecursiveDeref(dd, B);
                  Cudd_RecursiveDeref(dd, temp3);
                  B = temp;
              }
              Cudd_RecursiveDeref(dd, temp2);
          }

          if (B != Cudd_ReadLogicZero(dd)) {
              temp2 = bdd_not(dd, B);
              temp = bdd_and(dd, temp2, C);
              Cudd_RecursiveDeref(dd, temp2);
              Cudd_RecursiveDeref(dd, C);
              C = temp;
              removedsome = true;
              removedlastiter = true;

              aux = newsize(C, 1, classical);
              if (!classical) {
                  if (aux->size() == 0 || (aux->size() == 1 && selfloops.count((*aux)[0]) == 0)) {
                      Cudd_RecursiveDeref(dd, C);
                      done = true;
                  }
              }
              delete aux;
          }
          Cudd_RecursiveDeref(dd, B);

      } while (removedlastiter && badloop && !done);

      if (done)
          continue;

      if (removedsome) {
          reset_sccs();
          ImprovedSCCFind(C, bdd_zero(dd), bdd_zero(dd));
          for (int i = 0; i < sccs.size(); i++) {
              temp = sccs[i];
              Cudd_Ref(temp);
              candidates.push(temp);
          }
      } else {
          temp2 = post(C);
          temp3 = bdd_and(dd, temp2, C);
          if (temp3 != Cudd_ReadLogicZero(dd)) {
              temp = bdd_or(dd, C, good_components);
              Cudd_RecursiveDeref(dd, good_components);
              good_components = temp;
          }
          Cudd_RecursiveDeref(dd, temp2);
          Cudd_RecursiveDeref(dd, temp3);
          Cudd_RecursiveDeref(dd, C);
      }
      mainLoopCounter++;
    }

    //answer = reach(good_components, universe);
    return;
}

void streettgraph::ImprovedSTROBJ(int threshold, bool checkedge)
{
    presteps = poststeps = 0;
    case1 = case2 = case3 = 0;
    mainLoopCounter = 0;

    bdd_ptr temp, temp2, temp3, temp4, temp5;
    bdd_ptr B, S, H, T, C;

    std::queue <bdd_ptr> Q;
    for (int i = 0; i < sccs.size(); i++) {
      temp = sccs[i];
      Cudd_Ref(temp);
      Q.push(temp);
      temp2 = bdd_zero(dd);
      Q.push(temp2);
      temp3 = bdd_zero(dd);
      Q.push(temp3);
    }
    while (!Q.empty()) {
      S = Q.front();  Q.pop();
      H = Q.front();  Q.pop();
      T = Q.front();  Q.pop();

      bool classical2, classical;
      std::vector<int>* aux = newsize(S, 1, classical);
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
      for (int i = 0; i < k; i++) {
        temp2 = bdd_and(dd, Ustreett_pair[i], S);
        if (temp2 == Cudd_ReadLogicZero(dd)) {
          temp3 = bdd_and(dd, Lstreett_pair[i], S);
          temp = bdd_or(dd, B, temp3);
          Cudd_RecursiveDeref(dd, temp3);
          Cudd_RecursiveDeref(dd, B);
          B = temp;
        }
        Cudd_RecursiveDeref(dd, temp2);
      }
      bool done = false;
      while (B != Cudd_ReadLogicZero(dd)) {
            temp2 = bdd_not(dd, B);
            temp = bdd_and(dd, temp2, S);
            Cudd_RecursiveDeref(dd, temp2);
            Cudd_RecursiveDeref(dd, S);
            S = temp;

            aux = newsize(S, 1, classical);
            if (!classical) {
                if (aux->size() == 0 || (aux->size() == 1 && selfloops.count((*aux)[0]) == 0)) {
                    Cudd_RecursiveDeref(dd, S);
                    Cudd_RecursiveDeref(dd, H);
                    Cudd_RecursiveDeref(dd, T);
                    delete aux;
                    done = true;
                    break;
                }
            }
            delete aux;

            temp2 = post(B);
            temp3 = bdd_or(dd, temp2, H);
            Cudd_RecursiveDeref(dd, temp2);
            temp = bdd_and(dd, temp3, S);
            Cudd_RecursiveDeref(dd, temp3);
            Cudd_RecursiveDeref(dd, H);
            H = temp;

            temp2 = pre(B);
            temp3 = bdd_or(dd, temp2, T);
            Cudd_RecursiveDeref(dd, temp2);
            temp = bdd_and(dd, temp3, S);
            Cudd_RecursiveDeref(dd, temp3);
            Cudd_RecursiveDeref(dd, T);
            T = temp;

            Cudd_RecursiveDeref(dd, B);
            B = bdd_zero(dd);
            for (int i = 0; i < k; i++) {
              temp2 = bdd_and(dd, Ustreett_pair[i], S);
              if (temp2 == Cudd_ReadLogicZero(dd)) {
                temp3 = bdd_and(dd, Lstreett_pair[i], S);
                temp = bdd_or(dd, B, temp3);
                Cudd_RecursiveDeref(dd, temp3);
                Cudd_RecursiveDeref(dd, B);
                B = temp;
              }
              Cudd_RecursiveDeref(dd, temp2);
            }
      }
      Cudd_RecursiveDeref(dd, B);
      if (done) continue;

      // if |Hs| + |Ts| = 0 then S is good
      // In particular it has an edge (so no need
      // to check) since by Inv.1 strongly connected
      // and by above check not '1node-noselfloop'
      if (H == Cudd_ReadLogicZero(dd) && T == Cudd_ReadLogicZero(dd)) {
          case1++;
          temp = bdd_or(dd, S, good_components);
          Cudd_RecursiveDeref(dd, good_components);
          Cudd_RecursiveDeref(dd, S);
          Cudd_RecursiveDeref(dd, H);
          Cudd_RecursiveDeref(dd, T);
          good_components = temp;
      } else {

        bool containsedge = !checkedge;
        if (checkedge) {
            temp2 = post(S);
            temp3 = bdd_and(dd, temp2, S);
            Cudd_RecursiveDeref(dd, temp2);
            containsedge = temp3 != Cudd_ReadLogicZero(dd);
            Cudd_RecursiveDeref(dd, temp3);
        }

        if (containsedge) {
            std::vector<int>* pos, *pos2;
            pos = newsize(H, threshold, classical);
            pos2 = newsize(T, threshold, classical2);
          if (classical || classical2 || ((pos->size() + pos2->size()) > threshold)) {
              case2++;
              reset_sccs();
              Cudd_RecursiveDeref(dd, H);
              Cudd_RecursiveDeref(dd, T);
              ImprovedSCCFind(S,bdd_zero(dd),bdd_zero(dd));
              if (sccs.size() == 1) {
                temp = bdd_or(dd, S, good_components);
                Cudd_RecursiveDeref(dd, good_components);
                Cudd_RecursiveDeref(dd, S);
                good_components = temp;
              } else {
                for (int i = 0; i < sccs.size(); i++) {
                  temp = sccs[i];
                  Cudd_Ref(temp);
                  Q.push(temp);
                  temp = bdd_zero(dd);
                  Q.push(temp);
                  temp = bdd_zero(dd);
                  Q.push(temp);
                }
              }
          } else {
              case3++;
              bool bottomscc;
              C = lockStepSearch(S, H, T, pos, pos2, bottomscc);
              temp2 = bdd_not(dd, C);
              temp3 = bdd_and(dd, temp2, S);
              Cudd_RecursiveDeref(dd, temp2);
              bool onescc = (temp3 == Cudd_ReadLogicZero(dd));
              Cudd_RecursiveDeref(dd, temp3);
            if (onescc) {
                temp = bdd_or(dd, S, good_components);
                Cudd_RecursiveDeref(dd, good_components);
                Cudd_RecursiveDeref(dd, S);
                Cudd_RecursiveDeref(dd, H);
                Cudd_RecursiveDeref(dd, T);
                good_components = temp;
            } else {

              // may skip C-part

              bool skipC = false;

              aux = newsize(S, 1, classical);
              if (!classical && aux->size()==1) {
                  if (selfloops.count((*aux)[0]) == 0) {
                              skipC = true;
                  }
              }
              delete aux;

              if (!skipC) {
                  Q.push(C);

                  temp = bdd_zero(dd);
                  Q.push(temp); // Hc = \emptyset

                  temp = bdd_zero(dd);
                  Q.push(temp); // Tc = \emptyset

                  /*
                  if (bottomscc)
                      temp = post(S);
                  else
                      temp = bdd_zero(dd); // if C is top, then (post(S) and C) = empty
                  temp2 = bdd_or(dd, temp, H);
                  Cudd_RecursiveDeref(dd, temp);
                  temp = bdd_and(dd, temp2, C);
                  Cudd_RecursiveDeref(dd, temp2);
                  Q.push(temp);

                  if (!bottomscc)
                      temp = pre(S);
                  else
                      temp = bdd_zero(dd); // if C is bottom, then (pre(S) and C) = empty
                  temp2 = bdd_or(dd, temp, T);
                  Cudd_RecursiveDeref(dd, temp);
                  temp = bdd_and(dd, temp2, C);
                  Cudd_RecursiveDeref(dd, temp2);
                  Q.push(temp);
                  */
              }

              // never skip S-part

              temp2 = bdd_not(dd, C);
              temp = bdd_and(dd, temp2, S);
              Cudd_RecursiveDeref(dd, temp2);
              Cudd_RecursiveDeref(dd, S);
              S = temp;

              Q.push(S);

              if (!bottomscc)
                  temp = post(C);
              else
                  temp = bdd_zero(dd); // if C is bottom, then (post(C) and S) = empty
              temp2 = bdd_or(dd, temp, H);
              Cudd_RecursiveDeref(dd, temp);
              temp = bdd_and(dd, temp2, S);
              Cudd_RecursiveDeref(dd, temp2);
              Q.push(temp);

              if (bottomscc)
                  temp = pre(C);
              else
                  temp = bdd_zero(dd); // if C is top, then (pre(C) and S) = empty
              temp2 = bdd_or(dd, temp, T);
              Cudd_RecursiveDeref(dd, temp);
              temp = bdd_and(dd, temp2, S);
              Cudd_RecursiveDeref(dd, temp2);
              Q.push(temp);
            }
          }
        }
      }
      mainLoopCounter++;
    }

  //answer = reach(good_components, universe);
  return;
}
