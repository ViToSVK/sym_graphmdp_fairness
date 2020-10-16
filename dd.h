#ifndef _DD_H
#define _DD_H

#include <stdio.h>
#include "cudd.h"

typedef struct DdNode * add_ptr;
typedef struct DdNode * bdd_ptr;

void common_error(void * R, const char * s);

extern bdd_ptr  bdd_zero (DdManager *);
extern bdd_ptr  bdd_cube_union (DdManager *, bdd_ptr, bdd_ptr);
extern bdd_ptr  bdd_and (DdManager *, bdd_ptr, bdd_ptr);
extern bdd_ptr  bdd_or (DdManager *, bdd_ptr, bdd_ptr);
extern bdd_ptr  bdd_one (DdManager *);
extern bdd_ptr  bdd_forsome (DdManager *, bdd_ptr, bdd_ptr);
extern bdd_ptr  bdd_forall (DdManager *, bdd_ptr, bdd_ptr);
extern bdd_ptr  bdd_cube_intersection (DdManager *, bdd_ptr, bdd_ptr);
extern bdd_ptr  bdd_cube_diff (DdManager *, bdd_ptr, bdd_ptr);
extern bdd_ptr  bdd_dup (bdd_ptr);
extern bdd_ptr  bdd_support (DdManager *, bdd_ptr);
extern bdd_ptr  bdd_new_var_with_index (DdManager *, int);
extern bdd_ptr  bdd_vector_support (DdManager *, bdd_ptr*, int);
extern bdd_ptr  bdd_cofactor (DdManager *, bdd_ptr, bdd_ptr);
extern void     bdd_and_accumulate (DdManager *, bdd_ptr *, bdd_ptr);
extern int      bdd_get_lowest_index (DdManager *, bdd_ptr);
extern void     bdd_free (DdManager *, bdd_ptr);
extern bdd_ptr  bdd_not  (DdManager *, bdd_ptr);
extern int      bdd_is_one (DdManager *, add_ptr);
extern void     bdd_or_accumulate (DdManager *, bdd_ptr *, bdd_ptr);
extern int      bdd_size (DdManager *, bdd_ptr);
extern void     bdd_print_minterms(DdManager *dd, bdd_ptr f);
extern bdd_ptr	bdd_compose(DdManager *dd, bdd_ptr f, bdd_ptr var, int ind);
extern bdd_ptr	bdd_and_abstract(DdManager * dd, bdd_ptr f, bdd_ptr g, bdd_ptr cube);

extern bdd_ptr	bdd_existential_quantify (DdManager*,bdd_ptr,bdd_ptr);
extern bdd_ptr	bdd_universal_quantify (DdManager*,bdd_ptr,bdd_ptr);
extern bdd_ptr	bdd_then (bdd_ptr);
extern bdd_ptr	bdd_else (bdd_ptr);
extern int		bdd_value (bdd_ptr);
extern int		bdd_is_constant (bdd_ptr);

#endif /* _DD_H */
