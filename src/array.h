
 /*
   Author: Rouben Rostamian, <rostamian@umbc.edu>
   June 2000
   Revised Jul 2000
   Revised Sep 2002
   Revised Jun 2003 - macros don't call exit anymore
                      changed macro names to all-caps
   Revised Aug 2003 - changed reserved names to all caps
   Revised Dec 2003 - changed array index types to size_t (were int)
                    - added an "out of memory" message, printed to stderr
                    - most FREE* macros now come before MAKE* macros,
		      for possible improved efficiency in preprocessing
*/






/* array.h

   o This file defines the following macros:

     MAKE_1ARRAY(a,n)       make a 1D array of length "n"
     MAKE_2ARRAY(a,m,n)     make a 2D array of dimensions "m x n"
     MAKE_3ARRAY(a,l,m,n)   make a 3D array of dimensions "l x m x n"
     MAKE_4ARRAY(a,k,l,m,n) make a 4D array of dimensions "k x l x m x n"

     FREE_1ARRAY(a)         free memory allocated by MAKE_1ARRAY()
     FREE_2ARRAY(a)         free memory allocated by MAKE_2ARRAY()
     FREE_3ARRAY(a)         free memory allocated by MAKE_3ARRAY()
     FREE_4ARRAY(a)         free memory allocated by MAKE_4ARRAY()

   o Additionally, it defines the following convenience macros as synonyms:

     MAKE_VECTOR(a,n)       same as MAKE_1ARRAY(a,n)
     FREE_VECTOR(a)         same as FREE_1ARRAY(a)
     MAKE_MATRIX(a,m,n)     same as MAKE_2ARRAY(a,m,n)
     FREE_MATRIX(a)         same as FREE_2ARRAY(a)

   o Additionally, it declares and uses the identifiers:

     ARRAY_H2RESERVED
     ARRAY_H3RESERVED
     ARRAY_H4RESERVED

     within local blocks within the macros.  THESE IDENTIFIERS
     ARE RESERVED.  The user should not use identifiers with this
     names within the scope of these macros.

   o vector/matrix/array elements can be of _any_ type.

   o If malloc() fails during the execution of any of the MAKE_* macros:
         . An "out of memory" message is printed to stderr.
         . The macro's first argument is set to NULL.

   o After a call to any of the FREE_*() macros, the macro's argument
     is set to NULL.

   o The FREE_*() macros can be applied to previously FREE_*()-ed
     arrays with no ill effect.

   o Note that macro arguments are evaluated more than once.

   o Customization

      When malloc() returns NULL, MAKE_1ARRAY prints a message on stderr
      and the program continues.  The user can alter this behavior by
      redefining the MAKE_1ARRAY macro.  For instance, to call exit()
      whenever malloc() fails, the user can do:

      #undef MAKE_1ARRAY
      #define MAKE_1ARRAY(a,n) do {                                          \
          (a) = malloc((n) * sizeof *(a));                                   \
          if ((a)==NULL) {                                                   \
              fprintf(stderr, "*** in file %s, function %s(), line %d: "     \
                      "out of memory!\n",  __FILE__, __func__, __LINE__);    \
              exit(EXIT_FAILURE);                                            \
          }                                                                  \
      } while (0)

      Since only MAKE_1ARRAY calls malloc() explicitly, this change affects
      the behavior of not only MAKE_1ARRAY but all other MAKE_* macros as well.


---- SAMPLE USAGE -------------
#include "array.h"
int main(void)
{
    float ***a;            // can use any other type instead of "float"
    size_t p=3, q=4, r=5;  // will make a 3x4x5 3-D array of float

    MAKE_3ARRAY(a, p, q, r);
    if (a==NULL)
        return EXIT_FAILURE;

    a[2][0][1] = 3.14;
    printf("%g \n", a[2][0][1]);
    FREE_3ARRAY(a);
    return EXIT_SUCCESS;
}
---- END OF SAMPLE USAGE -------


   Author: Rouben Rostamian, <rostamian@umbc.edu>
   June 2000
   Revised Jul 2000
   Revised Sep 2002
   Revised Jun 2003 - macros don't call exit anymore
                      changed macro names to all-caps
   Revised Aug 2003 - changed reserved names to all caps
   Revised Dec 2003 - changed array index types to size_t (were int)
                    - added an "out of memory" message, printed to stderr
                    - most FREE* macros now come before MAKE* macros,
		      for possible improved efficiency in preprocessing
*/


#ifndef H_ARRAY_
#define H_ARRAY_

#include <stdio.h>
#include <stdlib.h>

/* ---------- 1D arrays ---------------------- */

/* WCC: R don't like fprintf(stderr, ...) */
#ifndef __HAVE_R_
#define MAKE_1ARRAY(a,n) do {                                                \
    (a) = malloc((n) * sizeof *(a));                                         \
    if ((a)==NULL)                                                           \
        REprintf("*** in file %s, function %s(), line %d: "           \
                "out of memory!\n",  __FILE__, __func__, __LINE__);          \
} while (0)                                                                  
#else
#include <R.h>
#define MAKE_1ARRAY(a,n) do {                                                \
    (a) = malloc((n) * sizeof *(a));                                         \
    if ((a)==NULL)                                                           \
        REprintf("*** in file %s, function %s(), line %d: "           \
                "out of memory!\n",  __FILE__, __func__, __LINE__);          \
} while (0)                                                                  
#endif

#define FREE_1ARRAY(a)  do {                                                 \
    free(a);                                                                 \
    a = NULL;                                                                \
} while (0)

/* ---------- 2D arrays ---------------------- */

#define FREE_2ARRAY(a) do {                                                  \
    size_t ARRAY_H2RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H2RESERVED=0; (a)[ARRAY_H2RESERVED]!=NULL; ARRAY_H2RESERVED++)\
        FREE_1ARRAY((a)[ARRAY_H2RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)

/* note: parenthesize first arg because it may be given as `*a' */
#define MAKE_2ARRAY(a,m,n) do {                                              \
    size_t ARRAY_H2RESERVED;                                                 \
    MAKE_1ARRAY(a,(m)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[m] = NULL;                                                           \
    for (ARRAY_H2RESERVED=0; ARRAY_H2RESERVED<(m); ARRAY_H2RESERVED++) {     \
        MAKE_1ARRAY((a)[ARRAY_H2RESERVED],(n));                              \
        if ((a)[ARRAY_H2RESERVED]==NULL) {                                   \
            FREE_2ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

/* ---------- 3D arrays ---------------------- */

#define FREE_3ARRAY(a) do {                                                  \
    size_t ARRAY_H3RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H3RESERVED=0; (a)[ARRAY_H3RESERVED]!=NULL; ARRAY_H3RESERVED++)\
        FREE_2ARRAY((a)[ARRAY_H3RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)

#define MAKE_3ARRAY(a,p,q,r) do {                                            \
    size_t ARRAY_H3RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H3RESERVED=0; ARRAY_H3RESERVED<(p); ARRAY_H3RESERVED++) {     \
        MAKE_2ARRAY((a)[ARRAY_H3RESERVED],(q),(r));                          \
        if ((a)[ARRAY_H3RESERVED]==NULL) {                                   \
            FREE_3ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

/* ---------- 3D arrays ---------------------- */

#define FREE_4ARRAY(a) do {                                                  \
    size_t ARRAY_H4RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H4RESERVED=0; (a)[ARRAY_H4RESERVED]!=NULL; ARRAY_H4RESERVED++)\
        FREE_3ARRAY((a)[ARRAY_H4RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)

#define MAKE_4ARRAY(a,p,q,r,s) do {                                          \
    size_t ARRAY_H4RESERVED;                                                 \
    MAKE_1ARRAY(a,(p)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[p] = NULL;                                                           \
    for (ARRAY_H4RESERVED=0; ARRAY_H4RESERVED<(p); ARRAY_H4RESERVED++) {     \
        MAKE_3ARRAY((a)[ARRAY_H4RESERVED],(q),(r),(s));                      \
        if ((a)[ARRAY_H4RESERVED]==NULL) {                                   \
            FREE_4ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

/* ---------- synonyms ---------------------- */

#define MAKE_VECTOR(a,n)    MAKE_1ARRAY(a,n)
#define MAKE_MATRIX(a,m,n)  MAKE_2ARRAY(a,m,n)

#define FREE_VECTOR(a)      FREE_1ARRAY(a)
#define FREE_MATRIX(a)      FREE_2ARRAY(a)

#endif /* H_ARRAY_ */
