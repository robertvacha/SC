/** @file macros.h*/

#ifndef MACROS2_H
#define MACROS2_H

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

//#define EXTRA_HYDROPHOBIC_ALL_BODY_ATTRACTION // eCpscCpsc, all CPSC
//#define E_ISO 4

//
//  for assert() -> if NDEBUG defined assert not compiled
//  #ifndef NDEBUG someting() #endif -> if NDEBUG defined -> not compiled
//
#define NDEBUG

// DEFAULT CUBOID
//#define WEDGE

//#define TESTING
#ifdef TESTING
    #define RAN2
#endif

// Default in e2ScaOr2Spa() - WEEK-CHANDLER-ANDERSEN
//#define LJ

#ifdef ENABLE_MPI
# include <mpi.h>
#endif

//#define OMP1 // energy functions openMP
#ifdef OMP1
#include "omp.h"
#endif



#ifdef DEBUGGING_INIT
#define DEBUG_INIT(...) fprintf(stderr, "DB in INIT: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); fflush(stderr);
#else
#define DEBUG_INIT(...)
#endif

#ifdef DEBUGGING_SIM
#define DEBUG_SIM(...) fprintf(stderr, "DB in SIM: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); fflush(stderr);
#else
#define DEBUG_SIM(...)
#endif

#ifdef DEBUGGING
#define DEBUG(...) fprintf(stderr, "DB: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); fflush(stderr);
#else
#define DEBUG(...)
#endif
/* End of DEBUG macros */

#define MAXN 14000           /* Maximum number of particles */
#define MAXCHL 20          /* Maximum length of chain */
#define ZEROTOL 1.0e-12     /* Dot products below ZEROTOL are deemed zero */
#define ZEROTOL2 1.0e-8     /* numbers below ZEROTOL are deemed zero */
#define PI 3.141592653589793238462643383279  /* pi */
#define PIH 1.57079632679489661923132169163975  /* pi half*/

/*Particle types*/
#define SC 10             /*spherocylinder*/
#define SCN SC+0          /*spherocylinder non-attractive*/
#define SCA SC+1          /*spherocylinder isotropicaly attractive*/
#define PSC SC+2          /*spherocylinder with patchy attraction*/
#define CPSC SC+3         /*spherocylinder with cylindrical patchy attraction*/
#define CHPSC SC+4        /* chiral psc */
#define CHCPSC SC+5       /* chiral cpsc */
#define TPSC SC+6          /*spherocylinder with two patches*/
#define TCPSC SC+7         /*spherocylinder with two cylindrical patches*/
#define TCHPSC SC+8        /* chiral 2psc */
#define TCHCPSC SC+9       /* chiral 2cpsc */

#define SP 30             /*sphere - should be over all apherocylinders*/ // changed problem?
#define SPN SP+0          /* sphere non-attractive*/
#define SPA SP+1          /* spherocylinder isotropicaly attractive*/

#define MAXT 40           /* Maximum number of types we have*/
#define MAXMT 100         /* Maximum number of molecular types */

/*Reading topology*/
#define SMSTR 64           /* Small string length */
#define STRLEN 400         /* maximum length of line*/
#define CONTINUE    '\\'   /* symbol for line continue*/
#define COMMENTSIGN '#'    /* symbol for comment*/
#define OPENKEY  '['       /* starting sign for keyword*/
#define CLOSEKEY ']'       /* ending sign for keyword*/
#define SEPARATOR ':'      /* sign for separator*/
#define OPENMOL '{'        /* starting sign for molecules*/
#define CLOSEMOL '}'       /* ending sign for molecules*/
#define BOXSEP 'x'       /* extraction of box*/

/* Wang Landau method */
#define WL_GERR 0.0001           /* Max roughnes in histogram */
#define WL_ALPHATOL 0.000001     /* Covergence crietria for detailed balance */
#define WL_MINHIST 1000         /* Minimum histogram sampling for considering roughness */
#define WL_ZERO 0.000000000000  /* Zero for histogram with current weights*/
#define WL_CONTACTS 36.0          /* Square distance under which are particles in contact */

/* Math */
#define DOT(a,b) ((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)   /* Dot product */
#define AVER(a,b) ((a+b)*0.5)                                      /* Arithmetic average*/
#define ROUND(a) (a > 0.0) ? floor(a + 0.5) : ceil(a - 0.5);       /* Round double*/
#define PMONE(a) (1 - 2 * a)                                       /* Takes 1 or 0, return +-1 */

/* Acceptance ratio */
#define RATIO(a) ( ((a).acc+(a).rej) > 0 ? 1.0*(a).acc/((a).acc+(a).rej) : 0.0 )

#define INBOX(a,b) ( a > 0 ? modf(a,&b) : modf(a,&b)+1 )

#endif // MACROS2_H
