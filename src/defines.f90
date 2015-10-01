!===================================================================================================================================
! Here, preprocessor variables for different equation systems and abbrevbiations for specific expressions are defined
!===================================================================================================================================

! Abbreviations
!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef SUN
#define __DATE__ '__TIME__ and __DATE__ not'
#define __TIME__ 'available for SUN COMPILER'
#define IEEE_ISNAN
#elif PGI
#define NO_ISNAN
#endif
#define __STAMP__ __FILE__,__LINE__,__DATE__,__TIME__

#define SDEALLOCATE(A) IF(ASSOCIATED(A)) DEALLOCATE(A)
#define ERRWRITE(a,b) WRITE(UNIT_errFile,b)
#define LOGWRITE(a,b) IF(Logging) WRITE(UNIT_logOut,b)

! General tolerance for mesh node comparisons and sfc (1.E.-13 should be ok too)
#define PP_MeshTolerance 1.E-12
! General tolerance for Real comparisons (double precision)
#define PP_RealTolerance 1.E-16

! Settings for different compilers
!-----------------------------------------------------------------------------------------------------------------------------------
! Colored output
#if defined(GNU)||defined(PGI)||defined(BLUEGENE)
! No colors available (still)
#define __FGCOLORRED__    ''
#define __FGCOLORGREEN__  ''
#define __FGCOLORYELLOW__ ''
#define __FGCOLORBLUE__   ''
#define __FGCOLOREND__    ''
#else
! Use default colors
#define __FGCOLORRED__    '\033[31m'
#define __FGCOLORGREEN__  '\033[32m'
#define __FGCOLORYELLOW__ '\033[33m'
#define __FGCOLORBLUE__   '\033[34m'
#define __FGCOLOREND__    '\033[0m'
#endif

#if(PP_CGNS_INT==64)
#  define PP_CGNS_INT_TYPE INTEGER(KIND=8)
#endif
#if(PP_CGNS_INT==32)
#  define PP_CGNS_INT_TYPE INTEGER
#endif
