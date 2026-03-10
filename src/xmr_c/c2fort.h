#ifndef _C2FORT_H
#define _C2FORT_H
//======================================================================================
// Defines a macro EXTFORTNAME(x) to facilitate declaring external stuff in C
// that is compiled in Fortran.
//
#ifdef EXTFORTNAMES_NO_
#define EXTFORTNAME(NAME) NAME
#else
#define EXTFORTNAME(NAME) NAME##_
#endif
//======================================================================================
#endif
