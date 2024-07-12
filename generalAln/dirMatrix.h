/*########################################################
# Name: dirMatrix
#  - holds functions for dealing with the dirMatrix
#    returned by water and needle
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - forward declerations and guards
'   o .h st01: alnMatrixStruct
'     - Holds the direction matrix and best score(s) for a
'       single aligment
'   o .h fun01: blank_dirMatrix
'     - blanks all score, index, and length variables in a
'       dirMatrix structure
'   o .h fun02: init_dirMatrix
'     - initializes a dirMatrix structure
'   o fun03: freeStack_dirMatrix
'     - frees heap allocated variables in a dirMatrix
'       structure
'   o fun04: freeHeap_dirMatrix
'     - frees a dirMatrix structure
'   o fun05: getAln_dirMatrix
'     - gets a sam file entry (alignment) from a direction
'       matrix (inside the dirMatrix structure)
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - forward declerations and guards
\-------------------------------------------------------*/

#ifndef DIRECTION_MATRIX_H
#define DIRECTION_MATRIX_H

typedef struct samEntry samEntry;
typedef struct seqST seqST;
typedef struct alnSet alnSet;

/*-------------------------------------------------------\
| St01: dirMatrix
|  - Holds the direction matrix and best score(s) for a
|    single aligment
\-------------------------------------------------------*/
typedef struct dirMatrix
{ /*alnStruct*/
  signed char *dirMatrixSC; /*directional matrix*/
  unsigned long lenMatrixUL;/*size of directional matrix*/

  signed long scoreSL;       /*best score in alignment*/
  unsigned long indexUL;     /*index for best score*/

  signed long *scoreArySL;   /*score array used to score*/
  unsigned long lenScoreUL;  /*length of score array*/

  unsigned long lenRefUL;    /*reference length*/
  unsigned long refOffsetUL; /*first ref base to align*/
  unsigned long refEndUL;    /*last ref base to align*/

  unsigned long lenQryUL;    /*length of query*/
  unsigned long qryOffsetUL; /*first query base to align*/
  unsigned long qryEndUL;    /*last query base to align*/
}dirMatrix;

/*-------------------------------------------------------\
| Fun01: blank_dirMatrix
|   - blanks all score, index, and length variables in a
|     dirMatrix structure
| Input:
|   - dirMatrixSTPtr:
|     o a pointer to a dirMatrix to blank
| Output:
|   - Sets:
|     o all score, length, and index variables are set to
|       0 (direction matrix is not touched)
\-------------------------------------------------------*/
#define \
blank_dirMatrix( \
  dirMatrixSCSTPtr \
){ \
   (dirMatrixSCSTPtr)->scoreSL = 0; \
   (dirMatrixSCSTPtr)->indexUL = 0; \
   \
   (dirMatrixSCSTPtr)->lenRefUL = 0; \
   (dirMatrixSCSTPtr)->refOffsetUL = 0; \
   (dirMatrixSCSTPtr)->refEndUL = 0; \
   \
   (dirMatrixSCSTPtr)->lenQryUL = 0; \
   (dirMatrixSCSTPtr)->qryOffsetUL = 0; \
   (dirMatrixSCSTPtr)->qryEndUL = 0; \
} /*blank_dirMatrixSC*/

/*-------------------------------------------------------\
| Fun02: init_dirMatrix
|   - initializes a dirMatrix structure
| Input:
|   - dirMatrixSTPtr:
|     o a pointer to a dirMatrix structure initialize
| Output:
|   - Sets:
|     o all variables in matrixSTPtr to 0
\-------------------------------------------------------*/
#define \
init_dirMatrix( \
  dirMatrixSTPtr \
){ \
   blank_dirMatrix((dirMatrixSTPtr)); \
   \
   (dirMatrixSTPtr)->dirMatrixSC = 0; \
   (dirMatrixSTPtr)->lenMatrixUL = 0; \
   \
   (dirMatrixSTPtr)->scoreArySL = 0; \
   (dirMatrixSTPtr)->lenScoreUL = 0; \
} /*init_dirMatrixSC*/

/*-------------------------------------------------------\
| Fun03: freeStack_dirMatrix
|   - frees heap allocated variables in a dirMatrix
|     structure
| Input:
|   - matrixSTPtr
|     o pointer to dirMatrix structure with variables to
|       free
| Output:
|   - Frees:
|     o dirMatrix->dirMatrixSC
|   - Sets:
|     o all non-freeded variables to 0
\-------------------------------------------------------*/
void
freeStack_dirMatrix(
   struct dirMatrix *matrixSTPtr
);

/*-------------------------------------------------------\
| Fun04: freeHeap_dirMatrix
|   - Frees a dirMatrix structure
| Input:
|   - matrixSTPtr
|     o pointer to a dirMatrix structure to free
| Output:
|   - Frees:
|     o matrixSTPtr
\-------------------------------------------------------*/
void
freeHeap_dirMatrix(
   struct dirMatrix *matrixSTPtr
);

/*-------------------------------------------------------\
| Fun05: getAln_dirMatrix
|   - gets a sam file entry (alignment) from a direction
|     matrix (inside the dirMatrix structure)
| Input:
|   - matrixSTPtr
|     o pointer to a dirMatrix structure to get alignment
|       from
|   - indexUL:
|     o index of last base in the alignment
|     o 0 to use index from matirxSTPtr
|   - revBl:
|     o 1: reverse alignment (sam flag is 16)
|     o 0: foward alignment (sam flag is 0)
|   - qrySTPtr:
|     o pointer to a seqST with the query sequence
|   - refSTPtr:
|     o pointer to a seqST with the reference sequence
|   - samSTPtr:
|     o pointer to a samEntry struct to hold the alignment
|   - numAnonUI:
|     o pointer to unsigned in to hold the number of
|       anonymous bases (matches only)
|   - alnSetSTPtr:
|     o pointer to alnSet structure with the match matrix
| Output:
|   - Modifies:
|     o samSTPtr to have the aligned sequence
|   - Returns:
|     o 0 for no errors
|     o 1 for memory error (only error possible)
\-------------------------------------------------------*/
signed char
getAln_dirMatrix(
   struct dirMatrix *matrixSTPtr,
   unsigned long indexUL,
   signed char revBl,
   struct seqST *qrySTPtr,
   struct seqST *refSTPtr,
   struct samEntry *samSTPtr,
   unsigned int *numAnonUI,
   struct alnSet *alnSetSTPtr
);

#endif
