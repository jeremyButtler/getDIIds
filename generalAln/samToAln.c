/*#######################################################\
# Name: samToAln
#   - has functions to convert a sam entry and the
#     reference sequence to an alignment
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../generalLib/samEntry.h"

#include "../generalAln/alnSet.h"

/*no .h files*/
#include "../generalLib/dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .h  #include "../generalLib/base10str.h"
!   o .h  #include "../generalLib/ulCp.h"
!   o .h  #include "../generalLib/numToStr.h"
!   o .h  #include "../generalLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: paln_samToAln
|   - prints an alignment in expanded cigar format
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct with aligned sequence
|   - refSTPtr:
|     o pointer to seqST struct with reference sequence
|   - scoreSL:
|     o score for the alignment
|   - numAnonUI:
|     o number of anonymous matches
|   - setSTPtr:
|     o pointer to alnSet struct with print settings
|   - outFILE:
|     o pointer to file to write alignment to
| Outputs:
|   - Prints:
|     o aligned sequences to outFILE
|   - Returns:
|     o 0 for success
|     o 1 for memory error
\-------------------------------------------------------*/
signed char
paln_samToAln(
   struct samEntry *samSTPtr,
   struct seqST *refSTPtr,
   signed long scoreSL,
   uint numAnonUI,
   struct alnSet *setSTPtr,
   void *outFILE
){
   ushort lenBuffUS = 1 << 10; /*1 << 10 ~ 1024*/
   schar qryBuffStr[lenBuffUS];
   schar refBuffStr[lenBuffUS];
   schar cigBuffStr[lenBuffUS];
   ushort numNtUS = 0;   /*number nucleotides in buffers*/
   ushort ntInLineUS = 0; /*number nucleotides in line*/

   uint uiNt = 0;      /*nucleotide on*/

   fprintf(
      (FILE *) outFILE,
      "##################################################"
   );
} /*paln_samToAln*/
