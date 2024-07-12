/*########################################################
# Name memwater
# Use:
#  o Holds functions doing a memory efficent Smith
#    Waterman pairwise alignments.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'   o header:
'     - Included libraries
'   o .c fun01: max_memwater
'     - Find the maximum value
'   o .c fun02: ifmax_memwater
'     - Set a value (ret) to a value based on which value
'       is greater.
'   o .c fun03: hiScore_memwater
'     - Selects the max score and direction selected for
'       the max score.
'   o .c fun04: scoreIndel_memwater
'     - Gets the indel score for a memwater alignment
'   o .c fun05: scoreGt0_memwater
'     - Checks to see if the score is greater then zero.
'       If not, this resets the input values
'   o fun06 memwater:
'     - Run a memory efficent Waterman Smith alignment on
'       input sequences
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "memwater.h"
#include "../generalLib/seqST.h"

#include "../generalAln/alnSet.h"

/*.h files only*/
#include "../generalLib/dataTypeShortHand.h"

#include "../generalAln/alnDefs.h"
#include "../generalAln/indexToCoord.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden files
!   o std #include <stdio.h>
!   o .h  #include "../generalLib/base10StrToNum.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: max_memwater
|  - Find the maximum value (branchless)
| Input:
|  o ret:
|    - Value to hold the maximum value
|  o x:
|    - First value to compare, ret is set to x if x >= y
|  o y:
|    - second value to compare, ret is set to y if y > x
| Output:
|  - Sets:
|    - Sets ret to x if x is greater than or equal to y
|    - Sets ret to y if y is greather than x
\-------------------------------------------------------*/
#define \
max_memwater(\
   ret,\
   x,\
   y\
)(\
   (ret) = (x) ^ (((x) ^ (y)) & (-((x) < (y))))\
   /*Logic:
   `  x < y:
   `    if x < y I get 0 (x > y); else I get 1
   `  -(x < y):
   `    If x < y I get -1 (111...)
   `    If x >= y I get 0
   `  x ^ ((x ^ y) & 111...): (only when x < y)
   `    This gives me x
   `  x ^ (x ^ y) & 0: (only when x >= y)
   `    This gives me y
   */\
)

/*-------------------------------------------------------\
| Fun02: ifmax_memwater
|  - Set a value (ret) to a value based on which value
|    is greater.
| Input:
|  o ret:
|    - This will hold the return value
|  o x:
|    - First value to compare, (if x >= y)
|  o y:
|    - second value to compare, (if y > x)
|  o xRet:
|    - Value to set ret of x is >= y
|  o yRet:
|    - Value to set ret of y is > x
| Output:
|  - Sets:
|    - ret to xRet if x is greater than or equal to y
|    - ret to yRet if y is greater than x
\-------------------------------------------------------*/
#define \
ifmax_memwater(\
   ret,\
   x,\
   y,\
   xRet,\
   yRet\
)(\
   (ret) = (xRet) ^ (((xRet) ^ (yRet)) & (-((x) < (y))))\
   /*This follows the same logic as max_memwater(ret,x,y),
   ' except instead of setting ret to the highest value, I
   ' set ret to xRet if x is >= to y or yRet if y > x.
   */ \
)

/*-------------------------------------------------------\
| Fun03: hiScore_memwater
|  - Selects the max score and direction selected for the
|    max score.
| Input:
|  - maxSc:
|    o This will hold the max score
|  - maxDir
|    o This will hold the direction of the max score
|  - maxPos:
|    o Index to store
|  - insSc:
|    o Score for having an insertion at this position
|  - snpSc
|    o Score for having an SNP/match at this position
|  - delSc:
|    o Score for having an deletion at this position
|  - snpPos:
|    o Index for and snp
|  - insPos:
|    o Index for a insertion
|  - delPos:
|    o Index for an deletion
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\-------------------------------------------------------*/

#define \
hiScore_memwater(\
   maxSc,\
   maxDir,\
   maxPos,\
   snpSc,\
   insSc,\
   delSc,\
   snpPos,\
   insPos,\
   delPos\
){\
   max_memwater((maxSc), (insSc), (snpSc));      /*5 Op*/\
   maxDir=((maxSc > delSc) << (snpSc < insSc))+1;/*4 Op*/\
   \
   /*Logic:
   `   - noDel: maxSC > delSc:
   `     o 1 if deletion not max score
   `     o 0 if deletion is max score
   `   - type: noDel << (snpSc < insSc):
   `     o 1 << 1 = 2 if insertion is maximum
   `     o 1 << 0 = 1 if snp is maximum
   `     o 0 << 0 = 0 if deletion is max, and snp > ins
   `     o 0 << 1 = 0 if deletion is max, but ins >= snp
   `   - dir: type + 1
   `     o adds 1 to change from stop to direction
   */\
   \
   ifmax_memwater(\
      (maxPos),\
      (insSc),\
      (snpSc),\
      (insPos),\
      (delPos)\
   );\
   \
   ifmax_memwater(\
      (maxPos),\
      (delSc),\
      (maxSc),\
      (maxPos),\
      (snpPos)\
   );\
   \
   max_memwater((maxSc), (delSc), (maxSc));  /*5 Op*/\
} /*hiScore_memwater*/

/*-------------------------------------------------------\
| Fun04: scoreIndel_memwater
|   - Gets the indel score for a memwater alignment
| Input:
|   - lastScoreSL:
|     o Previous score that I am using to find this score
|   - gapOpenSC:
|     o Gap opening penalty
|   - gapDiffSS:
|     o The gap_open_penalty - gap_extend_penalty. Used to
|       get the gap extension from the gap opening penalty
|   - dirSC:
|     o Direction moved (indel or snp) for last score
| Variations:
|   - fun04 var-a: nogapextend
|     o Does not use the gap extension penatly
|   - fun04 var-b: default
|     o Applies the gap extension penalty
| Output:
|   - Returns the calcualted score for an indel
\-------------------------------------------------------*/

/*_______________________________________________________\
@ Fun04 Var-A: nogapextend
@   - Does not use the gap extension penatly
\_______________________________________________________*/
#ifdef NOGAPEXTEND
   #define \
   scoreIndel_memwater(\
      lastScoreSL,\
      gapOpenSC,\
      gapDiffSS,\
      dirSC,\
      mvDirSC\
   )( (lastScoreSL) + (gapOpenSC) )

/*_______________________________________________________\
@ Fun04 Var-B: default
@   - Applies the gap extension penalty
\_______________________________________________________*/
#else
   #define \
   scoreIndel_memwater(\
      lastScoreSL,\
      gapOpenSC,\
      gapDiffSS,\
      dirSC,\
      mvDirSC\
   )(\
        (lastScoreSL) + (gapOpenSC)\
      + ((gapDiffSS) & (long) (-((dirSC) == (mvDirSC))))\
   )
   /*Logic:
   `   - extend: -(dirSC == mvDirSC)
   `     o Is -1 if direction is current mvDirSC direction
   `       x for ins mvDirSC should be def_mvIns_alnDefs
   `       x for del mvDirSC should be def_mvDel_alnDefs
   `     o Is 0 if direction is 
   `   - diff: gapDiffSS & extend
   `     o -1: using gap extension (gapDiffSS is kept)
   `     o 0: gap open (no gapDiffSS)
   `   - penalty: gapOpenSC + diff
   `     o diff = -1 makes gap exetension penalty
   `     o diff = 0 makes the gap opening penalty
   `   - lastScoreSL + penalty
   `     o score for insertion
   */
#endif

/*-------------------------------------------------------\
| Fun05: scoreGt0_memwater
|   - Checks to see if the score is greater then zero. If
|     not, this resets the input values
| Input:
|   - scoreSL:
|     o Score to check if greater than zero
|   - dirOnSC:
|     o Direction for the score
|   - indexUL:
|     o Index at in the memwater alignment
|   - curIndexUL:
|     o Current index (applied if score is <= zero)
| Output:
|   - Modifies:
|     o Sets scoreSL to 0 if scoreSL <= 0
|     o Sets dirOnSC to 0 if scoreSL <= 0
|     o sets indexUL to curIndexUL if scoreSL <= 0
\-------------------------------------------------------*/
#define \
scoreGt0_memwater(\
   scoreSL,\
   dirOnSC,\
   indexUL,\
   curIndexUL\
){\
   long keepDirSL = (long) -((scoreSL) > 0);\
   (dirOnSC) &= keepDirSL;\
   (scoreSL) &= keepDirSL;\
   (indexUL) &= keepDirSL;\
   (indexUL) |= ( (curIndexUL) & (~keepDirSL) );\
} /*scoreGt0_memwater*/

/*-------------------------------------------------------\
| Fun06: memwater
|   - Performs a memory efficent Smith Waterman alignment
|     on a pair of sequences
| Input;
|   - qrySeqSTVoidPtr:
|     o Point to an seqST with the query sequence and
|       index 0 coordinates to start (offsetUL)/end
|       (endAlnUL) the alignment.
|   - refSeqSTVoidPtr:
|     o Point to an seqST with the reference sequence
|       and index 0 coordinates to start (offsetUL)/end
|       (endAlnUL) the alignment.
|   - refStartUL:
|     o Pointer to unsigned long to hold the frist
|       reference base in the alignment
|   - refEndUL:
|     o Pointer to unsigned long to hold the last
|       reference base in the alignment
|   - qryStartUL:
|     o Pointer to unsigned long to hold the frist query 
|       base in the alignment
|   - qryEndUL:
|     o Pointer to unsigned long to hold the last query
|       base in the alignment
|   - alnSetVoidPtr:
|     o Pointer to an alnSet structure with the gap open,
|       gap extend, and scoring matrix for the alingment
| Output:
|  - Modifies:
|    o refStartUL to have 1st reference base in alignment
|    o refEndUL to have last reference base in alignment
|    o qryStartUL to have first query base in alignment
|    o qryEndUL to have last query base in alignment
|  - Returns:
|    o Score for aligment
|    o 0 for memory errors
\-------------------------------------------------------*/
long
memwater(
    void *qrySeqSTVoidPtr,
    void *refSeqSTVoidPtr,
    unsigned long *refStartUL,
    unsigned long *refEndUL,
    unsigned long *qryStartUL,
    unsigned long *qryEndUL,
    void *alnSetVoidPtr      /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun06 TOC: memwaterAln
   '  - Run a memory efficent Waterman Smith alignment on
   '    input sequences
   '  o fun06 sec01:
   '    - Variable declerations
   '  o fun06 sec02:
   '    - Allocate memory for alignment
   '  o fun06 sec03:
   '    - Fill in initial negatives for ref
   '  o fun0 sec04:
   '    - Fill the matrix with scores
   '  o fun06 sec05:
   '    - Set up for returing matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec01: Variable declerations
   ^  o fun06 sec01 sub01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun06 sec01 sub02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun06 sec01 sub03:
   ^    - Directinol matrix variables
   ^  o fun06 sec01 sub04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec01 Sub01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   /*Set up the void pointers*/
   struct seqST *qryST =
      (struct seqST *) qrySeqSTVoidPtr;

   struct seqST *refST =
      (struct seqST *) refSeqSTVoidPtr;

   struct alnSet *settings =
      (struct alnSet *) alnSetVoidPtr;

   long scoreSL = 0; /*Score to return*/
   ulong bestStartUL = 0; /*Records best starting index*/
   ulong bestEndUL = 0;   /*Records best ending index*/

   /*Get start & end of query and reference sequences*/
   char *refSeqStr = 0;
   char *qrySeqStr = 0;

   ulong lenRefUL = refST->endAlnUL - refST->offsetUL + 1;
   ulong lenQryUL = qryST->endAlnUL - qryST->offsetUL + 1;

   /*Iterators for loops*/
   ulong ulRef = 0;
   ulong ulQry = 0;

   /*****************************************************\
   * Fun06 Sec01 Sub02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   long snpScoreSL = 0;
   long insScoreSL = 0;
   long delScoreSL = 0;    /*Score for doing an deletion*/
   long nextSnpScoreSL = 0;/*Score for next match/snp*/
   long *scoreHeapSL = 0;  /*matrix to use in alignment*/

   #ifndef NOGAPEXTEND
      /*Used in finding if useing gap extension*/
      short gapDiffS =
         settings->extendSS - settings->gapSS;
   #endif

   /*****************************************************\
   * Fun06 Sec01 Sub03:
   *  - Directional matrix variables
   \*****************************************************/

   /*Direction matrix (1 cell holds a single direction)*/
   char *dirRowHeapSC = 0;  /*Holds directions*/

   /*Keeping track of alignment starting positions*/
   ulong indexUL = 0;      /*Index I am at in the matrix*/
   ulong *indexHeapUL=0;   /*Row of starting indexes*/
   ulong *oldIndexHeapUL=0;/*Last round starting indexes*/
   ulong *swapPtrUL = 0;   /*For swapping ulongs*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec02:
   ^  - Allocate memory for alignment
   ^  o fun06 sec02 sub01:
   ^    - Allocate memory for the alignment
   ^  o fun06 sec02 sub02:
   ^    - Allocate memory for keeping track of indexes
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec02 Sub01:
   *  - Allocate memory for the alignment
   \****************************************************/

   dirRowHeapSC = malloc((lenRefUL + 1) * sizeof(char));

   if(dirRowHeapSC == 0)
      goto memErr_fun06_sec05_sub01;

   scoreHeapSL = malloc((lenRefUL + 1) * sizeof(long));
   /*+ 1 is for the indel column*/

   if(scoreHeapSL == 0)
      goto memErr_fun06_sec05_sub01;

   /*****************************************************\
   * Fun06 Sec02 Sub02:
   *  - Get memory for keeping track of starting indexes
   \*****************************************************/

   /*Set up the first row of starting indexes*/
   indexHeapUL = malloc((lenRefUL + 1) * sizeof(ulong));

   if(indexHeapUL == 0)
      goto memErr_fun06_sec05_sub01;

   /*Set up the second row of indexs (so have two rows)*/
   oldIndexHeapUL = malloc((lenRefUL +1) * sizeof(ulong));

   if(oldIndexHeapUL == 0)
      goto memErr_fun06_sec05_sub01;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      indexUL = 0;
      indexUL <= lenRefUL;
      ++indexUL
   ){ /*loop; till have initalized the first row*/
      dirRowHeapSC[indexUL] = def_mvStop_alnDefs;
      indexHeapUL[indexUL] = indexUL;
      scoreHeapSL[indexUL] = 0;
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec04:
   ^  - Fill the matrix with scores
   ^  o fun06 sec04 sub01:
   ^    - Final set up before scoring the matrix
   ^  o fun06 sec04 sub02:
   ^    - Start loops and get each score
   ^  o fun06 sec04 sub03:
   ^    - Check if is an alternative base best score
   ^  o fun06 sec04 sub04:
   ^    - Find the best score for the last base
   ^  o fun06 sec04 sub05:
   ^    - Is last base in row an alternative alignment?
   ^  o fun06 sec04 sub06:
   ^    - Prepare to score the next row in the matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec04 Sub01:
   *  - Final set up before scoring the matrix
   \*****************************************************/

   /*Move the row of starting indexes to the last row*/
   swapPtrUL = indexHeapUL;
   indexHeapUL = oldIndexHeapUL;
   oldIndexHeapUL = swapPtrUL;

   indexHeapUL[0] = indexUL;

   nextSnpScoreSL = 0;

   /*These are always negative*/
   delScoreSL = 0;
   indexHeapUL[0] = indexUL;
   /*dirRowHeapSC is already set to stop*/
   /*scoreHeapSL is already set to 0*/

   /*Incurment to the frist base*/
   ++indexUL;
   qrySeqStr = qryST->seqStr + qryST->offsetUL;
   refSeqStr = refST->seqStr + refST->offsetUL - 1;
      /*Offseting reference by 1 to account for the gap
      `  column
      */

   /*****************************************************\
   * Fun06 Sec04 Sub02:
   *  - Start loops and get each score
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(
      ulQry = 0;
      ulQry < lenQryUL;
      ++ulQry
   ){ /*loop; compare query base against all ref bases*/

      for(
         ulRef = 1;
         ulRef <= lenRefUL;
         ++ulRef
      ){ /*loop; compare one query to one reference base*/

         snpScoreSL =
            getScore_alnSet(
               qrySeqStr[ulQry],
               refSeqStr[ulRef],
               settings
            ); /*Find the score for the base pairs*/

         snpScoreSL += nextSnpScoreSL;
         nextSnpScoreSL = scoreHeapSL[ulRef];

         insScoreSL =
            scoreIndel_memwater(
               scoreHeapSL[ulRef],
               settings->gapSS,
               gapDiffS,
               dirRowHeapSC[ulRef],
               def_mvIns_alnDefs
            );

         hiScore_memwater(
            scoreHeapSL[ulRef],
            dirRowHeapSC[ulRef],
            indexHeapUL[ulRef],
            snpScoreSL,
            insScoreSL,
            delScoreSL,
            oldIndexHeapUL[ulRef - 1], /*snp index*/
            oldIndexHeapUL[ulRef], /*insertion index*/
            indexHeapUL[ulRef - 1] /*Deletion index*/
         );

         scoreGt0_memwater(
            scoreHeapSL[ulRef],
            dirRowHeapSC[ulRef],
            indexHeapUL[ulRef],
            indexUL
         );

         delScoreSL =
            scoreIndel_memwater(
               scoreHeapSL[ulRef],
               settings->gapSS,
               gapDiffS,
               dirRowHeapSC[ulRef],
               def_mvDel_alnDefs
            );

         /***********************************************\
         * Fun06 Sec04 Sub03:
         *  - Determine if is best score (keep as primary)
         \***********************************************/

         if(scoreSL < scoreHeapSL[ulRef])
         { /*If: this is the best score*/
            scoreSL = scoreHeapSL[ulRef];
            bestStartUL = indexHeapUL[ulRef];
            bestEndUL = indexUL;
         } /*If: this was an snp or match*/

         ++indexUL;
      } /*loop; compare one query to one reference base*/

     /***************************************************\
     *  Fun06 Sec04 Sub06:
     *   - Prepare for the next round
     \***************************************************/

     /*Get scores set up for the gap column*/
	  nextSnpScoreSL = 0;
     delScoreSL = 0;

     /*Swap index arrays so the current is last*/
     swapPtrUL = indexHeapUL;
     indexHeapUL = oldIndexHeapUL;
     oldIndexHeapUL = swapPtrUL;

     indexHeapUL[0] = indexUL;

     ++indexUL; /*Set index for the next base pair*/
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   ^  o fun06 sec05 sub01:
   ^    - clean up
   ^  o fun06 sec05 sub02:
   ^    - find the best score
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec05 Sub01:
   *  - clean up
   \*****************************************************/

   indexToCoord(
      lenRefUL,
      bestStartUL,
      *refStartUL,
      *qryStartUL
   ); /*Convert the starting index to coordinates*/

   *refStartUL += refST->offsetUL;
   *qryStartUL += qryST->offsetUL;

   indexToCoord(
      lenRefUL,
      bestEndUL,
      *refEndUL,
      *qryEndUL
   ); /*Convert the ending index to coordinates*/

   *refEndUL += refST->offsetUL;
   *qryEndUL += qryST->offsetUL;

   free(dirRowHeapSC);
   dirRowHeapSC = 0;

   free(scoreHeapSL);
   scoreHeapSL = 0;

   free(indexHeapUL);
   indexHeapUL = 0;

   free(oldIndexHeapUL);
   oldIndexHeapUL = 0;

   return scoreSL;

   memErr_fun06_sec05_sub01:;

   free(dirRowHeapSC);
   dirRowHeapSC = 0;

   free(scoreHeapSL);
   scoreHeapSL = 0;

   free(indexHeapUL);
   indexHeapUL = 0;

   free(oldIndexHeapUL);
   oldIndexHeapUL = 0;

   return 0;
} /*memwater*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible, this code is under the
:   MIT license.
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy, modify, publish, use, compile,
:   sell, or distribute this software, either in source
:   code form or as a compiled binary, for any purpose,
:   commercial or non-commercial, and by any means.
: 
: In jurisdictions that recognize copyright laws, the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
:   CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information, please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2024 jeremyButtler
: 
: Permission is hereby granted, free of charge, to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"), to
:   deal in the Software without restriction, including
:   without limitation the rights to use, copy, modify,
:   merge, publish, distribute, sublicense, and/or sell
:   copies of the Software, and to permit persons to whom
:   the Software is furnished to do so, subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
:   AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
:   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
