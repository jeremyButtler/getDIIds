/*########################################################
# Name water
# Use:
#  o Holds functions to do a Smith Waterman pairwise
#    alignment
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'   o header:
'     - included libraries
'   o .c fun01: max_water
'     - find the maximum value
'   o .c fun02: ifmax_water
'     - set a value (ret) to a value based on which value
'       is greater.
'   o .c fun03: hiScore_water
'     - selects the max score and direction selected for
'       the max score.
'   o .c fun04: scoreIndel_water
'     - gets the indel score for a water alignment
'   o .c fun05: scoreGt0_water
'     - checks to see if the score is greater then zero.
'       If not, this resets the input values
'   o fun06 water:
'     - run a memory efficent Waterman Smith alignment on
'       input sequences
'   o license:
'     - licensing for this code (public domain / mit)
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

#include "water.h"

#include "../generalLib/seqST.h"

#include "../generalAln/alnSet.h"
#include "../generalAln/dirMatrix.h"

/*.h files only*/
#include "../generalLib/dataTypeShortHand.h"
#include "../generalAln/alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden files
!   o .c  #include <stdio.h>
!   o .h  #include "../generalLib/base10str.h"
!   o .h  #include "../generalLib/ulCp.h"
!   o .h  #include "../generalLib/numToStr.h"
!   o .h  #include "../generalLib/ntTo5Bit.h"
!   o .h  #include "../generalAln/indexToCoord.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: max_water
|  - find the maximum value (branchless)
| Input:
|  o ret:
|    - value to hold the maximum value
|  o x:
|    - first value to compare, ret is set to x if x >= y
|  o y:
|    - second value to compare, ret is set to y if y > x
| Output:
|  - Sets:
|    - sets ret to x if x is greater than or equal to y
|    - sets ret to y if y is greather than x
\-------------------------------------------------------*/
#define \
max_water(\
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
| Fun02: ifmax_water
|  - set a value (ret) to a value based on which value
|    is greater.
| Input:
|  o ret:
|    - this will hold the return value
|  o x:
|    - first value to compare, (if x >= y)
|  o y:
|    - second value to compare, (if y > x)
|  o xRet:
|    - value to set ret of x is >= y
|  o yRet:
|    - value to set ret of y is > x
| Output:
|  - Sets:
|    - ret to xRet if x is greater than or equal to y
|    - ret to yRet if y is greater than x
\-------------------------------------------------------*/
#define \
ifmax_water(\
   ret,\
   x,\
   y,\
   xRet,\
   yRet\
)(\
   (ret) = (xRet) ^ (((xRet) ^ (yRet)) & (-((x) < (y))))\
   /*This follows the same logic as max_water(ret, x, y),
   ' except instead of setting ret to the highest value, I
   ' set ret to xRet if x is >= to y or yRet if y > x.
   */ \
)

/*-------------------------------------------------------\
| Fun03: hiScore_water
|  - selects the max score and direction selected for the
|    max score. Priority is deltions, insertions, then
|    snps
| Input:
|  - maxSc:
|    o this will hold the max score
|  - maxDir
|    o this will hold the direction of the max score
|  - insSc:
|    o score for having an insertion at this position
|  - snpSc
|    o score for having an SNP/match at this position
|  - delSc:
|    o score for having an deletion at this position
| Output:
|  - Sets:
|    o sets maxDir to the direction of the max score
|    0 sets maxSc to the max score
\-------------------------------------------------------*/
#define \
hiScore_water(\
   maxSc,\
   maxDir,\
   snpSc,\
   insSc,\
   delSc\
){\
   max_water((maxSc), (insSc), (snpSc));         /*5 Op*/\
   maxDir=((maxSc > delSc) << (snpSc < insSc))+1;/*4 Op*/\
   max_water((maxSc), (delSc), (maxSc));         /*5 Op*/\
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
}

/*-------------------------------------------------------\
| Fun04: scoreIndel_water
|   - gets the indel score for a water alignment
| Input:
|   - lastScoreSL:
|     o previous score that I am using to find this score
|   - gapOpenSC:
|     o gap opening penalty
|   - gapDiffSSS:
|     o the gap_open_penalty - gap_extend_penalty. Used to
|       get the gap extension from the gap opening penalty
|   - dirSC:
|     o direction moved (indel or snp) for last score
| Variations:
|   - fun04 var-a: nogapextend
|     o does not use the gap extension penatly
|   - fun04 var-b: default
|     o applies the gap extension penalty
| Output:
|   - Returns:
|     o the calcualted score for an indel
\-------------------------------------------------------*/

/*_______________________________________________________\
@ Fun04 Var-A: nogapextend
@   - does not use the gap extension penatly
\_______________________________________________________*/
#ifdef NOGAPEXTEND
   #define \
   scoreIndel_water(\
      lastScoreSL,\
      gapOpenSC,\
      gapDiffSSS,\
      dirSC,\
      mvDirBl
   )( (lastScoreSL) + (gapOpenSC) )

/*_______________________________________________________\
@ Fun04 Var-B: default
@   - applies the gap extension penalty
\_______________________________________________________*/
#else
   #define \
   scoreIndel_water(\
      lastScoreSL,\
      gapOpenSC,\
      gapDiffSS,\
      dirSC,\
      mvDirSC\
   )(\
        (lastScoreSL) + (gapOpenSC) \
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
| Fun05: scoreGt0_water
|   - checks to see if the score is greater then zero. If
|     not, this resets the input values
| Input:
|   - scoreSL:
|     o score to check if greater than zero
|   - dirOnSC:
|     o direction for the score
| Output:
|   - Modifies:
|     o sets scoreSL to 0 if scoreSL <= 0
|     o sets dirOnSC to 0 if scoreSL <= 0
\-------------------------------------------------------*/
#define \
scoreGt0_water( \
   scoreSL, \
   dirOnSC \
){ \
   long keepDirSL = (long) -((scoreSL) > 0); \
   (dirOnSC) &= keepDirSL; \
   (scoreSL) &= keepDirSL; \
} /*scoreGt0_water*/

/*-------------------------------------------------------\
| Fun06: water
|   - run a Waterman Smith alignment on input sequences
| Input:
|   - qrySTPtr:
|     o pointer to seqST with query sequence 
|       - qrySTPtr->offsetUL; first query base to align
|       - qrySTPtr->endAlnUL; last query base to align
|   - refSTPtr:
|     o pointer to seqST with reference sequence 
|       - refSTPtr->offsetUL; 1st reference base to align
|       - refSTPtr->endAlnUL; last reference base to align
|   - matrixSTPtr:
|     o pointer to dirMatrix to use for the alingment
|   - alnSet:
|     o pointer to alnSet with alignment settings
| Output:
|  - Modifies:
|    o allocates memory for dirMatrixSC and scoreAryUL
|      if they are to small
|    o updates lenMatrixUL and lenScoreUL if dirMatrixSC
|      or scoreAryUL are resized
|  - Returns:
|    o 0 for memory error
|    o score for alignment
\-------------------------------------------------------*/
signed long
water(
   struct seqST *qrySTPtr, /*query sequence and data*/
   struct seqST *refSTPtr, /*ref sequence and data*/
   struct dirMatrix *matrixSTPtr, /*direction matrix*/
   struct alnSet *settings     /*settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun06 TOC: WatermanAln
   '  - Run a Waterman Smith alignment on input sequences
   '  o fun06 sec01:
   '    - variable declerations
   '  o fun06 sec02:
   '    - allocate memory for alignment
   '  o fun06 sec03:
   '    - fill in initial negatives for ref
   '  o fun06 sec04:
   '    - fill the matrix with scores
   '  o fun06 sec05:
   '    - set up for returing matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec01: Variable declerations
   ^  o fun06 sec01 sub01:
   ^    - variables dealing with the query and reference
   ^      starting positions
   ^  o fun06 sec01 sub02:
   ^    - variables holding the scores (only two rows)
   ^  o fun06 sec01 sub03:
   ^    - directinol matrix variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec01 Sub01:
   *  - variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   /*Get start & end of query and reference sequences*/
   char *refSeqStr =
      refSTPtr->seqStr + refSTPtr->offsetUL;

   char *qrySeqStr =
      qrySTPtr->seqStr + qrySTPtr->offsetUL;

   /*Find the length of the reference and query*/
   ulong lenQryUL =
      qrySTPtr->endAlnUL - qrySTPtr->offsetUL + 1;

   ulong lenRefUL =
      refSTPtr->endAlnUL - refSTPtr->offsetUL + 1;
     /*The + 1 is to account for index 0 of endAlnUL*/

   ulong lenMatrixUL = (lenRefUL + 1) * (lenQryUL + 1);
     /*+1 for the gap column and row*/

   ulong ulRef = 0;
   ulong ulQry = 0;

   /*Set up counters for the query and reference base
   `  index
   */
   /*****************************************************\
   * Fun06 Sec01 Sub02:
   *  - variables holding the scores (only two rows)
   \*****************************************************/

   slong snpScoreSL = 0;    /*Score for deletion*/
   slong nextSnpScoreSL = 0;/*Score for match/snp*/

   slong insScoreSL = 0;    /*Score for deletion*/
   slong delScoreSL = 0;    /*Score for deletion*/

   /*Marks when to reset score buffer (every second row)*/
   slong *scoreArySL = 0;/*scoring row for alignment*/

   /*Gap penalities*/
   sshort gapDiffSS= settings->extendSS - settings->gapSS;

   /*****************************************************\
   * Fun06 Sec01 Sub03:
   *  - directinol matrix variables
   \*****************************************************/

   /*Direction matrix (one cell holds one direction)*/
   schar *dirMatrixSC = 0;/*Direction matrix*/
   schar *insDir = 0;    /*Direction above cell*/
   ulong indexUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec02:
   ^   - allocate memory for alignment
   ^   o fun06 sec02 sub01:
   ^     - set up the directional matrix
   ^   o fun06 sec02 sub02:
   ^     - set up score array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec02 Sub01:
   *   - set up the directional matrix
   \*****************************************************/

   if(matrixSTPtr->lenMatrixUL < lenMatrixUL)
   { /*If: need to resize the matrix*/
      free(matrixSTPtr->dirMatrixSC);
      matrixSTPtr->dirMatrixSC = 0;
      matrixSTPtr->lenMatrixUL = 0;

      matrixSTPtr->dirMatrixSC =
         malloc((lenMatrixUL + 1) * sizeof(char));

      if(matrixSTPtr->dirMatrixSC == 0)
         goto memErr_fun06_sec05;

       matrixSTPtr->lenMatrixUL = lenMatrixUL;
   } /*If: need to resize the matrix*/

   dirMatrixSC = matrixSTPtr->dirMatrixSC;

   blank_dirMatrix(matrixSTPtr);

   matrixSTPtr->lenRefUL = lenRefUL;
   matrixSTPtr->refOffsetUL = refSTPtr->offsetUL;
   matrixSTPtr->refEndUL = refSTPtr->endAlnUL;

   matrixSTPtr->lenQryUL = lenQryUL;
   matrixSTPtr->qryOffsetUL = qrySTPtr->offsetUL;
   matrixSTPtr->qryEndUL = qrySTPtr->endAlnUL;

   /*****************************************************\
   * Fun06 Sec02 Sub02:
   *   - set up score array
   \*****************************************************/

   if(matrixSTPtr->lenScoreUL < lenRefUL + 1)
   { /*If: need to make a larger score array*/
      free(matrixSTPtr->scoreArySL);
      matrixSTPtr->scoreArySL = 0;
      matrixSTPtr->lenScoreUL = 0;
      
      matrixSTPtr->scoreArySL =
         calloc((lenRefUL + 1), sizeof(long));

      if(! matrixSTPtr->scoreArySL)
         goto memErr_fun06_sec05;

      matrixSTPtr->lenScoreUL = lenRefUL;
   } /*If: need to make a larger score array*/

   scoreArySL = matrixSTPtr->scoreArySL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec03:
   ^  - fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      indexUL = 0;
      indexUL <= lenRefUL;
      ++indexUL
   ){ /*Loop: initialize the first row*/
      dirMatrixSC[indexUL] = def_mvStop_alnDefs;
      scoreArySL[indexUL] = 0;
   } /*Loop: initialize the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec04:
   ^  - fill the matrix with scores
   ^  o fun06 sec04 sub01:
   ^    - get initial scores
   ^  o fun06 sec04 sub02:
   ^    - start loops
   ^  o fun06 sec04 sub03:
   ^    - find bests score for each base pair
   ^  o fun06 sec04 sub04:
   ^    - check if have the best score
   ^  o fun06 sec04 sub05:
   ^    - set up for scoring next row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec04 Sub01:
   *  - get initial scores
   \*****************************************************/

   /*set up pointers*/
   insDir = dirMatrixSC;

   /*set up scores*/
   delScoreSL = 0;
   nextSnpScoreSL = 0;
   dirMatrixSC[indexUL] = def_mvStop_alnDefs;

   /*These are always negative*/
   ++indexUL;
   refSeqStr = refSTPtr->seqStr + refSTPtr->offsetUL - 1;
   qrySeqStr = qrySTPtr->seqStr + qrySTPtr->offsetUL;

   /*****************************************************\
   * Fun06 Sec04 Sub02:
   *   - start loops
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
     ){ /* loop; compare one query to one reference base*/

        /************************************************\
        * Fun06 Sec04 Sub03:
        *   - find bests score for each base pair
        \************************************************/

        snpScoreSL =
           getScore_alnSet(
              qrySeqStr[ulQry],
              refSeqStr[ulRef],
              settings
           ); /*Find the score for the base pairs*/

        snpScoreSL += nextSnpScoreSL;
        nextSnpScoreSL = scoreArySL[ulRef];

        insScoreSL =
           scoreIndel_water(
              scoreArySL[ulRef],
              settings->gapSS,
              gapDiffSS,
              insDir[ulRef],
              def_mvIns_alnDefs
           );

        hiScore_water(
           scoreArySL[ulRef],
           dirMatrixSC[indexUL],
           snpScoreSL,
           insScoreSL,
           delScoreSL
        );

        scoreGt0_water(
           scoreArySL[ulRef],
           dirMatrixSC[indexUL]
        );

        delScoreSL =
           scoreIndel_water(
              scoreArySL[ulRef],
              settings->gapSS,
              gapDiffSS,
              dirMatrixSC[indexUL],
              def_mvDel_alnDefs
           );

       /*************************************************\
       * Fun06 Sec04 Sub04:
       *   - check if have the best score
       \*************************************************/

        /*This is faster than the branchless option.
        ` I am guessing to much is done in the if and 
        ` that the if if fired rarely.
        */
        if(
            matrixSTPtr->scoreSL
          < scoreArySL[ulRef]
        ){ /*if have a new best score*/
           matrixSTPtr->scoreSL = scoreArySL[ulRef];
           matrixSTPtr->indexUL = indexUL;
        } /*if have a new best score*/

        ++indexUL;
     } /*loop; compare one query to one reference base*/

      /**************************************************\
      *  Fun06 Sec04 Sub05:
      *    - set up for scoring next row
      \**************************************************/

       delScoreSL = 0;
       nextSnpScoreSL = 0; /*indel column is always 0*/

       insDir += ulRef;
       dirMatrixSC[indexUL] = def_mvStop_alnDefs;

       ++indexUL;
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec05:
   ^  - set up for returing the matrix (clean up/wrap up)
   ^  o fun06 sec05 sub01:
   ^    - no error clean up and return
   ^  o fun06 sec05 sub02:
   ^    - error clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec05 Sub01:
   *   - no error clean up and return
   \*****************************************************/

   /*Move back to the lower right conor cell
   ` This is not needed, but is nice.
   */

   --indexUL; /*get off last -1*/
   dirMatrixSC[indexUL] = def_mvStop_alnDefs;
   return matrixSTPtr->scoreSL; /*best score*/

   /*****************************************************\
   * Fun06 Sec05 Sub02:
   *   - error clean up and return
   \*****************************************************/

   memErr_fun06_sec05:;
   goto errCleanUp_fun06_sec05;
   errCleanUp_fun06_sec05:;

   return 0;
} /*water*/

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
