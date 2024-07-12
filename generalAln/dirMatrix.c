/*########################################################
# Name: dirMatrix
#  - holds functions for dealing with the dirMatrix
#    returned by water and needle
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
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
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "../generalLib/samEntry.h"
#include "../generalLib/seqST.h"

#include "dirMatrix.h"
#include "alnSet.h"

/*.h files only*/
#include "../generalLib/dataTypeShortHand.h"
#include "../generalLib/ulCp.h"
#include "../generalLib/charCp.h"
#include "../generalLib/numToStr.h"

#include "alnDefs.h"
#include "indexToCoord.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .h  #include "../generalLib/base10str.h"
!   o .h  #include "../generalLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
){
   if(matrixSTPtr)
   { /*If: not null*/
      free(matrixSTPtr->dirMatrixSC);
      free(matrixSTPtr->scoreArySL);

      init_dirMatrix(matrixSTPtr); /*sets to null*/
   } /*If: not null*/

   return;
} /*freeAlnMatrixStack*/

/*-------------------------------------------------------\
| Fun04: freeHeap_dirMatrix
|   - frees a dirMatrix structure
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
){
   freeStack_dirMatrix(matrixSTPtr);
   free(matrixSTPtr);
} /*freeAlnMatrix*/

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
|       - this means I had to reverse complement the
|         reference sequence
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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC:
   '   - gets an alignment from a dirMatrix structure
   '   o fun05 sec01:
   '     - variable declerations
   '   o fun05 sec02:
   '     - find start and ending positions
   '   o fun05 sec03:
   '     - allocate memroy and copy query
   '   o fun05 sec04:
   '     - get alignment form matrix
   '   o fun05 sec05:
   '     - add starting softmasked bases and invert cigar
   '   o fun05 sec06:
   '     - add tags (NM, AS, nn)
   '   o fun05 sec07:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint qryPosUI = 0;
   uint refPosUI = 0;
   uint lenRefUI = matrixSTPtr->lenRefUL;

   schar *qrySeqStr = (schar *) qrySTPtr->seqStr;
   schar *refSeqStr = (schar *) refSTPtr->seqStr;
 
   schar *tmpStr = 0;
   uchar tmpUC = 0;
   schar matchBl = 0; /*check if had match or snp*/

   schar *dirMatrixSC = matrixSTPtr->dirMatrixSC;

   struct seqST seqDoNotFreeST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - find sequence start and ending positions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! indexUL)
      indexUL = matrixSTPtr->indexUL; /*no index input*/

   indexToCoord(
      lenRefUI,
      indexUL,
      refPosUI,
      qryPosUI
   ); /*find last aligned base in query/reference*/

   /*bases not aligned by user*/
   refSeqStr += matrixSTPtr->refOffsetUL;
   qrySeqStr += matrixSTPtr->qryOffsetUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec03:
   ^   - allocate memroy and copy query
   ^   o fun05 sec03 sub01:
   ^     - set up memory for cigar entry and blank entry
   ^   o fun05 sec03 sub02:
   ^     - add ending soft masked bases to cigar
   ^   o fun05 sec03 sub03:
   ^     - copy read id
   ^   o fun05 sec03 sub04:
   ^     - copy the reference id
   ^   o fun05 sec03 sub05:
   ^     - copy query sequence
   ^   o fun05 sec03 sub06:
   ^     - copy query q-score entry
   ^   o fun05 sec03 sub07:
   ^     - set flag and reference end
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec03 Sub01:
   *   - set up memory for cigar entry and blank entry
   \*****************************************************/

   /*is inefficent, but works*/
   if(samSTPtr->lenCigBuffUI < qrySTPtr->lenSeqUL)
   { /*If: want more cigar memory*/
      free(samSTPtr->cigTypeStr);
      samSTPtr->cigTypeStr = 0;

      /*cigar entry types*/
      samSTPtr->cigTypeStr =
         (char *)
         malloc((qrySTPtr->lenSeqUL + 1) * sizeof(schar));

      if(! samSTPtr->cigTypeStr)
         goto memErr_fun05_sec07;

      /*count for number bases supporting entry*/
      free(samSTPtr->cigValAryI);
      samSTPtr->cigValAryI = 0;

      samSTPtr->cigValAryI =
         (sint *)
         malloc((qrySTPtr->lenSeqUL + 1) * sizeof(sint));

      if(! samSTPtr->cigValAryI)
         goto memErr_fun05_sec07;

      samSTPtr->lenCigBuffUI = qrySTPtr->lenSeqUL;
   } /*If: want more cigar memory*/

   blank_samEntry(samSTPtr);

   /*****************************************************\
   * Fun05 Sec03 Sub02:
   *   - add ending soft masked bases to cigar
   \*****************************************************/

   if(qryPosUI + 1 < qrySTPtr->lenSeqUL)
   { /*If: need to add ending softmasked bases*/
      samSTPtr->cigTypeStr[0] = 'S';

      samSTPtr->cigValAryI[0] =
         qrySTPtr->lenSeqUL - qryPosUI - 1;
         /*-1 to account for lenSeqUL being index 1 and
         `   qryPosUI being index 0
         */
      samSTPtr->numMaskUI =
         qrySTPtr->lenSeqUL - qryPosUI - 1;
   } /*If: need to add ending softmasked bases*/

   /*else is set to null*/

   /*****************************************************\
   * Fun05 Sec03 Sub03:
   *   - copy read id
   \*****************************************************/

   /*check for header marks*/
   tmpUC = (qrySTPtr->idStr[0] == '>');
   tmpUC += (qrySTPtr->idStr[0] == '@');

   tmpStr = (schar *) qrySTPtr->idStr + tmpUC;

   while(*tmpStr++ > 32) ;

   --tmpStr;
   matchBl = *tmpStr;
   *tmpStr = '\0';

   samSTPtr->lenQryIdUC = 
      cpDelim_ulCp(
         (char *) samSTPtr->qryIdStr, 
         (char *) qrySTPtr->idStr + tmpUC,
         0,
         0
      ); /*copy read id*/

   *tmpStr = matchBl;

   /*****************************************************\
   * Fun05 Sec03 Sub04:
   *   - copy the reference id
   \*****************************************************/

   /*check for header marks*/
   tmpUC = (refSTPtr->idStr[0] == '>');
   tmpUC += (refSTPtr->idStr[0] == '@');

   tmpStr = (schar *) refSTPtr->idStr;

   while(*tmpStr++ > 32) ;

   --tmpStr;
   matchBl = *tmpStr;
   *tmpStr = '\0';

   samSTPtr->lenRefIdUC = 
      cpDelim_ulCp(
         (char *) samSTPtr->refIdStr, 
         (char *) refSTPtr->idStr + tmpUC,
         0,
         0
      ); /*copy read id*/

   *tmpStr = matchBl;

   /*****************************************************\
   * Fun05 Sec03 Sub05:
   *   - copy query sequence
   \*****************************************************/

   if(samSTPtr->lenSeqBuffUI < qrySTPtr->lenSeqUL)
   { /*If: I need more memory for the sequence*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         (char *)
         malloc((qrySTPtr->lenSeqUL+1024) * sizeof(char));

      if(! samSTPtr->seqStr)
         goto memErr_fun05_sec07;

      samSTPtr->lenSeqBuffUI = qrySTPtr->lenSeqUL + 1023;
      tmpStr = 0;
   } /*If: I need more memory for the sequence*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      qrySTPtr->seqStr,
      qrySTPtr->lenSeqUL
   );

   samSTPtr->readLenUI = qrySTPtr->lenSeqUL;

   /*convert the index to a real sequence*/
   indexToSeq_alnSet(samSTPtr->seqStr);

   /*****************************************************\
   * Fun05 Sec03 Sub06:
   *   - copy query q-score entry
   \*****************************************************/

   if(qrySTPtr->lenQUL)
   { /*If: have q-score entry*/
      if(samSTPtr->lenQBuffUI < qrySTPtr->lenQUL)
      { /*If: I need more memory for the sequence*/
         free(samSTPtr->qStr);
         samSTPtr->qStr = 0;

         samSTPtr->qStr =
            (char *)
            malloc((qrySTPtr->lenQUL+1024) *sizeof(char));

         if(! samSTPtr->qStr)
            goto memErr_fun05_sec07;

         samSTPtr->lenQBuffUI = qrySTPtr->lenQUL + 1023;
      } /*If: I need more memory for the sequence*/

      cpQEntry_samEntry(
         samSTPtr,
         (char *) qrySTPtr->qStr,
         1 /*need to make sure histogram is blanked*/
      );

      samSTPtr->readLenUI = qrySTPtr->lenSeqUL;
   } /*If: have q-score entry*/

   else
   { /*Else: set the q-score entry to nothing*/
      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';
   } /*Else: set the q-score entry to nothing*/

   /*****************************************************\
   * Fun05 Sec03 Sub07:
   *   - set flag and reference end
   \*****************************************************/

   if(revBl)
   { /*If: read is reverse complment forward*/
      samSTPtr->flagUS = 16;
      samSTPtr->refStartUI = lenRefUI - refPosUI - 1;
   } /*If: read is reverse complment forward*/

   else
   { /*Else: read is forward*/
      samSTPtr->flagUS = 0;
      samSTPtr->refEndUI = refPosUI;
   } /*Else: read is forward*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec04:
   ^   - get alignment form matrix
   ^   o fun05 sec04 sub01:
   ^     - find alignment end and start loop
   ^   o fun05 sec04 sub02:
   ^     - insertion cases
   ^   o fun05 sec04 sub03:
   ^     - snp or match cases
   ^   o fun05 sec04 sub04:
   ^     - deletion cases
   ^   o fun05 sec04 sub05:
   ^     - note ending referencde position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec04 Sub01:
   *   - find alignment end and start loop
   \*****************************************************/

   dirMatrixSC += indexUL;

   while(*dirMatrixSC != def_mvStop_alnDefs)
   { /*Loop: trace alignment path*/

      switch(*dirMatrixSC)
      { /*Switch: find direction (snp/del/ins)*/
         case def_mvStop_alnDefs:
            break;

         /***********************************************\
         * Fun05 Sec04 Sub02:
         *   - insertion cases
         \***********************************************/

         case def_mvIns_alnDefs:
         /*Case: is a insertion*/
            if(
                  samSTPtr->cigTypeStr[samSTPtr->lenCigUI]
               != 'I'
            ){ /*If: need to make new cigar entry*/
               if(
                  samSTPtr->cigTypeStr[samSTPtr->lenCigUI]
               ) ++samSTPtr->lenCigUI; /*if not null*/

               samSTPtr->cigTypeStr[
                  samSTPtr->lenCigUI
               ] = 'I';

               samSTPtr->cigValAryI[
                  samSTPtr->lenCigUI
               ] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++samSTPtr->cigValAryI[samSTPtr->lenCigUI];

            ++samSTPtr->numInsUI;

            --qryPosUI;
            dirMatrixSC -= (lenRefUI + 1);
            break;
         /*Case: is a insertion*/

         /***********************************************\
         * Fun05 Sec04 Sub03:
         *   - snp or match cases
         \***********************************************/

         case def_mvSnp_alnDefs:
         /*Case: is a snp or match*/
            matchBl =
               getMatch_alnSet(
                  qrySeqStr[qryPosUI],
                  refSeqStr[refPosUI],
                  alnSetSTPtr
               ); /*find if had a match*/

            if(matchBl & def_ntEql_alnDefs)
            { /*If: had a match*/
               matchBl = '='; /*match*/
               ++samSTPtr->numMatchUI;

               /*count number of anonymous matches*/
               *numAnonUI +=
                  (matchBl & def_anonymous_alnDefs);
            } /*If: had a match*/

            else
            { /*Else: had a mismatch*/
               matchBl = 'X'; /*mismatch*/
               ++samSTPtr->numSnpUI;
            } /*Else: had a mismatch*/

            if(
                  samSTPtr->cigTypeStr[samSTPtr->lenCigUI]
               != matchBl
            ){ /*If: need to make new cigar entry*/
               if(
                  samSTPtr->cigTypeStr[samSTPtr->lenCigUI]
               ) ++samSTPtr->lenCigUI; /*if not null*/

               samSTPtr->cigTypeStr[
                  samSTPtr->lenCigUI
               ] = matchBl;

               samSTPtr->cigValAryI[
                  samSTPtr->lenCigUI
               ] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++samSTPtr->cigValAryI[samSTPtr->lenCigUI];

            --refPosUI;
            --qryPosUI;

            ++samSTPtr->alnReadLenUI;

            dirMatrixSC -= (lenRefUI + 2);
            break;
         /*Case: is a snp or match*/

         /***********************************************\
         * Fun05 Sec04 Sub04:
         *   - deletion cases
         \***********************************************/

         case def_mvDel_alnDefs:
         /*Case: is a deletion*/
            if(
                  samSTPtr->cigTypeStr[samSTPtr->lenCigUI]
               != 'D'
            ){ /*If: need to make new cigar entry*/
               if(
                  samSTPtr->cigTypeStr[samSTPtr->lenCigUI]
               ) ++samSTPtr->lenCigUI; /*if not null*/

               samSTPtr->cigTypeStr[
                  samSTPtr->lenCigUI
               ] = 'D';

               samSTPtr->cigValAryI[
                  samSTPtr->lenCigUI
               ] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++samSTPtr->cigValAryI[samSTPtr->lenCigUI];

            ++samSTPtr->numDelUI;
            ++samSTPtr->alnReadLenUI;

            --refPosUI;
            --dirMatrixSC;
            break;
         /*Case: is a deletion*/
      } /*Switch: find direction (snp/del/ins)*/
   } /*Loop: trace alignment path*/

   /*****************************************************\
   * Fun05 Sec04 Sub06:
   *   - note ending reference position
   \*****************************************************/

   /*account for ending on one index behind*/
   ++refPosUI;
   ++qryPosUI;

   if(revBl)
      samSTPtr->refEndUI = lenRefUI - refPosUI - 1;
   else
      samSTPtr->refStartUI = refPosUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec05:
   ^   - add starting softmasked bases and invert cigar
   ^   o fun05 sec05 sub01:
   ^     - add starting softmasked bases to cigar
   ^   o fun05 sec05 sub02:
   ^     - if forward sequence invert cigar; is backwards
   ^   o fun05 sec05 sub03:
   ^     - if reverse sequence; reverse complement
   ^   o fun05 sec05 sub04:
   ^     - make sure cigar length is index 1 (not 0)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec05 Sub01:
   *   - add starting softmasked bases to cigar
   \*****************************************************/

   if(qryPosUI > 0)
   { /*If: need to add starting softmasked bases*/
      ++samSTPtr->lenCigUI;
      samSTPtr->cigTypeStr[samSTPtr->lenCigUI] = 'S';

      samSTPtr->cigValAryI[samSTPtr->lenCigUI] = qryPosUI;
      samSTPtr->numMaskUI += qryPosUI;
   } /*If: need to add starting softmasked bases*/

   /*****************************************************\
   * Fun05 Sec05 Sub02:
   *   - if forward sequence invert cigar; is backwards
   \*****************************************************/

   if(! revBl)
   { /*If: need to reverse (ivert) cigar*/
      qryPosUI = 0;
      refPosUI = samSTPtr->lenCigUI;

      while(qryPosUI < refPosUI)
      { /*Loop: invert the cigar to forward direction*/
         /*swap was from standfords bithacking guide*/

         /*this allows ref ^ qry to give qry value*/
         samSTPtr->cigTypeStr[qryPosUI] ^=
            samSTPtr->cigTypeStr[refPosUI];

         samSTPtr->cigValAryI[qryPosUI] ^=
            samSTPtr->cigValAryI[refPosUI];

         /*this (ref ^ qry) sets ref position to qry value
         `  it also allows qry ^ ref to give rev
         */
         samSTPtr->cigTypeStr[refPosUI] ^=
            samSTPtr->cigTypeStr[qryPosUI];

         samSTPtr->cigValAryI[refPosUI] ^=
            samSTPtr->cigValAryI[qryPosUI];

         /*this (qry ^ ref) sets qry position to ref*/
         samSTPtr->cigTypeStr[qryPosUI] ^=
            samSTPtr->cigTypeStr[refPosUI];

         samSTPtr->cigValAryI[qryPosUI] ^=
            samSTPtr->cigValAryI[refPosUI];

         ++qryPosUI;
         --refPosUI;
      } /*Loop: invert the cigar to forward direction*/

      /*null/0 or softmask at start flipped to end*/
   } /*If: need to reverse (ivert) cigar*/

   /*****************************************************\
   * Fun05 Sec05 Sub03:
   *   - if reverse sequence; reverse complement
   \*****************************************************/

   else
   { /*Else: reverse complement sequence*/
      seqDoNotFreeST.seqStr = samSTPtr->seqStr;
      seqDoNotFreeST.lenSeqUL = samSTPtr->readLenUI;

      if(
            samSTPtr->qStr[0] != '*'
         && samSTPtr->qStr[1] != '\0'
         && samSTPtr->qStr[0] != '\0'
      ){ /*If: have q-score entry*/
         seqDoNotFreeST.qStr = samSTPtr->qStr;
         seqDoNotFreeST.lenQUL = samSTPtr->readLenUI;
      } /*If: have q-score entry*/

      else
      { /*Else: do not have q-score entry*/
         seqDoNotFreeST.qStr = 0;
         seqDoNotFreeST.lenQUL = 0;
      } /*Else: do not have q-score entry*/

      revComp_seqST(&seqDoNotFreeST);

      seqDoNotFreeST.seqStr = 0;
      seqDoNotFreeST.lenSeqUL = 0;

      seqDoNotFreeST.qStr = 0;
      seqDoNotFreeST.lenQUL = 0;
   } /*Else: reverse complement sequence*/


   /*****************************************************\
   * Fun05 Sec05 Sub04:
   *   - make sure cigar length is index 1 (not 0)
   \*****************************************************/

   if(samSTPtr->cigTypeStr[samSTPtr->lenCigUI] != '\0')
   { /*If: is index 0*/
      ++samSTPtr->lenCigUI; /*convert to index 1*/
      samSTPtr->cigTypeStr[samSTPtr->lenCigUI] = '\0';
   } /*If: is index 0*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec06:
   ^   - add tags (NM, AS, nn)
   ^   o fun05 sec06 sub01:
   ^     - add NM (edit distance/number differences) flag
   ^   o fun05 sec06 sub02:
   ^     - add AS (score) flag
   ^   o fun05 sec06 sub03:
   ^     - add nn (anonymoys bases) flag
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec06 Sub01:
   *   - add NM (edit distance/number differences) flag
   \*****************************************************/

   tmpStr = (schar *) samSTPtr->extraStr;

   samSTPtr->lenExtraUI +=
      cpDelim_charCp(
         (char *) &tmpStr[samSTPtr->lenExtraUI],
         "NM:i:",
         '\0'
      ); /*copy edit distance*/

   samSTPtr->lenExtraUI +=
      numToStr(
         (char *) &tmpStr[samSTPtr->lenExtraUI],
         samSTPtr->numSnpUI
            + samSTPtr->numDelUI
            + samSTPtr->numInsUI
      ); /*find edit distance*/

   /*****************************************************\
   * Fun05 Sec06 Sub02:
   *   - add AS (score) flag
   \*****************************************************/

   tmpStr[samSTPtr->lenExtraUI++] = '\t';

   samSTPtr->lenExtraUI +=
      cpDelim_charCp(
         (char *) &tmpStr[samSTPtr->lenExtraUI],
         "AS:i:",
          '\0'
      ); /*copy assembler score*/

   samSTPtr->lenExtraUI +=
      numToStr(
         (char *) &tmpStr[samSTPtr->lenExtraUI],
         matrixSTPtr->scoreSL /  def_scoreAdj_alnDefs
      ); /*get alignment score*/

   /*****************************************************\
   * Fun05 Sec06 Sub03:
   *   - add nn (anonymoys bases) flag
   \*****************************************************/

   tmpStr[samSTPtr->lenExtraUI++] = '\t';

   samSTPtr->lenExtraUI +=
      cpDelim_charCp(
         (char *) &tmpStr[samSTPtr->lenExtraUI],
         "nn:i:",
         '\0'
      ); /*copy number anonymous bases*/

   samSTPtr->lenExtraUI +=
      numToStr(
         (char *) &tmpStr[samSTPtr->lenExtraUI],
         *numAnonUI
      ); /*get number anonymous bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec07:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun05_sec07:;

   return 1;
} /*getAln_dirMatrix*/

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
