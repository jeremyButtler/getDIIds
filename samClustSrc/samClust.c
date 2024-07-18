#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

/*.h files only*/
#include "../generalLib/dataTypeShortHand.h"
#include "../generalLib/ulCp.h"

/*-------------------------------------------------------\
| Fun01: getNumReads_samClust
|   - gets number of reads in the sam file
| Input:
|   - samFILE:
|     o sam file to get read counts for
| Output:
|   - Returns:
|     o number of reads in sam file
\-------------------------------------------------------*/
unsigned int
getNumReads_samClust(
   void *samFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   o fun01 sec01:
   '     - variable declarations
   '   o fun01 sec02:
   '     - get number of reads in sam file
   '   o fun01 sec03:
   '     - return count
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ushort lenBuffUS = 4096;  /*size of buffer*/
   schar buffStr[lenBuffUS]; /*buffer for read counting*/
   uint posUI = 0;       /*position at in buffer*/
   ulong bytesUL = 0;    /*number of bytes in buffer*/
   uint readsUI = 0;     /*number of reads in file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - get number of reads in sam file
   ^   o fun01 sec02 sub01:
   ^     - read in first chunk of file
   ^   o fun01 sec02 sub02:
   ^     - find number of lines in sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *   - read in first chunk of file
   \*****************************************************/

   bytesUL =
      fread(
         (char *) buffStr,
         sizeof(schar),
         lenBuffUS - 1,
         (FILE *) samFILE
      );

   buffStr[bytesUL] = '\0';
   posUI = 0;

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *   - find number of lines in sam file
   \*****************************************************/

   while(bytesUL)
   { /*Loop: find number of reads*/
      if(buffStr[posUI] != '@')   /*true most of time*/
         ++readsUI; /*count number of reads*/

       posUI +=
          endLine_ulCp((char *) &buffStr[posUI]);

       while(buffStr[posUI] == '\0')
       { /*Loop: find next new line*/
          bytesUL =
             fread(
                (char *) buffStr,
                sizeof(schar),
                lenBuffUS - 1,
                (FILE *) samFILE
             ); /*at null, so need to read in more file*/

          if(bytesUL == 0)
             break; /*hit end of file*/

          buffStr[bytesUL] = '\0';
          posUI = 0;

          posUI += 
             endLine_ulCp((char *) &buffStr[posUI]);

          if(bytesUL == 0)
             break; /*at end of file*/
       } /*Loop: find next new line*/

       ++posUI; /*move off new line*/

       if(buffStr[posUI] == '\0')
       { /*If: new line was last character in buffer*/
          bytesUL =
             fread(
                (char *) buffStr,
                sizeof(schar),
                lenBuffUS - 1,
                (FILE *) samFILE
             );

          buffStr[bytesUL] = '\0';
          posUI = 0; /*is true until EOF*/
       } /*If: new line was last character in buffer*/
   } /*Loop: find number of reads*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - return count
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fseek(samFILE, 0, SEEK_SET);
   return readsUI;
} /*getNumReads_samClust*/

/*-------------------------------------------------------\
| Fun02: getDist_samClust
|   - gets edit distances between the query and reference
|   - deletions and insertions are only counted if they
|     execed a minimum length.
|   - length of deletions and insertions are reduced to
|     minimum length for adding. the idea is to put
|     priority on snps.
| Input:
|   - qrySTPtr:
|     o pointer to samEntry structure with read (query) to
|       find the edit distance for
|   - refSTPtr:
|     o pointer to samEntry structure with reference to
|       compare query (qrySTPtr) to
|   - numReadsUI:
|     o number of reads in the sam file
|   - delLenUI:
|     o minimum length for a deletion to count as an event
|     o for edit distanace all kept deletions are reduced
|       to this
|   - insLenUI:
|     o minimum length for a insertion to count as a event
|     o for edit distanace all kept insertions are reduced
|       to this
|   - minQUC:
|     o minimum q-score to keep an snp or count insertion
| Output:
|   - Returns:
|     o edit distance between query and ref
|     o -1 if reads to not overlap
\-------------------------------------------------------*/
signed int
getDist_samClust(
   struct samEntry *qrySTPtr, /*read for edit distance*/
   struct samEntry *refSTPtr, /*ref to compare*/
   samReadSTPtr,            /*reads to get edit distance*/
   unsigned int numReadsUI, /*number of reads in samFILE*/
   unsigned int delLenUI,   /*min deletion length*/
   unsigned int insLenUI,   /*min insertion length*/
   unsigned int maxIndeUI,  /*reduce all indels to this*/
   unsigned char minQUC     /*min Q-score*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - gets edit distances between the query & reference
   '   o fun02 sec01:
   '     - variable declerations
   '   o fun02 sec02:
   '     - find start of overlap between query & reference
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint distSI = 0;    /*edit distance*/

   sint siQryCig = 0;  /*query cigar position*/
   sint siRefCig = 0;  /*reference cigar position*/
   sint refValSI = 0;
   sint qryValSI = 0;

   sint tmpSI = 0;

   uint uiQry = 0;     /*query nucleotide on*/
   uint uiRef = 0;     /*reference nucleotide on*/

   uint refKeptUI = 0;  /*kept reference snps/ins*/
   uint qryKeptUI = 0;  /*kept query snp/ins*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - find start of overlap between query & reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(qrySTPtr->refEndUI < refSTPtr->refStartUI)
      goto noOverlap_fun02_sec0x;

   if(qrySTPtr->refStartUI > refSTPtr->refEndUI)
      goto noOverlap_fun02_sec0x;

   refValSI = refSTPtr->cigValAryI[0];
   qryValSI = refSTPtr->cigValAryI[0];

   if(qrySTPtr->refStartUI > refSTPtr->refStartUI)
   { /*If: move the reference foward*/
      uiQry = refSTPtr->refStartUI;

      findRefPos_samEntry(
         refSTPtr,
         siRefCig,
         refValSI,
         qrySTPtr->refStartUI,
         uiQry,                /*reference pos (discard)*/
         siRefNt
      ); /*set reference to first query base*/

      uiQry = 0;
   } /*If: move reference foward*/

   else if(qrySTPtr->refStartUI < refSTPtr->refStartUI)
   { /*Else If: move query foward*/
      uiRef = refSTPtr->refStartUI;

      findRefPos_samEntry(
         qrySTPtr,
         siQryCig,
         qryValSI,
         refSTPtr->refStartUI, /*end query position*/
         uiRef,                /*query pos (discard)*/
         siQryNt               /*end query position*/
      ); /*set reference to first query base*/

      uiRef = 0;
   } /*Else If: move query foward*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - find edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(refSTPtr->cigTypeStr[siRefCig] == 'S')
   { /*If: have softmasked region*/
      uiRef += refValSI;
      ++siRefCig;
      refValSI = refSTPtr->cigValAryI[siRefCig];
   } /*If: have softmasked region*/

   if(qrySTPtr->cigTypeStr[siQryCig] == 'S')
   { /*If: have softmasked region*/
      uiQry += qryValSI;
      ++siQryCig;
      qryValSI = qrySTPtr->cigValAryI[siQryCig];
   } /*If: have softmasked region*/

   if(refSTPtr->cigTypeStr[siRefCig] == 'H')
   { /*If: have hardmasked region*/
      ++siRefCig;
      refValSI = refSTPtr->cigValAryI[siRefCig];
   } /*If: have hardmasked region*/

   if(qrySTPtr->cigTypeStr[siQryCig] == 'H')
   { /*If: have hardmaksed region*/
      ++siQryCig;
      qryValSI = qrySTPtr->cigValAryI[siQryCig];
   } /*If: have hardmasked region*/

   while(
         siQryCig < (sint) qrySTPtr->lenCigUI
      && siRefCig < (sint) qrySTPtr->lenCigUI
   ){ /*Loop: get edit distance*/
      if(refSTPtr->cigTypeStr[siRefCig == 'S')
         break; /*soft masking only at ends (finished)*/

      if(qrySTPtr->cigTypeStr[siRefCig == 'S')
         break; /*soft masking only at ends (finished)*/

      if(refSTPtr->cigTypeStr[siRefCig == 'H')
         break; /*hard masking only at ends (finished)*/

      if(qrySTPtr->cigTypeStr[siRefCig == 'H')
         break; /*hard masking only at ends (finished)*/

      /*Values:
      ` =: 61 (match)
      ` X: 88 (snp)
      ` M: 77 (match/snp)
      ` D: 68 (deletion)
      ` I: 73 (insertion)
      ` M + M: 154
      ` = + =: 122
      ` X + X: 176 (need to check if snps agree/disagree)
      ` I + I: 146 (check size; rm low q-scores)
      ` D + D: 136 (check size)
      `
      ` X + M: 165 (need to check if agree/disagree)
      ` = + M: 138 (need to check if agree/disagree)
      ` = + X: 149 (disagree, check q-scores)
      `
      ` D + M: 145 (disagree, check size)
      ` D + =: 129 (disagree, check size)
      ` D + X: 156 (disagree, check siz)
      `
      ` I + M: 150
      ` I + =: 134
      ` I + X: 161
      ` I + D: 141
      */

      switch(refSTPtr->cigTypeStr[siRefCig])
      { /*Switch: find mutation combination*/
         /*cases were I have to compare sequences*/
         case 154: /*M/M*/
         case 176: /*X/X*/
         case 165: /*X/M*/
         case 138: /*=/M*/
         /*Case: M / (M, X, =) count differences*/
            break;
         /*Case: M / (M, X, =) count differences*/

         /*need to check q-scores*/
         case 138: /*=/X*/
         /*Case: snp and match, find number pass diff*/
            break;
         /*Case: snp and match, find number pass diff*/

         /*cases I treat as same*/
         case 122: /*=/=*/
         case 146: /*I/I*/
         case 136: /*D/D*/
         /*Case: treat as same*/
            tmpSI =
               min_genMath(
                  refValSI,
                  qryValSI
               );

            refValSI -= tmpSI;
            qryValSI -= tmpSI;
            break;
         /*Case: treat as same*/
 
         /*deletion cases*/
         case 145: /*D/M*/
         case 129: /*D/=*/
         case 156: /*D/X*/
         /*Case: deletion (only one sequence)*/
            if(refSTPtr->cigTypeStr[siRefCig] == 'D')
               tmpSI = refSTPtr->cigValAryI[siRefCig];
            else
               tmpSI = qrySTPtr->cigValAryI[siRefCig];

            if(tmpSI >= delLenUI)
            { /*If: keeping deletion (>= min size)*/
                  min_genMath(
                     tmpSI,
                     (sint) maxIndelUI
                  );

               /*reducing deletion to min length to reduce
               `   impact on edit distance
               */
            } /*If: keeping deletion (>= min size)*/


            while(tmpSI > 0)
            { /*Loop: move past deletion entry*/
               --refValSI;
               --qryValSI;

               if(refValSI == 0)
               { /*If: on next cigar entry*/
                  ++siRefCig;
                  refValSI=refSTPtr->cigValAryI[siRefCig];
               } /*If: on next cigar entry*/

               if(qryValSI == 0)
               { /*If: on next cigar entry*/
                  ++siQryCig;
                  qryValSI=qrySTPtr->cigValAryI[siQryCig];
               } /*If: on next cigar entry*/
            } /*Loop: move past deletion entry*/

            refValSI -= tmpSI;
            qryValSI -= tmpSI;
            break;
         /*Case: deletion (only one sequence)*/

         /*insertion cases*/
         case 150: /*I/M*/
         case 134: /*I/=*/
         case 161: /*I/X*/
         case 141: /*I/D*/

         case 'X':
         case 'D':
      } /*Switch: find mutation combination*/

      if(refValSI == 0)
      { /*If: on next cigar entry*/
         ++siRefCig;
         refValSI = refSTPtr->cigValAryI[siRefCig];
      } /*If: on next cigar entry*/

      if(qryValSI == 0)
      { /*If: on next cigar entry*/
         ++siQryCig;
         qryValSI = qrySTPtr->cigValAryI[siQryCig];
      } /*If: on next cigar entry*/
   } /*Loop: get edit distance*/

   return distSI;

   noOverlap_fun02_sec0x:;

   return -1;
} /*getDist_samClust*/
