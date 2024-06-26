/*########################################################
# Name: primFind
#   - finds sequences with input primers
########################################################*/

/*-------------------------------------------------------\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o main:
'     - driver function
'   o license:
'     - Licensing for this code (public domain / mit)
\-------------------------------------------------------*/

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

#include "kmerFind.h"
#include "memwater/seqST.h"
#include "memwater/alnSetST.h"
#include "inputPrimFind.h"

/*.h files only*/
#include "generalLib/dataTypeShortHand.h"
#include "generalLib/charCp.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c   #include "memwater/memwater.h"
!   - .h   #include "generalLib/ulCp.h"
!   - .h   #include "generalLib/genMath.h"
!   - .h   #include "generalLib/shellSort.h"
!   - .h   #include "generalLib/base10str.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Main:
|   - main driver function for primFind
| Input:
|   - numArgsSI:
|     o number of arguments user input in argAryStr
|   - argAryStr:
|     o array of c-strings with user arguments
| Output:
|   - Prints:
|     o
\-------------------------------------------------------*/
#ifdef PLAN9
void
#else
int
#endif
main(
   int numArgsSI,
   char *argAryStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '   o main sec01:
   '     - variable declerations
   ^   - get user input and initialize structures
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - variable declerations
   ^   o main sec01 sub01:
   ^     - output/input variables (one tempory)
   ^   o main sec01 sub02:
   ^     - variables for filtering or storing results
   ^   o main sec01 sub03:
   ^     - variables unique to kmer search
   ^   o main sec01 sub04:
   ^     - structure and file variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec01 Sub01:
   *   - output/input variables (one tempory)
   \*****************************************************/

   schar errSC = 0;      /*returned error*/

   schar *seqFileStr = 0; /*fasta/q (x) with sequences*/
   schar fqBl = 1;       /*1 is a fastq file; 0 is fasta*/

   schar *primFileStr = 0; /*file with primer sequences*/
   schar primFileTypeSC = def_tsvPrimFile_inputPrimFind;
      /*file type with primer sequences*/

   schar *outFileStr = 0; /*output file*/

   ulong entryOnUL = 0;
   schar *tmpStr = 0;

   /*****************************************************\
   * Main Sec01 Sub02:
   *   - variables for filtering or storing results
   \*****************************************************/

   /*variables for filtering*/
   uint minLenUI = def_minLen_inputPrimFind;
   uint maxLenUI = def_maxLen_inputPrimFind;

   /*variables for holding search output*/
   uint *codeHeapAryUI = 0; /*matchs/primer*/

   slong *scoreHeapArySL = 0;
   schar *dirHeapArySC = 0;

   ulong *seqStartHeapAryUL = 0;
   ulong *seqEndHeapAryUL = 0;

   ulong *primStartHeapAryUL = 0;
   ulong *primEndHeapAryUL = 0;

   /*****************************************************\
   * Main Sec01 Sub03:
   *   - variables unique to kmer search
   \*****************************************************/

   schar fastBl = def_fastSearch_inputPrimFind;
      /*1 = kmer search, 0 = waterman*/

   float minPercScoreF = def_minPercScore_kmerFind;
      /*min score to keep mapping*/

   uchar lenKmerUC = def_lenKmer_kmerFind; /*kmer length*/

   float frameShiftF = def_percShift_kmerFind;
      /*percentage of bases to move when shifting a win*/

   float minPercKmerF = def_minKmerPerc_kmerFind;
      /*minimum percentage of kmers to do a waterman*/

   float extraNtInWinF = def_extraNtInWin_kmerFind;
     /*number of extra (beyond primer) bases in window*/

   /*****************************************************\
   * Main Sec01 Sub04:
   *   - structure and file variables
   \*****************************************************/

   /*structures for primer search*/
   struct refST_kmerFind *refHeapAryST = 0;
   sint numRefsSI = 0; /*references in refHeapAryST*/

   struct tblST_kmerFind kmerTblStackST;
   struct seqStruct seqStackST;
   struct alnSet alnSetStackST;

   FILE *seqFILE = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - get user input and initialize structures
   ^   o main sec02 sub01:
   ^     - initialize structures
   ^   o main sec02 sub02:
   ^     - get input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize structures
   \*****************************************************/

   init_tblST_kmerFind(
      &kmerTblStackST,
      lenKmerUC
   );

   init_seqST(&seqStackST);
   init_alnSetST(&alnSetStackST);

   alnSetStackST.gapOpenC = -4;

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get input
   \*****************************************************/

   errSC =
      getInput_inputPrimFind(
         numArgsSI,
         argAryStr,
         &primFileStr,
         &primFileTypeSC,
         &seqFileStr,
         &fqBl,
         &minLenUI,
         &maxLenUI,
         &minPercScoreF,
         &fastBl,
         &lenKmerUC,
         &minPercKmerF
      );

   if(errSC)
   { /*If: had error*/
      --errSC; /*remove help message error*/
      errSC <<= 5; /*avoid conflicts with other errors*/
      goto cleanUp_main_sec05_sub04;
   } /*If: had error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - open files
   *   o main sec02 sub03 cat01:
   *     - open output file
   *   o main sec02 sub03 cat02:
   *     - open read file
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat01:
   +   - open output file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(!outFileStr || *outFileStr == '-')
      outFILE = stdout;

   else
   { /*Else: user supplied an output file*/
      outFILE =
         fopen(
            (char *) outFileStr,
            "w"
         );

      if(! outFILE)
      { /*If: could not open output file*/
         fprintf(
            stderr,
            "could not open -out %s\n",
            outFileStr
         );
      } /*If: could not open output file*/

      goto fileErr_main_sec05_sub03;
   } /*Else: user supplied an output file*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat02:
   +   - open read file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(! seqFileStr || *seqFileStr == '-')
      seqFILE = stdin;

   else
   { /*Else: user supplied an read file*/
      seqFILE =
         fopen(
            (char *) seqFileStr,
            "r"
         );

      if(! seqFILE)
      { /*If: could not open readFxput file*/
         if(fqBl)
            tmpStr = (schar *) "-fq";
         else
            tmpStr = (schar *) "-fa";

         fprintf(
            stderr,
            "could not open %s %s\n",
            tmpStr,
            seqFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: could not open readFxput file*/
   } /*Else: user supplied an read file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - read in references and allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub02:
   *   - read in references
   \*****************************************************/

   if(primFileTypeSC == def_tsvPrimFile_inputPrimFind)
   { /*If: input was from tsv file*/
      refHeapAryST =
         tsvToAry_refST_kmerFind(
            primFileStr,
            lenKmerUC,
            &numRefsSI,
            minPercKmerF,
            &kmerTblStackST,
            extraNtInWinF,
            frameShiftF,
            &alnSetStackST,
            (uchar *) &errSC
         );
   } /*If: input was from tsv file*/

   else
   { /*Else: input was from fasta file*/
      refHeapAryST =
         faToAry_refST_kmerFind(
            primFileStr,
            lenKmerUC,
            &numRefsSI,
            minPercKmerF,
            &kmerTblStackST,
            extraNtInWinF,
            frameShiftF,
            &alnSetStackST,
            (uchar *) &errSC
         );
   } /*Else: input was from fasta file*/

   if(errSC)
   { /*If: I had an error*/
      if(errSC == def_memErr_kmerFind)
      { /*If: I had a memory error*/
         fprintf(
            stderr,
            "Memory error when reading primer sequences\n"
         );

         goto memErr_main_sec05_sub02;
      } /*If: I had a memory error*/

      if(errSC == def_fileErr_kmerFind)
      { /*If: had a file error*/
         fprintf(
            stderr,
            "Problem reading -prim %s\n",
            primFileStr
         );
      } /*If: had a file error*/

      goto fileErr_main_sec05_sub03;
   } /*If: I had an error*/

   /*****************************************************\
   * Main Sec03 Sub02:
   *   - allocate arrays for reference results
   \*****************************************************/

   dirHeapArySC = malloc(numRefsSI * sizeof(schar));

   if(!dirHeapArySC)
      goto memErr_main_sec05_sub02;
  
   codeHeapAryUI = malloc(numRefsSI * sizeof(uint));

   if(! codeHeapAryUI)
      goto memErr_main_sec05_sub02;

   scoreHeapArySL = malloc(numRefsSI * sizeof(slong));

   if(! scoreHeapArySL)
      goto memErr_main_sec05_sub02;

   seqStartHeapAryUL = malloc(numRefsSI * sizeof(ulong));

   if(! seqStartHeapAryUL)
      goto memErr_main_sec05_sub02;

   seqEndHeapAryUL = malloc(numRefsSI * sizeof(ulong));

   if(! seqEndHeapAryUL)
      goto memErr_main_sec05_sub02;

   primStartHeapAryUL = malloc(numRefsSI * sizeof(ulong));

   if(! primStartHeapAryUL)
      goto memErr_main_sec05_sub02;

   primEndHeapAryUL = malloc(numRefsSI * sizeof(ulong));

   if(! primEndHeapAryUL)
      goto memErr_main_sec05_sub02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - find and print primer positions
   ^   o main sec04 sub01:
   ^     - read in first sequence and start loop
   ^   o main sec04 sub02:
   ^     - check if sequence passes filters
   ^   o main sec04 sub03:
   ^     - find and print primer positions
   ^   o main sec04 sub04:
   ^     - get the next sequence
   ^   o main sec04 sub05:
   ^     - check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - read in first sequence and start loop
   \*****************************************************/

   pHeaderHit_kmerFind(outFILE);

   if(fqBl)
   { /*If: i have a fastq file*/
      errSC =
         (schar)
         getFqSeq_seqST(
            seqFILE,
            &seqStackST
         );
   } /*If: i have a fastq file*/

   else
   { /*Else: i have a fasta file*/
      errSC =
         (schar)
         getFaSeq_seqST(
            seqFILE,
            &seqStackST
         );
   } /*Else: i have a fasta file*/

   while(! errSC)
   { /*Loop: find primers in sequences*/

      /**************************************************\
      * Main Sec04 Sub02:
      *   - check if sequence passes filters
      \**************************************************/

      ++entryOnUL;

      if(seqStackST.lenSeqUL < minLenUI)
         goto nextSeq_main_sec04_sub04;

      if(seqStackST.lenSeqUL > maxLenUI)
         goto nextSeq_main_sec04_sub04;

      /**************************************************\
      * Main Sec04 Sub03:
      *   - find and print primer positions
      \**************************************************/

      if(fastBl)
      { /*If: i am searching with kmers*/
         errSC =
            fxFindPrims_kmerFind(
               &kmerTblStackST,    /*table structure*/
               refHeapAryST,       /*primers looking for*/
               numRefsSI,          /*number of primers*/
               &seqStackST,        /*sequence to look at*/
               minPercScoreF,      /*min % for waterman*/
               codeHeapAryUI,      /*# times found prim*/
               dirHeapArySC,       /*direction of primer*/
               scoreHeapArySL,     /*scores for prims*/
               seqStartHeapAryUL,  /*1st align seq base*/
               seqEndHeapAryUL,    /*last align seq base*/
               primStartHeapAryUL, /*1st align prim base*/
               primEndHeapAryUL,   /*last align prim bas*/
               &alnSetStackST      /*settings*/
            ); 
      } /*If: i am searching with kmers*/

      else
      { /*Else: slower waterman search*/
         errSC =
            waterFindPrims_kmerFind(
               refHeapAryST,       /*primers looking for*/
               numRefsSI,          /*number of primers*/
               &seqStackST,        /*sequence to look at*/
               minPercScoreF,      /*min % for waterman*/
               codeHeapAryUI,      /*# times found prim*/
               dirHeapArySC,       /*direction of primer*/
               scoreHeapArySL,     /*scores for prims*/
               seqStartHeapAryUL,  /*1st align seq base*/
               seqEndHeapAryUL,    /*last align seq base*/
               primStartHeapAryUL, /*1st align prim base*/
               primEndHeapAryUL,   /*last align prim bas*/
               &alnSetStackST      /*settings*/
            ); 
      } /*Else: slower waterman search*/

      if(! errSC)
      { /*If: I found at least on primer*/
         phit_kmerFind(
            refHeapAryST,       /*primers looking for*/
            numRefsSI,          /*number of primers*/
            &seqStackST,        /*sequence to look at*/
            codeHeapAryUI,      /*# times found a primer*/
            dirHeapArySC,       /*direction of primer*/
            scoreHeapArySL,     /*scores for each primer*/
            seqStartHeapAryUL,  /*first aligned seq base*/
            seqEndHeapAryUL,    /*last aligned seq base*/
            primStartHeapAryUL, /*1st aligned prim base*/
            primEndHeapAryUL,   /*last alinged prim base*/
            outFILE
         );
      } /*If: I found at least on primer*/

      /**************************************************\
      * Main Sec04 Sub04:
      *   - get the next sequence
      \**************************************************/

      nextSeq_main_sec04_sub04:;

      if(fqBl)
      { /*If: i have a fastq file*/
         errSC =
            (schar)
            getFqSeq_seqST(
               seqFILE,
               &seqStackST
            );
      } /*If: i have a fastq file*/

      else
      { /*Else: i have a fasta file*/
         errSC =
            (schar)
            getFaSeq_seqST(
               seqFILE,
               &seqStackST
            );
      } /*Else: i have a fasta file*/
   } /*Loop: find primers in sequences*/

   /*****************************************************\
   * Main Sec04 Sub05:
   *   - check for errors
   \*****************************************************/

   if(errSC != def_EOF_seqST)
   { /*If: i had an error of some kind*/
      if(fqBl)
         tmpStr = (schar *) "-fq";
      else
         tmpStr = (schar *) "-fa";

      if(errSC & def_memErr_seqST)
      { /*If: i had a memory error*/
         fprintf(
            stderr,
            "Memory error when %s %s\n",
            tmpStr,
            seqFileStr
         );

         goto memErr_main_sec05_sub02;
      } /*If: i had a memory error*/

      fprintf(
         stderr,
         "Error on entry %lu of %s %s\n",
         entryOnUL,
         tmpStr,
         seqFileStr
      );

      goto fileErr_main_sec05_sub03;
   } /*If: i had an error of some kind*/

   else
      errSC = 0; /*no errors*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05:
   ^   - clean up
   ^   o main sec05 sub01:
   ^     - clean up after no errors
   ^   o main sec05 sub02:
   ^     - deal with memory errors
   ^   o main sec05 sub03:
   ^     - deal with file errors
   ^   o main sec05 sub04:
   ^     - clean up and exit
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec05 Sub01:
   *   - clean up after no errors
   \*****************************************************/

   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub02:
   *   - deal with memory errors
   \*****************************************************/

   memErr_main_sec05_sub02:;

   errSC = def_memErr_kmerFind;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub03:
   *   - deal with file errors
   \*****************************************************/

   fileErr_main_sec05_sub03:;

   errSC = def_fileErr_kmerFind;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub04:
   *   - clean up and exit
   \*****************************************************/

   cleanUp_main_sec05_sub04:;

   /*waterman smith alignment arrays*/
   free(dirHeapArySC);
   dirHeapArySC = 0;

   free(codeHeapAryUI);
   codeHeapAryUI = 0;

   free(scoreHeapArySL);
   scoreHeapArySL = 0;

   free(seqStartHeapAryUL);
   seqStartHeapAryUL = 0;

   free(seqEndHeapAryUL);
   seqEndHeapAryUL = 0;

   free(primStartHeapAryUL);
   primStartHeapAryUL = 0;

   free(primEndHeapAryUL);
   primEndHeapAryUL = 0;

   /*alignment settings*/
   freeStack_alnSetST(&alnSetStackST);

   /*kmer find structures*/
   freeStack_seqST(&seqStackST);
   freeStack_tblST_kmerFind(&kmerTblStackST);

   freeHeapAry_refST_kmerFind(
      refHeapAryST,
      numRefsSI
   );

   refHeapAryST = 0;

   if(outFILE && outFILE != stdout)
      fclose(outFILE);

   outFILE = 0;

   if(seqFILE && seqFILE != stdin)
      fclose(seqFILE);

   seqFILE = 0;

   exit(errSC);
} /*main*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconveint / not possible, this code is under the
:   MIT license
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
