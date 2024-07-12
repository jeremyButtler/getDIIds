/*########################################################
# Name: getDiDis
#   - finds diRNA and mvRNA sequences in reads
########################################################*/

/*-------------------------------------------------------\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o main:
'     - driver function
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

#include "../generalLib/seqST.h"
#include "../generalAln/alnSet.h"

#include "kmerFind.h"
#include "inputGetDIIds.h"
#include "fluST.h"

/*.h files only*/
#include "../generalLib/dataTypeShortHand.h"
#include "../generalLib/charCp.h"
#include "../generalLib/ulCp.h"
#include "../generalLib/genMath.h"

#include "../generalAln/alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../memwater/memwater.h"
!   - .h  #include "../generalAln/memGetCoord.h"
!   - .h  #include "../generalLib/shellSort.h"
!   - .h  #include "../generalLib/base10str.h"
!   - .h  #include "../generalLib/ntToTwoBit.h"
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
   '   o main sec02:
   '     - get user input and initialize structures
   '   o main sec03:
   '     - set up primer sequences and kmer tables
   '   o main sec04:
   '     - find and print primer positions
   '   o main sec05:
   '     - clean up
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

   schar *forSeqStr = forPrimStr_inputGetIds;
   schar *revSeqStr = revPrimStr_inputGetIds;

   /*settings for printing out ids*/
   schar diRnaBl = def_pDIRna_inputGetDIIds;
   schar vRnaBl = def_pVRna_inputGetDIIds;
   schar partBl = def_partSeq_inputGetDIIds;

   schar keepBl = 0;  /*holds if keep di/v RNA read*/
   schar diFlagSC = 0;/*holds result from detectDI_fluST*/

   schar *outFileStr = 0; /*output file*/

   ulong entryOnUL = 0;
   schar *tmpStr = 0;

   /*****************************************************\
   * Main Sec01 Sub02:
   *   - variables for filtering or storing results
   \*****************************************************/

   /*variables for filtering*/
   uint minLenUI = def_minLen_inputGetDIIds;
   uint maxLenUI = def_maxLen_inputGetDIIds;

   /*variables for holding search output*/
   /*usign 3 to avoid overflow errors*/
   uint codeAryUI[3]; /*matchs/primer*/

   slong scoreArySL[3];
   schar dirArySC[3];

   ulong seqStartAryUL[3];
   ulong seqEndAryUL[3];

   ulong primStartAryUL[3];
   ulong primEndAryUL[3];

   float minPercLenF = def_minPercLen_inputGetDIIds;
      /*minimum % of read length between primers*/

   float maxPercLenF = def_maxPercLen_inputGetDIIds;
      /*length between primers / expected segment length*/

   /*variables holding return values*/
   float percLenDiffF = 0; /*primers map len versus read*/

   schar segNumSC = 0; /*holds segment number of read*/

   ulong mapLenUL = 0; /*length between primers*/

   /*****************************************************\
   * Main Sec01 Sub03:
   *   - variables unique to kmer search
   \*****************************************************/

   schar fastBl = def_fastSearch_inputGetDIIds;
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
   struct refST_kmerFind refStackAryST[2];
   unsigned int longestSeqUI = 0;

   struct tblST_kmerFind kmerTblStackST;
   struct seqST seqStackST;
   struct alnSet alnSetStackST;
   struct fluST fluStackST;

   FILE *seqFILE = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - get user input and initialize structures
   ^   o main sec02 sub01:
   ^     - initialize structures
   ^   o main sec02 sub02:
   ^     - get input
   ^   o main sec02 sub03:
   ^     - initialize structures needing user input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize structures
   \*****************************************************/

   init_seqST(&seqStackST);
   init_alnSet(&alnSetStackST);
   errSC = init_fluST(&fluStackST);

   if(errSC)
   { /*If: had a memory error*/
      fprintf(
         stderr,
         "memory error setting up fluST\n"
      );

      goto memErr_main_sec05_sub02;
   } /*If: had a memory error*/

   alnSetStackST.gapSS = -4 * def_scoreAdj_alnDefs;
   alnSetStackST.extendSS = -1 * def_scoreAdj_alnDefs;

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get input
   \*****************************************************/

   errSC =
      getInput_inputGetDIIds(
         numArgsSI,
         argAryStr,
         &fluStackST,
         &seqFileStr,
         &fqBl,
         &forSeqStr,
         &revSeqStr,
         &minLenUI,
         &maxLenUI,
         &minPercLenF,
         &maxPercLenF,
         &minPercScoreF,
         &diRnaBl,
         &vRnaBl,
         &partBl,
         &fastBl,
         &lenKmerUC,
         &minPercKmerF
      );

   if(errSC)
   { /*If: had error*/
      --errSC; /*remove help message error*/
      errSC <<= 5; /*avoid conflicts with other errors*/

      /*I need to initialize these for the free step*/
      init_tblST_kmerFind(
         &kmerTblStackST,
         lenKmerUC
      );

      init_refST_kmerFind(
         &refStackAryST[0],
         lenKmerUC
      ); /*foward primer*/

      init_refST_kmerFind(
         &refStackAryST[1],
         lenKmerUC
      ); /*reverse primer*/

      goto cleanUp_main_sec05_sub04;
   } /*If: had error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - initialize structures needing user input
   \*****************************************************/

   init_tblST_kmerFind(
      &kmerTblStackST,
      lenKmerUC
   );

   init_refST_kmerFind(
      &refStackAryST[0],
      lenKmerUC
   ); /*foward primer*/

   init_refST_kmerFind(
      &refStackAryST[1],
      lenKmerUC
   ); /*reverse primer*/

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
   ^   - set up primer sequences and kmer tables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   seqStackST.seqStr = (char *) forSeqStr;

   seqStackST.lenSeqUL =
     lenStr_ulCp(
        forSeqStr,
        0,
        0
      );
  
   seqStackST.idStr = "for";

   longestSeqUI = 
      addSeqToRefST_kmerFind(
         &kmerTblStackST,
         &refStackAryST[0],
         &seqStackST,
         minPercKmerF,
         longestSeqUI,
         &alnSetStackST
      ); /*set up the foward primer sequence*/

   seqStackST.seqStr = 0;
   seqStackST.lenSeqUL = 0;
   seqStackST.idStr = 0;

   if(longestSeqUI == 0)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error when reading primer sequences\n"
      );

      goto memErr_main_sec05_sub02;
   } /*If: I had an error*/


   seqStackST.seqStr = (char *) revSeqStr;

   seqStackST.lenSeqUL =
     lenStr_ulCp(
        revSeqStr,
        0,
        0
      );
  
   seqStackST.idStr = "rev";

   longestSeqUI = 
      addSeqToRefST_kmerFind(
         &kmerTblStackST,
         &refStackAryST[1],
         &seqStackST,
         minPercKmerF,
         longestSeqUI,
         &alnSetStackST
      ); /*set up the reverse primer sequence*/

   seqStackST.seqStr = 0;
   seqStackST.lenSeqUL = 0;
   seqStackST.idStr = 0;

   if(longestSeqUI == 0)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error when reading primer sequences\n"
      );

      goto memErr_main_sec05_sub02;
   } /*If: I had an error*/

   errSC =
      (schar)
      setUpTblST_kmerFind(
         &kmerTblStackST,
         extraNtInWinF,
         frameShiftF,
         longestSeqUI
      ); /*set up the kmer table*/

   if(errSC)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error when reading primer sequences\n"
      );

      goto memErr_main_sec05_sub02;
   } /*If: I had an error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - find and print primer positions
   ^   o main sec04 sub01:
   ^     - read in first sequence and start loop
   ^   o main sec04 sub02:
   ^     - check if sequence passes filters
   ^   o main sec04 sub03:
   ^     - find primer positions
   ^   o main sec04 sub04:
   ^     - check if has primers and is di/mvRNA
   ^   o main sec04 sub05:
   ^     - check if read is within print settings
   ^   o main sec04 sub06:
   ^     - print read (passed checks)
   ^   o main sec04 sub07:
   ^     - get the next sequence
   ^   o main sec04 sub08:
   ^     - check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - read in first sequence and start loop
   \*****************************************************/

   pidHeader_fluST(outFILE);

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
         goto nextSeq_main_sec04_sub05;

      if(seqStackST.lenSeqUL > maxLenUI)
         goto nextSeq_main_sec04_sub05;

      /**************************************************\
      * Main Sec04 Sub03:
      *   - find primer positions
      \**************************************************/

      if(fastBl)
      { /*If: i am searching with kmers*/
         errSC =
            fxFindPrims_kmerFind(
               &kmerTblStackST,    /*table structure*/
               refStackAryST,      /*primers looking for*/
               2,                  /*number of primers*/
               &seqStackST,        /*sequence to look at*/
               minPercScoreF,      /*min % for waterman*/
               codeAryUI,      /*# times found prim*/
               dirArySC,       /*direction of primer*/
               scoreArySL,     /*scores for prims*/
               seqStartAryUL,  /*1st align seq base*/
               seqEndAryUL,    /*last align seq base*/
               primStartAryUL, /*1st align prim base*/
               primEndAryUL,   /*last align prim bas*/
               &alnSetStackST      /*settings*/
            ); 
      } /*If: i am searching with kmers*/

      else
      { /*Else: slower waterman search*/
         errSC =
            waterFindPrims_kmerFind(
               refStackAryST,      /*primers looking for*/
               2,                  /*number of primers*/
               &seqStackST,        /*sequence to look at*/
               minPercScoreF,      /*min % for waterman*/
               codeAryUI,      /*# times found prim*/
               dirArySC,       /*direction of primer*/
               scoreArySL,     /*scores for prims*/
               seqStartAryUL,  /*1st align seq base*/
               seqEndAryUL,    /*last align seq base*/
               primStartAryUL, /*1st align prim base*/
               primEndAryUL,   /*last align prim bas*/
               &alnSetStackST      /*settings*/
            ); 
      } /*Else: slower waterman search*/

      /**************************************************\
      * Main Sec04 Sub04:
      *   - check if has primers and is di/mvRNA
      \**************************************************/

      if(! errSC)
      { /*If: I found at least on primer*/
         diFlagSC =
            detectDI_fluST(
               &fluStackST,
               (schar *) seqStackST.seqStr,  /*sequence*/
               codeAryUI,      /*# times found a primer*/
               dirArySC,       /*direction of primer*/
               seqStartAryUL,  /*first aligned seq base*/
               seqEndAryUL,    /*last aligned seq base*/
               &segNumSC,      /*gets segment number*/
               &mapLenUL
         );

         /***********************************************\
         * Main Sec04 Sub05:
         *   - check if read is within print settings
         \***********************************************/

         /*check to see if printing out*/
         keepBl |=
            ((!!(diFlagSC & def_fullFound_fluST))&vRnaBl);

         keepBl |=
            ((!!(diFlagSC & def_diFound_fluST)) &diRnaBl);

         keepBl &=
              ((!!(diFlagSC & def_partSeg_fluST)) &partBl)
            | (!!(diFlagSC & def_segFound_fluST));

         /*find percent of read is mapped region*/
         percLenDiffF = (float) mapLenUL;
         percLenDiffF /= (float) seqStackST.lenSeqUL;
         keepBl &= (percLenDiffF >= minPercLenF);

         /*see if meets maximum length for segment*/
         percLenDiffF = (float) mapLenUL;
         percLenDiffF /=
            (float) fluStackST.lenSegArySI[segNumSC];

         keepBl &= (percLenDiffF <= maxPercLenF);

         /*check if could id the segment (found result)*/
         if(keepBl)
         { /*If: print segment (passes filters)*/

            /********************************************\
            * Main Sec04 Sub06:
            *   - print read (passed checks)
            \********************************************/

            pid_fluST(
                &fluStackST,
                (schar *) seqStackST.idStr, /*read id*/
                segNumSC,      /*segment number*/
                diFlagSC,
                dirArySC,      /*primer direction*/
                scoreArySL,    /*primer scores*/
                refStackAryST[0].maxForScoreF, /*for*/
                refStackAryST[1].maxForScoreF, /*rev*/
                seqStartAryUL, /*1st seq base*/
                seqEndAryUL,   /*last seq base*/
                primStartAryUL,/*1st primer base*/
                primEndAryUL,  /*last primer base*/
                seqStackST.lenSeqUL,
                mapLenUL,
                outFILE
             ); /*print out read ids*/
         } /*If: print segment (passes filters)*/
      } /*If: I found at least on primer*/

      /**************************************************\
      * Main Sec04 Sub07:
      *   - get the next sequence
      \**************************************************/

      nextSeq_main_sec04_sub05:;

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
   * Main Sec04 Sub08:
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

   /*fluST for segment searching*/
   freeStack_fluST(&fluStackST);

   /*alignment settings*/
   freeStack_alnSet(&alnSetStackST);

   /*kmer find structures*/
   freeStack_seqST(&seqStackST);
   freeStack_tblST_kmerFind(&kmerTblStackST);

   freeStack_refST_kmerFind(&refStackAryST[0]);
   freeStack_refST_kmerFind(&refStackAryST[1]);

   if(outFILE && outFILE != stdout)
      fclose(outFILE);

   outFILE = 0;

   if(seqFILE && seqFILE != stdin)
      fclose(seqFILE);

   seqFILE = 0;

   exit(errSC);
} /*main*/
