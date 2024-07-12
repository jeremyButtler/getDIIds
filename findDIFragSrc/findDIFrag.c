/*########################################################
# Name: findDIFrag
#   - find defective influenza segments using a waterman
#     alignment
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries and defined variables
'   o fun01: pversion_findDIFrag
'     - prints the version number for findDIFrag
'   o fun02: phelp_findDIFrag
'     - prints the help message for findDIFrag
'   o fun03: input_findDIFrag
'     - gets user input
'   o main:
'     - driver function for findDIFrag
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and defined variables
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../generalLib/samEntry.h"
#include "../generalLib/seqST.h"
#include "../generalLib/kmerCnt.h"

#include "../generalAln/dirMatrix.h"
#include "../generalAln/alnSet.h"
#include "../water/water.h"

#include "diScan.h"

/*.h only*/
#include "../generalLib/dataTypeShortHand.h"
#include "../generalLib/base10str.h"
#include "../generalLib/charCp.h"
#include "../generalLib/ulCp.h"

#include "fragSeg.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .h  #include "../generalLib/ntTo2Bit.h"
!   - .h  #include "../generalLib/ntTo5Bit.h"
!   - .h  #include "../generalLib/shellSort.h"
!   - .h  #include "../generalLib/genMath.h"
!   - .h  #include "../generalAln/alnDefs.h"
!   - .h  #include "../generalAln/indexToCoord.h"
!   - .h  #include "../getDICoordsSrc/diCoords.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_year_findDIFrag 2024
#define def_month_findDIFrag 7
#define def_day_findDIFrag 12

#define def_lenKmer_findDIFrag 7 /*length of one kmer*/
#define def_minPercScore_findDIFrag 0.5f /*90% min socre*/
#define def_minDels_findDIFrag 20 /*min del size to be DI*/
#define def_minPadNt_findDIFrag 12

#define def_pDI_findDIFrag 1 /*1 = print DI reads*/
#define def_pVRna_findDIFrag 1 /*1 = print vRNA reads*/

/*-------------------------------------------------------\
| Fun01: pversion_findDIFrag
|   - prints the version number for findDIFrag
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints:
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_findDIFrag(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "findDIFrag version: %i-%02i-%02i\n",
      def_year_findDIFrag,
      def_month_findDIFrag,
      def_day_findDIFrag
   );
} /*pversion_findDIFrag*/

/*-------------------------------------------------------\
| Fun02: phelp_findDIFrag
|   - prints the help message for findDIFrag
| Input:
|   - outFILE:
|     o file to print help message to
| Output:
|   - Prints:
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_findDIFrag(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints the help message for findDIFrag
   '   o fun02 sec01:
   '     - usage
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - usage
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
    (FILE *) outFILE,
    "findDIFrag -fq reads.fq -out-sam out.sam > out.tsv\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - finds the number of DI events in reads\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print the input options
   ^   o fun02 sec02 sub01:
   ^     - print input block header
   ^   o fun02 sec02 sub02:
   ^     - print file io
   ^   o fun02 sec02 sub03:
   ^     - print filter for DI events
   ^   o fun02 sec02 sub04:
   ^     - filters for printing
   ^   o fun02 sec02 sub05:
   ^     - help and version flag
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec02 Sub01:
   *   - print input block header
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "Input:\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub02:
   *   - print file io
   \*****************************************************/

   /*input fastq*/
   fprintf(
      (FILE *) outFILE,
      "  -fq reads.fq: [Required]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o fastq file with reads to find DI events in\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o use \"-\" for stdin\n"
   );

   /*output tsv*/
   fprintf(
      (FILE *) outFILE,
      "  -out out.tsv: [Optinal; stdout]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o file to save tsv to \n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o use \"-\" for stdout\n"
   );

   /*output alignments (sam file)*/
   fprintf(
      (FILE *) outFILE,
      "  -out-sam out.sam: [Optinal; No]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o sam file to save aligned reads to\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub03:
   *   - print filter for DI events
   \*****************************************************/

   /*print min del size*/
   fprintf(
      (FILE *) outFILE,
      "  -min-del %i: [Optinal; %i]\n",
      def_minDels_findDIFrag,
      def_minDels_findDIFrag
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum number of deletions to have a DI\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      event\n"
   );

   /*min nucleotides at start/end of sequence*/
   fprintf(
      (FILE *) outFILE,
      "  -min-nt-pad %u: [Optinal; %u]\n",
      def_minPadNt_findDIFrag,
      def_minPadNt_findDIFrag
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum number of matches, SNPs, and\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      insertions needed at start/end of read\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      before counting large deltions as DI\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub04:
   *   - filters for printing
   \*****************************************************/

   /*diRNA flag*/
   if(def_pDI_findDIFrag)
      fprintf(
         (FILE *) outFILE,
         "  -di: [Optinal; Yes]\n"
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -di: [Optinal; No]\n"
      );

   fprintf(
      (FILE *) outFILE,
      "    o print out reads with DI events\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o disable with \"-no-di\"\n"
   );

   /*vRNA flag*/
   if(def_pVRna_findDIFrag)
      fprintf(
         (FILE *) outFILE,
         "  -vrna: [Optinal; Yes]\n"
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -vrna: [Optinal; No]\n"
      );

   fprintf(
      (FILE *) outFILE,
      "    o print out reads with no DI events\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o disable with \"-no-vrna\"\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub05:
   *   - help and version flag
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  -h: print this help message and exit\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -v: print version number and exit\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - outpub block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Output:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -prints: tsv with read ids, mapped segment,\n"
   );

   fprintf(
      (FILE *) outFILE,
      " and number of DI events to stdout/file (-out)\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -prints: alignments to -out-sam if file input\n"
   );
} /*phelp_findDIFrag*/


/*-------------------------------------------------------\
| Fun03: input_findDIFrag
|   - gets user input
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o c-string array with user input
|   - fqFileStrPtr:
|     o pointer to c-sting to set to fastq file name
|   - outFileStrPtr:
|     o pointer to c-sting to set to output tsv file name
|   - samFileStrPtr:
|     o pointer to c-sting to set to output sam file name
|   - minDIDelUIPtr:
|     o pointer to unsigned int to hold min deletion size
|       to call DI
|   - minPadNtUIPtr:
|     o pointer to unsigned int to hold the minimum
|       number of nucleotides before/after a large DI
|       event
|   - pDIBlPtr:
|     o pointer to signed char, set to 1 if printing DI
|       reads
|   - pVRnaBlPtr:
|     o pointer to signed char, set to 1 if printing vRNA
|       reads
| Output:
|   - Modifies:
|     o all input except numArgsSI and argAryStr to hold
|       user input
|   - Returns:
|     o 0 for no errors
|     o 1 if printed help message or version number
|     o 2 for an error
\-------------------------------------------------------*/
signed char
input_findDIFrag(
   int numArgsSI,
   char *argAryStr[],
   signed char **fqFileStrPtr,
   signed char **outFileStrPtr,
   signed char **samFileStrPtr,
   unsigned int *minDIDelUIPtr,
   unsigned int *minPadNtUIPtr,  /*ends non-del length*/
   signed char *pDIBlPtr,
   signed char *pVRnaBlPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - gets user input
   '   o fun03 sec01:
   '     - variable declerations
   '   o fun03 sec02:
   '     - check if something was input
   '   o fun03 sec03:
   '     - get user input
   '   o fun03 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siArg = 1;
   schar *tmpStr = 0;
   schar errSC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - check if something was input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numArgsSI < 1)
      goto phelp_fun03_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - get user input
   ^   o fun03 sec03 sub01:
   ^     - start loop and check for file io
   ^   o fun03 sec03 sub02:
   ^     - get deletion counts
   ^   o fun03 sec03 sub03:
   ^     - print settings
   ^   o fun03 sec03 sub04:
   ^     - check for help message requests
   ^   o fun03 sec03 sub04:
   ^     - check for version number requests
   ^   o fun03 sec03 sub05:
   ^     - check for invalid input
   ^   o fun03 sec03 sub06:
   ^     - move to the next argument
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - start loop and check for file io
   \*****************************************************/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/

      if(! eql_charCp("-fq", argAryStr[siArg], 0))
      { /*If: a fastq file was input*/
         ++siArg;
         *fqFileStrPtr = (schar *) argAryStr[siArg];
      } /*If: a fastq file was input*/

      else if(! eql_charCp("-out", argAryStr[siArg], 0))
      { /*Else If: an output file was input*/
         ++siArg;
         *outFileStrPtr = (schar *) argAryStr[siArg];
      } /*Else If: an ouput file was input*/

      else if(! eql_charCp("-out-sam",argAryStr[siArg],0))
      { /*Else If: an sam output file was input*/
         ++siArg;
         *samFileStrPtr = (schar *) argAryStr[siArg];
      } /*Else If: an sam output file was input*/

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - get deletion counts
      \**************************************************/

      else if(! eql_charCp("-min-del",argAryStr[siArg],0))
      { /*Else If: number minimum deletions input*/
         ++siArg;

         tmpStr =
            (schar *)
            strToUI_base10str(
               argAryStr[siArg],
               *minDIDelUIPtr
            );

         if(*tmpStr != '\0')
         { /*If: had non-numeric input or overflow*/
            fprintf(
               stderr,
               "-min-del %s is non-numeric or to large\n",
               argAryStr[siArg]
            );

            goto err_fun03_sec04;
         } /*If: had non-numeric input or overflow*/
      } /*Else If: number minimum deletions input*/

      else if(
         ! eql_charCp("-min-nt-pad",argAryStr[siArg],0)
      ){ /*Else If: is the minimum DI deletion length*/
         ++siArg;
         tmpStr =
            (schar *)
            strToUI_base10str(
               argAryStr[siArg],
               *minPadNtUIPtr
            );

         if(*tmpStr != '\0')
         { /*If: had non-numeric entry*/
            fprintf(
               stderr,
               "-min-nt-pad %s; non-numeric or to big\n",
               argAryStr[siArg]
            );
 
            goto err_fun03_sec04;
         } /*If: had non-numeric entry*/
      } /*Else If: is the minimum DI deletion length*/

      /**************************************************\
      * Fun03 Sec03 Sub03:
      *   - print settings
      \**************************************************/

      else if(! eql_charCp("-di",argAryStr[siArg],0))
         *pDIBlPtr = 1;

      else if(! eql_charCp("-no-di",argAryStr[siArg],0))
         *pDIBlPtr = 0;

      else if(! eql_charCp("-vrna",argAryStr[siArg],0))
         *pVRnaBlPtr = 1;

      else if(! eql_charCp("-no-vrna",argAryStr[siArg],0))
         *pVRnaBlPtr = 0;

      /**************************************************\
      * Fun03 Sec03 Sub04:
      *   - check for help message requests
      \**************************************************/

      else if(! eql_charCp("-h", argAryStr[siArg], 0))
         goto phelp_fun03_sec04;

      else if(! eql_charCp("--h", argAryStr[siArg], 0))
         goto phelp_fun03_sec04;

      else if(! eql_charCp("help", argAryStr[siArg], 0))
         goto phelp_fun03_sec04;

      else if(! eql_charCp("-help", argAryStr[siArg], 0))
         goto phelp_fun03_sec04;

      else if(! eql_charCp("--help", argAryStr[siArg], 0))
         goto phelp_fun03_sec04;

      /**************************************************\
      * Fun03 Sec03 Sub04:
      *   - check for version number requests
      \**************************************************/

      else if(! eql_charCp("-v", argAryStr[siArg], 0))
         goto pversion_fun03_sec04;

      else if(! eql_charCp("--v", argAryStr[siArg], 0))
         goto pversion_fun03_sec04;

      else if(! eql_charCp("version",argAryStr[siArg],0))
         goto pversion_fun03_sec04;

      else if(! eql_charCp("-version",argAryStr[siArg],0))
         goto pversion_fun03_sec04;

      else if(!eql_charCp("--version",argAryStr[siArg],0))
         goto pversion_fun03_sec04;

      /**************************************************\
      * Fun03 Sec03 Sub05:
      *   - check for invalid input
      \**************************************************/

      else
      { /*Else: input not recongnized*/
         fprintf(
            stderr,
            "%s is not recongized\n",
            argAryStr[siArg]
         );

         goto err_fun03_sec04;
      } /*Else: input not recongnized*/

      /**************************************************\
      * Fun03 Sec03 Sub06:
      *   - move to the next argument
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto cleanUp_fun03_sec04;

   phelp_fun03_sec04:;
   phelp_findDIFrag(stdout);
   errSC = 1;
   goto cleanUp_fun03_sec04;

   pversion_fun03_sec04:;
   pversion_findDIFrag(stdout);
   errSC = 1;
   goto cleanUp_fun03_sec04;

   err_fun03_sec04:;
   errSC = 2;
   goto cleanUp_fun03_sec04;

   cleanUp_fun03_sec04:;
   return(errSC);
} /*input_findDIFrag*/

/*-------------------------------------------------------\
| Main:
|   - driver function for findDIFrag
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o array of c-strings with user input
| Output:
|   - Prints:
|     o prints tsv with ids, scores, reference mapped to
|       and other data of detected DI fragment to stdout
\-------------------------------------------------------*/
#ifdef PLAN9
void
#else
int 
#endif
main(
   int numArgsSI,         /*number input arguments*/
   char *argAryStr[]     /*arguments input*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main:
   '   - driver function for findDIFrag
   '   o main sec01:
   '     - variable declerations
   '   o main sec02:
   '     - initialize and get/check user input
   '   o main sec03:
   '     - open files
   '   o main sec04:
   '     - check for DI sequences
   '   o main sec05:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar *fqFileStr = 0;
   schar *outFileStr = 0;
   schar *samFileStr = 0; /*save alignments to*/

   schar pDIRnaBl = def_pDI_findDIFrag;
   schar pVRnaBl = def_pVRna_findDIFrag;

   schar errSC = 0;

   sint segSI = 0; /*index of segment read mapped to*/
   sint numDIEventsSI = 0; /*number DI events in seq*/
   uint minDIDelUI = def_minDels_findDIFrag;
   uint minPadNtUI = def_minPadNt_findDIFrag;

   /*scoring variables for waterman*/
   float minPercScoreF = def_minPercScore_findDIFrag;

   uchar lenKmerUC = def_lenKmer_findDIFrag;/*kmer length*/
   sint *kmerHeapArySI = 0; /*for get_kmerCnt*/;
   sint *cntHeapArySI = 0;  /*for get_kmerCnt*/
   uint maxKmerUI = 0;      /*for allocating arrays*/

   struct samEntry samStackST;
   schar *buffHeapStr = 0;
   ulong lenBuffUL = 0;

   struct seqST seqStackST;
   struct alnSet alnSetStackST;
   struct dirMatrix matrixStackST;
   struct kmerCnt kmerStackAryST[def_NSNum_fragSeg + 1];

   uint seqEntryUI = 0; /*entry error happend on*/
   FILE *fqFILE = 0;
   FILE *outFILE = 0;
   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - initialize and get/check user input
   ^   o main sec02 sub01:
   ^     - initialize non-error causing structures
   ^   o main sec02 sub02:
   ^     - get user input
   ^   o main sec02 sub03:
   ^     - allocate memory for samEntry structure
   ^   o main sec02 sub04:
   ^     - setup and add references to kmerCnt structures
   ^   o main sec02 sub05:
   ^     - allocate memory for kmer arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize non-error causing structures
   \*****************************************************/

   init_seqST(&seqStackST);
   init_alnSet(&alnSetStackST);
   init_dirMatrix(&matrixStackST);
   init_samEntry(&samStackST);

   for(
     segSI = 0;
     segSI < def_NSNum_fragSeg + 1;
     ++segSI
   ) init_kmerCnt(&kmerStackAryST[segSI]);

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get user input
   \*****************************************************/

   errSC =
      input_findDIFrag(
         numArgsSI,
         argAryStr,
         &fqFileStr,
         &outFileStr,
         &samFileStr,
         &minDIDelUI,
         &minPadNtUI,   /*min pad length to count DI*/
         &pDIRnaBl,
         &pVRnaBl
      ); /*get user input*/

   if(errSC)
   { /*If: had input error*/
      --errSC; /*convert help/version error to no error*/
      goto cleanUp_main_sec05_sub04; /*help/version number*/
   } /*If: had input error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - allocate memory for samEntry structure
   \*****************************************************/

   errSC = setup_samEntry(&samStackST);

   if(errSC)
      goto memErr_main_sec05_sub02;

   /*****************************************************\
   * Main Sec02 Sub04:
   *   - setup and add references to kmerCnt structures
   \*****************************************************/

   for(
     segSI = 0;
     segSI <= def_NSNum_fragSeg;
     ++segSI
   ){ /*Loop: initialize kmerCnt array*/
      errSC =
         setup_kmerCnt(
            &kmerStackAryST[segSI],
            lenKmerUC
         ); /*allocate memory for kmerCnt structures*/

      if(errSC)
         goto memErr_main_sec05_sub02;
         
      /*add in sequence*/
      seqStackST.seqStr =
          (char *)
          segSeqAryStr_fragSeg[segSI];

      seqStackST.lenSeqUL =
         lenStr_ulCp(
            segSeqAryStr_fragSeg[segSI],
            0,
            0
         ); /*get sequence length*/

      /*add in segment id*/
      seqStackST.idStr =
         (char *)
         segIdAryStr_fragSeg[segSI];

      seqStackST.lenIdUL =
         lenStr_charCp(
            (char *) segIdAryStr_fragSeg[segSI],
            0
         ); /*get sequence length*/

      errSC =
         addSeq_kmerCnt(
            &kmerStackAryST[segSI],
            &seqStackST
         ); /*add flu segment to array*/

       seqStackST.seqStr = 0;
       seqStackST.idStr = 0;

      if(errSC)
         goto memErr_main_sec05_sub02;
   } /*Loop: initialize kmerCnt array*/

   if(errSC)
      goto memErr_main_sec05_sub02;

   /*****************************************************\
   * Main Sec02 Sub05:
   *   - allocate memory for kmer arrays
   \*****************************************************/

   /*these variables are set to null, so safe to free
   `  when no memory was allocated
   */
   maxKmerUI = 1;

   for(
     segSI = 0;
     segSI < lenKmerUC;
     ++segSI
   ) maxKmerUI <<= 2;

   kmerHeapArySI = malloc((maxKmerUI + 1) * sizeof(sint));

   if(! kmerHeapArySI)
      goto memErr_main_sec05_sub02;

   cntHeapArySI = malloc((maxKmerUI + 1) * sizeof(sint));

   if(! cntHeapArySI)
      goto memErr_main_sec05_sub02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - open files
   ^   o main sec03 sub01:
   ^     - open fastq file
   ^   o main sec03 sub02:
   ^     - open output file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub01:
   *   - open fastq file
   \*****************************************************/

   if(! fqFileStr || *fqFileStr == '-')
      fqFILE = stdin;

   else
   { /*Else: user input fastq file*/
      fqFILE = 
         fopen(
            (char *) fqFileStr,
            "r"
         ); /*open the fastq file*/

      if(! fqFILE)
      { /*If: could not open file*/
         fprintf(
            stderr,
            "could not open -fq %s\n",
            fqFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: could not open file*/
   } /*Else: user input fastq file*/

   /*****************************************************\
   * Main Sec03 Sub02:
   *   - open output file
   \*****************************************************/

   if(! outFileStr || *outFileStr == '-')
      outFILE = stdout;

   else
   { /*Else: user input an output file*/
      outFILE = 
         fopen(
            (char *) outFileStr,
            "w"
         ); /*open output file*/

      if(! outFILE)
      { /*If: could not open file*/
         fprintf(
            stderr,
            "could not open -out %s\n",
            outFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: could not open file*/
   } /*Else: user input an output file*/


   /*****************************************************\
   * Main Sec03 Sub03:
   *   - open sam output file
   \*****************************************************/

   if(samFileStr)
   { /*If: user input an sam output file*/
      samFILE = 
         fopen(
            (char *) samFileStr,
            "w"
         ); /*open output file*/

      if(! samFILE)
      { /*If: could not open file*/
         fprintf(
            stderr,
            "could not open -out-sam %s\n",
            samFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: could not open file*/
   } /*If: user input an sam output file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - check for DI sequences
   ^   o main sec04 sub01:
   ^     - print out header(s)
   ^   o main sec04 sub02:
   ^     - read in first sequence and start loop
   ^   o main sec04 sub03:
   ^     - scan for DI sequences
   ^   o main sec04 sub04:
   ^     - print out results
   ^   o main sec04 sub05:
   ^     - move to next sequence
   ^   o main sec04 sub05:
   ^     - check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - print out header(s)
   \*****************************************************/

   phead_diScan(outFILE); /*print header*/

   if(samFILE)
   { /*If: printing to sam file*/
      fprintf(
         samFILE,
         "@HD\tVN:1.6\tSO:unsorted\tGO:none\n"
      );

      for(
         segSI = 0;
         segSI <= def_NSNum_fragSeg;
         ++segSI
      ){ /*Loop: print out sequence headers*/
         fprintf(
            samFILE,
            "@SQ\tSN:%s\tLN:%lu\n",
            kmerStackAryST[segSI].forSeqST->idStr,
            kmerStackAryST[segSI].forSeqST->lenSeqUL
         );
      } /*Loop: print out sequence headers*/

      fprintf(
        samFILE,
        "@PG\tID:%s\tVN:%i-%02i-%02i\t%s",
        "findDIFrag",
        def_year_findDIFrag,
        def_month_findDIFrag,
        def_day_findDIFrag,
        "findDIFrag"
     ); /*print out first part of program id tag*/

      if(pDIRnaBl)
        fprintf( samFILE, " -di");
      else
        fprintf( samFILE, " -no-di");

      if(pVRnaBl)
        fprintf( samFILE, " -vrna");
      else
        fprintf( samFILE, " -no-vrna");

     if(fqFileStr)
        fprintf( samFILE, " -fq %s", fqFileStr);
      else
        fprintf( samFILE, " -fq -");

      if(outFileStr)
        fprintf( samFILE, " -out %s", outFileStr);
      else
        fprintf( samFILE, " -out -");

     fprintf(samFILE, " -out-sam %s\n", samFileStr);
   } /*If: printing to sam file*/

   /*****************************************************\
   * Main Sec04 Sub02:
   *   - read in first sequence and start loop
   \*****************************************************/
   errSC =
     (schar)
     getFqSeq_seqST(
        fqFILE,
        &seqStackST
     );

   while(! errSC)
   { /*Loop: find DI sequences*/
      ++seqEntryUI;

      /**************************************************\
      * Main Sec04 Sub03:
      *   - scan for DI sequences
      \**************************************************/

      /*TODO: FIX MEMORY CRASH WHEN RELLOCATING Q-SCORE
      `   ARRAY FOR SAMSTACKST
      */
      numDIEventsSI =
         waterScan_diScan(
            &seqStackST,    /*sequence*/
            kmerStackAryST,/*references*/
            def_NSNum_fragSeg + 1, /*last segment in FLU*/
            kmerHeapArySI,  /*gets kmers in sequence*/
            cntHeapArySI,   /*gets sequenced kmer counts*/
            minPercScoreF,  /*min align score to keep*/
            lenKmerUC,      /*length of one kmer*/
            minDIDelUI,     /*min deletion size to be DI*/
            minPadNtUI,     /*min nt at ends before DI*/
            &samStackST,    /*will have aligned sequence*/
            &segSI,         /*returned segment number*/
            &alnSetStackST, /*alignment settings*/
            &matrixStackST  /*matrix for waterman to use*/
         );

      if(numDIEventsSI < 0)
      { /*If: had an error*/
         if(numDIEventsSI == def_noMatch_diScan)
            goto getNextSeq_main_sec04_sub04;
            /*likely not a flu segment*/

         if(numDIEventsSI == def_memErr_diScan)
            goto memErr_main_sec05_sub02;
            /*def_memErr_diScan is always negative*/
      } /*If: had an error*/

      /**************************************************\
      * Main Sec04 Sub04:
      *   - print out results
      \**************************************************/

      samStackST.lenRefIdUC =
         cpDelim_charCp(
            (char *) samStackST.refIdStr,
            (char *) segIdAryStr_fragSeg[segSI],
            '\0'
         ); /*copy the mapped segment id*/

      if(numDIEventsSI)
      { /*If: DI sequence*/
         if(pDIRnaBl)
            pfrag_diScan(
               &samStackST,
               numDIEventsSI,
               segLenArySS_fragSeg[segSI],
               matrixStackST.scoreSL,
               outFILE
            );

         if(samFILE)
            p_samEntry(
               &samStackST,
               (char **) &buffHeapStr,
               &lenBuffUL,
               0,
               samFILE
            ); /*print alignment*/
      } /*If: DI sequence*/

      else if(pVRnaBl)
      { /*Else: not DI sequence*/
         pfrag_diScan(
            &samStackST,
            numDIEventsSI,
            segLenArySS_fragSeg[segSI],
            matrixStackST.scoreSL,
            outFILE
         );

         if(samFILE)
            p_samEntry(
               &samStackST,
               (char **) &buffHeapStr,
               &lenBuffUL,
               0,
               samFILE
            ); /*print alignment*/
      } /*Else: not DI sequence*/

      /**************************************************\
      * Main Sec04 Sub05:
      *   - move to next sequence
      \**************************************************/

      getNextSeq_main_sec04_sub04:;

      errSC =
         (schar)
         getFqSeq_seqST(
            fqFILE,
            &seqStackST
         );
   } /*Loop: find DI sequences*/

   /*****************************************************\
   * Main Sec04 Sub06:
   *   - check for errors
   \*****************************************************/

   if(errSC != def_EOF_seqST)
   { /*If: had an error*/
      if(errSC == def_memErr_seqST)
      { /*If: memory error*/
         fprintf(
            stderr,
            "memory error\n"
         );

         goto memErr_main_sec05_sub02;
      } /*If: memory error*/
      
      fprintf(
         stderr,
         "problem with entry %u in -fq %s\n",
          seqEntryUI,
          fqFileStr
      );

      goto memErr_main_sec05_sub02;
   } /*If: had an error*/

   if(! seqEntryUI)
   { /*If: no sequences in file*/
      fprintf(
         stderr,
         "-fq %s has no reads\n",
          fqFileStr
      );

      goto memErr_main_sec05_sub02;
   } /*If: no sequences in file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05:
   ^   - clean up
   ^   o main sec05 sub01:
   ^     - clean up after no errors
   ^   o main sec05 sub02:
   ^     - memory error; set return value
   ^   o main sec05 sub03:
   ^     - file error; set return value
   ^   o main sec05 sub04:
   ^     - clean up after no error or after errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec05 Sub01:
   *   - clean up after no errors
   \*****************************************************/

   goto cleanUp_main_sec05_sub04; /*no error*/

   /*****************************************************\
   * Main Sec05 Sub02:
   *   - memory error; set return value
   \*****************************************************/

   memErr_main_sec05_sub02:;
   errSC = 1;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub03:
   *   - file error; set return value
   \*****************************************************/

   fileErr_main_sec05_sub03:;
   errSC = 2;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub04:
   *   - clean up after no error or after errors
   \*****************************************************/

   cleanUp_main_sec05_sub04:;

   if(
         fqFILE
      && fqFILE != stdin
      && fqFILE != stdout
   ) fclose(fqFILE);

   fqFILE = 0;

   if(
         outFILE
      && outFILE != stdin
      && outFILE != stdout
   ) fclose(outFILE);

   outFILE = 0;

   if(
         samFILE
      && samFILE != stdin
      && samFILE != stdout
   ) fclose(samFILE);

   samFILE = 0;

   free(kmerHeapArySI);
   free(cntHeapArySI);

   freeStack_seqST(&seqStackST);
   freeStack_alnSet(&alnSetStackST);
   freeStack_dirMatrix(&matrixStackST);

   freeStackAry_kmerCnt(
      kmerStackAryST,
      def_NSNum_fragSeg
   );

   free(buffHeapStr);
   buffHeapStr = 0;

   exit(errSC);
} /*main*/
