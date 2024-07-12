/*########################################################
# Name: inputGetDIIds
#   - holds user input, help message, and version
#     functions for primFind
#   - Also has some default variables and versio numbers
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o fun01: pversion_inputGetDIIds
'     - prints the version for primFind to the input file
'   o fun02: phelp_inputGetDIIds
'     - prints the help message for primFind to the input
'       file
'   o .c fun03: getSeq_inputGetDIIds
'     - gets a sequence from a sequence user input flag
'   o fun04: getInput_inputGetDIIds
'     - gets user input for primFind
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
# else
   #include <stdlib.h> /*not really needed*/
#endif

#include <stdio.h>

#include "inputGetDIIds.h"
#include "fluST.h"

/*.h files only*/
#include "../generalLib/dataTypeShortHand.h"
#include "../generalLib/charCp.h"
#include "../generalLib/base10str.h"
#include "../generalLib/ntTo2Bit.h"
#include "kmerFind.h" /*only for default variables*/
#include "fluSeg.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: pversion_inputGetDIIds
|   - prints the version for primFind to the input file
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints;
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_inputGetDIIds(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "primFind verssion: %i-%02i-%02i\n",
       def_year_inputGetDIIds,
       def_month_inputGetDIIds,
       def_day_inputGetDIIds
   );
} /*pverson_inputGetDIIds*/

/*-------------------------------------------------------\
| Fun02: phelp_inputGetDIIds
|   - prints the help message for primFind to the input
|     file
| Input:
|   - pSegHelpBl:
|     o 1 print the segment help message
|     o 0 print the shorter help message
|   - outFILE:
|     o file to print the help message to
| Output:
|   - Prints;
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_inputGetDIIds(
   signed char pSegHelpBl,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints help message
   '   o fun02 sec01:
   '     - print out usage
   '   o fun02 sec02:
   '     - print out input
   '   o fun02 sec03:
   '     - print output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - print out usage
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siSeg = 0;

   fprintf(
      (FILE *) outFILE,
      "getDIIds -fq reads.fasta\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - Finds full length flu segments and checks if\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    segment is diRNA or vRNA (full length)\n"
   );

   if(pSegHelpBl)
       goto segInput_fun02_sec02_sub03_cat02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print out input
   ^   o fun02 sec02 sub01:
   ^      - input block header
   ^   o fun02 sec02 sub02:
   ^     - file IO
   ^   o fun02 sec02 sub03:
   ^     - changing sequence input
   ^   o fun02 sec02 sub04:
   ^     - filtering input
   ^   o fun02 sec02 sub05:
   ^     - kmer search options
   ^   o fun02 sec02 sub06:
   ^     - help message and version number options
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec02 Sub01:
   *   - input block header
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "Input:\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub02:
   *   - file IO
   *   o fun02 sec02 sub02 cat01:
   *     - read file (fastx)
   *   o fun02 sec02 sub02 cat02:
   *     - print di reads
   *   o fun02 sec02 sub02 cat03:
   *     - print full length reads
   *   o fun02 sec02 sub02 cat04:
   *     - print genomes with one primer id for segment
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat01:
   +   - read file (fastx)
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    fprintf(
      (FILE *) outFILE,
       "  -fq reads.fastq: [Required]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o fastq file with reads to search\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -fa reads.fasta: [Replaces -fq]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o fasta file with reads to search\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat02:
   +   - print di reads
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    if(def_pDIRna_inputGetDIIds)
       fprintf(
         (FILE *) outFILE,
          "  -di: [Yes]\n"
       );

    else
       fprintf(
         (FILE *) outFILE,
          "  -di: [No]\n"
       );

    fprintf(
      (FILE *) outFILE,
       "    o print out diRNA read ids\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o disable with \"-no-di\"\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat03:
   +   - print full length reads
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    if(def_pVRna_inputGetDIIds)
       fprintf(
         (FILE *) outFILE,
          "  -v: [Yes]\n"
       );

    else
       fprintf(
         (FILE *) outFILE,
          "  -v: [No]\n"
       );

    fprintf(
      (FILE *) outFILE,
       "    o print out read ids for vRNA (full length)\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o disable with \"-no-v\"\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat04:
   +   - print genoems with on solid primer
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    if(def_partSeq_inputGetDIIds)
       fprintf(
         (FILE *) outFILE,
          "  -part: [Yes]\n"
       );

    else
       fprintf(
         (FILE *) outFILE,
          "  -part: [No]\n"
       );

    fprintf(
      (FILE *) outFILE,
       "    o print out read ids with only one primer\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      having a segment id\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o disable with \"-no-part\"\n"
    );

   /*****************************************************\
   * Fun02 Sec02 Sub03:
   *   - changing sequence input
   *   o fun02 sec02 sub03 cat01:
   *     - primer sequence input
   *   o fun02 sec02 sub03 cat02:
   *     - segment input
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub03 Cat01:
   +   - primer sequence input
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    fprintf(
      (FILE *) outFILE,
       "  -prims %s,%s: [Optional]\n",
       forPrimStr_inputGetIds,
       revPrimStr_inputGetIds
    );

    fprintf(
      (FILE *) outFILE,
       "    o primer pair (forward,reverse) to find\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub03 Cat02:
   +   - segment input
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    fprintf(
      (FILE *) outFILE,
       "  -segment forward-seq,reverse-seq: [Optional]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o sequences used to id a segment\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      HA would be -HA GG,GTGT\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -segment-len length: [Optional]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o expected length of a segment\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      HA would be -HA 1778\n"
    );

   if(pSegHelpBl)
   { /*If: printing segment help message*/
      segInput_fun02_sec02_sub03_cat02:;

      for(
         siSeg = 0;
         siSeg < def_NSNum_fluSeg + 1;
         ++siSeg
      ){ /*Loop: print out segment sequences and length*/
         fprintf(
            (FILE *) outFILE,
            "    o -%s-len %i\t-%s %s,%s\n",
            segIdAryStr_fluSeg[siSeg],
            segLenArySS_fluSeg[siSeg],
            segIdAryStr_fluSeg[siSeg],
            forSeqAryStr_fluSeg[siSeg],
            revSeqAryStr_fluSeg[siSeg]
         );
      } /*Loop: print out segment sequences and length*/

      goto helpInput_fun02_sec02_sub06;
   } /*If: printing segment help message*/

   /*****************************************************\
   * Fun02 Sec02 Sub04:
   *   - filtering input
   \*****************************************************/

    fprintf(
      (FILE *) outFILE,
       "  -min-len %i: [Optional]\n",
       def_minLen_inputGetDIIds
    );

    fprintf(
      (FILE *) outFILE,
       "    o minimum read length to keep a read\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -max-len %i: [Optional]\n",
       def_maxLen_inputGetDIIds
    );

    fprintf(
      (FILE *) outFILE,
       "    o maximum read length to keep a read\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -min-perc-len %0.2f: [Optional]\n",
       def_minPercLen_inputGetDIIds
    );

    fprintf(
      (FILE *) outFILE,
       "    o minimum percent length to keep a read id\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o (length between primers) / (read length)\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -max-perc-len %0.2f: [Optional]\n",
       def_maxPercLen_inputGetDIIds
    );

    fprintf(
      (FILE *) outFILE,
       "    o maximum percent length to keep a read id\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o (mapped read length) / (expected length)\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -min-perc-score %0.2f: [Optional]\n",
       def_minPercScore_kmerFind
    );

    fprintf(
      (FILE *) outFILE,
       "    o minimum percent score (alignment) to\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      count a primer as mapped\n"
    );

   /*****************************************************\
   * Fun02 Sec02 Sub05:
   *   - kmer search options
   *   o fun02 sec02 sub05 cat01:
   *     - fast/slow options
   *   o fun02 sec02 sub05 cat02:
   *     - kmer length options
   *   o fun02 sec02 sub05 cat03:
   *     - min percent kmers setting
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat01:
   +   - fast/slow options
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_fastSearch_inputGetDIIds)
      fprintf(
        (FILE *) outFILE,
         "  -fast: [True]\n"
      );

   else
      fprintf(
        (FILE *) outFILE,
         "  -fast: [False]\n"
      );

    fprintf(
      (FILE *) outFILE,
       "    o do faster, but less percise kmer search\n"
    );

   fprintf(
     (FILE *) outFILE,
      "    o use \"-slow\" to search by Waterman\n"
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat02:
   +   - kmer length options
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
     (FILE *) outFILE,
      "  -len-kmer %i: [Optional]\n",
      def_lenKmer_kmerFind
   );

   fprintf(
     (FILE *) outFILE,
      "    o kmer size for -fast; bigger is faster, but\n"
   );

   fprintf(
     (FILE *) outFILE,
      "      also less sensitive (misses more)\n"
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat03:
   +   - min percent kmers setting
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
     (FILE *) outFILE,
      "  -min-perc-kmer %0.2f: [Optional]\n",
      def_minKmerPerc_kmerFind
   );

   fprintf(
     (FILE *) outFILE,
      "    o minimum percent of total kmers needed to\n"
   );

   fprintf(
     (FILE *) outFILE,
      "      do a Waterman on a window\n"
   );

   fprintf(
     (FILE *) outFILE,
      "    o higher is faster, but also less sensitive\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub06:
   *   - help message and version number options
   \*****************************************************/

   helpInput_fun02_sec02_sub06:;

   if(! pSegHelpBl)
   { /*If: printing the normal shorter help message*/
      fprintf(
        (FILE *) outFILE,
         "  -h: print this help message and exit\n"
      );

      fprintf(
        (FILE *) outFILE,
         "  -h-seg: print segment defaults and exit\n"
      );
   } /*If: printing the normal shorter help message*/

   else
   { /*Else: this is the segment help message*/
      fprintf(
        (FILE *) outFILE,
         "  -h: print non-segmet help message and exit \n"
      );

      fprintf(
        (FILE *) outFILE,
         "  -h-seg: print this help message and exit\n"
      );
   } /*Else: this is the segment help message*/

   fprintf(
     (FILE *) outFILE,
      "  -v: print version number and exit\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - print output block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Output:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - prints reads with detected primers to stdout\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o column 1 is read ID; see output for rest\n"
   );
} /*phelp_inputGetDIIds*/

/*-------------------------------------------------------\
| Fun03: getSeq_inputGetDIIds
|   - gets a sequence from a sequence user input flag
| Input:
|   - firstSeqStr:
|     o c-string to hold the first sequence read in
|   - secSeqStr:
|     o c-string to hold the second sequence read in
|   - inStr:
|     o c-string with user input
| Output:
|   - Modifies:
|     o firstSeqStr to point to frist sequence in user
|       input
|     o secSeqStr to point to second sequence in user
|       input
|     o inStr to have a null separating the two sequences
|   - Returns:
|     o 0: for no problems
|     o 1: for only one sequence
|     o 2: for more than two sequences
\-------------------------------------------------------*/
signed char
getSeq_inputGetDIIds(
   signed char **firstSeqStr,
   signed char **secSeqStr,
   signed char *inStr
){
   *firstSeqStr = inStr; 

   while(
      ntTo2Bit[(uchar) *inStr++] != def_err3rdBit_ntTo2Bit
   ) ;

   --inStr;

   if(*inStr == '\0')
      return 1; /*only one sequence*/

   *inStr = '\0';
   ++inStr;

   *secSeqStr = inStr;

   while(
      ntTo2Bit[(uchar) *inStr++] != def_err3rdBit_ntTo2Bit
   ) ;

   --inStr;

   if(*inStr != '\0')
      return 2; /*more than two sequences*/

   return 0;
} /*getSeq_inputGetDIIds*/

/*-------------------------------------------------------\
| Fun04: getInput_inputGetDIIds
|   - gets user input for primFind
| Input:
|   - numArgsSI:
|     o number of arguments/parameters in argAryStr
|   - argAryStr:
|     o array of c-strings with user input
|   - fluSTPtr:
|     o pointer to a fluST structure with settings to
|       change
|   - readFileStr:
|     o pointer to c-string to point to file with reads
|   - fqBl:
|     o pointer to char to be set to
|       - 1 = readFileStr is a fastq file
|       - 0 = readFileStr is a fasta file
|   - forPrimStr:
|     o pointer to c-string to hold forward primer seq
|   - revPrimStr:
|     o pointer to c-string to hold reverse primer seq
|   - minLenUI:
|     o pointer to unsigned long to hold minimum primer
|       length
|   - maxlenUI:
|     o pointer to unsigned long to hold maximum primer
|       length
|   - minPercLenF:
|     o pointer to float to hold the minimum percent
|       length to keep a read
|   - maxPercLenF:
|     o pointer to float to hold the maximum percent
|       length (of expected segment length) to keep a read
|   - minPerScoreF:
|     o pointer to float to hold the minimum percent score
|       to count a primer as mapped (found by waterman)
|   - pDiRnaBl:
|     o pointer to set to 1 (print diRNA) or 0 (no print)
|   - pVRnaBl:
|     o pointer to set to 1 (print vRNA) or 0 (no print)
|   - pPartBl:
|     o pointer to set to 1 to print out segments with
|       only one primer having segment id
|     o pointer to set to 0 to not print out 
|   - fastBl:
|     o ponter to char to be set to
|       - 1 = use faster, but less senistive kmer method
|       - 0 = use slower, but more percise waterman method
|   - lenKmerUC:
|     o pointer to unsigned char to hold kmer length for
|       the faster kmer search
|   - minPercKmerF:
|     o pointer to float to hold min percentage of kmers
|       needed to do waterman alignment (kmer search only)
| Output:
|   - Modifies:
|     o all input, except numArgsSI and argAryStr to hold
|       the users input
|   - Prints:
|     o help message and version number to stdout if
|       requested
|     o prints any errors to stderr
|   - Returns:
|     o 0 for no errors
|     o 1 for printed help message or version number
|     o 2 for an error
\-------------------------------------------------------*/
signed char
getInput_inputGetDIIds(
   int numArgsSI,
   char *argAryStr[],
   struct fluST *fluSTPtr,
   signed char **readFileStr,
   signed char *fqBl,
   signed char **forPrimStr,
   signed char **revPrimStr,
   unsigned int *minLenUI,
   unsigned int *maxLenUI,
   float *minPercLenF,
   float *maxPercLenF,
   float *minPercScoreF,
   signed char *pDiRnaBl,
   signed char *pVRnaBl,
   signed char *pPartBl,
   signed char *fastBl,
   unsigned char *lenKmerUC,
   float *minPercKmerF
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun04 TOC:
   '   - gets user input for primFind
   '   o fun04 sec01:
   '     - variable declerations
   '   o fun04 sec01:
   '     - get input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar *tmpStr = 0;
   sint siArg = 1;
   schar errSC = 0;
   sint eqlSI = 0;

   schar segFlagsStr[(def_NSNum_fluSeg + 1) << 1][10];
   sint siSeg = 0;

   schar *forSeqStr = 0;
   schar *revSeqStr = 0;

   for(
      siSeg = 0;
      siSeg < def_NSNum_fluSeg + 1;
      siSeg += 2
   ){ /*Loop: build segment id flags*/
      /*set up sequence flag*/
      tmpStr = segFlagsStr[siSeg];

      *tmpStr++ = '-';

      tmpStr +=
         cpDelim_charCp(
            (char *) tmpStr,
            (char *) segIdAryStr_fluSeg[siSeg >> 1],
            '\0'
         );

      /*set up length flag*/
      tmpStr = segFlagsStr[siSeg + 1];

      *tmpStr++ = '-';

      tmpStr +=
         cpDelim_charCp(
            (char *) tmpStr,
            (char *) segIdAryStr_fluSeg[siSeg >> 1],
            '\0'
         );

      tmpStr +=
         cpDelim_charCp(
            (char *) tmpStr,
            (char *) "-len",
            '\0'
         );
   } /*Loop: build segment id flags*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec02:
   ^   - get input
   ^   o fun04 sec02 sub01:
   ^     - check if have input and start loop
   ^   o fun04 sec02 sub02:
   ^     - file IO
   ^   o fun04 sec02 sub03:
   ^     - check if primer sequences
   ^   o fun04 sec02 sub04:
   ^     - filtering input (length/score)
   ^   o fun04 sec02 sub05:
   ^     - fast alignment input
   ^   o fun04 sec02 sub06:
   ^     - help messages
   ^   o fun04 sec02 sub07:
   ^     - check segment input or if invalid
   ^   o fun04 sec02 sub08:
   ^     - move to next argument
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec02 Sub01:
   *   - check if have input and start loop
   \*****************************************************/

   if(numArgsSI < 2)
   { /*If: nothing was input; assume help message wanted*/
      phelp_inputGetDIIds(
         0,      /*print short help message*/
         stdout
      );

      return 1;
   } /*If: nothing was input; assume help message wanted*/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/

      /**************************************************\
      * Fun04 Sec02 Sub02:
      *   - file IO
      \**************************************************/

      if(! eql_charCp("-fa",argAryStr[siArg],0)) 
      { /*If: a fasta sequence file was input*/
         ++siArg; /*move to argument*/
         *fqBl = 0;
         *readFileStr = (schar *) argAryStr[siArg];
      } /*If: a fasta sequence file was input*/

      else if(! eql_charCp("-fq",argAryStr[siArg],0)) 
      { /*Else If: a fastq sequence file was input*/
         ++siArg; /*move to argument*/
         *fqBl = 1;
         *readFileStr = (schar *) argAryStr[siArg];
      } /*Else If: a fastq sequence file was input*/

      else if(! eql_charCp("-di",argAryStr[siArg],0)) 
         *pDiRnaBl = 1;

      else if(! eql_charCp("-no-di",argAryStr[siArg],0)) 
         *pDiRnaBl = 0;

      else if(! eql_charCp("-v",argAryStr[siArg],0)) 
         *pVRnaBl = 1;

      else if(! eql_charCp("-no-v",argAryStr[siArg],0)) 
         *pVRnaBl = 0;

      else if(! eql_charCp("-part",argAryStr[siArg],0)) 
         *pPartBl = 1;

      else if(! eql_charCp("-no-part",argAryStr[siArg],0)) 
         *pPartBl = 0;

      /**************************************************\
      * Fun04 Sec02 Sub03:
      *   - check if primer sequences
      \**************************************************/

      else if(! eql_charCp("-prims",argAryStr[siArg],0)) 
      { /*Else If: primer sequences input*/
         ++siArg; /*move to argument*/

         errSC = 
            getSeq_inputGetDIIds(
               forPrimStr,
               revPrimStr,
               (schar *) argAryStr[siArg]
            ); /*set up pointers for primer sequence*/

         if(errSC)
         { /*If: there was a error*/
            if(errSC == 1)
               fprintf(
                  stderr,
                  "-prims %s has only one sequence\n",
                  argAryStr[siArg]
               );
            else
               fprintf(
                  stderr,
                  "-prims %s has 3 or more sequences\n",
                  argAryStr[siArg]
               );

            return 2;
         } /*If: there was a error*/
      } /*Else If: primer sequences input*/

      /**************************************************\
      * Fun04 Sec02 Sub04:
      *   - filtering input (length/score)
      \**************************************************/

      else if(! eql_charCp("-min-len",argAryStr[siArg],0)) 
      { /*Else If: the min length was input*/
         ++siArg; /*move to argument*/

         tmpStr =
            (schar *)
            strToUI_base10str(
               argAryStr[siArg],
               *minLenUI
            );

         if(*tmpStr != '\0')
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-min-len %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: the min length was input*/

      else if(! eql_charCp("-max-len",argAryStr[siArg],0)) 
      { /*Else If: the max length was input*/
         ++siArg; /*move to argument*/

         tmpStr =
            (schar *)
            strToUI_base10str(
               argAryStr[siArg],
               *maxLenUI
            );

         if(*tmpStr != '\0')
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-max-len %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: the max length was input*/

      else if(!
         eql_charCp("-min-perc-len",argAryStr[siArg],0)
      ){ /*Else If: the min percent length was input*/
         ++siArg; /*move to argument*/

          *minPercLenF = atof(argAryStr[siArg]);

         if(minPercLenF == 0)
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-min-perc-len %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: the min percent length was input*/

      else if(!
         eql_charCp("-max-perc-len",argAryStr[siArg],0)
      ){ /*Else If: the max percent length was input*/
         ++siArg; /*move to argument*/

          *maxPercLenF = atof(argAryStr[siArg]);

         if(minPercLenF == 0)
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-max-perc-len %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: the max percent length was input*/

      else if(
          ! eql_charCp(
               "-min-perc-score",
               argAryStr[siArg],
               0
            )
      ){ /*Else If: the min percent score was input*/
         ++siArg; /*move to argument*/
         *minPercScoreF = atof(argAryStr[siArg]);

         if(*minPercScoreF == 0)
         { /*If: invalid input*/
            fprintf(
               stderr,
               "-min-perc-score %s is non-numeric or 0\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: invalid input*/
      } /*Else If: the min percent score was input*/

      /**************************************************\
      * Fun04 Sec02 Sub05:
      *   - fast alignment input
      \**************************************************/

      /*check if using a fast (kmer) or slow alignment*/
      else if(! eql_charCp( "-fast", argAryStr[siArg], 0))
         *fastBl = 1;

      else if(! eql_charCp( "-slow", argAryStr[siArg], 0))
         *fastBl = 0;

      else if(!eql_charCp("-len-kmer",argAryStr[siArg],0))
      { /*Else If: kmer length was input*/
         ++siArg; /*move to argument*/

         tmpStr =
            (schar *)
            strToUC_base10str(
               argAryStr[siArg],
               *lenKmerUC
            );

         if(*tmpStr != '\0')
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-len-kmer %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: kmer length was input*/

      else if(
          ! eql_charCp(
               "-min-perc-kmer",
               argAryStr[siArg],
               0
            )
      ){ /*Else If: min percent kmers was input*/
         ++siArg; /*move to argument*/
         *minPercKmerF = atof(argAryStr[siArg]);

         if(*minPercScoreF == 0)
         { /*If: invalid input*/
            fprintf(
               stderr,
               "-min-perc-kmer %s is non-numeric or 0\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: invalid input*/
      } /*Else If: min percent kmers was input*/

      /**************************************************\
      * Fun04 Sec02 Sub06:
      *   - help and version messages
      \**************************************************/

      else if(
            ! eql_charCp("-h",argAryStr[siArg],0)
         || ! eql_charCp("--h",argAryStr[siArg],0)
         || ! eql_charCp("help",argAryStr[siArg],0)
         || ! eql_charCp("-help",argAryStr[siArg],0)
         || ! eql_charCp("--help",argAryStr[siArg],0)
      ){ /*Else If: help message*/
         phelp_inputGetDIIds(
            0,     /*print shorter help message*/
            stdout
         );

         return 1;
      } /*Else If: help message*/

      else if(
            ! eql_charCp("-h-seg",argAryStr[siArg],0)
         || ! eql_charCp("--h-seg",argAryStr[siArg],0)
         || ! eql_charCp("help-seg",argAryStr[siArg],0)
         || ! eql_charCp("-help-seg",argAryStr[siArg],0)
         || ! eql_charCp("--help-seg",argAryStr[siArg],0)
      ){ /*Else If: help message*/
         phelp_inputGetDIIds(
            1,    /*print help message with segment parm*/
            stdout
         );

         return 1;
      } /*Else If: help message*/

      else if(
            ! eql_charCp("-v",argAryStr[siArg],0)
         || ! eql_charCp("--v",argAryStr[siArg],0)
         || ! eql_charCp("version",argAryStr[siArg],0)
         || ! eql_charCp("-version",argAryStr[siArg],0)
         || ! eql_charCp("--version",argAryStr[siArg],0)
      ){ /*Else If: help message*/
         pversion_inputGetDIIds(stdout);
         return 1;
      } /*Else If: help message*/

      /**************************************************\
      * Fun04 Sec02 Sub07:
      *   - check segment input or if invalid
      *   o fun04 sec02 sub07 cat01:
      *     - check if segment length (and start loop)
      *   o fun04 sec02 sub07 cat02:
      *     - check if segment sequence
      *   o fun04 sec02 sub07 cat03:
      *     - check if valid input
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun04 Sec02 Sub07 Cat01:
      +   - check if segment length (and start loop)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else
      { /*Else: segment input or not recognized*/

         for(
            siSeg = 0;
            siSeg < def_NSNum_fluSeg + 1;
            siSeg += 2
         ){ /*Loop: check if segment input*/
            eqlSI =
               eql_charCp(
                  (char *) segFlagsStr[siSeg + 1],
                  (char *) argAryStr[siArg],
                  '\0'
               ); /*check if is segment length*/

            if(! eqlSI)
            { /*If: input is length of segment*/
               ++siArg;

               tmpStr =
                  (schar *)
                  strToSI_base10str(
                     argAryStr[siArg],
                     fluSTPtr->lenSegArySI[siSeg >> 1]
                  );

               if(*tmpStr != '\0')
               { /*If: I had a non-numeric entry*/
                  fprintf(
                     stderr,
                     "%s %s is non-numeric\n",
                     segFlagsStr[siSeg + 1],
                     argAryStr[siArg]
                  );

                  return 2;
               } /*If: I had a non-numeric entry*/

               break;
            } /*If: input is length of segment*/

            eqlSI =
               eql_charCp(
                  (char *) segFlagsStr[siSeg],
                  (char *) argAryStr[siArg],
                  '\0'
               ); /*check if is segement sequence*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun04 Sec02 Sub07 Cat02:
            +   - check if segment sequence
            \+++++++++++++++++++++++++++++++++++++++++++*/

            if(! eqlSI)
            { /*Else If: segment sequences input*/
               ++siArg; /*move to argument*/

               errSC = 
                  getSeq_inputGetDIIds(
                     &forSeqStr,
                     &revSeqStr,
                     (schar *) argAryStr[siArg]
                  ); /*set up pointers for sequences*/


               errSC =
                  rmSegSeqFrom_fluST(
                     fluSTPtr,
                     forSeqAryStr_fluSeg[siSeg >> 1],
                     0,
                     segLenArySS_fluSeg[siSeg >> 1]
                  );

               errSC =
                  addSegTo_fluST(
                     fluSTPtr,
                     forSeqStr,
                     segLenArySS_fluSeg[siSeg >> 1],
                     0,
                     siSeg >> 1
                  );

               if(errSC)
                  goto segErr_fun04_sec0x_sub0y;

               errSC =
                  rmSegSeqFrom_fluST(
                     fluSTPtr,
                     revSeqAryStr_fluSeg[siSeg >> 1],
                     1,
                     segLenArySS_fluSeg[siSeg >> 1]
                  );

               errSC =
                  addSegTo_fluST(
                     fluSTPtr,
                     revSeqStr,
                     segLenArySS_fluSeg[siSeg >> 1],
                     1,
                     siSeg >> 1
                  );

               if(errSC)
               { /*If: there was a error*/
                  segErr_fun04_sec0x_sub0y:;

                  if(errSC == 1)
                     fprintf(
                        stderr,
                        "%s %s has only one sequence\n",
                        segFlagsStr[siSeg],
                        argAryStr[siArg]
                     );
                  else
                     fprintf(
                        stderr,
                        "%s %s has 3 or more sequences\n",
                        segFlagsStr[siSeg],
                        argAryStr[siArg]
                     );

                  return 2;
               } /*If: there was a error*/
            } /*Else If: segment sequences input*/
         } /*Loop: check if segment input*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun04 Sec02 Sub07 Cat03:
         +   - check if valid input
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(siSeg >= def_NSNum_fluSeg + 1)
         { /*If: input no recongized*/
            fprintf(
               stderr,
               "%s is not recognized\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: input no recongized*/
      } /*Else: segment input or not recognized*/

      /**************************************************\
      * Fun04 Sec02 Sub08:
      *   - move to next argument
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   return 0;
} /*getInput_inputGetDIIds*/

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
