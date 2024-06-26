/*########################################################
# Name: inputPrimFind
#   - holds user input, help message, and version
#     functions for primFind
#   - Also has some default variables and versio numbers
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o fun01: pversion_inputPrimFind
'     - prints the version for primFind to the input file
'   o fun02: phelp_inputPrimFind
'     - prints the help message for primFind to the input
'       file
'   o fun03: getInput_inputPrimFind
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

#include "inputPrimFind.h"

/*.h files only*/
#include "generalLib/dataTypeShortHand.h"
#include "generalLib/charCp.h"
#include "generalLib/base10str.h"
#include "kmerFind.h" /*only for default variables*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: pversion_inputPrimFind
|   - prints the version for primFind to the input file
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints;
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_inputPrimFind(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "primFind verssion: %i-%02i-%02i\n",
       def_year_inputPrimFind,
       def_month_inputPrimFind,
       def_day_inputPrimFind
   );
} /*pverson_inputPrimFind*/

/*-------------------------------------------------------\
| Fun02: phelp_inputPrimFind
|   - prints the help message for primFind to the input
|     file
| Input:
|   - outFILE:
|     o file to print the help message to
| Output:
|   - Prints;
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_inputPrimFind(
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

   fprintf(
      (FILE *) outFILE,
      "primFind -prim-tsv primers.tsv -fq reads.fasta\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - Finds reads with input primers\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print out input
   ^   o fun02 sec02 sub01:
   ^      - input block header
   ^   o fun02 sec02 sub02:
   ^     - file input
   ^   o fun02 sec02 sub03:
   ^     - filtering input
   ^   o fun02 sec02 sub04:
   ^     - kmer search options
   ^   o fun02 sec02 sub05:
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
   *   - file input
   *   o fun02 sec02 sub02 cat01:
   *     - primer tsv file input
   *   o fun02 sec02 sub02 cat02:
   *     - primer fasta file input
   *   o fun02 sec02 sub02 cat03:
   *     - read input
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat01:
   +   - primer tsv file input
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    fprintf(
      (FILE *) outFILE,
       "  -prim-tsv primers.tsv: [Required]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o tsv file with primers to search for\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      - id\tpaired\tfor_seq\trev_seq\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      - column 1: primer name\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      - column 2: True = primers paired\n"
    );

    fprintf(
      (FILE *) outFILE,
       "        - both present in opposite directions\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      - column 3: forward primer sequence or NA\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      - column 4: reverse primer sequence or NA\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat02:
   +   - primer fasta file input
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    fprintf(
      (FILE *) outFILE,
       "  -prim-fa primers.fasta: [Replaces -prim-tsv]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o fasta file with primers to search for\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat03:
   +   - read input
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

   /*****************************************************\
   * Fun02 Sec02 Sub03:
   *   - filtering input
   \*****************************************************/

    fprintf(
      (FILE *) outFILE,
       "  -min-len %i: [Optional]\n",
       def_minLen_inputPrimFind
    );

    fprintf(
      (FILE *) outFILE,
       "    o minimum read length to keep a read\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -max-len %i: [Optional]\n",
       def_maxLen_inputPrimFind
    );

    fprintf(
      (FILE *) outFILE,
       "    o maximum read length to keep a read\n"
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
   * Fun02 Sec02 Sub04:
   *   - kmer search options
   *   o fun02 sec02 sub04 cat01:
   *     - fast/slow options
   *   o fun02 sec02 sub04 cat02:
   *     - kmer length options
   *   o fun02 sec02 sub04 cat03:
   *     - min percent kmers setting
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub04 Cat01:
   +   - fast/slow options
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_fastSearch_inputPrimFind)
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
   + Fun02 Sec02 Sub04 Cat02:
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
   + Fun02 Sec02 Sub04 Cat03:
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
   * Fun02 Sec02 Sub05:
   *   - help message and version number options
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
} /*phelp_inputPrimFind*/

/*-------------------------------------------------------\
| Fun03: getInput_inputPrimFind
|   - gets user input for primFind
| Input:
|   - numArgsSI:
|     o number of arguments/parameters in argAryStr
|   - argAryStr:
|     o array of c-strings with user input
|   - primFileStr:
|     o pointer to c-string to point to the primer file
|       argument
|   - primTypeSC:
|     o pointer to char to hold the primer file type
|       (tsv or fasta)
|   - readFileStr:
|     o pointer to c-string to point to file with reads
|   - fqBl:
|     o ponter to char to be set to
|       - 1 = readFileStr is a fastq file
|       - 0 = readFileStr is a fasta file
|   - minLenUI:
|     o pointer to unsigned long to hold minimum primer
|       length
|   - maxlenUI:
|     o pointer to unsigned long to hold maximum primer
|       length
|   - minPerScoreF:
|     o pointer to float to hold the minimum percent score
|       to count a primer as mapped (found by waterman)
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
getInput_inputPrimFind(
   int numArgsSI,
   char *argAryStr[],
   signed char **primFileStr,
   signed char *primTypeSC,
   signed char **readFileStr,
   signed char *fqBl,
   unsigned int *minLenUI,
   unsigned int *maxLenUI,
   float *minPercScoreF,
   signed char *fastBl,
   unsigned char *lenKmerUC,
   float *minPercKmerF
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - gets user input for primFind
   '   o fun03 sec01:
   '     - variable declerations
   '   o fun03 sec01:
   '     - get input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *tmpStr = 0;
   sint siArg = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - get input
   ^   o fun03 sec02 sub01:
   ^     - check if have input and start loop
   ^   o fun03 sec02 sub02:
   ^     - file input (primer sequences and reads)
   ^   o fun03 sec02 sub03:
   ^     - filtering input (length/score)
   ^   o fun03 sec02 sub04:
   ^     - fast alignment input
   ^   o fun03 sec02 sub05:
   ^     - help messages
   ^   o fun03 sec02 sub06:
   ^     - check invalid input and move on
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec02 Sub01:
   *   - check if have input and start loop
   \*****************************************************/

   if(numArgsSI < 2)
   { /*If: nothing was input; assume help message wanted*/
      phelp_inputPrimFind(stdout);
      return 1;
   } /*If: nothing was input; assume help message wanted*/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/

      /**************************************************\
      * Fun03 Sec02 Sub02:
      *   - file input (primer sequences and reads)
      \**************************************************/

      if(! eql_charCp("-prim-tsv", argAryStr[siArg], 0) ) 
      { /*If: a tsv primer file was input*/
         ++siArg; /*move to argument*/
         *primTypeSC = def_tsvPrimFile_inputPrimFind;
         *primFileStr = (schar *) argAryStr[siArg];
      } /*If: a tsv primer file was input*/

      else if(! eql_charCp("-prim-fa",argAryStr[siArg],0)) 
      { /*Else If: a fasta primer file was input*/
         ++siArg; /*move to argument*/
         *primTypeSC = def_faPrimFile_inputPrimFind;
         *primFileStr = (schar *) argAryStr[siArg];
      } /*Else If: a fasta primer file was input*/

      else if(! eql_charCp("-fa",argAryStr[siArg],0)) 
      { /*Else If: a fasta sequence file was input*/
         ++siArg; /*move to argument*/
         *fqBl = 0;
         *readFileStr = (schar *) argAryStr[siArg];
      } /*Else If: a fasta sequence file was input*/

      else if(! eql_charCp("-fq",argAryStr[siArg],0)) 
      { /*Else If: a fastq sequence file was input*/
         ++siArg; /*move to argument*/
         *fqBl = 1;
         *readFileStr = (schar *) argAryStr[siArg];
      } /*Else If: a fastq sequence file was input*/

      /**************************************************\
      * Fun03 Sec02 Sub03:
      *   - filtering input (length/score)
      \**************************************************/

      else if(! eql_charCp("-min-len",argAryStr[siArg],0)) 
      { /*Else If: the min length was input*/
         ++siArg; /*move to argument*/

         tmpStr =
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
      * Fun03 Sec02 Sub04:
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
      * Fun03 Sec02 Sub05:
      *   - help and version messages
      \**************************************************/

      else if(
            ! eql_charCp("-h",argAryStr[siArg],0)
         || ! eql_charCp("--h",argAryStr[siArg],0)
         || ! eql_charCp("help",argAryStr[siArg],0)
         || ! eql_charCp("-help",argAryStr[siArg],0)
         || ! eql_charCp("--help",argAryStr[siArg],0)
      ){ /*Else If: help message*/
         phelp_inputPrimFind(stdout);
         return 1;
      } /*Else If: help message*/

      else if(
            ! eql_charCp("-v",argAryStr[siArg],0)
         || ! eql_charCp("--v",argAryStr[siArg],0)
         || ! eql_charCp("version",argAryStr[siArg],0)
         || ! eql_charCp("-version",argAryStr[siArg],0)
         || ! eql_charCp("--version",argAryStr[siArg],0)
      ){ /*Else If: help message*/
         pversion_inputPrimFind(stdout);
         return 1;
      } /*Else If: help message*/

      /**************************************************\
      * Fun03 Sec02 Sub06:
      *   - check invalid input and move on
      \**************************************************/

      else
      { /*Else: input not recognized*/
         fprintf(
            stderr,
            "%s is not recognized\n",
            argAryStr[siArg]
         );

         return 2;
      } /*Else: input not recognized*/

      ++siArg;
   } /*Loop: get user input*/

   return 0;
} /*getInput_inputPrimFind*/

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
