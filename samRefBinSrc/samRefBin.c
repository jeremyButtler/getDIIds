/*#######################################################\
# Name: samRefBin
#   - bins reads in a sam file by reference
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries and default variables
'   o fun01: pversion_samRefBin
'     - prints the version number for samRefBin
'   o fun02: phelp_samRefBin
'     - prints the help message for samRefBin
'   o fun03: input_samRefBin
'     - gets user input
'   o main:
'     - driver fuction to bin sam file reads by reference
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and default variables
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../generalLib/samEntry.h"
#include "samBin.h"

/*.h files only*/
#include "../generalLib/dataTypeShortHand.h"
#include "../generalLib/ulCp.h"
#include "../generalLib/charCp.h"
#include "../generalLib/numToStr.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .h  #include "../generalLib/dataTypeShortHand.h"
!   o .h  #include "../generalLib/base10str.h"
!   o .h  #include "../generalLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_year_samRefBin 2024
#define def_month_samRefBin 7
#define def_day_samRefBin 16

#define def_phelp_samRefBin 1
#define def_pversion_samRefBin 1
#define def_inputErr_samRefBin 2

schar *defPrefixStr = (schar *) "out"; /*default prefix*/

/*-------------------------------------------------------\
| Fun01: pversion_samRefBin
|   - prints the version number for samRefBin
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints:
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_samRefBin(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "samRefBin version: %i-%02i-%02i\n",
      def_year_samRefBin,
      def_month_samRefBin,
      def_day_samRefBin
   );
} /*pversion_samRefBin*/

/*-------------------------------------------------------\
| Fun02: phelp_samRefBin
|   - prints the help message for samRefBin
| Input:
|   - outFILE:
|     o file to print help message to
| Output:
|   - Prints:
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_samRefBin(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints the help message for samRefBin
   '   o fun02 sec01:
   '     - usage block
   '   o fun02 sec02:
   '     - input block
   '   o fun02 sec03:
   '     - output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - usage block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "samRefBin -sam file.sam -prefix out\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - splits (bins) sam file reads by reference\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - input block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Input:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -sam file.sam: [Required; stdin]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o sam file to bin\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o use \"-\" for stdin input\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -prefix %s: [Optional; %s]\n",
      defPrefixStr,
      defPrefixStr
   );

   fprintf(
      (FILE *) outFILE,
      "    o prefix to add to binned file names\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - output block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Output:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - prints reads to prefix-mapped_reference.sam\n"
   );
} /*phelp_samRefBin*/

/*-------------------------------------------------------\
| Fun03: input_samRefBin
|   - gets user input
| Input:
|   - numArgsSI:
|     o number of arguments input
|   - argAryStr:
|     o c-string array with user input
|   - samFileStrPtr:
|     o pointer to c-string to point to the name of the
|       sam file to bin
|   - prefixStrPtr:
|     o pointer to c-string to point to prefix
| Output:
|   - Modifies:
|     o samFileStrPtr and prefixStrPtr to hold user input
|   - Returns:
|     o 0 for no errors
|     o def_phelp_samRefBin if printed the help message
|     o def_pversion_samRefBin if printed the version #
|     o def_inputErr_samRefBin if had an error
\-------------------------------------------------------*/
signed char
input_samRefBin(
   int numArgsSI,         /*# argumnets input*/
   char *argAryStr[],     /*arguments input*/
   schar **samFileStrPtr, /*will hold path to sam file*/
   schar **prefixStrPtr   /*will hold prefix*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - gets user input
   '   o fun03 sec01:
   '     - variable declaration
   '   o fun03 sec02:
   '     - check number arguments
   '   o fun03 sec03:
   '     - get user input
   '   o fun03 sec04:
   '     - set error type and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declaration
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0; /*return error value*/
   sint siArg = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - check number arguments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! numArgsSI)
   { /*If: no argmuments input*/
      phelp_samRefBin(stdout);
      goto phelp_fun03;
   } /*If: no argmuments input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - get user input
   ^   o fun03 sec03 sub01:
   ^     - start loop and get settings input
   ^   o fun03 sec03 sub02:
   ^     - check for help message requests
   ^   o fun03 sec03 sub03:
   ^     - check for version number requests
   ^   o fun03 sec03 sub04:
   ^     - unkown input
   ^   o fun03 sec03 sub05:
   ^     - move to next argument
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - start loop and get settings input
   \*****************************************************/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/
      if(! eql_charCp("-sam", argAryStr[siArg], 0))
      { /*If: sam file input*/
         ++siArg;
         *samFileStrPtr = (schar *) argAryStr[siArg];
      } /*If: sam file input*/

      else if(! eql_charCp("-prefix",argAryStr[siArg],0))
      { /*Else If: prefix input*/
         ++siArg;
         *prefixStrPtr = (schar *) argAryStr[siArg];
      } /*Else If: prefix input*/

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - check for help message requests
      \**************************************************/

      else if(! eql_charCp("-h",argAryStr[siArg],0))
      { /*Else If: help message requested*/
         phelp_samRefBin(stdout);
         goto phelp_fun03;
      } /*Else If: help message requested*/

      else if(! eql_charCp("--h",argAryStr[siArg],0))
      { /*Else If: help message requested*/
         phelp_samRefBin(stdout);
         goto phelp_fun03;
      } /*Else If: help message requested*/

      else if(! eql_charCp("help",argAryStr[siArg],0))
      { /*Else If: help message requested*/
         phelp_samRefBin(stdout);
         goto phelp_fun03;
      } /*Else If: help message requested*/

      else if(! eql_charCp("-help",argAryStr[siArg],0))
      { /*Else If: help message requested*/
         phelp_samRefBin(stdout);
         goto phelp_fun03;
      } /*Else If: help message requested*/

      else if(! eql_charCp("--help",argAryStr[siArg],0))
      { /*Else If: help message requested*/
         phelp_samRefBin(stdout);
         goto phelp_fun03;
      } /*Else If: help message requested*/

      /**************************************************\
      * Fun03 Sec03 Sub03:
      *   - check for version number requests
      \**************************************************/

      else if(! eql_charCp("-v",argAryStr[siArg],0))
      { /*Else If: version number requested*/
         pversion_samRefBin(stdout);
         goto pversion_fun03;
      } /*Else If: version number requested*/

      else if(! eql_charCp("--v",argAryStr[siArg],0))
      { /*Else If: version number requested*/
         pversion_samRefBin(stdout);
         goto pversion_fun03;
      } /*Else If: version number requested*/

      else if(! eql_charCp("version",argAryStr[siArg],0))
      { /*Else If: version number requested*/
         pversion_samRefBin(stdout);
         goto pversion_fun03;
      } /*Else If: version number requested*/

      else if(! eql_charCp("-version",argAryStr[siArg],0))
      { /*Else If: version number requested*/
         pversion_samRefBin(stdout);
         goto pversion_fun03;
      } /*Else If: version number requested*/

      else if(!eql_charCp("--version",argAryStr[siArg],0))
      { /*Else If: version number requested*/
         pversion_samRefBin(stdout);
         goto pversion_fun03;
      } /*Else If: version number requested*/

      /**************************************************\
      * Fun03 Sec03 Sub04:
      *   - unkown input
      \**************************************************/

      else
         goto err_fun03;

      /**************************************************\
      * Fun03 Sec03 Sub05:
      *   - move to next argument
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - set error type and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun03;

   phelp_fun03:;
   errSC = def_phelp_samRefBin;
   goto ret_fun03;

   pversion_fun03:;
   errSC = def_pversion_samRefBin;
   goto ret_fun03;

   err_fun03:;
   errSC = def_inputErr_samRefBin;
   goto ret_fun03;

   ret_fun03:;

   return errSC;
} /*input_samRefBin*/

/*-------------------------------------------------------\
| Main:
|   - driver fuction to bin sam file reads by reference
| Input:
|   - numArgsSI:
|     o number of arguments input
|   - argAryStr:
|     o c-string array with user input
| Output:
|   - Prints:
|     o reads to their mapped references sam file
|       (prefix-reference.sam)
|     o errors to stderr
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
   '     - initialize, get input, and open files
   '   o main sec03:
   '     - build program header entry
   '   o main sec04:
   '     - print headers and bin reads
   '   o main sec05:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   schar *samFileStr = 0;
   schar *prefixStr = defPrefixStr;

   schar *tmpStr = 0;
   schar pgEntryStr[2048]; /*sam file program line*/

   struct samEntry samStackST;
   schar *buffHeapStr = 0;
   ulong lenBuffUL = 0;
   ulong entryUL = 0;

   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - initialize, get input, and open files
   ^   o main sec02 sub01:
   ^     - initialize and get input
   ^   o main sec02 sub02:
   ^     - open sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize and get input
   \*****************************************************/

   init_samEntry(&samStackST);

   errSC =
      input_samRefBin(
         numArgsSI,
         argAryStr,
         &samFileStr,
         &prefixStr
      );

   if(errSC)
   { /*If: had an error*/
      --errSC; /*convert help/version requests to 0*/
      goto cleanUp_main_sec05_sub04;
   } /*If: had an error*/

   errSC = setup_samEntry(&samStackST);

   if(errSC)
      goto memErr_main_sec05_sub02;

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - open sam file
   \*****************************************************/

   if(
         ! samFileStr
      || *samFileStr == '\0'
   ) samFILE = stdin;

   else
   { /*Else: user suppplied a sam file*/
      samFILE =
         fopen(
            (char *) samFileStr,
            "r"
         );

      if(! samFILE)
      { /*If: could not open sam file*/
         fprintf(
            stderr,
            "could not open -sam %s\n",
            samFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: could not open sam file*/
   } /*Else: user suppplied a sam file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - build program header entry
   ^   o main sec03 sub01:
   ^     - copy first general part
   ^   o main sec03 sub02:
   ^     - copy the version number
   ^   o main sec03 sub03:
   ^     - copy the users command
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub01:
   *   - copy first general part
   \*****************************************************/

   tmpStr = pgEntryStr;

   tmpStr +=
      cpDelim_ulCp(
         (char *) tmpStr,
         "@PG\tID:samRefBin\tPN:samRefBin\tVN:",
         0,
         '\0'
      );
         
   /*****************************************************\
   * Main Sec03 Sub02:
   *   - copy the version number
   \*****************************************************/

   tmpStr +=
      numToStr(
         (char *) tmpStr,
         def_year_samRefBin
      );

   *tmpStr++ = '-';

   tmpStr +=
      numToStr(
         (char *) tmpStr,
         def_month_samRefBin
      );

   if(*(tmpStr - 2) == '-')
   { /*If: single digit*/
      *tmpStr = *(tmpStr - 1);
      *(tmpStr - 1) = '0';
      ++tmpStr;
      *tmpStr = '\0';
   } /*If: single digit*/

   *tmpStr++ = '-';

   tmpStr +=
      numToStr(
         (char *) tmpStr,
         def_day_samRefBin
      );

   if(*(tmpStr - 2) == '-')
   { /*If: single digit*/
      *tmpStr = *(tmpStr - 1);
      *(tmpStr - 1) = '0';
      ++tmpStr;
      *tmpStr = '\0';
   } /*If: single digit*/

   /*****************************************************\
   * Main Sec03 Sub03:
   *   - copy the users command
   \*****************************************************/

   tmpStr +=
      cpDelim_ulCp(
         (char *) tmpStr,
         "\tCL:samRefBin -prefix ",
         0,
         '\0'
      );

   tmpStr +=
      cpDelim_ulCp(
         (char *) tmpStr,
         (char *) prefixStr,
         0,
         '\0'
      );

   *tmpStr++ = ' ';
   *tmpStr++ = '-';
   *tmpStr++ = 's';
   *tmpStr++ = 'a';
   *tmpStr++ = 'm';
   *tmpStr++ = ' ';

   if(! samFileStr)
      *tmpStr++ = '-';

   else
   { /*Else: sam file was input*/
      tmpStr +=
         cpDelim_ulCp(
            (char *) tmpStr,
            (char *) samFileStr,
            0,
            '\0'
         );
   } /*Else: sam file was input*/

   *tmpStr++ = '\n';
   *tmpStr++ = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - print headers and bin reads
   ^   o main sec04 sub01:
   ^     - print headers to bins (reference sam files)
   ^   o main sec04 sub02:
   ^     - print reads to bins (reference sam files)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - print headers to bins (reference sam files)
   \*****************************************************/

   errSC =
      refPHead_samBin(
         samFILE,
         pgEntryStr,
         prefixStr,
         &samStackST,
         &buffHeapStr,
         &lenBuffUL
      ); /*print out headers to each references sam file*/

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_samBin)
      { /*If: memory error*/
         fprintf(
           stderr,
           "MEMORY ERROR while printing headers to bins\n"
         );

         goto memErr_main_sec05_sub02;
      } /*If: memory error*/

      fprintf(
         stderr,
         "could not make bins for -sam %s\n",
         samFileStr
      );

      goto fileErr_main_sec05_sub03;
   } /*If: had error*/

   /*****************************************************\
   * Main Sec04 Sub02:
   *   - print reads to bins (reference sam files)
   \*****************************************************/

   while(! errSC)
   { /*Loop: bin reads*/
      ++entryUL;

      errSC =
         refPRead_samBin(
            &samStackST,
            prefixStr,
            &buffHeapStr,
            &lenBuffUL,
            0            /*not adding extra entries*/
         ); /*bin the reads*/

      if(errSC)
      { /*If: had error*/
         if(errSC == def_memErr_samBin)
         { /*If: memory error*/
            fprintf(
               stderr,
               "MEMORY ERROR while binning reads\n"
            );

            goto memErr_main_sec05_sub02;
         } /*If: memory error*/

         fprintf(
           stderr,
           "could not open reference file for read %lu\n",
           entryUL
         );

         fprintf(
            stderr,
            "  in -sam %s\n",
            samFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: had error*/

      errSC =
         (schar)
         get_samLine(
            &samStackST,
            (char **) &buffHeapStr,
            &lenBuffUL,
            samFILE
         ); /*get the next samfile entry*/
   } /*Loop: bin reads*/

   /*****************************************************\
   * Main Sec04 Sub03:
   *   - handle errors from reading sam file
   \*****************************************************/

   if(errSC == 64)
   { /*If: had a memory error*/
      fprintf(
         stderr,
         "MEMORY error getting entry %lu in -sam %s\n",
         entryUL,
         samFileStr
      );

      goto memErr_main_sec05_sub02;
   } /*If: had a memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05:
   ^   - clean up
   ^   o main sec05 sub01:
   ^     - no error clean up
   ^   o main sec05 sub02:
   ^     - memory error clean up
   ^   o main sec05 sub03:
   ^     - file error clean up
   ^   o main sec05 sub04:
   ^     - general clean up (always called)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec05 Sub01:
   *   - no error clean up
   \*****************************************************/

   errSC = 0;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub02:
   *   - memory error clean up
   \*****************************************************/

   memErr_main_sec05_sub02:;
   errSC = 1;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub03:
   *   - file error clean up
   \*****************************************************/

   fileErr_main_sec05_sub03:;
   errSC = 2;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub04:
   *   - general clean up (always called)
   \*****************************************************/

   cleanUp_main_sec05_sub04:;

   if(
         samFILE
      && samFILE != stdin
      && samFILE != stdout
   ) fclose(samFILE);

   samFILE = 0;

   free(buffHeapStr);
   buffHeapStr = 0;

   freeStack_samEntry(&samStackST);

   exit(errSC);
} /*main*/
