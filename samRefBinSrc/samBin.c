/*#######################################################\
# Name: samBin
#   - holds functions for binning sam file reads by ref
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o .c st01: sqST_samBin
'     - holds one @SQ entry (reference) from a sam file
'     - hear to support refPHead_samBin (fun01) only, so
'       no support functions
'   o fun01: refPHead_samBin
'     - prints the header for the sam file to dedicate to
'       each reference in the sam file binning step
'   o fun02: refPRead_samBin
'     - prints the read into the correct references sam
'       file
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

#include <stdio.h>

#include "samBin.h"

#include "../generalLib/samEntry.h"

/*These have no .c files*/
#include "../generalLib/dataTypeShortHand.h"
#include "../generalLib/ulCp.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .h  #include "../generalLib/base10str.h"
!   o .h  #include "../generalLib/numToStr.h"
!   o .h  #include "../generalLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| ST01: sqST_samBin
|   - holds a single @SQ entry (reference) from a sam file
|   - this is only hear to support refPHead_samBin (fun01)
|     so there are no support functions
\-------------------------------------------------------*/
typedef struct sqST_samBin
{
   signed char sqStr[256]; /*entire @SQ entry*/
   signed char *refStr;    /*ref position in sqStr*/
   uint lenSQUI;           /*length of c-string in sqStr*/
   struct sqST_samBin *nextRef;/*next reference in list*/
}sqST_samBin;

/*-------------------------------------------------------\
| Fun01: refPHead_samBin
|   - prints the header for the sam file to dedicate to
|     each reference in the sam file binning step
| Input:
|   - samFILE:
|     o sam file to get headers from
|   - pgEntryStr:
|     o c-string with the program entry to add to the
|       header
|   - prefixStr:
|     o c-string with prefix for file name
|   - samSTPtr:
|     o pointer to samEntry structure to use in reading
|       samFILE
|   - buffStrPtr:
|     o pointer to c-string to use in getting each sam
|       file entry
|   - lenBuffULPtr:
|     o pointer to unsigned long with/to hold current size
|       of buffStrPtr
| Output:
|   - Prints:
|     o headers to each references sam file
|   - Modifies:
|     o samFILE to be on the second read entry
|     o samSTPtr to hold the first read entry
|     o buffStrPtr to resized (bigger) if needed
|     o lenBuffULPtr to have the updated size of
|       buffStrPtr if buffStrPtr is resized
|   - Returns:
|     o 0 for no errors
|     o def_memErr_samBin for memory errors
|     o def_fileErr_samBin for file errors
\-------------------------------------------------------*/
signed char
refPHead_samBin(
   void *samFILE,              /*sam to get header from*/
   signed char *pgEntryStr,    /*program entry to print*/
   signed char *prefixStr,     /*prefix for sam file*/
   struct samEntry *samSTPtr,  /*for getting headers*/
   signed char **buffStrPtr,   /*printing/read sam file*/
   unsigned long *lenBuffULPtr /*size of buffStrPtr*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   o fun01 sec01:
   '     - variable declarations
   '   o fun01 sec02:
   '     - get the pre @SQ (ref) entry comments
   '   o fun01 sec03:
   '     - get the references (@SQ entries)
   '   o fun01 sec04:
   '     - get the post @SQ (ref) entry comments
   '   o fun01 sec05:
   '     - print out headers to each reference file
   '   o fun01 sec06:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   schar outFileStr[255]; /*file to save to*/
   schar preSQStr[2048];      /*entries before @SQ*/
   schar postSQStr[4096];     /*entries after @SQ*/
   schar *tmpStr = 0;

   ulong tabUL = mkDelim_ulCp('\t');

   struct sqST_samBin *refHeapListST = 0;
   struct sqST_samBin *nextRefST = 0;

   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - get the pre @SQ (ref) entry comments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC =
      (schar)
      get_samLine(
         samSTPtr,
         (char **) buffStrPtr,
         lenBuffULPtr,
         (FILE *) samFILE
      ); /*get the first line from the file*/

   if(errSC == 64)
      goto memErr_fun01_sec06_sub02;

   if(errSC)
      goto fileErr_fun01_sec06_sub03;

   tmpStr = preSQStr;
   *tmpStr = '\0';

   while(! errSC)
   { /*Loop: get pre reference data*/
      if(samSTPtr->extraStr[0] != '@')
         break; /*finished*/

      if(
            samSTPtr->extraStr[1] == 'S'
         && samSTPtr->extraStr[2] == 'Q'
      ) break; /*found seq entry*/

      cpLen_ulCp(
         (char *) tmpStr,
         (char *) samSTPtr->extraStr,
         samSTPtr->lenExtraUI
      );

      tmpStr += samSTPtr->lenExtraUI;
      *tmpStr++ = '\n';
      *tmpStr = '\0';

      errSC =
         (schar)
         get_samLine(
            samSTPtr,
            (char **) buffStrPtr,
            lenBuffULPtr,
            (FILE *) samFILE
         ); /*get the first line from the file*/

   } /*Loop: get pre reference data*/

   if(errSC == 64)
      goto memErr_fun01_sec06_sub02;

   if(errSC)
      goto fileErr_fun01_sec06_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - get the references (@SQ entries)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   refHeapListST = malloc(sizeof(struct sqST_samBin));

   if(! refHeapListST)
      goto memErr_fun01_sec06_sub02;

   refHeapListST->nextRef = 0;

   while(! errSC)
   { /*Loop: get reference entries*/
      if(samSTPtr->extraStr[0] != '@')
         break; /*finished*/

      if(
            samSTPtr->extraStr[1] != 'S'
         || samSTPtr->extraStr[2] != 'Q'
      ) break; /*finished seq entry*/

      if(! nextRefST)
         nextRefST = refHeapListST; /*first round*/
      else
      { /*Else: need a new structure*/
         nextRefST->nextRef =
            malloc(sizeof(struct sqST_samBin));

         if(! nextRefST->nextRef)
            goto memErr_fun01_sec06_sub02;

         nextRefST = nextRefST->nextRef;
         nextRefST->nextRef = 0;
      } /*Else: need a new structure*/

      nextRefST->lenSQUI =
         (uint)
         cpDelim_ulCp(
            (char *) nextRefST->sqStr,
            (char *) samSTPtr->extraStr,
            0,
            '\0'
         );

      nextRefST->refStr = nextRefST->sqStr;

      while(
            nextRefST->refStr[0] != 'S'
         || nextRefST->refStr[1] != 'N'
         || nextRefST->refStr[2] != ':'
      ) ++nextRefST->refStr; /*find reference id start*/

      nextRefST->refStr += 3; /*get of SN:*/

      errSC =
         (schar)
         get_samLine(
            samSTPtr,
            (char **) buffStrPtr,
            lenBuffULPtr,
            (FILE *) samFILE
         ); /*get the first line from the file*/
   } /*Loop: get reference entries*/

   if(errSC == 64)
      goto memErr_fun01_sec06_sub02;

   if(errSC)
      goto fileErr_fun01_sec06_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - get the post @SQ (ref) entry comments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tmpStr = postSQStr;
   *tmpStr = '\0';

   while(! errSC)
   { /*Loop: get post reference data*/
      if(samSTPtr->extraStr[0] != '@')
         break; /*finished*/

      cpLen_ulCp(
         (char *) tmpStr,
         (char *) samSTPtr->extraStr,
         samSTPtr->lenExtraUI
      );

      tmpStr += samSTPtr->lenExtraUI;
      *tmpStr++ = '\n';
      *tmpStr = '\0';

      errSC =
         (schar)
         get_samLine(
            samSTPtr,
            (char **) buffStrPtr,
            lenBuffULPtr,
            (FILE *) samFILE
         ); /*get the first line from the file*/
   } /*Loop: get post reference data*/

   if(errSC == 64)
      goto memErr_fun01_sec06_sub02;

   if(errSC)
      goto fileErr_fun01_sec06_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec05:
   ^   - print out headers to each reference file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   nextRefST = refHeapListST; /*first round*/

   while(nextRefST)
   { /*Loop: print out headers for all references*/
      tmpStr = outFileStr;
      *tmpStr = '\0';

      tmpStr +=
         cpDelim_ulCp(
            (char *) tmpStr,        
            prefixStr,
            0,
            '\0'
         );

      tmpStr +=
         cpDelim_ulCp(
            (char *) tmpStr,        
            (char *) nextRefST->refStr,
            tabUL,
            '\t'
         );

      *tmpStr++ = '.';
      *tmpStr++ = 's';
      *tmpStr++ = 'a';
      *tmpStr++ = 'm';
      *tmpStr = '\0';

      outFILE =
         fopen(
            (char *) outFileStr,
            "w"
         );

      if(! outFILE)
         goto fileErr_fun01_sec06_sub03;

      fprintf(
         outFILE,
         "%s%s\n%s",
         preSQStr,
         nextRefST->sqStr,
         postSQStr
      );

      if(pgEntryStr)
         fprintf(
            outFILE,
            "%s",
            pgEntryStr
         ); /*if user supplied a program line*/

      fclose(outFILE);
      outFILE = 0;

      nextRefST = nextRefST->nextRef;
   } /*Loop: print out headers for all references*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec06:
   ^   - clean up and return
   ^   o fun01 sec06 sub01:
   ^     - clean up after no errors
   ^   o fun01 sec06 sub02:
   ^     - memory error clean up
   ^   o fun01 sec06 sub03:
   ^     - file error clean up
   ^   o fun01 sec06 sub04:
   ^     - general (no error/error) clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec06 Sub01:
   *   - clean up after no errors
   \*****************************************************/

   errSC = 0;
   goto cleanUp_fun01_sec06_sub04;

   /*****************************************************\
   * Fun01 Sec06 Sub02:
   *   - memory error clean up
   \*****************************************************/

   memErr_fun01_sec06_sub02:;
   errSC = def_memErr_samBin;
   goto cleanUp_fun01_sec06_sub04;

   /*****************************************************\
   * Fun01 Sec06 Sub03:
   *   - file error clean up
   \*****************************************************/

   fileErr_fun01_sec06_sub03:;
   errSC = def_fileErr_samBin;
   goto cleanUp_fun01_sec06_sub04;

   /*****************************************************\
   * Fun01 Sec06 Sub04:
   *   - general (no error/error) clean up
   \*****************************************************/

   cleanUp_fun01_sec06_sub04:;
   
   if(
         outFILE
      && outFILE != stdin
      && outFILE != stdout
   ) fclose(outFILE);

   outFILE = 0;

   while(refHeapListST)
   { /*Loop: free all sqST structures*/
      nextRefST = refHeapListST->nextRef;
      free(refHeapListST);
      refHeapListST = nextRefST;
   } /*Loop: free all sqST structures*/

   return errSC;
} /*refPHead_samBin*/

/*-------------------------------------------------------\
| Fun02: refPRead_samBin
|   - prints the read into the correct references sam
|     file
| Input:
|   - samSTPtr:
|     o pointer to samEntry structure with read to bin
|   - prefixStr:
|     o c-string with prefix for file name
|   - buffStrPtr:
|     o pointer to c-string to use as buffer in printing
|   - lenBuffULPtr:
|     o current length of buffer
|   - pNoNewLineBl:
|     o 1: do not print a new line (do not end entry)
|     o 0: print a new line (end entry)
| Output:
|   - Prints:
|     o read to its references sam file (prefixRef.sam)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_samBin for memory errors
|     o def_fileErr_samBin for file errors
\-------------------------------------------------------*/
signed char
refPRead_samBin(
   struct samEntry *samSTPtr,   /*read to bin*/
   signed char *prefixStr,      /*prefix for file name*/
   signed char **buffStrPtr,    /*buffer for printing*/
   unsigned long *lenBuffULPtr, /*size of buffer*/
   signed char pNoNewLineBl     /*1: do not end entry*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   '   o fun02 sec01:
   '     - variable declerations
   '   o fun02 sec02:
   '     - build file name
   '   o fun02 sec03:
   '     - open file and print read
   '   o fun02 sec04:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   schar outFileStr[255];
   schar *tmpStr = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - build file name
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tmpStr = outFileStr;

   tmpStr +=
      cpDelim_ulCp(
         (char *) tmpStr,        
         prefixStr,
         0,
         '\0'
      );

   tmpStr +=
      cpDelim_ulCp(
         (char *) tmpStr,
         (char *) samSTPtr->refIdStr,
         0,
         '\0'
      );

   *tmpStr++ = '.';
   *tmpStr++ = 's';
   *tmpStr++ = 'a';
   *tmpStr++ = 'm';
   *tmpStr++ = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - open file and print read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   outFILE =
      fopen(
         (char *) outFileStr,
         "a"
      );

   if(! outFILE)
      goto fileErr_fun02;

   errSC =
      p_samEntry(
         samSTPtr,
         (char **) buffStrPtr,
         lenBuffULPtr,
         (char) pNoNewLineBl,
         outFILE
      );

   if(errSC)
      goto memErr_fun02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun02;

   memErr_fun02:;
   errSC = def_memErr_samBin;
   goto ret_fun02;

   fileErr_fun02:;
   errSC = def_fileErr_samBin;
   goto ret_fun02;

   ret_fun02:;

   if(
         outFILE
      && outFILE != stdin
      && outFILE != stdout
   ) fclose(outFILE);

   outFILE = 0;

   return errSC;
} /*refPRead_samBin*/
   
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
