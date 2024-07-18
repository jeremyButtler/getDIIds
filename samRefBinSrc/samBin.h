/*#######################################################\
# Name: samBin
#   - holds functions for binning sam file reads by ref
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - guards, forward declerations, & defined variables
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
|   - guards, forward declerations, and defined variables
\-------------------------------------------------------*/

#ifndef SAM_FILE_BIN_READS_H
#define SAM_FILE_BIN_READS_H

typedef struct samEntry samEntry;

#define def_memErr_samBin 2
#define def_fileErr_samBin 4

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
);

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
);

#endif

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
