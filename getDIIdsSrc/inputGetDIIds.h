/*########################################################
# Name: inputGetDIIds
#   - holds user input, help message, and version
#     functions for primFind
#   - Also has some default variables and versio numbers
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - defined variables and guards
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
|   - defined variables and guards
\-------------------------------------------------------*/

#ifndef INPUT_PRIMER_FIND_H
#define INPUT_PRIMER_FIND_H

typedef struct fluST fluST;

#define def_year_inputGetDIIds 2024
#define def_month_inputGetDIIds 7
#define def_day_inputGetDIIds 11

#define def_faPrimFile_inputGetDIIds 1
#define def_tsvPrimFile_inputGetDIIds 2

#define def_minLen_inputGetDIIds 40
#define def_maxLen_inputGetDIIds 3000
#define def_minPercLen_inputGetDIIds 0.85f
#define def_maxPercLen_inputGetDIIds 1.1f /*110%*/

#define def_fastSearch_inputGetDIIds 1 /*1 = fast search*/

static signed char
   *forPrimStr_inputGetIds =
      (signed char *) "AGCGAAAGCAGG";
static signed char
   *revPrimStr_inputGetIds =
      (signed char *) "AGTAGAAACAAGG";

#define def_pDIRna_inputGetDIIds 1
#define def_pMVRna_inputGetDIIds 1
#define def_pVRna_inputGetDIIds 1
#define def_partSeq_inputGetDIIds 1

#define def_noAgree_inputGetDIIds 0
   /*primers in read not agreeing for segment*/

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
);

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
);

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
);

#endif

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
