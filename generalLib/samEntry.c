/*########################################################
# Name: samEntryStruct
#  - Holds structer to hold a sam file entry. This also
#    includes the functions needed to support this
#    structure.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - Included libraries
'  o .h st01 samEntry:
'    - Holds a single samfile entry
'  o .h fun01 blank_samEntry:
'    - Sets all non-alloacted variables in samEntryST to 0
'  o .h fun02 init_samEntry:
'    - Initalize a samEntry struct to 0's
'  o fun03: setup_samEntry
'    - allocates memory for a samEntry structure (call
'      after init)
'  o fun04 freeStack_samEntry:
'    - Frees heap allocations in a stack allocated
'      samEntry struct
'  o fun05 freeHeap_samEntry:
'    - Frees a samEntry structer (and sets to null)
'  o fun06: mk_samEntry
'    - Makes an heap allocated samEntry structure
'  o .h fun07: qhistToMed_samEntry
'    - Gets the median q-score for an histogram of
'      q-scores in a samStruct
'  o fun08: findQScores_samEntry
'     - Gets the median and mean q-scores from a samEntry
'       Structure.
'  o fun09: cpQEntry_samEntry
'    - Copies q-scores from a string into a samEntry
'      structure
'  o fun10: get_samLine
'    - Reads in a single line from a sam file
'  o .h fun11: findRefPos_samEntry
'    - Find an reference coordinate in an sequence in
'      an sam entry structure
'  o fun12: p_samEntry
'    - Prints the sam file entry to a file. This does not
'      print any extra stats that were found.
'  o fun13: pfq_samEntry
'    - Prints the sam entry as a fastq entry to a fastq
'      file
'  o fun14: pfa_samEntry
'    - Prints the sam entry as a fasta entry to a fasta
'      file
'  o fun15: pstats_samEntry
'    - Prints out the stats in a samEntry struct to a file
'  o .h note01:
'     - Notes about the sam file format from the sam file
'       pdf
'  o license:
'    - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "samEntry.h"

#include <stdio.h>

/*These have no .c files*/
#include "dataTypeShortHand.h"
#include "base10str.h"
#include "ulCp.h"
#include "numToStr.h"
#include "ntTo5Bit.h"
   /*look up table to see if have anonymous bases*/

#define def_newLine_samEntry mkDelim_ulCp('\n')
#define def_tab_samEntry mkDelim_ulCp('\t')
#define def_one_samEntry mkDelim_ulCp(0x01)
#define def_highBit_samEntry mkDelim_ulCp(0x80)

/*-------------------------------------------------------\
| Fun03: setup_samEntry
|  - allocates memory for a samEntry structure (call after
|    init)
| Input:
|  - samSTPtr:
|    o pointer to samEntry struct to allocate memory to
| Output:
|  - Allocates:
|    o memory for seqStr, qStr, cigTypStr, cigValAryI,
|      and extraStr
|  - Returns:
|    o 0 for success
|    o 1 for memory error
\-------------------------------------------------------*/
unsigned char
setup_samEntry(
   struct samEntry *samSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03: setup_samEntry
   '   o fun03 sec01:
   '     - allocate sequence memory
   '   o fun03 sec02:
   '     - allocate q-score memory
   '   o fun03 sec03:
   '     - allocate cigar type (snp/match/ins/del) memory
   '   o fun03 sec04:
   '     - allocate cigar count array
   '   o fun03 sec05:
   '     - allocate extra array
   '   o fun03 sec06:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec01:
    ^   - allocate sequence memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samSTPtr->seqStr)
       free(samSTPtr->seqStr);
  
    samSTPtr->seqStr = 0;

    samSTPtr->seqStr =
       calloc(
          1024,
          sizeof(char)
       ); /*allocate some memory for sequence*/

    if(! samSTPtr->seqStr)
       goto memErr_fun03_sec06;

    samSTPtr->lenSeqBuffUI = 1024;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec02:
    ^   - allocate q-score memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->qStr)
       free(samSTPtr->qStr);
  
    samSTPtr->qStr = 0;

    samSTPtr->qStr =
       calloc(
          1024,
          sizeof(char)
       ); /*allocate some memory for q-score*/

    if(! samSTPtr->qStr)
       goto memErr_fun03_sec06;

    samSTPtr->lenQBuffUI = 1024;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec03:
    ^   - allocate cigar type (snp/match/ins/del) memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->cigTypeStr)
       free(samSTPtr->cigTypeStr);
  
    samSTPtr->cigTypeStr = 0;

    samSTPtr->cigTypeStr =
       calloc(
          256,
          sizeof(char)
       ); /*allocate memory for cigary type*/

    if(! samSTPtr->cigTypeStr)
       goto memErr_fun03_sec06;

    samSTPtr->lenCigBuffUI = 256;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec04:
    ^   - allocate cigar count array
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->cigValAryI)
       free(samSTPtr->cigValAryI);
  
    samSTPtr->cigValAryI = 0;

    samSTPtr->cigValAryI =
       calloc(
          256,
          sizeof(sint)
       ); /*allocate memory for cigar counts*/

    if(! samSTPtr->cigValAryI)
       goto memErr_fun03_sec06;
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec05:
    ^   - allocate extra array
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samSTPtr->extraStr)
       free(samSTPtr->extraStr);
  
    samSTPtr->extraStr = 0;

    samSTPtr->extraStr =
       calloc(
          1024,
          sizeof(char)
       ); /*allocate memory for extra entry*/

    if(! samSTPtr->extraStr)
       goto memErr_fun03_sec06;

    samSTPtr->lenExtraBuffUI = 1024;
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec06:
    ^   - return
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    return 0;

    memErr_fun03_sec06:;
    return 1;
} /*setup_samEntry*/

/*-------------------------------------------------------\
| Fun04: freeStack_samEntry
| Use:
|  - Frees all variables in samEntry, but not samEntry
| Input:
|  - samSTPtr
|    o Pointer to samEntry struct to free the interal
|      memory of
| Output:
|  - Frees:
|    o allocated memory in samSTPtr
\-------------------------------------------------------*/
void
freeStack_samEntry(
   struct samEntry *samSTPtr
){
    blank_samEntry((samSTPtr));
    
    free((samSTPtr)->seqStr);
    (samSTPtr)->seqStr = 0;
    (samSTPtr)->lenSeqBuffUI = 0;
    
    free((samSTPtr)->qStr);
    (samSTPtr)->qStr = 0;
    (samSTPtr)->lenQBuffUI = 0;
    
    free((samSTPtr)->cigTypeStr);
    (samSTPtr)->cigTypeStr = 0;
    (samSTPtr)->lenCigBuffUI = 0;
    
    free((samSTPtr)->cigValAryI);
    (samSTPtr)->cigValAryI = 0;
    
    free((samSTPtr)->extraStr);
    (samSTPtr)->extraStr = 0;
    (samSTPtr)->lenExtraBuffUI = 0;
} /*freeStack_samEntry*/

/*-------------------------------------------------------\
| Fun05: freeHeap_samEntry
|  - Frees a samEntry struct
| Input:
|  - samSTPtr
|    o Pointer to Sam entry struct to free
| Output:
|  - Frees:
|    o samSTPtr and its variables with allocated memory
\-------------------------------------------------------*/
void freeHeap_samEntry(
   struct samEntry **samSTPtr
){
    freeStack_samEntry((*samSTPtr));
    free((*samSTPtr));
    (*samSTPtr) = 0;
} /*freeHeap_samEntry*/

/*-------------------------------------------------------\
| Fun06: mk_samEntry
|  - Makes an heap allocated samEntry structure
| Input:
| Output:
|  - Returns:
|    o Pointer to a samEntry structure
|    o 0 if had an memory error
\-------------------------------------------------------*/
struct samEntry *
mk_samEntry(
){
  struct samEntry *samST=malloc(sizeof(struct samEntry));
  uchar errUC = 0;
  
  if(samST)
  { /*If: I did not have a memory error*/
     init_samEntry(samST);
     errUC = setup_samEntry(samST);

     if(errUC)
        freeHeap_samEntry(&samST);
  } /*If: I did not have a memory error*/
  
  return samST;
} /*makeSamEntry*/

/*-------------------------------------------------------\
| Fun08: findQScores_samEntry
|   - Gets the median and mean q-scores from a samEntry
|     Structure.
| Input:
|   - samSTPTr:
|     o Pointer to samEntry struct to find the median and
|       mean q-scores for
| Output:
|   - Modifies:
|     o samSTPtr->medianQF to have the median q-score
|     o samSTPtr->meanQF to have the mean q-score
|     o samSTPtr->qHistUI to have a histogram of q-scores
|     o samSTPtr->sumQUL to have the sum of all q-scores
\-------------------------------------------------------*/
void
findQScores_samEntry(
   struct samEntry *samSTPtr
){
    ulong qAdjustUL = mkDelim_ulCp(def_ajdQ_samEntry);
    ulong qScoresUL = 0;
    ulong *qPtrUL = (ulong *) (samSTPtr)->qStr;

    uchar *scoreAryUC = 0;
    uint uiQScore = 0;
    uint uiChar = 0;
    
    /*Find the number of q-score characters in buffer*/
    for(
       uiQScore = 0;
       uiQScore <
          ((samSTPtr)->readLenUI >> def_shiftULBy_ulCp);
       ++uiQScore
    ) { /*Loop: Update the q-score historgram and sum*/
       qScoresUL = qPtrUL[uiQScore] - qAdjustUL;

       scoreAryUC =
          (uchar *)
          qScoresUL;
       
       for(
          uiChar = 0;
          uiChar < def_charInUL_ulCp;
          ++uiChar
       ){ /*Loop: Get the q-score entries*/
         ++(samSTPtr)->qHistUI[scoreAryUC[uiChar]];
         (samSTPtr)->sumQUL += (uchar) scoreAryUC[uiChar];
       } /*Loop: Get the q-score entries*/
    } /*Loop: Update the q-score historgram and sum*/
    
    uiQScore = (samSTPtr)->readLenUI;
    scoreAryUC = (uchar *) (samSTPtr)->qStr;
    
    for(
       uiQScore -=
          ((samSTPtr)->readLenUI & def_modUL_ulCp);
       uiQScore < (samSTPtr)->readLenUI;
       ++uiQScore
    ) { /*Loop: Copy the q-score entries*/
       (samSTPtr)->qStr[uiQScore] = scoreAryUC[uiQScore];

        uiChar = scoreAryUC[uiQScore] - def_ajdQ_samEntry;
       
       ++(samSTPtr)->qHistUI[uiChar];
       (samSTPtr)->sumQUL += uiChar;
    } /*Loop: Copy the q-score entries*/
    
    (samSTPtr)->meanQF =
       (float) (samSTPtr)->sumQUL /(samSTPtr)->readLenUI;
    
    qhistToMed_samEntry((samSTPtr));
} /*findQScores_samEntry*/

/*-------------------------------------------------------\
| Fun09: cpQEntry_samEntry
|   - Copies q-scores from a string into a samEntry
|     structure
| Input:
|   - samSTPtr:
|     o Pionter to sam entry struct to copys q-scores to
|   _ cpQStr:
|     o C-string with q-scores to copy to samSTPtr
|   - blankQHistBl:
|     o 1: Blank q-score vars (histogram/sum/mean/median)
|     o 0: do not blank the q-score variables
| Output:
|   - Mofidies:
|     o qStr in samSTPtry to have the q-scores
|     o medianQF in samSTPtr to have the median q-score
|     o meanQF in samSTPtr to have the mean q-score
|     o qHistUI in samSTPtr to have histogram of q-scores
|     o samQUL in samSTPtr to have sum off all q-scores
|   - Returns
|     o The value in samSTPtr->readLenUI
\-------------------------------------------------------*/
int
cpQEntry_samEntry(
   struct samEntry *samSTPtr, /*Copy q-scores to*/
   char *cpQStr,              /*q-scores to copy*/
   char blankQHistBl          /*1: to blank q-score hist*/
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ' Fun09 TOC: cpQEntry_samEntry
  '   o fun09 sec01:
  '     - Variable declerations
  '   o fun09 sec02:
  '     - Check and if asked blank the q-score values
  '   o fun09 sec03:
  '     - Copy q-scores using unsigned longs
  '   o fun09 sec04:
  '     - Copy last q-scores I could not copy with
  '       unsigned longs
  '   o fun09 sec05:
  '     - Find the median and mean q-scores
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec01:
  ^   - Variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  char *tmpStr = 0;
  uint uiQ = 0;
  uint uiChar = 0;
  ulong qAdjustUL = mkDelim_ulCp(def_ajdQ_samEntry);
  ulong *cpPtrUL = (ulong *) (cpQStr);
  ulong *dupPtrUL = (ulong *) samSTPtr->qStr;
  ulong qScoreUL = 0;
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec02:
  ^   - Check and if asked blank the q-score values
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  if(blankQHistBl)
  { /*If: I need to blank the q-score histograms*/
    for(uiChar = 0; uiChar < def_maxQ_samEntry; ++uiChar)
       samSTPtr->qHistUI[uiChar] = 0;

    samSTPtr->sumQUL = 0;
    samSTPtr->medianQF = 0;
    samSTPtr->meanQF = 0;
  } /*If: I need to blank the q-score histograms*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec03:
  ^   - Copy q-scores using unsigned longs
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  for(
     uiQ = 0;
     uiQ < ((samSTPtr)->readLenUI >> def_shiftULBy_ulCp);
     ++uiQ
  ) { /*Loop: Copy the q-score entries*/
     dupPtrUL[uiQ] = cpPtrUL[uiQ];
     qScoreUL = dupPtrUL[uiQ] - qAdjustUL;
     tmpStr = (char *) &qScoreUL;
     
     for(
        uiChar = 0;
        uiChar < def_charInUL_ulCp;
        ++uiChar
     ) { /*Loop: Get the q-score entries*/
        ++(samSTPtr)->qHistUI[ (uchar) tmpStr[uiChar] ];
        (samSTPtr)->sumQUL += (uchar) tmpStr[uiChar];
     } /*Loop: Get the q-score entries*/
  } /*Loop: Copy the q-score entries*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec04:
  ^   - Copy last q-scores I could not copy with unsgined
  ^     longs
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  uiQ = (samSTPtr)->readLenUI;
  
  for(
     uiQ -= ((samSTPtr)->readLenUI & def_modUL_ulCp);
     uiQ < (samSTPtr)->readLenUI;
     ++uiQ
  ) { /*Loop: Copy the q-score entries*/
     (samSTPtr)->qStr[uiQ] = (cpQStr)[uiQ];
     qScoreUL = (uchar) (cpQStr)[uiQ] - def_ajdQ_samEntry;

     ++(samSTPtr)->qHistUI[qScoreUL];
     (samSTPtr)->sumQUL += qScoreUL;
  } /*Loop: Copy the q-score entries*/
  
  (samSTPtr)->qStr[uiQ] = '\0';
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec05:
  ^   - Find the median and mean q-scores
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  (samSTPtr)->meanQF =
     (float)(samSTPtr)->sumQUL/(samSTPtr)->readLenUI;
  
  qhistToMed_samEntry((samSTPtr));
  return uiQ;
} /*cpQEntry_samEntry*/

/*-------------------------------------------------------\
| Fun10: get_samLine
|  - Reads in a single line from a sam file
| Input:
|  - samSTPtr:
|    o Pionter to samEntry structure to store the sam file
|      line in. This structure should be initialized
|  - buffStr:
|    o Buffer to read the sam file line temporarly into.
|      This is resized if needed. You can input NULL to
|      create a new buffer.
|  - lenBuffUL:
|    o Length of buffStr (updated if buffStr is resized)
|  - samVoidFILE:
|    o Sam file to read a line from. The void is so that
|      I can use samFILE in the function.
| Output:
|  - Modifies:
|    o samSTPtr to have the next line
|      - Comments are in extraStr
|    o samFILE to be on the next line
|    o buffStr to hold a sam file line (resized if needed)
|    o lenBuffUL to hold the resized length of buffStr
|  - Returns:
|    o 0 for success
|    o 1 for EOF (End Of File)
|    o 64 for memory errors
\-------------------------------------------------------*/
char
get_samLine(
   struct samEntry *samSTPtr,
   char **buffStr,
   ulong *lenBuffUL,
   void *samVoidFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC: get_samLine
   '   - Reads a single line from a sam file into samSTPtr
   '   o fun10 sec01:
   '     - Variable declerations
   '   o fun10 sec02:
   '     - Read in one line from the sam file
   '   o fun10 sec03:
   '     - Get the query id from the buffer
   '   o fun10 sec04:
   '     - Get the flag
   '   o fun10 sec05:
   '     - REad in the reference name/id
   '   o fun10 sec06:
   '     - Get the reference position
   '   o fun10 sec07:
   '     - Get the mapping quality
   '   o fun10 sec08:
   '     - Get the cigar entry
   '   o fun10 sec09:
   '     - Get the RNEXT entry
   '   o fun10 sec10:
   '     - Get the PNEXT entry
   '   o fun10 sec11:
   '     - Get the TLEN entry
   '   o fun10 sec12:
   '     - Get the sequence entry
   '   o fun10 sec13:
   '     - Get the q-score entry
   '   o fun10 sec14:
   '     - Copy the extra entry; after strict sam entries
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01:
   ^   - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ushort extraBuffUS = 4096;

   ulong oldLenUL = 0;
   ulong *ulStr = 0;

   char *tmpStr = 0;
   char *iterStr = 0;

   FILE *samFILE = (FILE *) samVoidFILE;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^   - Read in one line from the sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_samEntry(samSTPtr);

   if(*buffStr == 0)
   { /*If: I need to create a buffer*/
      *lenBuffUL = 0;
      *buffStr = malloc(extraBuffUS * sizeof(char));

      if(*buffStr == 0) return 64;

      *lenBuffUL = extraBuffUS;
   } /*If: I need to create a buffer*/

   /*I originally tried reallocs with fgets by setting
   `   the *lenBuffUL - 2 charater to null. But this did
   `   not always work. Aparatenly fgets is a lttle buggy,
   `   at least on void linux. However, it was faster to 
   `   use
   */
   
   tmpStr = fgets(*buffStr, *lenBuffUL, samFILE);

   if(! tmpStr ) return 1; /*EOF*/

   goto samEntry_Fun_12_Sec_02_checkEOL;

   while(*tmpStr != '\n')
   { /*Loop: Find the length of  the line*/
      *lenBuffUL <<= 1;
        /*This is a little agressive for memory usage, but
        `   it is fast*/

      tmpStr = malloc((*lenBuffUL + 1) * sizeof(char));

      /*Check for memory errors; let user handle
      `   freeing buffStr when have memory errors
      */
      if(! tmpStr) return 64;

      /*This avoids odd memory issues with realloc*/
      oldLenUL =
         cpDelim_ulCp(
            tmpStr,
            (char *) *buffStr,
            0,
            '\0'
         ); /*copy old buffer*/

      free(*buffStr);
      *buffStr = tmpStr;
      
      tmpStr = *buffStr + oldLenUL;
      tmpStr = fgets(tmpStr, *lenBuffUL >> 1, samFILE);

      if(! tmpStr) break; /*End of file*/

      samEntry_Fun_12_Sec_02_checkEOL:;

      ulStr = (ulong *) tmpStr;

      while(!
         ((   (*ulStr & ( ~ def_newLine_samEntry) )
            - def_one_samEntry
          ) & def_highBit_samEntry
        )
      ) ++ulStr;

      /*Logic:
      ` This is a little faster then cpLen_ulCp's if check
      `   This is due to this not having to worry about
      `   unintended results.
      `   - *ulStr & ~(0x0a0a0a...0a):
      `     o Converts new lines, start of text
      `       (STX or 2) and backspace (BS or 8) to 0
      `     o So I only have to worry about values not
      `       in a text file (STX and BS)
      `     o def_newLine_samEntry 0x0a0a0a...0a, were
      `       0x0a = '\n'. So, it is  long with every
      `       byte being a new line
      `   - (*ulStr & ~(0x0a0a0a...0a)) - 0x010101...01
      `     o Converts 0 (newline or null) to -1
      `     o def_one_samEntry = 0x010101...01
      `   - (-1 or postive number) & 0x808080...80
      `     o Converts -1 to -127 (smallest char)
      `     o Converts any positive number to 0
      `     o def_highBit_samEntry = 0x808080...80
      */

      tmpStr = (char *) ulStr;

      while( *tmpStr & ~'\n' ) ++tmpStr;
   } /*Loop: Find the length of  the line*/
   
   /*Add null in for line end (not end of file)*/
   *tmpStr = '\0';
   iterStr = *buffStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^   - Get the query id from the buffer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*This is a comment*/
   if(iterStr[0] == '@') goto extraEntry;
 
   /*Query id, single byte copy may be better here,
    ` hard to tell
   */
   samSTPtr->lenQryIdUC = (uchar)
      cpDelim_ulCp(
         samSTPtr->qryIdStr,
         iterStr,
         def_tab_samEntry,
         '\t'
      ); /*Copy the reference id/name*/

   iterStr += samSTPtr->lenQryIdUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec04:
   ^   - Get the flag for the alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr =
      strToUS_base10str(
         iterStr,
         samSTPtr->flagUS
      ); /*Get the flag for the alignment*/

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec05:
   ^   - REad in the reference name/id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Not clear wich version is faster here*/
   samSTPtr->lenRefIdUC = (uchar)
      cpDelim_ulCp(
         samSTPtr->refIdStr,
         iterStr,
         def_tab_samEntry,
         '\t'
      ); /*Copy the reference id/name*/
   
   iterStr += samSTPtr->lenRefIdUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec06:
   ^   - Get the reference position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr =
      strToUI_base10str(
         iterStr,
         samSTPtr->refStartUI
      ); /*Get the starting base in the reference*/

    /*Convert the starting positionto index 0*/
   samSTPtr->refStartUI -= (samSTPtr->refStartUI > 0);

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec07:
   ^   - Get the mapping quality
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr =
      strToUC_base10str(
         iterStr,
         samSTPtr->mapqUC
      ); /*Get the mapping quality of the alignment*/

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec08:
   ^   - Get the cigar entry
   ^   o fun10 sec08 sub01:
   ^     - Check if there is a cigar entry
   ^   o fun10 sec08 sub02:
   ^     - Read in the cigar entry
   ^   o fun10 sec08 sub03:
   ^     - Count number of matchs/snps/dels/inss/masks in
   ^       the cigar entry
   ^   o fun10 sec08 sub04:
   ^     - Get the read lengths from the cigar entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun10 Sec08 Sub01:
   *   - Check if there is a cigar entry
   \*****************************************************/

   if(iterStr[0] == '*')
   { /*If: the cigar entry was not present*/
      samSTPtr->cigTypeStr[0] = '*';
      samSTPtr->cigValAryI[0] = 0;
      iterStr += 2;
      (samSTPtr)->lenCigUI = 1;

      goto rNextEntry;
   } /*If: the cigar entry was not present*/

   /*****************************************************\
   * Fun10 Sec08 Sub02:
   *   - Read in the cigar entry
   \*****************************************************/
   
   while(iterStr[0] > 32)
   { /*Loop: Read in the cigar entry*/

      /*Using -1 to account for the null I will add at the
      `end
      */
      if((samSTPtr)->lenCigUI >= samSTPtr->lenCigBuffUI-1)
      { /*If: I need to increase the cigar buff size*/
         samSTPtr->lenCigBuffUI <<= 1;

         tmpStr =
            realloc(
               samSTPtr->cigTypeStr,
               samSTPtr->lenCigBuffUI * sizeof(char)
            ); /*Resize the type cigar buffer*/

          if(tmpStr == 0) return 64;
          samSTPtr->cigTypeStr = tmpStr;

         tmpStr = (char *)
            realloc(
               samSTPtr->cigValAryI,
               samSTPtr->lenCigBuffUI * sizeof(int)
            ); /*Resize the value cigar buffer*/

          if(tmpStr == 0) return 64;
          samSTPtr->cigValAryI = (int *) tmpStr;
      } /*If: I need to increase the cigar buff size*/

      /*Record the cigar entry*/
      iterStr =
         strToSI_base10str(
            iterStr,
            samSTPtr->cigValAryI[samSTPtr->lenCigUI]
          ); /*Get the number of bases for this type*/

      samSTPtr->cigTypeStr[samSTPtr->lenCigUI]=iterStr[0];

      /**************************************************\
      * Fun10 Sec08 Sub03:
      *   - Count number of matchs/snps/dels/inss/masks in
      *     the cigar entry
      \**************************************************/

      switch(iterStr[0])
      { /*Switch: Check the cigar entry type*/
         case '=':
         /*Case: This was a match*/
             samSTPtr->numMatchUI +=
                samSTPtr->cigValAryI[samSTPtr->lenCigUI];

             break;
         /*Case: This was a match*/

         case 'M':
         /*Case: This is an snp or match*/
             samSTPtr->numMatchUI +=
                samSTPtr->cigValAryI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an snp or match*/

         case 'X':
         /*Case: This is an snp*/
             samSTPtr->numSnpUI +=
                samSTPtr->cigValAryI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an snp*/

         case 'I':
         /*Case: This is an insertion*/
             samSTPtr->numInsUI +=
                samSTPtr->cigValAryI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an insertion*/

         case 'D':
         /*Case: This is an deletion*/
             samSTPtr->numDelUI +=
                samSTPtr->cigValAryI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an deletion*/

         case 'S':
         /*Case: This is an softmasked region*/
             samSTPtr->numMaskUI +=
                samSTPtr->cigValAryI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an softmasked region*/
      } /*Switch: Check the cigar entry type*/

      ++iterStr; /*Get off the tab*/
      ++(samSTPtr)->lenCigUI;
   } /*Loop: Read in the cigar entry*/

   /*****************************************************\
   * Fun10 Sec08 Sub04:
   *   - Get the read lengths from the cigar entries
   \*****************************************************/

    samSTPtr->cigValAryI[samSTPtr->lenCigUI] = 0;
    samSTPtr->cigTypeStr[samSTPtr->lenCigUI] = '\0';

   samSTPtr->readLenUI =
        samSTPtr->numMatchUI
      + samSTPtr->numSnpUI
      + samSTPtr->numInsUI
      + samSTPtr->numMaskUI;

   samSTPtr->alnReadLenUI =
        samSTPtr->numMatchUI
      + samSTPtr->numSnpUI
      + samSTPtr->numDelUI;

   samSTPtr->refEndUI = samSTPtr->refStartUI;
   samSTPtr->refEndUI += samSTPtr->alnReadLenUI;

   samSTPtr->refEndUI -= (samSTPtr->alnReadLenUI > 0);
      /*-1 from (alnReadLen > 0) converts to index 0*/

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec09:
   ^   - Get the RNEXT entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   rNextEntry:

   /*Not sure which is better here*/
   samSTPtr->lenRNextUC = (uchar)
      cpDelim_ulCp(
         samSTPtr->rNextStr,
         iterStr,
         def_tab_samEntry,
         '\t'
      ); /*Copy the query id/name*/

   iterStr += samSTPtr->lenRNextUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec10:
   ^   - Get the PNEXT entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr =
      strToSI_base10str(
         iterStr,
         samSTPtr->pNextI
       ); /*Get the number of bases for this type*/

   samSTPtr->pNextI -= (samSTPtr->pNextI > 0);

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec11:
   ^   - Get the TLEN entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr =
      strToSI_base10str(
         iterStr,
         samSTPtr->tLenI
       ); /*Get the number of bases for this type*/

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec12:
   ^   - Get the sequence entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(samSTPtr->readLenUI == 0 && iterStr[0] != '*')
      samSTPtr->readLenUI =
         lenStr_ulCp(iterStr, def_tab_samEntry, '\t');

   else if(iterStr[0] == '*')
   { /*Else If: There  is no sequence entry*/
      samSTPtr->seqStr[0] = '*';
      samSTPtr->seqStr[1] = '\0';
      iterStr += 2;

      goto noQEntry;
   } /*Else If: There  is no sequence entry*/

   if(samSTPtr->readLenUI > samSTPtr->lenSeqBuffUI)
   { /*If: I need to resize sequence & q-score buffers*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         malloc((samSTPtr->readLenUI + 1) * sizeof(char));

      if(!samSTPtr->seqStr) return 64;

      free(samSTPtr->qStr);
      samSTPtr->qStr = 0;

      samSTPtr->qStr =
         malloc((samSTPtr->readLenUI + 1) * sizeof(char));

      if(!samSTPtr->qStr) return 64;

      samSTPtr->lenSeqBuffUI = samSTPtr->readLenUI;
      samSTPtr->lenQBuffUI = samSTPtr->readLenUI;
   } /*If: I need to resize sequence & q-score buffers*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      iterStr,
      samSTPtr->readLenUI
   );

   iterStr += samSTPtr->readLenUI + 1; /*+1 gets off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec13:
   ^   - Get the q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(iterStr[0] == '*' && iterStr[1] == '\t')
   { /*If: there is no q-score entry*/
      noQEntry:

      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';
      iterStr += 2;
   } /*If: there is no q-score entry*/

   else
   { /*Else: is a q-score entry*/
      iterStr +=
         cpQEntry_samEntry(
            samSTPtr,
            iterStr,
            0
         ) + 1;
   } /*Else: is a q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec14:
   ^   - Copy the extra entry; after strict sam entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   extraEntry:

   /*not sure if char or ul copy better here*/
   samSTPtr->lenExtraUI = lenStr_ulCp(iterStr, 0, '\0');

   if(samSTPtr->lenExtraUI > samSTPtr->lenExtraBuffUI)
   { /*If: I need to resize the buffer*/
      free(samSTPtr->extraStr);
      samSTPtr->extraStr = 0;

      samSTPtr->extraStr =
         malloc((samSTPtr->lenExtraUI +1)* sizeof(char));
      
      if(samSTPtr->extraStr == 0) return 64;
   } /*If: I need to resize the buffer*/

   cpLen_ulCp(
      samSTPtr->extraStr,
      iterStr,
      samSTPtr->lenExtraUI
   ); /*Copy the extra entry*/

   return 0;
} /*get_samLine*/

/*-------------------------------------------------------\
| Fun12: p_samEntry
|  - Prints the sam file entry to a file. This does not
|    print any extra stats that were found.
| Input:
|  - samST
|    o Pointer to samEntry struct with sam entry to print
|  - buffStr:
|    o Pointer to c-string buffer to temporarly hold the
|      cigar entry (speeds things up)
|  - lenBuffUL:
|    o Current length of buffer, adjusted if buffStr is
|      expanded
|  - pNoNewLineBl:
|    o 1: do not print a new line after (you will do this)
|         this is for when you want to add in extra
|         entires
|    o 0: end of sam entry, print a new line
|  - outFILE:
|    o File to print the sam entry to
| Output:
|  - Prints:
|    o Sam file entry in samST to outFILE.
|  - Modifies:
|    o inceases the size of buffStr, if is not 8x the
|      cigar length
|    o Sets lenBuffUL to the new buffStr size when buffStr
|      is resized
|  - Returns:
|    o 0 for no problems
|    o 64 for memory errors
\-------------------------------------------------------*/
char
p_samEntry(
   struct samEntry *samSTPtr,
   char **buffStr,
   ulong *lenBuffUL,
   char pNoNewLineBl,
   void *outFILE
){
   uint uiCig = 0;
   char *tmpStr = *buffStr;
   ulong maxLenUL = 0; 

   maxLenUL =
        (samSTPtr)->lenQryIdUC
      + (samSTPtr)->lenRefIdUC
      + ((samSTPtr)->lenCigUI << 3)
      + (samSTPtr)->lenRNextUC
      + ((samSTPtr)->readLenUI << 1)
      + (samSTPtr)->lenExtraBuffUI
      + 128; /*128 to account for numbers/tabs*/

   if(*(lenBuffUL) < maxLenUL)
   { /*If: I need to add more buffer*/
      free(*buffStr);
      *(lenBuffUL) = maxLenUL;
      *buffStr = malloc((*lenBuffUL + 1) * sizeof(char));
   } /*If: I need to add more buffer*/

   tmpStr = *buffStr;
   
   if(! (*buffStr)) return 64; /*If I had a memory error*/
   
   else if(
         (samSTPtr)->extraStr[0] == '@'
      && (samSTPtr)->qryIdStr[0] == '\0'
   ){
      fprintf(
         (FILE *) (outFILE),
         "%s\n",
         (samSTPtr)->extraStr
      ); /*Print out the header (in extra since commnet)*/

      return 0;
   } /*Else If: this was a header*/
   
   /*Copy the query id to the buffer*/
   cpLen_ulCp(
      tmpStr,
      (samSTPtr)->qryIdStr,
      (samSTPtr)->lenQryIdUC
   );

   tmpStr += (samSTPtr)->lenQryIdUC;
   *tmpStr++ = '\t';

   /*Copy the flag over*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->flagUS);
   *tmpStr++ = '\t';

   /*Copy the referenced id to the buffer*/
   cpLen_ulCp(
      tmpStr,
      (samSTPtr)->refIdStr,
      (samSTPtr)->lenRefIdUC
   );

   tmpStr += (samSTPtr)->lenRefIdUC;
   *tmpStr++ = '\t';

   /*Reference position*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->refStartUI + 1);
   *tmpStr++ = '\t';

   /*mapq*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->mapqUC);
   *tmpStr++ = '\t';

   /*Check if there is a cigar entry*/
   if((samSTPtr)->cigTypeStr[0]=='*') *tmpStr++ ='*';

   else
   { /*Else: convert the cigar entry*/
      for(
         uiCig=0;
         uiCig < (samSTPtr)->lenCigUI;
         ++uiCig
      ){ /*Loop: Convert cigar to string*/
         tmpStr +=
            numToStr(
               tmpStr,
               (samSTPtr)->cigValAryI[uiCig]
            );
         *tmpStr++ = (samSTPtr)->cigTypeStr[uiCig];
      } /*Loop: Convert cigar to string*/
   } /*Else: convert the cigar entry*/

   *tmpStr++ = '\t';

   /*RNEXT*/
   cpLen_ulCp(
      tmpStr,
      (samSTPtr)->rNextStr,
      (samSTPtr)->lenRNextUC
   );

   tmpStr += (samSTPtr)->lenRNextUC;
   *tmpStr++ = '\t';

   /*PNEXT*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->mapqUC);
   *tmpStr++ = '\t';

   /*TLEN*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->mapqUC);
   *tmpStr++ = '\t';

   /*Copy the sequence to the buffer*/
   cpLen_ulCp(
      tmpStr,
      (samSTPtr)->seqStr,
      (samSTPtr)->readLenUI
   );

   tmpStr += (samSTPtr)->readLenUI;
   *tmpStr++ = '\t';

   /*Copy the q-score entry to the buffer*/
   if((samSTPtr)->qStr[1] == '\0') *tmpStr++ = '*';

   else
   { /*Else: there is an q-score entry*/
      cpLen_ulCp(
         tmpStr,
         (samSTPtr)->qStr,
         (samSTPtr)->readLenUI
      );

      tmpStr += (samSTPtr)->readLenUI;
   } /*Else: there is an q-score entry*/

   if((samSTPtr)->lenExtraUI)
   { /*If: have extra items*/
      *tmpStr++ = '\t';

      /*Copy the extra entry*/
      cpLen_ulCp(
         tmpStr,
         (samSTPtr)->extraStr,
         (samSTPtr)->lenExtraUI
      );

      tmpStr += (samSTPtr)->lenExtraUI;
   } /*If: have extra items*/

   if(! pNoNewLineBl)
      *tmpStr++ = '\n';

   fwrite(
      *buffStr,
      sizeof(char),
      (tmpStr - *buffStr),
      (FILE *) (outFILE)
   );
   
   return 0;
} /*p_samEntry*/

/*-------------------------------------------------------\
| Fun13: pfq_samEntry
|  - Prints the sam entry as a fastq entry to a fastq file
| Input:
|  - samST:
|    o Pointer to samEntry structure with fastq entry to
|      print out
|  - outFILE:
|    o Fastq file to print the new fastq entry to
| Output:
|  - Prints:
|    o fastq entry from samST to outFILE
\-------------------------------------------------------*/
void pfq_samEntry(
   struct samEntry *samSTPtr,
   void *outFILE
){
   if(   (samSTPtr)->seqStr[1] != '\0'
      && (samSTPtr)->qStr[1] != '\0'
   ){ /*If: This entry can be printed out*/
      fprintf(
        (FILE *) (outFILE),
        "@%s ref=%s start=%u len=%u refAlnLen=%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  flag=%u mapq=%u match=%u snp=%u ins=%u",
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->numMatchUI,
        (samSTPtr)->numSnpUI,
        (samSTPtr)->numInsUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  del=%u softMasked=%u meanQ=%f medianQ=%f\n",
        (samSTPtr)->numDelUI,
        (samSTPtr)->numMaskUI,
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF
      );
      fprintf(
        (FILE *) (outFILE),
        "%s\n+\n%s\n",
        (samSTPtr)->seqStr,
        (samSTPtr)->qStr
      );
   } /*If: This entry can be printed out*/
} /*p_samEntryAsFq*/

/*-------------------------------------------------------\
| Fun14: pfa_samEntry
|  - Prints the sam entry as a fasta entry to a fasta file
| Input:
|  - samST:
|    o Pointer to samEntry structure with fastq entry to
|      print out
|  - outFILE:
|    o Fasta file to print the new fasta entry to
| Output:
|  - Prints:
|    o fasta entry from samST to outFILE
\-------------------------------------------------------*/
void pfa_samEntry(
   struct samEntry *samSTPtr,
   void *outFILE
){
   if((samSTPtr)->seqStr[1] != '\0')
   { /*If: This entry can be printed out*/
      fprintf(
        (FILE *) (outFILE),
        ">%s ref=%s start=%u len=%u refAlnLen=%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  flag=%u mapq=%u match=%u snp=%u ins=%u",
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->numMatchUI,
        (samSTPtr)->numSnpUI,
        (samSTPtr)->numInsUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  del=%u softMasked=%u meanQ=%f medianQ=%f\n",
        (samSTPtr)->numDelUI,
        (samSTPtr)->numMaskUI,
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF
      );
      fprintf((outFILE), "%s\n", (samSTPtr)->seqStr);
   } /*If: This entry can be printed out*/
} /*pfa_samEntry*/

/*-------------------------------------------------------\
| Fun15: pstats_samEntry
|  - Prints out the stats in a samEntry struct to a file
| Input:
|  - samEntryST:
|    o Pointer to samEntry struct to print stats for
|  - pHeadBl:
|    o 1: Print the header for the stats tsv file
|    o 0: Do not print the header
|  - pNsBl:
|    o 1: find and print out the anonymous base counts
|    o 0: do not print out anonymous base counts
|  - outFILE:
|    o TSV (tab deliminated) file to print stats to
| Output:
|  - Prints:
|    o line with stats from samEntryStruct to file (tsv)
|  - Modifies:
|    o printHeaderChar to be 0 if set to 1
\-------------------------------------------------------*/
void pstats_samEntry(
   struct samEntry *samSTPtr,
   char *pHeadBl,
   char pNsBl,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun15 TOC:
   '   - Prints stats in a samEntry struct to file
   '   o fun15 sec01:
   '     - variable declerations
   '   o fun15 sec02:
   '     - check if comment, if not check if print header
   '   o fun15 sec03:
   '     - print out general stats
   '   o fun15 sec04:
   '     - print matches, snps, ins, dels, and masking.
   '     - if asked also include anonymous bases
   '   o fun15 sec05:
   '     - print out the accuracy stats
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec01:
   ^   - variable declerations 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint matchNCntUI = 0;
   uint snpNCntUI = 0;
   uint insNCntUI = 0;
   uint maskNCntUI = 0;
   uint *cntUI = 0;

   sint siCig = 0;
   sint numNtSI = 0;
   schar *tmpStr = 0;
   uchar ntUC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec02:
   ^   - check if comment, if not check if print header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(   (samSTPtr)->extraStr[0] != '@'
      && (samSTPtr)->qryIdStr[0] != '\0'
   ){ /*If: This is not a comment*/

      if(*(pHeadBl))
      { /*If: I need to print the header*/
        fprintf((FILE *) (outFILE), "Read\tRef\tFlag");
        fprintf((FILE *) (outFILE), "\tMapQ\tRefPos"); 
        fprintf((FILE *) (outFILE), "\tReadLength");
        fprintf((FILE *) (outFILE), "\tRefAlnLength");

        if((pNsBl))
        { /*If: I am printing the anonymous base counts*/
           fprintf((FILE *) (outFILE), "\tmatch_total");
           fprintf((FILE *) (outFILE), "\tmatch");
           fprintf((FILE *)(outFILE),"\tmatch_anonymous");

           fprintf((FILE *) (outFILE), "\tsnp_total");
           fprintf((FILE *) (outFILE), "\tsnp");
           fprintf((FILE *) (outFILE), "\tsnp_anonymous");

           fprintf((FILE *) (outFILE), "\tins_total");
           fprintf((FILE *) (outFILE), "\tins");
           fprintf((FILE *) (outFILE), "\tins_anonymous");

           fprintf((FILE *) (outFILE), "\tdel");

           fprintf((FILE *) (outFILE), "\tmask_total");
           fprintf((FILE *) (outFILE), "\tmask");
           fprintf((FILE *) (outFILE),"\tmask_anonymous");
        } /*If: I am printing the anonymous base counts*/

        else
        { /*Else: Print out the normal header*/
           fprintf((FILE *) (outFILE), "\tMatches\tSnps");
           fprintf((FILE *) (outFILE), "\tInss\tDels");
           fprintf((FILE *) (outFILE), "\tmask");
        } /*Else: Print out the normal header*/

        fprintf((FILE *) (outFILE), "\tMedianQ\tMeanQ\n");

        *(pHeadBl) = 0;
      } /*If: I need to print the header*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun15 Sec03:
      ^   - print out general stats
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      fprintf(
        (FILE *) (outFILE),
        "%s\t%s\t%u\t%u\t%u\t%u\t%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun15 Sec04:
      ^   - print matches, snps, ins, dels, and masking.
      ^   - if asked also include anonymous bases
      ^   o fun15 sec04 sub01:
      ^     - find the number of anonymous bases
      ^   o fun15 sec04 sub02:
      ^     - print out anonymous base and other counts
      ^   o fun15 sec04 sub03:
      ^     - not printing anoynous bases do regular print
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      /**************************************************\
      * Fun15 Sec04 Sub01:
      *   - find the number of anonymous bases
      \**************************************************/

      if( (pNsBl) )
      { /*If: I am finding the anonymous base counts*/
         tmpStr = (schar *) (samSTPtr)->seqStr;

         siCig = 0;

         while(siCig < (sint) (samSTPtr)->lenCigUI)
         { /*Loop: count number of anonymous bases*/
            numNtSI = (samSTPtr)->cigValAryI[siCig];

            switch((samSTPtr)->cigTypeStr[siCig])
            { /*Switch: check what type of cigar entry*/
               case 'I':
               /*Case: Insertions*/
                  cntUI = &insNCntUI;
                  break;
               /*Case: Insertions*/

               case 'X':
               /*Case: Mismatches*/
                  cntUI = &snpNCntUI;
                  break;
               /*Case: Mismatches*/

               case '=':
               case 'M':
               /*Case: Matches or non-eqx entry*/
                  cntUI = &matchNCntUI;
                  break;
               /*Case: Matches or non-eqx entry*/

               case 'S':
               /*Case: softmasking*/
                  cntUI = &maskNCntUI;
                  break;
               /*Case: softmasking*/

               /*Handle missing information cases*/
               default:
               /*Case: deletions or hard masks*/
                  cntUI = 0;
                  numNtSI = 0;
               /*Case: deletions or hard masks*/
            } /*Switch: check what type of cigar entry*/

            while(numNtSI > 0)
            { /*Loop: Check each nucleotide*/
               ntUC = ntTo5Bit[(uchar) *tmpStr];
               (*cntUI) +=
                  (!!(ntUC & def_n_fithBit_ntTo5Bit));

               ++tmpStr;
               --numNtSI;
            } /*Loop: Check each nucleotide*/

            ++siCig;
         } /*Loop: count number of anonymous bases*/

         /***********************************************\
         * Fun15 Sec04 Sub02:
         *   - print anonymous base and other counts
         \***********************************************/

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->numMatchUI,
           (samSTPtr)->numMatchUI - matchNCntUI,
           matchNCntUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->numSnpUI,
           (samSTPtr)->numSnpUI - snpNCntUI,
           snpNCntUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u\t%u",
           (samSTPtr)->numInsUI,
           (samSTPtr)->numInsUI - insNCntUI,
           insNCntUI,
           (samSTPtr)->numDelUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->numMaskUI,
           (samSTPtr)->numMaskUI - maskNCntUI,
           maskNCntUI
         );
      } /*If: I am finding the anonymous base counts*/

      /**************************************************\
      * Fun15 Sec04 Sub03:
      *   - not printing anoynous bases do regular print
      \**************************************************/

      else
      { /*Else: I am not counting anonymous bases*/
         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u\t%u\t%u",
           (samSTPtr)->numMatchUI,
           (samSTPtr)->numSnpUI,
           (samSTPtr)->numInsUI,
           (samSTPtr)->numDelUI,
           (samSTPtr)->numMaskUI
         );
      } /*Else: I am not counting anonymous bases*/
      
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun15 Sec05:
      ^   - print out the q-scores
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      fprintf(
        (FILE *) (outFILE),
        "\t%f\t%f\n",
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF
      );

   } /*If: This is not a comment*/
} /*pstats_samEntry*/

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
