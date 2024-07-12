/*########################################################
# Name alnSetStruct
# Use:
#  o Holds the settings structures and supporting
#    functions for setting structures for alnSeq's
#    pairwise aligners.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - included libraries
'  o .h st01 alnSet:
'     o Holds settings for my alignment program
'  o .h fun01 setScore_alnSet:
'    - Sets the score for a base pair (reference/query)
'  o .h fun02 setMatch_alnSet:
'    - Sets if two bases are a match or mismtach
'  o fun03 freeStack_alnSet:
'    o Frees variables inside alnSet
'  o fun04 freeHeap_alnSet:
'    o Frees an alnSet structure (and sets to 0)
'  o .h fun05 getScore_alnSet:
'    - Get the score for a pair of bases from an alignment
'  o .h fun06 getMatch_alnSet:
'    - Check if two bases were a match or mismatch
'  o fun07 readScoreFile_alnSet
'     - Reads in a file of scores for a scoring matrix
'  o fun08 readMatchFile_alnSet:
'    - Reads in a file of matches
'  o .h fun09 seqToIndex_alnSet:
'    - Converts a sequence to a look up table index
'  o .h fun10 indexToSeq_alnSet:
'    - Converts a sequence of lookup indexs back into
'      uppercase characters (a-z)
'  o fun11 init_alnSet:
'    - Set all values in altSet (alingment settings)
'      structure to defaults
'  o license:
'    - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "alnSet.h"

#include <stdio.h>

#include "alnDefs.h"
#include "../generalLib/base10str.h"
#include "../generalLib/dataTypeShortHand.h"

/*-------------------------------------------------------\
| Fun03: freeStack_alnSet
|  - Does a stack free of an alnSet structer
| Input:
|  - alnSetSTPtr:
|    o alnSetSTPtr to free internal variables in
| Output:
|  - Free:
|    o Nothing, there are no heap allocated variables.
|      Is here in case I decide to have heap allocated
|      variables on day.
\-------------------------------------------------------*/
void
freeStack_alnSet(
   struct alnSet *alnSetSTPtr
){alnSetSTPtr = alnSetSTPtr; /*quites error message*/}

/*-------------------------------------------------------\
| Fun04: freeHeap_alnSet
|  - Frees and alnSet (alignment settings) structure
| Input:
|  - alnSetSTPtr:
|    o Pionter to alnSetST to free
| Output:
|  - Free:
|    o alnSetSTPtr
|  - Set:
|    o alnSetSTPtr to 0 (NULL)
\-------------------------------------------------------*/
void
freeHeap_alnSet(
   struct alnSet *alnSetSTPtr
){
   freeStack_alnSet(alnSetSTPtr);
   free(alnSetSTPtr);
} /*freeHeap_alnSet*/

/*-------------------------------------------------------\
| Fun07: readScoreFile_alnSet
|  - Reads in a file of scores for a scoring matrix
| Input:
|  - alnSetSTPtr:
|    o pointer to an alnSetST (settings) structure with
|      the score matrix to modify
|  - scoreFILE:
|    o File to get scores from
| Output:
|  - Modifies:
|    o Score matrix in alngSetST to hold the scores from
|      the file (scoreFILE)
|    o scoreFILE to point to the end of the file
|  - Returns:
|    o 0 for no errors
|    o position of error in file
\-------------------------------------------------------*/
unsigned long
readScoreFile_alnSet(
    struct alnSet *alnSetSTPtr, /*score matrix to change*/
    void *scoreFILE  /*File scoring matrix scores*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC: readScoreFile_alnSet
   '  o fun07 sec01:
   '    - Variable declerations & set up
   '  o fun07 sec02:
   '    - Blank the scoring matrix
   '  o fun07 sec03:
   '    - Read in line and check if comment
   '  o fun07 sec04:
   '    - Convert score & add to matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^  - Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ushort lenBuffUS = 1024;
   char buffCStr[lenBuffUS];
   char *tmpCStr = 0;
   sshort scoreSS = 0;

   uchar colUC = 0;
   uchar rowUC = 0;

   buffCStr[lenBuffUS - 1] = '\0';
   buffCStr[lenBuffUS - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^  - Blank the scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < defMatrixCol; ++colUC)
   { /*Loop: blank all values in the scoring matrix*/
       for(rowUC = 0; rowUC < defMatrixCol; ++rowUC)
           alnSetSTPtr->scoreMatrixSS[colUC][rowUC] = 0;
   } /*Loop: blank all values in the scoring matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^  - Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(fgets(buffCStr, 1024, scoreFILE))
   { /*While I have scores to read in*/
       
       if(buffCStr[0] == '/' && buffCStr[1] == '/')
       { /*On a comment, move onto the next line*/
           while(
               buffCStr[lenBuffUS - 2] != '\0' &&
               buffCStr[lenBuffUS - 2] != '\n'
           ) { /*While have more buffer to read in*/
               buffCStr[lenBuffUS - 2] = '\0';
               fgets(buffCStr, 1024, (FILE *) scoreFILE);
           } /*While have more buffer to read in*/

           /*Reset the buffer*/
           buffCStr[lenBuffUS - 2] = '\0';

           continue;
       } /*On a comment, move onto the next line*/

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun07 Sec04:
       ^  - Convert score & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffCStr[0] == '\n')
           continue;                        /*Blank line*/

       if(buffCStr[0] < 64 && buffCStr[2] < 64)
           return ftell(scoreFILE);  /*Invalid character*/
       
       tmpCStr = strToSS_base10str(&buffCStr[4], scoreSS);

       setScore_alnSet(
         buffCStr[0],
         buffCStr[2],
         scoreSS,
         alnSetSTPtr
       ); /*Add the score to the matrix*/

       if(tmpCStr == &buffCStr[3])
           return ftell(scoreFILE);         /*No score*/

       while(
           buffCStr[lenBuffUS - 2] != '\0' &&
           buffCStr[lenBuffUS - 2] != '\n'
       ){ /*While have more buffer to read in*/
           buffCStr[lenBuffUS - 2] = '\0';
           fgets(buffCStr, 1024, (FILE *) scoreFILE);
       } /*While have more buffer to read in*/

       /*Reset the buffer*/
       buffCStr[lenBuffUS - 2] = '\0';
   } /*While I have scores to read in*/

   return 0;
} /*readScoreFile_alnSet*/

/*-------------------------------------------------------\
| Fun08: readMatchFile_alnSet
|  - Reads in a file of matches
| Input:
|  - alnSetSTPtr:
|    o pointer to an alnSetST (settings) structure with
|      the match matrix to modify
|  - matchFILE:
|    o File to get matchs from
| Output:
|  - Modifies:
|    o Match matrix in alngSetST to hold the matchs from
|      the file (matchFILE)
|    o matchFILE to point to the end of the file
|  - Returns:
|    o 0 for no errors
|    o position of error in file
\-------------------------------------------------------*/
unsigned long
readMatchFile_alnSet(
    struct alnSet *alnSetSTPtr,
    void *matchFILE  /*File of matchs for scoring matrix*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun08 TOC: readMatchFile_alnSet
   '  o fun08 sec01:
   '    - Variable declerations & set up
   '  o fun08 sec02:
   '    - Blank the match matrix
   '  o fun08 sec03:
   '    - Read in line and check if comment
   '  o fun08 sec04:
   '    - Add match/snp (mismatch) to match matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec01:
   ^  - Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ushort lenBuffUS = 1024;
   char buffCStr[lenBuffUS];

   uchar colUC = 0;
   uchar rowUC = 0;

   buffCStr[lenBuffUS - 1] = '\0';
   buffCStr[lenBuffUS - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec02:
   ^  - Blank the scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < defMatrixCol; ++colUC)
   { /*Loop: blank all values in the scoring matrix*/
       for(rowUC = 0; rowUC < defMatrixCol; ++rowUC)
           alnSetSTPtr->matchMatrixSC[colUC][rowUC] = 0;
   } /*Loop: blank all values in the scoring matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec03:
   ^  - Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(fgets(buffCStr, 1024, (FILE *) matchFILE))
   { /*While I have matchs to read in*/
       
       if(buffCStr[0] == '/' && buffCStr[1] == '/')
       { /*On a comment, move onto the next line*/
           while(
               buffCStr[lenBuffUS - 2] != '\0' &&
               buffCStr[lenBuffUS - 2] != '\n'
           ) { /*While have more buffer to read in*/
               buffCStr[lenBuffUS - 2] = '\0';
               fgets(buffCStr, 1024, (FILE *) matchFILE);
           } /*While have more buffer to read in*/

           /*Reset the buffer*/
           buffCStr[lenBuffUS - 2] = '\0';

           continue;
       } /*On a comment, move onto the next line*/

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun08 Sec04:
       ^  - Convert match & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffCStr[0] == '\n')
           continue;                        /*Blank line*/

       if(buffCStr[4] != '1' && buffCStr[4] != '0')
           return ftell((FILE *) matchFILE);
           /*This error means I do not know if match/snp*/

       setMatch_alnSet(
         buffCStr[0],
         buffCStr[2],
         buffCStr[4],
         alnSetSTPtr
       ); /*Add the match to the matrix*/

       while(
           buffCStr[lenBuffUS - 2] != '\0' &&
           buffCStr[lenBuffUS - 2] != '\n'
       ){ /*While have more buffer to read in*/
           buffCStr[lenBuffUS - 2] = '\0';
           fgets(buffCStr, 1024, (FILE *) matchFILE);
       } /*While have more buffer to read in*/

       /*Reset the buffer*/
       buffCStr[lenBuffUS - 2] = '\0';
   } /*While I have matchs to read in*/

   return 0;
} /*REAdMatchFile_alnSet*/

/*-------------------------------------------------------\
| Fun11: init_alnSet
|  - Set values in altSet (alingment settings) structure
|    to default values
| Input:
|  - alnSetSTPtr:
|    o poineter to an alnSet (settings) structure to
|      initialize
| Output:
|  o Modifies:
|    - alnSetST to have default alignment settings values
\-------------------------------------------------------*/
void
init_alnSet(
    struct alnSet *alnSetST /*Has settings to initialize*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC: init_alnSet
   '  - Set values in altSet (alingment settings)
   '    structure to defaults
   '  o fun11 sec01:
   '    - Set non-matrix variables
   '  o fun11 sec02:
   '    - Initialize scoring matrix
   '  o fun11 sec03:
   '    - Initialize match matrix
   '  o fun11 sec04:
   '    - set up scoring matrix for nucleotides
   '  o fun11 sec05:
   '    - set up matching matrix for nucleotides
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^  - Set non-matrix variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Variables for my for loop*/
   uchar colUC = 0;
   uchar rowUC = 0;

   alnSetST->gapSS = def_gapOpen_alnDefs;
   alnSetST->extendSS = def_gapExtend_alnDefs;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^  - Initialize scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < defMatrixCol; ++colUC)
   { /*loop for all columns in the comparison matrix*/
       for(rowUC = 0; rowUC < defMatrixCol; ++rowUC)
           alnSetST->scoreMatrixSS[colUC][rowUC] = 0;
           /*Most of these cells will never be used*/
           /*But are needed to build the table*/
   } /*loop for all columns in the comparison matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec03:
   ^  - Initialize match matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Both words, DNA, and AA all are the same when both
    ` characters/bases/amino acids are the same.
    */
    for(colUC = 0; colUC < defMatrixCol; ++colUC)
    { /*loop for all columns in the comparison matrix*/
        for(rowUC = 0; rowUC < defMatrixCol; ++rowUC)
        { /*Loop: Fill in the entire matrix*/ 
           if(colUC == rowUC)
             alnSetST->matchMatrixSC[colUC][rowUC] =
                defBaseMatch;
           else
             alnSetST->matchMatrixSC[colUC][rowUC] =
                defBaseSnp;
        } /*Loop: Fill in the entire matrix*/ 
    } /*loop for all columns in the comparison matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec04:
   ^  - set up scoring matrix for nucleotides
   ^  o fun11 sec04 sub01:
   ^    - a as first base
   ^  o fun11 sec04 sub02:
   ^    - t as first base
   ^  o fun11 sec04 sub03:
   ^    - u (t) as first base
   ^  o fun11 sec04 sub04:
   ^    - g as first base
   ^  o fun11 sec04 sub05:
   ^    - c as first base
   ^  o fun11 sec04 sub06:
   ^    - w (anonymous) as first base
   ^  o fun11 sec04 sub07:
   ^    - s (anonymous) as first base
   ^  o fun11 sec04 sub08:
   ^    - m (anonymous) as first base
   ^  o fun11 sec04 sub09:
   ^    - k (anonymous) as first base
   ^  o fun11 sec04 sub10:
   ^    - r (anonymous) as first base
   ^  o fun11 sec04 sub11:
   ^    - y (anonymous) as first base
   ^  o fun11 sec04 sub12:
   ^    - b (anonymous) as first base
   ^  o fun11 sec04 sub13:
   ^    - d (anonymous) as first base
   ^  o fun11 sec04 sub14:
   ^    - h (anonymous) as first base
   ^  o fun11 sec04 sub15:
   ^    - v (anonymous) as first base
   ^  o fun11 sec04 sub16:
   ^    - n (anonymous) as first base
   ^  o fun11 sec04 sub17:
   ^    - x (anonymous) as first base (technically aa)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun11 Sec04 Sub01:
   *   - a as first base
   \*****************************************************/

   setScore_alnSet('a','a',def_AToA_alnDefs,alnSetST);
   setScore_alnSet('a','t',def_AToT_alnDefs,alnSetST);
   setScore_alnSet('a','u',def_AToU_alnDefs,alnSetST);
   setScore_alnSet('a','g',def_AToG_alnDefs,alnSetST);
   setScore_alnSet('a','c',def_AToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('a','w',def_AToW_alnDefs,alnSetST);
   setScore_alnSet('a','s',def_AToS_alnDefs,alnSetST);
   setScore_alnSet('a','m',def_AToM_alnDefs,alnSetST);
   setScore_alnSet('a','k',def_AToK_alnDefs,alnSetST);
   setScore_alnSet('a','r',def_AToR_alnDefs,alnSetST);
   setScore_alnSet('a','y',def_AToY_alnDefs,alnSetST);
   setScore_alnSet('a','b',def_AToB_alnDefs,alnSetST);
   setScore_alnSet('a','d',def_AToD_alnDefs,alnSetST);
   setScore_alnSet('a','h',def_AToH_alnDefs,alnSetST);
   setScore_alnSet('a','v',def_AToV_alnDefs,alnSetST);
   setScore_alnSet('a','n',def_AToN_alnDefs,alnSetST);
   setScore_alnSet('a','x',def_AToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub02:
   *   - t as first base
   \*****************************************************/

   setScore_alnSet('t','a',def_TToA_alnDefs,alnSetST);
   setScore_alnSet('t','t',def_TToT_alnDefs,alnSetST);
   setScore_alnSet('t','u',def_TToU_alnDefs,alnSetST);
   setScore_alnSet('t','g',def_TToG_alnDefs,alnSetST);
   setScore_alnSet('t','c',def_TToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('t','w',def_TToW_alnDefs,alnSetST);
   setScore_alnSet('t','s',def_TToS_alnDefs,alnSetST);
   setScore_alnSet('t','m',def_TToM_alnDefs,alnSetST);
   setScore_alnSet('t','k',def_TToK_alnDefs,alnSetST);
   setScore_alnSet('t','r',def_TToR_alnDefs,alnSetST);
   setScore_alnSet('t','y',def_TToY_alnDefs,alnSetST);
   setScore_alnSet('t','b',def_TToB_alnDefs,alnSetST);
   setScore_alnSet('t','d',def_TToD_alnDefs,alnSetST);
   setScore_alnSet('t','h',def_TToH_alnDefs,alnSetST);
   setScore_alnSet('t','v',def_TToV_alnDefs,alnSetST);
   setScore_alnSet('t','n',def_TToN_alnDefs,alnSetST);
   setScore_alnSet('t','x',def_TToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub03:
   *   - u (t) as first base
   \*****************************************************/

   setScore_alnSet('u','a',def_UToA_alnDefs,alnSetST);
   setScore_alnSet('u','t',def_UToT_alnDefs,alnSetST);
   setScore_alnSet('u','u',def_UToU_alnDefs,alnSetST);
   setScore_alnSet('u','g',def_UToG_alnDefs,alnSetST);
   setScore_alnSet('u','c',def_UToC_alnDefs,alnSetST);

   /*Set u & t to same scores (U is RNA version of T)*/
   setScore_alnSet('u','w',def_UToW_alnDefs,alnSetST);
   setScore_alnSet('u','s',def_UToS_alnDefs,alnSetST);
   setScore_alnSet('u','m',def_UToM_alnDefs,alnSetST);
   setScore_alnSet('u','k',def_UToK_alnDefs,alnSetST);
   setScore_alnSet('u','r',def_UToR_alnDefs,alnSetST);
   setScore_alnSet('u','y',def_UToY_alnDefs,alnSetST);
   setScore_alnSet('u','b',def_UToB_alnDefs,alnSetST);
   setScore_alnSet('u','d',def_UToD_alnDefs,alnSetST);
   setScore_alnSet('u','h',def_UToH_alnDefs,alnSetST);
   setScore_alnSet('u','v',def_UToV_alnDefs,alnSetST);
   setScore_alnSet('u','n',def_UToN_alnDefs,alnSetST);
   setScore_alnSet('u','x',def_UToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub04:
   *   - g as first base
   \*****************************************************/

   setScore_alnSet('g','a',def_GToA_alnDefs,alnSetST);
   setScore_alnSet('g','t',def_GToT_alnDefs,alnSetST);
   setScore_alnSet('g','u',def_GToU_alnDefs,alnSetST);
   setScore_alnSet('g','g',def_GToG_alnDefs,alnSetST);
   setScore_alnSet('g','c',def_GToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('g','w',def_GToW_alnDefs,alnSetST);
   setScore_alnSet('g','s',def_GToS_alnDefs,alnSetST);
   setScore_alnSet('g','m',def_GToM_alnDefs,alnSetST);
   setScore_alnSet('g','k',def_GToK_alnDefs,alnSetST);
   setScore_alnSet('g','r',def_GToR_alnDefs,alnSetST);
   setScore_alnSet('g','y',def_GToY_alnDefs,alnSetST);
   setScore_alnSet('g','b',def_GToB_alnDefs,alnSetST);
   setScore_alnSet('g','d',def_GToD_alnDefs,alnSetST);
   setScore_alnSet('g','h',def_GToH_alnDefs,alnSetST);
   setScore_alnSet('g','v',def_GToV_alnDefs,alnSetST);
   setScore_alnSet('g','n',def_GToN_alnDefs,alnSetST);
   setScore_alnSet('g','x',def_GToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub05:
   *   - c as first base
   \*****************************************************/

   setScore_alnSet('c','a',def_CToA_alnDefs,alnSetST);
   setScore_alnSet('c','t',def_CToT_alnDefs,alnSetST);
   setScore_alnSet('c','u',def_CToT_alnDefs,alnSetST);
   setScore_alnSet('c','g',def_CToG_alnDefs,alnSetST);
   setScore_alnSet('c','c',def_CToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('c','w',def_CToW_alnDefs,alnSetST);
   setScore_alnSet('c','s',def_CToS_alnDefs,alnSetST);
   setScore_alnSet('c','m',def_CToM_alnDefs,alnSetST);
   setScore_alnSet('c','k',def_CToK_alnDefs,alnSetST);
   setScore_alnSet('c','r',def_CToR_alnDefs,alnSetST);
   setScore_alnSet('c','y',def_CToY_alnDefs,alnSetST);
   setScore_alnSet('c','b',def_CToB_alnDefs,alnSetST);
   setScore_alnSet('c','d',def_CToD_alnDefs,alnSetST);
   setScore_alnSet('c','h',def_CToH_alnDefs,alnSetST);
   setScore_alnSet('c','v',def_CToV_alnDefs,alnSetST);
   setScore_alnSet('c','n',def_CToN_alnDefs,alnSetST);
   setScore_alnSet('c','x',def_CToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub06:
   *   - w (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('w','a',def_WToA_alnDefs,alnSetST);
   setScore_alnSet('w','c',def_WToC_alnDefs,alnSetST);
   setScore_alnSet('w','g',def_WToG_alnDefs,alnSetST);
   setScore_alnSet('w','t',def_WToT_alnDefs,alnSetST);
   setScore_alnSet('w','u',def_WToT_alnDefs,alnSetST);

   setScore_alnSet('w','w',def_WToW_alnDefs,alnSetST);
   setScore_alnSet('w','s',def_WToS_alnDefs,alnSetST);
   setScore_alnSet('w','m',def_WToM_alnDefs,alnSetST);
   setScore_alnSet('w','k',def_WToK_alnDefs,alnSetST);
   setScore_alnSet('w','r',def_WToR_alnDefs,alnSetST);
   setScore_alnSet('w','y',def_WToY_alnDefs,alnSetST);
   setScore_alnSet('w','b',def_WToB_alnDefs,alnSetST);
   setScore_alnSet('w','d',def_WToD_alnDefs,alnSetST);
   setScore_alnSet('w','h',def_WToH_alnDefs,alnSetST);
   setScore_alnSet('w','v',def_WToV_alnDefs,alnSetST);
   setScore_alnSet('w','n',def_WToN_alnDefs,alnSetST);
   setScore_alnSet('w','x',def_WToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub07:
   *   - s (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('s','a',def_SToA_alnDefs,alnSetST);
   setScore_alnSet('s','c',def_SToC_alnDefs,alnSetST);
   setScore_alnSet('s','g',def_SToG_alnDefs,alnSetST);
   setScore_alnSet('s','t',def_SToT_alnDefs,alnSetST);
   setScore_alnSet('s','u',def_SToT_alnDefs,alnSetST);

   setScore_alnSet('s','w',def_SToW_alnDefs,alnSetST);
   setScore_alnSet('s','s',def_SToS_alnDefs,alnSetST);
   setScore_alnSet('s','m',def_SToM_alnDefs,alnSetST);
   setScore_alnSet('s','k',def_SToK_alnDefs,alnSetST);
   setScore_alnSet('s','r',def_SToR_alnDefs,alnSetST);
   setScore_alnSet('s','y',def_SToY_alnDefs,alnSetST);
   setScore_alnSet('s','b',def_SToB_alnDefs,alnSetST);
   setScore_alnSet('s','d',def_SToD_alnDefs,alnSetST);
   setScore_alnSet('s','h',def_SToH_alnDefs,alnSetST);
   setScore_alnSet('s','v',def_SToV_alnDefs,alnSetST);
   setScore_alnSet('s','n',def_SToN_alnDefs,alnSetST);
   setScore_alnSet('s','x',def_SToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub08:
   *   - m (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('m','a',def_MToA_alnDefs,alnSetST);
   setScore_alnSet('m','c',def_MToC_alnDefs,alnSetST);
   setScore_alnSet('m','g',def_MToG_alnDefs,alnSetST);
   setScore_alnSet('m','t',def_MToT_alnDefs,alnSetST);
   setScore_alnSet('m','u',def_MToT_alnDefs,alnSetST);

   setScore_alnSet('m','w',def_MToW_alnDefs,alnSetST);
   setScore_alnSet('m','s',def_MToS_alnDefs,alnSetST);
   setScore_alnSet('m','m',def_MToM_alnDefs,alnSetST);
   setScore_alnSet('m','k',def_MToK_alnDefs,alnSetST);
   setScore_alnSet('m','r',def_MToR_alnDefs,alnSetST);
   setScore_alnSet('m','y',def_MToY_alnDefs,alnSetST);
   setScore_alnSet('m','b',def_MToB_alnDefs,alnSetST);
   setScore_alnSet('m','d',def_MToD_alnDefs,alnSetST);
   setScore_alnSet('m','h',def_MToH_alnDefs,alnSetST);
   setScore_alnSet('m','v',def_MToV_alnDefs,alnSetST);
   setScore_alnSet('m','n',def_MToN_alnDefs,alnSetST);
   setScore_alnSet('m','x',def_MToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub09:
   *   - k (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('k','a',def_KToA_alnDefs,alnSetST);
   setScore_alnSet('k','c',def_KToC_alnDefs,alnSetST);
   setScore_alnSet('k','g',def_KToG_alnDefs,alnSetST);
   setScore_alnSet('k','t',def_KToT_alnDefs,alnSetST);
   setScore_alnSet('k','u',def_KToT_alnDefs,alnSetST);

   setScore_alnSet('k','w',def_KToW_alnDefs,alnSetST);
   setScore_alnSet('k','s',def_KToS_alnDefs,alnSetST);
   setScore_alnSet('k','m',def_KToM_alnDefs,alnSetST);
   setScore_alnSet('k','k',def_KToK_alnDefs,alnSetST);
   setScore_alnSet('k','r',def_KToR_alnDefs,alnSetST);
   setScore_alnSet('k','y',def_KToY_alnDefs,alnSetST);
   setScore_alnSet('k','b',def_KToB_alnDefs,alnSetST);
   setScore_alnSet('k','d',def_KToD_alnDefs,alnSetST);
   setScore_alnSet('k','h',def_KToH_alnDefs,alnSetST);
   setScore_alnSet('k','v',def_KToV_alnDefs,alnSetST);
   setScore_alnSet('k','n',def_KToN_alnDefs,alnSetST);
   setScore_alnSet('k','x',def_KToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub10:
   *   - r (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('r','a',def_RToA_alnDefs,alnSetST);
   setScore_alnSet('r','c',def_RToC_alnDefs,alnSetST);
   setScore_alnSet('r','g',def_RToG_alnDefs,alnSetST);
   setScore_alnSet('r','t',def_RToT_alnDefs,alnSetST);
   setScore_alnSet('r','u',def_RToT_alnDefs,alnSetST);

   setScore_alnSet('r','w',def_RToW_alnDefs,alnSetST);
   setScore_alnSet('r','s',def_RToS_alnDefs,alnSetST);
   setScore_alnSet('r','m',def_RToM_alnDefs,alnSetST);
   setScore_alnSet('r','k',def_RToK_alnDefs,alnSetST);
   setScore_alnSet('r','r',def_RToR_alnDefs,alnSetST);
   setScore_alnSet('r','y',def_RToY_alnDefs,alnSetST);
   setScore_alnSet('r','b',def_RToB_alnDefs,alnSetST);
   setScore_alnSet('r','d',def_RToD_alnDefs,alnSetST);
   setScore_alnSet('r','h',def_RToH_alnDefs,alnSetST);
   setScore_alnSet('r','v',def_RToV_alnDefs,alnSetST);
   setScore_alnSet('r','n',def_RToN_alnDefs,alnSetST);
   setScore_alnSet('r','x',def_RToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub11:
   *   - y (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('y','a',def_YToA_alnDefs,alnSetST);
   setScore_alnSet('y','c',def_YToC_alnDefs,alnSetST);
   setScore_alnSet('y','g',def_YToG_alnDefs,alnSetST);
   setScore_alnSet('y','t',def_YToT_alnDefs,alnSetST);
   setScore_alnSet('y','u',def_YToT_alnDefs,alnSetST);

   setScore_alnSet('y','w',def_YToW_alnDefs,alnSetST);
   setScore_alnSet('y','s',def_YToS_alnDefs,alnSetST);
   setScore_alnSet('y','m',def_YToM_alnDefs,alnSetST);
   setScore_alnSet('y','k',def_YToK_alnDefs,alnSetST);
   setScore_alnSet('y','r',def_YToR_alnDefs,alnSetST);
   setScore_alnSet('y','y',def_YToY_alnDefs,alnSetST);
   setScore_alnSet('y','b',def_YToB_alnDefs,alnSetST);
   setScore_alnSet('y','d',def_YToD_alnDefs,alnSetST);
   setScore_alnSet('y','h',def_YToH_alnDefs,alnSetST);
   setScore_alnSet('y','v',def_YToV_alnDefs,alnSetST);
   setScore_alnSet('y','n',def_YToN_alnDefs,alnSetST);
   setScore_alnSet('y','x',def_YToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub12:
   *   - b (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('b','a',def_BToA_alnDefs,alnSetST);
   setScore_alnSet('b','c',def_BToC_alnDefs,alnSetST);
   setScore_alnSet('b','g',def_BToG_alnDefs,alnSetST);
   setScore_alnSet('b','t',def_BToT_alnDefs,alnSetST);
   setScore_alnSet('b','u',def_BToT_alnDefs,alnSetST);

   setScore_alnSet('b','w',def_BToW_alnDefs,alnSetST);
   setScore_alnSet('b','s',def_BToS_alnDefs,alnSetST);
   setScore_alnSet('b','m',def_BToM_alnDefs,alnSetST);
   setScore_alnSet('b','k',def_BToK_alnDefs,alnSetST);
   setScore_alnSet('b','r',def_BToR_alnDefs,alnSetST);
   setScore_alnSet('b','y',def_BToY_alnDefs,alnSetST);
   setScore_alnSet('b','b',def_BToB_alnDefs,alnSetST);
   setScore_alnSet('b','d',def_BToD_alnDefs,alnSetST);
   setScore_alnSet('b','h',def_BToH_alnDefs,alnSetST);
   setScore_alnSet('b','v',def_BToV_alnDefs,alnSetST);
   setScore_alnSet('b','n',def_BToN_alnDefs,alnSetST);
   setScore_alnSet('b','x',def_BToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub13:
   *   - d (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('d','a',def_DToA_alnDefs,alnSetST);
   setScore_alnSet('d','c',def_DToC_alnDefs,alnSetST);
   setScore_alnSet('d','g',def_DToG_alnDefs,alnSetST);
   setScore_alnSet('d','t',def_DToT_alnDefs,alnSetST);
   setScore_alnSet('d','u',def_DToT_alnDefs,alnSetST);

   setScore_alnSet('d','w',def_DToW_alnDefs,alnSetST);
   setScore_alnSet('d','s',def_DToS_alnDefs,alnSetST);
   setScore_alnSet('d','m',def_DToM_alnDefs,alnSetST);
   setScore_alnSet('d','k',def_DToK_alnDefs,alnSetST);
   setScore_alnSet('d','r',def_DToR_alnDefs,alnSetST);
   setScore_alnSet('d','y',def_DToY_alnDefs,alnSetST);
   setScore_alnSet('d','b',def_DToB_alnDefs,alnSetST);
   setScore_alnSet('d','d',def_DToD_alnDefs,alnSetST);
   setScore_alnSet('d','h',def_DToH_alnDefs,alnSetST);
   setScore_alnSet('d','v',def_DToV_alnDefs,alnSetST);
   setScore_alnSet('d','n',def_DToN_alnDefs,alnSetST);
   setScore_alnSet('d','x',def_DToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub14:
   *   - h (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('h','a',def_HToA_alnDefs,alnSetST);
   setScore_alnSet('h','c',def_HToC_alnDefs,alnSetST);
   setScore_alnSet('h','g',def_HToG_alnDefs,alnSetST);
   setScore_alnSet('h','t',def_HToT_alnDefs,alnSetST);
   setScore_alnSet('h','u',def_HToT_alnDefs,alnSetST);

   setScore_alnSet('h','w',def_HToW_alnDefs,alnSetST);
   setScore_alnSet('h','s',def_HToS_alnDefs,alnSetST);
   setScore_alnSet('h','m',def_HToM_alnDefs,alnSetST);
   setScore_alnSet('h','k',def_HToK_alnDefs,alnSetST);
   setScore_alnSet('h','r',def_HToR_alnDefs,alnSetST);
   setScore_alnSet('h','y',def_HToY_alnDefs,alnSetST);
   setScore_alnSet('h','b',def_HToB_alnDefs,alnSetST);
   setScore_alnSet('h','d',def_HToD_alnDefs,alnSetST);
   setScore_alnSet('h','h',def_HToH_alnDefs,alnSetST);
   setScore_alnSet('h','v',def_HToV_alnDefs,alnSetST);
   setScore_alnSet('h','n',def_HToN_alnDefs,alnSetST);
   setScore_alnSet('h','x',def_HToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub15:
   *   - v (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('v','a',def_VToA_alnDefs,alnSetST);
   setScore_alnSet('v','c',def_VToC_alnDefs,alnSetST);
   setScore_alnSet('v','g',def_VToG_alnDefs,alnSetST);
   setScore_alnSet('v','t',def_VToT_alnDefs,alnSetST);
   setScore_alnSet('v','u',def_VToT_alnDefs,alnSetST);

   setScore_alnSet('v','w',def_VToW_alnDefs,alnSetST);
   setScore_alnSet('v','s',def_VToS_alnDefs,alnSetST);
   setScore_alnSet('v','m',def_VToM_alnDefs,alnSetST);
   setScore_alnSet('v','k',def_VToK_alnDefs,alnSetST);
   setScore_alnSet('v','r',def_VToR_alnDefs,alnSetST);
   setScore_alnSet('v','y',def_VToY_alnDefs,alnSetST);
   setScore_alnSet('v','b',def_VToB_alnDefs,alnSetST);
   setScore_alnSet('v','d',def_VToD_alnDefs,alnSetST);
   setScore_alnSet('v','h',def_VToH_alnDefs,alnSetST);
   setScore_alnSet('v','v',def_VToV_alnDefs,alnSetST);
   setScore_alnSet('v','n',def_VToN_alnDefs,alnSetST);
   setScore_alnSet('v','x',def_VToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub16:
   *   - n (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('n','a',def_NToA_alnDefs,alnSetST);
   setScore_alnSet('n','c',def_NToC_alnDefs,alnSetST);
   setScore_alnSet('n','g',def_NToG_alnDefs,alnSetST);
   setScore_alnSet('n','t',def_NToT_alnDefs,alnSetST);
   setScore_alnSet('n','u',def_NToT_alnDefs,alnSetST);

   setScore_alnSet('n','w',def_NToW_alnDefs,alnSetST);
   setScore_alnSet('n','s',def_NToS_alnDefs,alnSetST);
   setScore_alnSet('n','m',def_NToM_alnDefs,alnSetST);
   setScore_alnSet('n','k',def_NToK_alnDefs,alnSetST);
   setScore_alnSet('n','r',def_NToR_alnDefs,alnSetST);
   setScore_alnSet('n','y',def_NToY_alnDefs,alnSetST);
   setScore_alnSet('n','b',def_NToB_alnDefs,alnSetST);
   setScore_alnSet('n','d',def_NToD_alnDefs,alnSetST);
   setScore_alnSet('n','h',def_NToH_alnDefs,alnSetST);
   setScore_alnSet('n','v',def_NToV_alnDefs,alnSetST);
   setScore_alnSet('n','n',def_NToN_alnDefs,alnSetST);
   setScore_alnSet('n','x',def_NToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec04 Sub17:
   *   - x (anonymous) as first base (technically aa)
   \*****************************************************/

   setScore_alnSet('x','a',def_XToA_alnDefs,alnSetST);
   setScore_alnSet('x','c',def_XToC_alnDefs,alnSetST);
   setScore_alnSet('x','g',def_XToG_alnDefs,alnSetST);
   setScore_alnSet('x','t',def_XToT_alnDefs,alnSetST);
   setScore_alnSet('x','u',def_XToT_alnDefs,alnSetST);

   setScore_alnSet('x','w',def_XToW_alnDefs,alnSetST);
   setScore_alnSet('x','s',def_XToS_alnDefs,alnSetST);
   setScore_alnSet('x','m',def_XToM_alnDefs,alnSetST);
   setScore_alnSet('x','k',def_XToK_alnDefs,alnSetST);
   setScore_alnSet('x','r',def_XToR_alnDefs,alnSetST);
   setScore_alnSet('x','y',def_XToY_alnDefs,alnSetST);
   setScore_alnSet('x','b',def_XToB_alnDefs,alnSetST);
   setScore_alnSet('x','d',def_XToD_alnDefs,alnSetST);
   setScore_alnSet('x','h',def_XToH_alnDefs,alnSetST);
   setScore_alnSet('x','v',def_XToV_alnDefs,alnSetST);
   setScore_alnSet('x','n',def_XToN_alnDefs,alnSetST);
   setScore_alnSet('x','x',def_XToX_alnDefs,alnSetST);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec05:
   ^  - set up matching matrix for nucleotides
   ^  o fun11 sec05 sub01:
   ^    - a as first base
   ^  o fun11 sec05 sub02:
   ^    - t as first base
   ^  o fun11 sec05 sub03:
   ^    - u (t) as first base
   ^  o fun11 sec05 sub04:
   ^    - g as first base
   ^  o fun11 sec05 sub05:
   ^    - c as first base
   ^  o fun11 sec05 sub06:
   ^    - w (anonymous) as first base
   ^  o fun11 sec05 sub07:
   ^    - s (anonymous) as first base
   ^  o fun11 sec05 sub08:
   ^    - m (anonymous) as first base
   ^  o fun11 sec05 sub09:
   ^    - k (anonymous) as first base
   ^  o fun11 sec05 sub10:
   ^    - r (anonymous) as first base
   ^  o fun11 sec05 sub11:
   ^    - y (anonymous) as first base
   ^  o fun11 sec05 sub12:
   ^    - b (anonymous) as first base
   ^  o fun11 sec05 sub13:
   ^    - d (anonymous) as first base
   ^  o fun11 sec05 sub14:
   ^    - h (anonymous) as first base
   ^  o fun11 sec05 sub15:
   ^    - v (anonymous) as first base
   ^  o fun11 sec05 sub16:
   ^    - n (anonymous) as first base
   ^  o fun11 sec05 sub17:
   ^    - x (anonymous) as first base (technically aa)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun11 Sec05 Sub01:
   *   - a as first base
   \*****************************************************/

   setMatch_alnSet('a','a',def_AEqlA_alnDefs,alnSetST);
   setMatch_alnSet('a','t',def_AEqlT_alnDefs,alnSetST);
   setMatch_alnSet('a','u',def_AEqlU_alnDefs,alnSetST);
   setMatch_alnSet('a','g',def_AEqlG_alnDefs,alnSetST);
   setMatch_alnSet('a','c',def_AEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('a','w',def_AEqlW_alnDefs,alnSetST);
   setMatch_alnSet('a','s',def_AEqlS_alnDefs,alnSetST);
   setMatch_alnSet('a','m',def_AEqlM_alnDefs,alnSetST);
   setMatch_alnSet('a','k',def_AEqlK_alnDefs,alnSetST);
   setMatch_alnSet('a','r',def_AEqlR_alnDefs,alnSetST);
   setMatch_alnSet('a','y',def_AEqlY_alnDefs,alnSetST);
   setMatch_alnSet('a','b',def_AEqlB_alnDefs,alnSetST);
   setMatch_alnSet('a','d',def_AEqlD_alnDefs,alnSetST);
   setMatch_alnSet('a','h',def_AEqlH_alnDefs,alnSetST);
   setMatch_alnSet('a','v',def_AEqlV_alnDefs,alnSetST);
   setMatch_alnSet('a','n',def_AEqlN_alnDefs,alnSetST);
   setMatch_alnSet('a','x',def_AEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub02:
   *   - t as first base
   \*****************************************************/

   setMatch_alnSet('t','a',def_TEqlA_alnDefs,alnSetST);
   setMatch_alnSet('t','t',def_TEqlT_alnDefs,alnSetST);
   setMatch_alnSet('t','u',def_TEqlU_alnDefs,alnSetST);
   setMatch_alnSet('t','g',def_TEqlG_alnDefs,alnSetST);
   setMatch_alnSet('t','c',def_TEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('t','w',def_TEqlW_alnDefs,alnSetST);
   setMatch_alnSet('t','s',def_TEqlS_alnDefs,alnSetST);
   setMatch_alnSet('t','m',def_TEqlM_alnDefs,alnSetST);
   setMatch_alnSet('t','k',def_TEqlK_alnDefs,alnSetST);
   setMatch_alnSet('t','r',def_TEqlR_alnDefs,alnSetST);
   setMatch_alnSet('t','y',def_TEqlY_alnDefs,alnSetST);
   setMatch_alnSet('t','b',def_TEqlB_alnDefs,alnSetST);
   setMatch_alnSet('t','d',def_TEqlD_alnDefs,alnSetST);
   setMatch_alnSet('t','h',def_TEqlH_alnDefs,alnSetST);
   setMatch_alnSet('t','v',def_TEqlV_alnDefs,alnSetST);
   setMatch_alnSet('t','n',def_TEqlN_alnDefs,alnSetST);
   setMatch_alnSet('t','x',def_TEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub03:
   *   - u (t) as first base
   \*****************************************************/

   setMatch_alnSet('u','a',def_UEqlA_alnDefs,alnSetST);
   setMatch_alnSet('u','g',def_UEqlG_alnDefs,alnSetST);
   setMatch_alnSet('u','c',def_UEqlC_alnDefs,alnSetST);
   setMatch_alnSet('u','t',def_UEqlT_alnDefs,alnSetST);
   setMatch_alnSet('u','u',def_UEqlU_alnDefs,alnSetST);

   /*Set u & t to same scores (U is RNA version of T)*/
   setMatch_alnSet('u','w',def_UEqlW_alnDefs,alnSetST);
   setMatch_alnSet('u','s',def_UEqlS_alnDefs,alnSetST);
   setMatch_alnSet('u','m',def_UEqlM_alnDefs,alnSetST);
   setMatch_alnSet('u','k',def_UEqlK_alnDefs,alnSetST);
   setMatch_alnSet('u','r',def_UEqlR_alnDefs,alnSetST);
   setMatch_alnSet('u','y',def_UEqlY_alnDefs,alnSetST);
   setMatch_alnSet('u','b',def_UEqlB_alnDefs,alnSetST);
   setMatch_alnSet('u','d',def_UEqlD_alnDefs,alnSetST);
   setMatch_alnSet('u','h',def_UEqlH_alnDefs,alnSetST);
   setMatch_alnSet('u','v',def_UEqlV_alnDefs,alnSetST);
   setMatch_alnSet('u','n',def_UEqlN_alnDefs,alnSetST);
   setMatch_alnSet('u','x',def_UEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub04:
   *   - g as first base
   \*****************************************************/

   setMatch_alnSet('g','a',def_GEqlA_alnDefs,alnSetST);
   setMatch_alnSet('g','t',def_GEqlT_alnDefs,alnSetST);
   setMatch_alnSet('g','u',def_GEqlU_alnDefs,alnSetST);
   setMatch_alnSet('g','g',def_GEqlG_alnDefs,alnSetST);
   setMatch_alnSet('g','c',def_GEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('g','w',def_GEqlW_alnDefs,alnSetST);
   setMatch_alnSet('g','s',def_GEqlS_alnDefs,alnSetST);
   setMatch_alnSet('g','m',def_GEqlM_alnDefs,alnSetST);
   setMatch_alnSet('g','k',def_GEqlK_alnDefs,alnSetST);
   setMatch_alnSet('g','r',def_GEqlR_alnDefs,alnSetST);
   setMatch_alnSet('g','y',def_GEqlY_alnDefs,alnSetST);
   setMatch_alnSet('g','b',def_GEqlB_alnDefs,alnSetST);
   setMatch_alnSet('g','d',def_GEqlD_alnDefs,alnSetST);
   setMatch_alnSet('g','h',def_GEqlH_alnDefs,alnSetST);
   setMatch_alnSet('g','v',def_GEqlV_alnDefs,alnSetST);
   setMatch_alnSet('g','n',def_GEqlN_alnDefs,alnSetST);
   setMatch_alnSet('g','x',def_GEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub05:
   *   - c as first base
   \*****************************************************/

   setMatch_alnSet('c','a',def_CEqlA_alnDefs,alnSetST);
   setMatch_alnSet('c','t',def_CEqlT_alnDefs,alnSetST);
   setMatch_alnSet('c','u',def_CEqlT_alnDefs,alnSetST);
   setMatch_alnSet('c','g',def_CEqlG_alnDefs,alnSetST);
   setMatch_alnSet('c','c',def_CEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('c','w',def_CEqlW_alnDefs,alnSetST);
   setMatch_alnSet('c','s',def_CEqlS_alnDefs,alnSetST);
   setMatch_alnSet('c','m',def_CEqlM_alnDefs,alnSetST);
   setMatch_alnSet('c','k',def_CEqlK_alnDefs,alnSetST);
   setMatch_alnSet('c','r',def_CEqlR_alnDefs,alnSetST);
   setMatch_alnSet('c','y',def_CEqlY_alnDefs,alnSetST);
   setMatch_alnSet('c','b',def_CEqlB_alnDefs,alnSetST);
   setMatch_alnSet('c','d',def_CEqlD_alnDefs,alnSetST);
   setMatch_alnSet('c','h',def_CEqlH_alnDefs,alnSetST);
   setMatch_alnSet('c','v',def_CEqlV_alnDefs,alnSetST);
   setMatch_alnSet('c','n',def_CEqlN_alnDefs,alnSetST);
   setMatch_alnSet('c','x',def_CEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub06:
   *   - w (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('w','a',def_WEqlA_alnDefs,alnSetST);
   setMatch_alnSet('w','c',def_WEqlC_alnDefs,alnSetST);
   setMatch_alnSet('w','g',def_WEqlG_alnDefs,alnSetST);
   setMatch_alnSet('w','t',def_WEqlT_alnDefs,alnSetST);
   setMatch_alnSet('w','u',def_WEqlT_alnDefs,alnSetST);

   setMatch_alnSet('w','w',def_WEqlW_alnDefs,alnSetST);
   setMatch_alnSet('w','s',def_WEqlS_alnDefs,alnSetST);
   setMatch_alnSet('w','m',def_WEqlM_alnDefs,alnSetST);
   setMatch_alnSet('w','k',def_WEqlK_alnDefs,alnSetST);
   setMatch_alnSet('w','r',def_WEqlR_alnDefs,alnSetST);
   setMatch_alnSet('w','y',def_WEqlY_alnDefs,alnSetST);
   setMatch_alnSet('w','b',def_WEqlB_alnDefs,alnSetST);
   setMatch_alnSet('w','d',def_WEqlD_alnDefs,alnSetST);
   setMatch_alnSet('w','h',def_WEqlH_alnDefs,alnSetST);
   setMatch_alnSet('w','v',def_WEqlV_alnDefs,alnSetST);
   setMatch_alnSet('w','n',def_WEqlN_alnDefs,alnSetST);
   setMatch_alnSet('w','x',def_WEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub07:
   *   - s (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('s','a',def_SEqlA_alnDefs,alnSetST);
   setMatch_alnSet('s','c',def_SEqlC_alnDefs,alnSetST);
   setMatch_alnSet('s','g',def_SEqlG_alnDefs,alnSetST);
   setMatch_alnSet('s','t',def_SEqlT_alnDefs,alnSetST);
   setMatch_alnSet('s','u',def_SEqlT_alnDefs,alnSetST);

   setMatch_alnSet('s','w',def_SEqlW_alnDefs,alnSetST);
   setMatch_alnSet('s','s',def_SEqlS_alnDefs,alnSetST);
   setMatch_alnSet('s','m',def_SEqlM_alnDefs,alnSetST);
   setMatch_alnSet('s','k',def_SEqlK_alnDefs,alnSetST);
   setMatch_alnSet('s','r',def_SEqlR_alnDefs,alnSetST);
   setMatch_alnSet('s','y',def_SEqlY_alnDefs,alnSetST);
   setMatch_alnSet('s','b',def_SEqlB_alnDefs,alnSetST);
   setMatch_alnSet('s','d',def_SEqlD_alnDefs,alnSetST);
   setMatch_alnSet('s','h',def_SEqlH_alnDefs,alnSetST);
   setMatch_alnSet('s','v',def_SEqlV_alnDefs,alnSetST);
   setMatch_alnSet('s','n',def_SEqlN_alnDefs,alnSetST);
   setMatch_alnSet('s','x',def_SEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub08:
   *   - m (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('m','a',def_MEqlA_alnDefs,alnSetST);
   setMatch_alnSet('m','c',def_MEqlC_alnDefs,alnSetST);
   setMatch_alnSet('m','g',def_MEqlG_alnDefs,alnSetST);
   setMatch_alnSet('m','t',def_MEqlT_alnDefs,alnSetST);
   setMatch_alnSet('m','u',def_MEqlT_alnDefs,alnSetST);

   setMatch_alnSet('m','w',def_MEqlW_alnDefs,alnSetST);
   setMatch_alnSet('m','s',def_MEqlS_alnDefs,alnSetST);
   setMatch_alnSet('m','m',def_MEqlM_alnDefs,alnSetST);
   setMatch_alnSet('m','k',def_MEqlK_alnDefs,alnSetST);
   setMatch_alnSet('m','r',def_MEqlR_alnDefs,alnSetST);
   setMatch_alnSet('m','y',def_MEqlY_alnDefs,alnSetST);
   setMatch_alnSet('m','b',def_MEqlB_alnDefs,alnSetST);
   setMatch_alnSet('m','d',def_MEqlD_alnDefs,alnSetST);
   setMatch_alnSet('m','h',def_MEqlH_alnDefs,alnSetST);
   setMatch_alnSet('m','v',def_MEqlV_alnDefs,alnSetST);
   setMatch_alnSet('m','n',def_MEqlN_alnDefs,alnSetST);
   setMatch_alnSet('m','x',def_MEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub09:
   *   - k (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('k','a',def_KEqlA_alnDefs,alnSetST);
   setMatch_alnSet('k','c',def_KEqlC_alnDefs,alnSetST);
   setMatch_alnSet('k','g',def_KEqlG_alnDefs,alnSetST);
   setMatch_alnSet('k','t',def_KEqlT_alnDefs,alnSetST);
   setMatch_alnSet('k','u',def_KEqlT_alnDefs,alnSetST);

   setMatch_alnSet('k','w',def_KEqlW_alnDefs,alnSetST);
   setMatch_alnSet('k','s',def_KEqlS_alnDefs,alnSetST);
   setMatch_alnSet('k','m',def_KEqlM_alnDefs,alnSetST);
   setMatch_alnSet('k','k',def_KEqlK_alnDefs,alnSetST);
   setMatch_alnSet('k','r',def_KEqlR_alnDefs,alnSetST);
   setMatch_alnSet('k','y',def_KEqlY_alnDefs,alnSetST);
   setMatch_alnSet('k','b',def_KEqlB_alnDefs,alnSetST);
   setMatch_alnSet('k','d',def_KEqlD_alnDefs,alnSetST);
   setMatch_alnSet('k','h',def_KEqlH_alnDefs,alnSetST);
   setMatch_alnSet('k','v',def_KEqlV_alnDefs,alnSetST);
   setMatch_alnSet('k','n',def_KEqlN_alnDefs,alnSetST);
   setMatch_alnSet('k','x',def_KEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub10:
   *   - r (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('r','a',def_REqlA_alnDefs,alnSetST);
   setMatch_alnSet('r','c',def_REqlC_alnDefs,alnSetST);
   setMatch_alnSet('r','g',def_REqlG_alnDefs,alnSetST);
   setMatch_alnSet('r','t',def_REqlT_alnDefs,alnSetST);
   setMatch_alnSet('r','u',def_REqlT_alnDefs,alnSetST);

   setMatch_alnSet('r','w',def_REqlW_alnDefs,alnSetST);
   setMatch_alnSet('r','s',def_REqlS_alnDefs,alnSetST);
   setMatch_alnSet('r','m',def_REqlM_alnDefs,alnSetST);
   setMatch_alnSet('r','k',def_REqlK_alnDefs,alnSetST);
   setMatch_alnSet('r','r',def_REqlR_alnDefs,alnSetST);
   setMatch_alnSet('r','y',def_REqlY_alnDefs,alnSetST);
   setMatch_alnSet('r','b',def_REqlB_alnDefs,alnSetST);
   setMatch_alnSet('r','d',def_REqlD_alnDefs,alnSetST);
   setMatch_alnSet('r','h',def_REqlH_alnDefs,alnSetST);
   setMatch_alnSet('r','v',def_REqlV_alnDefs,alnSetST);
   setMatch_alnSet('r','n',def_REqlN_alnDefs,alnSetST);
   setMatch_alnSet('r','x',def_REqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub11:
   *   - y (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('y','a',def_YEqlA_alnDefs,alnSetST);
   setMatch_alnSet('y','c',def_YEqlC_alnDefs,alnSetST);
   setMatch_alnSet('y','g',def_YEqlG_alnDefs,alnSetST);
   setMatch_alnSet('y','t',def_YEqlT_alnDefs,alnSetST);
   setMatch_alnSet('y','u',def_YEqlT_alnDefs,alnSetST);

   setMatch_alnSet('y','w',def_YEqlW_alnDefs,alnSetST);
   setMatch_alnSet('y','s',def_YEqlS_alnDefs,alnSetST);
   setMatch_alnSet('y','m',def_YEqlM_alnDefs,alnSetST);
   setMatch_alnSet('y','k',def_YEqlK_alnDefs,alnSetST);
   setMatch_alnSet('y','r',def_YEqlR_alnDefs,alnSetST);
   setMatch_alnSet('y','y',def_YEqlY_alnDefs,alnSetST);
   setMatch_alnSet('y','b',def_YEqlB_alnDefs,alnSetST);
   setMatch_alnSet('y','d',def_YEqlD_alnDefs,alnSetST);
   setMatch_alnSet('y','h',def_YEqlH_alnDefs,alnSetST);
   setMatch_alnSet('y','v',def_YEqlV_alnDefs,alnSetST);
   setMatch_alnSet('y','n',def_YEqlN_alnDefs,alnSetST);
   setMatch_alnSet('y','x',def_YEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub12:
   *   - b (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('b','a',def_BEqlA_alnDefs,alnSetST);
   setMatch_alnSet('b','c',def_BEqlC_alnDefs,alnSetST);
   setMatch_alnSet('b','g',def_BEqlG_alnDefs,alnSetST);
   setMatch_alnSet('b','t',def_BEqlT_alnDefs,alnSetST);
   setMatch_alnSet('b','u',def_BEqlT_alnDefs,alnSetST);

   setMatch_alnSet('b','w',def_BEqlW_alnDefs,alnSetST);
   setMatch_alnSet('b','s',def_BEqlS_alnDefs,alnSetST);
   setMatch_alnSet('b','m',def_BEqlM_alnDefs,alnSetST);
   setMatch_alnSet('b','k',def_BEqlK_alnDefs,alnSetST);
   setMatch_alnSet('b','r',def_BEqlR_alnDefs,alnSetST);
   setMatch_alnSet('b','y',def_BEqlY_alnDefs,alnSetST);
   setMatch_alnSet('b','b',def_BEqlB_alnDefs,alnSetST);
   setMatch_alnSet('b','d',def_BEqlD_alnDefs,alnSetST);
   setMatch_alnSet('b','h',def_BEqlH_alnDefs,alnSetST);
   setMatch_alnSet('b','v',def_BEqlV_alnDefs,alnSetST);
   setMatch_alnSet('b','n',def_BEqlN_alnDefs,alnSetST);
   setMatch_alnSet('b','x',def_BEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub13:
   *   - d (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('d','a',def_DEqlA_alnDefs,alnSetST);
   setMatch_alnSet('d','c',def_DEqlC_alnDefs,alnSetST);
   setMatch_alnSet('d','g',def_DEqlG_alnDefs,alnSetST);
   setMatch_alnSet('d','t',def_DEqlT_alnDefs,alnSetST);
   setMatch_alnSet('d','u',def_DEqlT_alnDefs,alnSetST);

   setMatch_alnSet('d','w',def_DEqlW_alnDefs,alnSetST);
   setMatch_alnSet('d','s',def_DEqlS_alnDefs,alnSetST);
   setMatch_alnSet('d','m',def_DEqlM_alnDefs,alnSetST);
   setMatch_alnSet('d','k',def_DEqlK_alnDefs,alnSetST);
   setMatch_alnSet('d','r',def_DEqlR_alnDefs,alnSetST);
   setMatch_alnSet('d','y',def_DEqlY_alnDefs,alnSetST);
   setMatch_alnSet('d','b',def_DEqlB_alnDefs,alnSetST);
   setMatch_alnSet('d','d',def_DEqlD_alnDefs,alnSetST);
   setMatch_alnSet('d','h',def_DEqlH_alnDefs,alnSetST);
   setMatch_alnSet('d','v',def_DEqlV_alnDefs,alnSetST);
   setMatch_alnSet('d','n',def_DEqlN_alnDefs,alnSetST);
   setMatch_alnSet('d','x',def_DEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub14:
   *   - h (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('h','a',def_HEqlA_alnDefs,alnSetST);
   setMatch_alnSet('h','c',def_HEqlC_alnDefs,alnSetST);
   setMatch_alnSet('h','g',def_HEqlG_alnDefs,alnSetST);
   setMatch_alnSet('h','t',def_HEqlT_alnDefs,alnSetST);
   setMatch_alnSet('h','u',def_HEqlT_alnDefs,alnSetST);

   setMatch_alnSet('h','w',def_HEqlW_alnDefs,alnSetST);
   setMatch_alnSet('h','s',def_HEqlS_alnDefs,alnSetST);
   setMatch_alnSet('h','m',def_HEqlM_alnDefs,alnSetST);
   setMatch_alnSet('h','k',def_HEqlK_alnDefs,alnSetST);
   setMatch_alnSet('h','r',def_HEqlR_alnDefs,alnSetST);
   setMatch_alnSet('h','y',def_HEqlY_alnDefs,alnSetST);
   setMatch_alnSet('h','b',def_HEqlB_alnDefs,alnSetST);
   setMatch_alnSet('h','d',def_HEqlD_alnDefs,alnSetST);
   setMatch_alnSet('h','h',def_HEqlH_alnDefs,alnSetST);
   setMatch_alnSet('h','v',def_HEqlV_alnDefs,alnSetST);
   setMatch_alnSet('h','n',def_HEqlN_alnDefs,alnSetST);
   setMatch_alnSet('h','x',def_HEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub15:
   *   - v (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('v','a',def_VEqlA_alnDefs,alnSetST);
   setMatch_alnSet('v','c',def_VEqlC_alnDefs,alnSetST);
   setMatch_alnSet('v','g',def_VEqlG_alnDefs,alnSetST);
   setMatch_alnSet('v','t',def_VEqlT_alnDefs,alnSetST);
   setMatch_alnSet('v','u',def_VEqlT_alnDefs,alnSetST);

   setMatch_alnSet('v','w',def_VEqlW_alnDefs,alnSetST);
   setMatch_alnSet('v','s',def_VEqlS_alnDefs,alnSetST);
   setMatch_alnSet('v','m',def_VEqlM_alnDefs,alnSetST);
   setMatch_alnSet('v','k',def_VEqlK_alnDefs,alnSetST);
   setMatch_alnSet('v','r',def_VEqlR_alnDefs,alnSetST);
   setMatch_alnSet('v','y',def_VEqlY_alnDefs,alnSetST);
   setMatch_alnSet('v','b',def_VEqlB_alnDefs,alnSetST);
   setMatch_alnSet('v','d',def_VEqlD_alnDefs,alnSetST);
   setMatch_alnSet('v','h',def_VEqlH_alnDefs,alnSetST);
   setMatch_alnSet('v','v',def_VEqlV_alnDefs,alnSetST);
   setMatch_alnSet('v','n',def_VEqlN_alnDefs,alnSetST);
   setMatch_alnSet('v','x',def_VEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub16:
   *   - n (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('n','a',def_NEqlA_alnDefs,alnSetST);
   setMatch_alnSet('n','c',def_NEqlC_alnDefs,alnSetST);
   setMatch_alnSet('n','g',def_NEqlG_alnDefs,alnSetST);
   setMatch_alnSet('n','t',def_NEqlT_alnDefs,alnSetST);
   setMatch_alnSet('n','u',def_NEqlT_alnDefs,alnSetST);

   setMatch_alnSet('n','w',def_NEqlW_alnDefs,alnSetST);
   setMatch_alnSet('n','s',def_NEqlS_alnDefs,alnSetST);
   setMatch_alnSet('n','m',def_NEqlM_alnDefs,alnSetST);
   setMatch_alnSet('n','k',def_NEqlK_alnDefs,alnSetST);
   setMatch_alnSet('n','r',def_NEqlR_alnDefs,alnSetST);
   setMatch_alnSet('n','y',def_NEqlY_alnDefs,alnSetST);
   setMatch_alnSet('n','b',def_NEqlB_alnDefs,alnSetST);
   setMatch_alnSet('n','d',def_NEqlD_alnDefs,alnSetST);
   setMatch_alnSet('n','h',def_NEqlH_alnDefs,alnSetST);
   setMatch_alnSet('n','v',def_NEqlV_alnDefs,alnSetST);
   setMatch_alnSet('n','n',def_NEqlN_alnDefs,alnSetST);
   setMatch_alnSet('n','x',def_NEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun11 Sec05 Sub17:
   *   - x (anonymous) as first base (technically aa)
   \*****************************************************/

   setMatch_alnSet('x','a',def_XEqlA_alnDefs,alnSetST);
   setMatch_alnSet('x','c',def_XEqlC_alnDefs,alnSetST);
   setMatch_alnSet('x','g',def_XEqlG_alnDefs,alnSetST);
   setMatch_alnSet('x','t',def_XEqlT_alnDefs,alnSetST);
   setMatch_alnSet('x','u',def_XEqlT_alnDefs,alnSetST);

   setMatch_alnSet('x','w',def_XEqlW_alnDefs,alnSetST);
   setMatch_alnSet('x','s',def_XEqlS_alnDefs,alnSetST);
   setMatch_alnSet('x','m',def_XEqlM_alnDefs,alnSetST);
   setMatch_alnSet('x','k',def_XEqlK_alnDefs,alnSetST);
   setMatch_alnSet('x','r',def_XEqlR_alnDefs,alnSetST);
   setMatch_alnSet('x','y',def_XEqlY_alnDefs,alnSetST);
   setMatch_alnSet('x','b',def_XEqlB_alnDefs,alnSetST);
   setMatch_alnSet('x','d',def_XEqlD_alnDefs,alnSetST);
   setMatch_alnSet('x','h',def_XEqlH_alnDefs,alnSetST);
   setMatch_alnSet('x','v',def_XEqlV_alnDefs,alnSetST);
   setMatch_alnSet('x','n',def_XEqlN_alnDefs,alnSetST);
   setMatch_alnSet('x','x',def_XEqlX_alnDefs,alnSetST);

   return;
} /*init_alnSet*/

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
