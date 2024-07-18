/*########################################################
# Name: kmerCnt
#   - holds functions to count number of matching kmers in
#     sequences
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o .c tbl01: alnNtToBit_kmerCnt
'     - converts an nucleotide alignment code form
'       alnSetStruct.h to an two bit value, with an 3rd
'       bit being used for anonymous bases and the 4th bit
'       for errors
'   o .h fun01: blank_kmerCnt
'     - blanks all variables in a kmerCnt structure
'   o fun02: init_kmerCnt
'     - initialize a kmerCnt structure
'   o fun03: setup_kmerCnt
'     - sets up kmerCnt structure for use
'       (memory allocations)
'   o fun04: freeStack_kmerCnt
'     - frees variables inside a kmerCnt structure
'   o fun05: freeHeap_kmerCnt
'     - frees a kmerCnt structure (you set to 0)
'   o fun06: freeStackAry_kmerCnt
'     - frees internal memory form all kmerCnt structure
'       in an array
'   o fun07: freeHeapAry_kmerCnt
'     - frees an array of kmerCnt structure (you set to 0)
'   o .c fun08: mkKmerMask_kmerCnt
'     - makes a kmer mask for removing extra bases from
'       kmer
'   o .c fun09: addNtToKmer_kmerCnt
'     - adds an nucleotide to an kmer
'   o fun10: addSeq_kmerCnt
'     - adds a sequence to a kmerCnt structure
'   o fun11: ntToKmerAry_kmerCnt
'     - converts a nucleotide sequence to a array of kmer
'       counts
'   o fun12: get_kmerCnt
'     - gets number of matching kmers between a kmerCnt
'       structure and the kmer arrays
'   o fun13: faToKmerCnt_kmerCnt
'     - converts a the sequences in a fasta file to a
'       kmerCnt array
'   o license:
'     - licensing for this code (public domain / mit)
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

#include "kmerCnt.h"

#include "seqST.h"

/*.h files only*/
#include "dataTypeShortHand.h"
#include "genMath.h"
#include "shellSort.h"
#include "ntTo2Bit.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .h  #include "ulCp.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun02: init_kmerCnt
|   - initialize a kmerCnt structure
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure to initialize
| Output:
|   - Modifies:
|     o all values in kmerCntSTPtr to be 0 or null
\-------------------------------------------------------*/
void
init_kmerCnt(
   struct kmerCnt *kmerCntSTPtr
){
   kmerCntSTPtr->forKmerArySI = 0;
   kmerCntSTPtr->forKmersUI = 0;

   kmerCntSTPtr->revKmerArySI = 0;
   kmerCntSTPtr->revKmersUI = 0;

   kmerCntSTPtr->forSeqST = 0;
   kmerCntSTPtr->revSeqST = 0;

   kmerCntSTPtr->lenKmerUC = 0;
   kmerCntSTPtr->maxKmersUI = 0;
} /*init_kmerCnt*/

/*-------------------------------------------------------\
| Fun03: setup_kmerCnt
|   - sets up kmerCnt structure for use
|     (memory allocations)
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure to initialize
|   - lenKmerUC:
|     o length of one kmer
| Output:
|   - Modifies:
|     o lenKmerUC to be one kmer length
|     o maxKmersUI to be the maximum number of kmers
|     o totalKmersUI to be 0 (by blank_kmerCnt)
|   - Allocates:
|     o maxKmersUI uints for forKmerArySI and revKMerArySI
|     o seqST for forSeqST and revSeqST
|   - Returns:
|     o 0 for no problems
|     o def_memErr_kmerCnt for memory errors
\-------------------------------------------------------*/
signed char
setup_kmerCnt(
   struct kmerCnt *kmerCntSTPtr,
   unsigned char lenKmerUC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   o fun03 sec01:
   '     - variable declerations
   '   o fun03 sec02:
   '     - find the kmer array length
   '   o fun03 sec03:
   '     - allocate memory for kmer arrays
   '   o fun03 sec04:
   '     - allocate memory for seqSTs
   '   o fun03 sec05:
   '     - blank and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declerations
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   uchar ucKmer = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - find the kmer array length
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   kmerCntSTPtr->lenKmerUC = lenKmerUC;
   kmerCntSTPtr->maxKmersUI = 1;

   for(
      ucKmer = 0;
      ucKmer < lenKmerUC;
      ++ucKmer
   ) kmerCntSTPtr->maxKmersUI <<= 2; /*multiply by 4*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - allocate memory for kmer arrays
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   kmerCntSTPtr->forKmerArySI =
      malloc(sizeof(sint) * (kmerCntSTPtr->maxKmersUI+1));

   if(! kmerCntSTPtr->forKmerArySI)
      goto memErr_fun03;

   /*reverse kmer array*/

   kmerCntSTPtr->revKmerArySI =
      malloc(sizeof(sint) * (kmerCntSTPtr->maxKmersUI+1));

   if(! kmerCntSTPtr->revKmerArySI)
      goto memErr_fun03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - allocate memory for seqSTs
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   kmerCntSTPtr->forSeqST =
      malloc(sizeof(struct seqST));

   if(! kmerCntSTPtr->forSeqST)
      goto memErr_fun03;

   init_seqST(kmerCntSTPtr->forSeqST);

   /*reverse seqStuct*/

   kmerCntSTPtr->revSeqST =
      malloc(sizeof(struct seqST));

   if(! kmerCntSTPtr->revSeqST)
      goto memErr_fun03;

   init_seqST(kmerCntSTPtr->revSeqST);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec05:
   ^   - blank and return
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   blank_kmerCnt(kmerCntSTPtr);

   return 0;

   memErr_fun03:;

   blank_kmerCnt(kmerCntSTPtr);
   return def_memErr_kmerCnt;
} /*setup_kmerCnt*/

/*-------------------------------------------------------\
| Fun04: freeStack_kmerCnt
|   - frees variables inside a kmerCnt structure
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure with variables to
|       free
| Output:
|   - Frees:
|     o forKmerArySI, revKmerArySI, forSeqST, and revSeqST
|   - Sets:
|     o blanks structure and sets pointers to 0
|       - you will needed to call init_kmerCnt again
|       - only lenKmerUC and maxKMerUC are not set to 0
\-------------------------------------------------------*/
void
freeStack_kmerCnt(
   struct kmerCnt *kmerCntSTPtr
){
   if(kmerCntSTPtr)
   { /*If: input structure is not null*/
      free(kmerCntSTPtr->forKmerArySI);
      kmerCntSTPtr->forKmerArySI = 0;

      free(kmerCntSTPtr->revKmerArySI);
      kmerCntSTPtr->revKmerArySI = 0;

      freeHeap_seqST(kmerCntSTPtr->forSeqST);
      kmerCntSTPtr->forSeqST = 0;

      freeHeap_seqST(kmerCntSTPtr->revSeqST);
      kmerCntSTPtr->revSeqST = 0;

      blank_kmerCnt(kmerCntSTPtr);
   } /*If: input structure is not null*/
} /*freeStack_kmerCnt*/

/*-------------------------------------------------------\
| Fun05: freeHeap_kmerCnt
|   - frees a kmerCnt structure (you set to 0)
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure to free
| Output:
|   - Frees:
|     o kmerCntSTPtr
\-------------------------------------------------------*/
void
freeHeap_kmerCnt(
   struct kmerCnt *kmerCntSTPtr
){
   freeStack_kmerCnt(kmerCntSTPtr);
   free(kmerCntSTPtr);
} /*freeHeap_kmerCnt*/

/*-------------------------------------------------------\
| Fun06: freeStackAry_kmerCnt
|   - frees internal memory form all kmerCnt structure
|     in an array
| Input:
|   - kmerCntAryST:
|     o pointer to an array of kmerCnt structures to
|       all variables for (goes to uninitialzed state)
|   - numSTUI:
|     o number of initialized kmerCnt structures in
|       kmerCntAryST (actually have internal memory)
| Output:
|   - Frees:
|     o forKmerArySI, revKmerArySI, forSeqST, and revSeqST
|       for numSTUI kmerCnt structures in the array
|   - Sets:
|     o blanks structure and sets pointers to 0
|       - you will needed to call init_kmerCnt again
\-------------------------------------------------------*/
void
freeStackAry_kmerCnt(
   struct kmerCnt *kmerCntAryST,
   unsigned int numSTUI
){
   unsigned int uiST = 0;

   for(
      uiST = 0;
      uiST < numSTUI;
      ++uiST
   ) freeStack_kmerCnt(&kmerCntAryST[uiST]);
} /*freeStackAry_kmerCnt*/

/*-------------------------------------------------------\
| Fun07: freeHeapAry_kmerCnt
|   - frees an array of kmerCnt structure (you set to 0)
| Input:
|   - kmerCntAryST:
|     o pointer to an array of kmerCnt structures to free
|   - numSTUI:
|     o number of initialized kmerCnt structures in
|       kmerCntAryST (actually have internal memory)
| Output:
|   - Frees:
|     o kmerCntAryST
\-------------------------------------------------------*/
void
freeHeapAry_kmerCnt(
   struct kmerCnt *kmerCntAryST,
   unsigned int numSTUI
){
   freeStackAry_kmerCnt(
      kmerCntAryST,
      numSTUI
   );

   free(kmerCntAryST);
} /*freeHeapAry_kmerCnt*/

/*-------------------------------------------------------\
| Fun08: mkKmerMask_kmerCnt
|   - makes a kmer mask for removing extra bases from kmer
| Input:
|   - lenKmerUC:
|     o length of one kmer
| Output:
|   - Returns:
|     o unsigned long with mask
\-------------------------------------------------------*/
#define \
mkKmerMask_kmerCnt( \
   lenKmerUC \
)({ \
   unsigned long retMaskMacUL = \
      ((lenKmerUC) * def_bitsPerKmer_kmerCnt); \
   \
   retMaskMacUL = \
      (sizeof(ulong) << 3) - retMaskMacUL; \
   \
   retMaskMacUL = \
      ((ulong) -1) >> retMaskMacUL; \
   \
   retMaskMacUL; \
}) /*mkKmerMask_kmerCnt*/

/*-------------------------------------------------------\
| Fun09: addNtToKmer_kmerCnt
|   - adds an nucleotide to an kmer
| Input:
|   - ntUC:
|     o character with nucleotide to add
|   - kmerUL:
|     o unsigned long with current kmer to add an
|       nucleotide to
|   - lenKmerUL:
|     o unsigned long with the number of nuceotdies since
|       the last anonymous nucleotide
|   - kmerMaskUL:
|     o unsigned long to mask out bits not used in the
|       kmer (each nucleotide takes two bits)
|       - make with mkKMerMask_kmerCnt (fun08)
| Output:
|   - Modifies:
|     o kmerUL to
|       - have the new nucleotide and masked by kmerMaskUL
|       - be def_noKmer_kmerCnt if ntUC is an anonymous
|         base (nucleotide)
|     o lenKmerUL
|       - be incurmented by 1 (+1)
|       - 0 if ntUC is an anonymous base (nucleotide)
\-------------------------------------------------------*/
#define \
addNtToKmer_kmerCnt( \
   ntUC,     /*nucleotide to add*/ \
   kmerUL,   /*kmer to add nucleotide to*/ \
   lenKmerUL,/*number bases since last anonymous base*/ \
   kmerMaskUL/*bases to mask at end*/ \
){ \
   unsigned long ntMacUL = 0; \
   \
   ntMacUL = ntTo2Bit[(unsigned char) (ntUC)];\
   (kmerUL) <<= def_bitsPerKmer_kmerCnt; \
   (kmerUL) |= ntMacUL; \
   \
   if(ntMacUL != def_err3rdBit_ntTo2Bit) \
   { /*If: I have no anonymous baseses*/ \
      (kmerUL) &= kmerMaskUL; \
      ++(lenKmerUL); \
   } /*If: I have no anonymous baseses*/ \
   \
   else \
   { /*Else: I have an anonymous base*/ \
      (kmerUL) = def_noKmer_kmerCnt; \
      (lenKmerUL) = 0; \
   } /*Else: I have an anonymous base*/ \
} /*addNtToKmer_kmerCnt*/

/*-------------------------------------------------------\
| Fun10: addSeq_kmerCnt
|   - adds a sequence to a kmerCnt structure
| Input:
|   - kmerCntSTPtr:
|     o pointer to kmerCnt structure to add sequence to
|   - seqSTPtr:
|     o pointer to a seqST with sequence to add
| Output:
|   - Modifies:
|     o kmerCntSTPtr to have the new sequence
|   - Returns:
|     o 0 for no errors
|     o def_memErr_kmerCnt for memory errors
\-------------------------------------------------------*/
signed char
addSeq_kmerCnt(
   struct kmerCnt *kmerCntSTPtr,/*will hold new sequence*/
   struct seqST *seqSTPtr   /*sequence to copy*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC:
   '   - adds a sequence to a kmerCnt structure
   '   o fun10 sec01:
   '     - variable declerations
   '   o fun10 sec02:
   '     - initialize and add sequences to structure
   '   o fun10 sec03:
   '     - add kmers to table
   '   o fun10 sec04:
   '     - return after success or error
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uchar errUC = 0;
   ulong ulNt = 0;      /*nucleotide on*/

   /*for getting kmers*/
   ulong maskUL = 0;    /*for removing extra nt*/

   uchar forNtUC = 0;   /*nucleotide to add to kmer*/
   ulong forKmerUL = 0; /*holds forward kmer*/
   ulong forLenUL = 0;  /*# nt since last anonymous nt*/

   uchar revNtUC = 0;   /*nucleotide to add to kmer*/
   ulong revKmerUL = 0; /*holds reverse kmer*/
   ulong revLenUL = 0;  /*# nt since last anonymous nt*/

   struct seqST *tmpSeqSTPtr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^   - initialize and add sequences to structure
   ^   o fun10 sec02 sub01:
   ^     - set up mask and blank structure
   ^   o fun10 sec02 sub02:
   ^     - copy sequence (foward and reverse)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun10 Sec02 Sub01:
   *    - set up mask and blank structure
   \*****************************************************/

   maskUL = mkKmerMask_kmerCnt(kmerCntSTPtr->lenKmerUC);

   if(seqSTPtr == kmerCntSTPtr->forSeqST)
   { /*If: copying sequence for forward seqST*/ 
      tmpSeqSTPtr = kmerCntSTPtr->forSeqST;
      kmerCntSTPtr->forSeqST = 0;
   } /*If: copying sequence for forward seqST*/ 

   blank_kmerCnt(kmerCntSTPtr);

   if(tmpSeqSTPtr)
   { /*If: copying sequence for forward seqST*/ 
      kmerCntSTPtr->forSeqST = tmpSeqSTPtr;
      tmpSeqSTPtr = 0;
   } /*If: copying sequence for forward seqST*/ 

   /*****************************************************\
   * Fun10 Sec02 Sub02:
   *    - copy sequence (foward and reverse)
   \*****************************************************/

   if(! (seqSTPtr == kmerCntSTPtr->forSeqST))
   { /*If: seqST structure is not in kmerCntSTPtr*/
      errUC =
         cp_seqST(
            kmerCntSTPtr->forSeqST,
            seqSTPtr
         ); /*copy forward sequence*/

      if(errUC)
        goto memErr_fun08_sec0x;
   } /*If: seqST structure is not in kmerCntSTPtr*/

   errUC =
      cp_seqST(
         kmerCntSTPtr->revSeqST,
         seqSTPtr
      ); /*copy forward sequence*/

   if(errUC)
     goto memErr_fun08_sec0x;

   revComp_seqST(kmerCntSTPtr->revSeqST); 

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^   - add kmers to table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      ulNt = 0;
      ulNt < seqSTPtr->lenSeqUL;
      ++ulNt
   ){ /*Loop: build kmer table*/
      forNtUC =
         (uchar)
         kmerCntSTPtr->forSeqST->seqStr[ulNt];

      revNtUC =
         (uchar)
         kmerCntSTPtr->revSeqST->seqStr[ulNt];

      addNtToKmer_kmerCnt(
         forNtUC,    /*nucleotide to add*/
         forKmerUL,  /*kmer adding nucleotide to*/
         forLenUL,   /*length since last non-nucleotide*/
         maskUL      /*mask to remove old nucleotide*/
      );

      addNtToKmer_kmerCnt(
         revNtUC,    /*nucleotide to add*/
         revKmerUL,  /*kmer adding nucleotide to*/
         revLenUL,   /*length since last non-nucleotide*/
         maskUL      /*mask to remove old nucleotide*/
      );

      if(forLenUL >= kmerCntSTPtr->lenKmerUC)
      { /*If: had a forward kmer*/
         ++(kmerCntSTPtr->forKmerArySI[forKmerUL]);
         ++(kmerCntSTPtr->forKmersUI);
      } /*If: had a forward kmer*/

      if(revLenUL >= kmerCntSTPtr->lenKmerUC)
      { /*If: had a reverse kmer*/
         ++(kmerCntSTPtr->revKmerArySI[revKmerUL]);
         ++(kmerCntSTPtr->revKmersUI);
      } /*If: had a reverse kmer*/
   } /*Loop: build kmer table*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec04:
   ^   - return after success or error
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun08_sec0x:;

   return def_memErr_kmerCnt;
} /*addSeq_kmerCnt*/

/*-------------------------------------------------------\
| Fun11: ntToKmerAry_kmerCnt
|   - converts a nucleotide sequence to a array of kmer
|     counts
| Input:
|   - seqSTPtr:
|     o pointer to a seqST with sequence to get kmer
|       counts for
|   - lenKmerUC:
|     o length of one kmer
|   - kmerArySI:
|     o pointer to a signed int array to add kmer so
|     o needs to be size of ((lenKmerUC ^ 4) + 1)
|   - cntArySI:
|     o pointer to a signed int array to add kmer counts
|       to
|     o needs to be size of ((lenKmerUC ^ 4) + 1)
| Output:
|   - Modifies:
|     o kmerArySI to hold kmers in seqSTPtr
|     o cntArySI to hold number times each kmer happened
|     o sorts kmerArySI and cntArySI by kmer
|       - this converts the hash table to a list of kmers
|       - end will be marked with a -2
|   - Returns:
|     o number of kmers in sequence
\-------------------------------------------------------*/
signed int
ntToKmerAry_kmerCnt(
   struct seqST *seqSTPtr, /*sequence to convert*/
   unsigned char lenKmerUC,    /*length of one kmer*/
   signed int *kmerArySI,      /*will hold uniqe kmers*/
   signed int *cntArySI        /*gets # kmer duplicates*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC:
   '   - converts a nucleotide sequence to a array of kmer
   '     counts
   '   o fun11 sec01:
   '     - variable declerations
   '   o fun11 sec02:
   '     - initialize (blank) the array
   '   o fun11 sec03:
   '     - get kmer counts
   '   o fun11 sec04:
   '     - sort and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ulong ulNt = 0;      /*nucleotide on*/
   unsigned int maxKmersUI = 1;

   /*for getting kmers*/
   ulong maskUL = 0;    /*for removing extra nt*/
   ulong ulSeq = 0;     /*nucleotide on*/

   uchar ntUC = 0;   /*nucleotide to add to kmer*/
   ulong kmerUL = 0; /*holds forward kmer*/
   ulong lenUL = 0;  /*# nt since last anonymous nt*/

   sint retNumKmersSI = 0; /*number kmers in sequence*/
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^   - initialize (blank) the array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   maskUL = mkKmerMask_kmerCnt(lenKmerUC);

   for(
      ulNt = 0;
      ulNt < lenKmerUC;
      ++ulNt
   ) maxKmersUI <<= 2;

   for(
      ulSeq  = 0;
      ulSeq <= maxKmersUI;
      ++ulSeq
   ){ /*Loop: blank kmer array*/
      kmerArySI[ulSeq] = def_noKmer_kmerCnt;
      cntArySI[ulSeq] = 0;
   } /*Loop: blank kmer array*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec03:
   ^   - get kmer counts
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      ulSeq  = 0;
      ulSeq < seqSTPtr->lenSeqUL;
      ++ulSeq
   ){ /*Loop: convert kmers*/
      ntUC =
         (uchar)
         seqSTPtr->seqStr[ulSeq];

      addNtToKmer_kmerCnt(
         ntUC,       /*nucleotide to add*/
         kmerUL,     /*kmer adding nucleotide to*/
         lenUL,      /*length since last non-nucleotide*/
         maskUL      /*mask to remove old nucleotide*/
      );

      if(lenUL >= lenKmerUC)
      { /*If: have kmer*/
         ++retNumKmersSI;
         kmerArySI[kmerUL] = kmerUL; /*for sorting*/
         ++cntArySI[kmerUL];
      } /*If: have kmer*/
   } /*Loop: convert kmers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec04:
   ^   - sort and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   twoAryNum_shellSort(
       (uint *) kmerArySI,
       cntArySI,
       0,
       maxKmersUI
   ); /*sort kmers by kmer id*/

    return retNumKmersSI;
} /*ntToKmerAry_kmerCnt*/

/*-------------------------------------------------------\
| Fun12: get_kmerCnt
|   - gets number of matching kmers between a kmerCnt
|     structure and the kmer arrays
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCntSTPtr structure with the kmer
|       table to compare to
|   - kmerArySI:
|     o pointer to a singed in array with the unique kmers
|       in the sequence (use fun11; ntToKmerAry_kmerCnt)
|   - cntArySI:
|     o pointer to a singed in array with the number of
|       duplicates per unique kmers in the sequence (use
|       fun11; ntToKmerAry_kmerCnt)
| Output:
|   - Returns:
|     o highes number of kmers found
|       - + signed int if forward sequence more kmers
|       - - signed int if reverse sequence more kmers
\-------------------------------------------------------*/
signed int
get_kmerCnt(
   struct kmerCnt *kmerCntSTPtr, /*table to get counts*/
   signed int *kmerArySI,        /*sequence uniqe kmers*/
   signed int *cntArySI          /*sequence kmer counts*/
){
   sint forKmersSI = 0;
   sint revKmersSI = 0;

   uint uiKmer = 0;
   uint kmerUI = 0;

   while(kmerArySI[uiKmer] >= 0)
   { /*Loop: count kmer numbers*/
      kmerUI =
         (uint)
         kmerArySI[uiKmer];

      forKmersSI +=
         min_genMath(
            cntArySI[uiKmer],
            kmerCntSTPtr->forKmerArySI[kmerUI]
         );

      revKmersSI +=
         min_genMath(
            cntArySI[uiKmer],
            kmerCntSTPtr->revKmerArySI[kmerUI]
         );

      ++uiKmer;
   } /*Loop: count kmer numbers*/

   /*find highest count*/
   if(revKmersSI > forKmersSI)
      forKmersSI = -revKmersSI;

   return forKmersSI;
} /*get_kmerCnt*/

/*-------------------------------------------------------\
| Fun13: faToKmerCnt_kmerCnt
|   - converts a the sequences in a fasta file to a
|     kmerCnt array
| Input:
|   - faFileStr:
|     o path to and name of fasta file with sequnces
|   - lenKmerUC:
|     o length of one kmer
|   - numSeqUI:
|     o pointer to unsigned int to hold the number of
|       sequences in the kmerCnt array
|   - errSCPtr:
|     o pointer to signed char to hold the errors
| Output:
|   - Modfies:
|     o numSeqUI to have the number of sequences in the
|       returned kmerCnt array
|     o errSCPtr to hold errors
|       - 0 for no erors
|       - def_fileErr_kmerCnt if could not open faFileStr
|       - def_memErr_kmerCnt for memory errors
\-------------------------------------------------------*/
struct kmerCnt *
faToKmerCnt_kmerCnt(
   signed char *faFileStr,/*path to fasta file*/
   unsigned char lenKmerUC, /*length of one kmer*/
   unsigned int *numSeqUI,/*will hold number of sequnces*/
   signed char *errSCPtr  /*holds error message*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun13 TOC:
   '   - converts a the sequences in a fasta file to a
   '     kmerCnt array
   '   o fun13 sec01:
   '     - variable declerations
   '   o fun13 sec02:
   '     - find number of sequences in fasta file
   '   o fun13 sec03:
   '     - allocate and initialize kmerCnt array
   '   o fun13 sec04:
   '     - add sequences to kmerCnt array
   '   o fun13 sec05:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;

   uint uiSeq = 0;
   struct kmerCnt *retHeapAryST = 0;

   struct seqST seqStackST;
   FILE *faFILE = 0;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec02:
   ^   - find number of sequences in fasta file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_seqST(&seqStackST);

   faFILE =
      fopen(
         (char *) faFileStr,
         "r"
      );

   if(! faFILE)
      goto fileErr_fun13_sec05_sub03;

   *numSeqUI = 0;
   uiSeq = 0;

   errSC =
      (schar)
      getFaSeq_seqST(
         faFILE,
         &seqStackST
      );   

   while(! errSC)
   { /*Loop: Count number of sequences in file*/
      ++(uiSeq);

      errSC =
         (schar)
         getFaSeq_seqST(
            (char *) faFILE,
            &seqStackST
         );   
   } /*Loop: Count number of sequences in file*/

   if(uiSeq < 1)
      goto fileErr_fun13_sec05_sub03;

   if(errSC == def_memErr_seqST)
      goto memErr_fun13_sec05_sub02;

   if(errSC != def_EOF_seqST)
      goto fileErr_fun13_sec05_sub03;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec03:
   ^   - allocate and initialize kmerCnt array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retHeapAryST = malloc(uiSeq * sizeof(struct kmerCnt));

   if(! retHeapAryST)
      goto memErr_fun13_sec05_sub02;

   for(
      *numSeqUI = 0;
      *numSeqUI < uiSeq;
      ++(*numSeqUI)
   ){ /*Loop: initialize and set up structures*/
      init_kmerCnt(&retHeapAryST[*numSeqUI]);

      errSC =
         setup_kmerCnt(
            &retHeapAryST[*numSeqUI],
            lenKmerUC
         );

      if(errSC)
         goto memErr_fun13_sec05_sub02;
   } /*Loop: initialize and set up structures*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec04:
   ^   - add sequences to kmerCnt array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uiSeq = 0;

   fseek(faFILE, 0, SEEK_SET);

   errSC =
      (schar)
      getFaSeq_seqST(
         faFILE,
         &seqStackST
      );   

   while(! errSC)
   { /*Loop: Count number of sequences in file*/
      errSC =
         addSeq_kmerCnt(
            &retHeapAryST[uiSeq],
            &seqStackST
         );

      if(errSC)
         goto memErr_fun13_sec05_sub02;

      errSC =
         (schar)
         getFaSeq_seqST(
            faFILE,
            &seqStackST
         ); /*already checked for errors*/

      ++uiSeq;
   } /*Loop: Count number of sequences in file*/

   /*already checked for file errors*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec05:
   ^   - clean up and return
   ^   o fun13 sec05 sub01:
   ^     - clean up after no errors
   ^   o fun13 sec05 sub02:
   ^     - memory error clean
   ^   o fun13 sec05 sub03:
   ^     - file error clean
   ^   o fun13 sec05 sub04:
   ^     - clean up return varible (errors only)
   ^   o fun13 sec05 sub05:
   ^     - clean up variables always cleaned up
   ^       (error/no error)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
   /*****************************************************\
   * Fun13 Sec05 Sub01:
   *   - clean up after no errors
   \*****************************************************/

   goto cleanUp_fun13_sec05_sub05;

   /*****************************************************\
   * Fun13 Sec05 Sub02:
   *   - memory error clean
   \*****************************************************/

   memErr_fun13_sec05_sub02:;

   *errSCPtr = def_memErr_kmerCnt;

   goto errCleanUp_fun13_sec05_sub04;

   /*****************************************************\
   * Fun13 Sec05 Sub03:
   *   - file error clean
   \*****************************************************/

   fileErr_fun13_sec05_sub03:;

   *errSCPtr = def_fileErr_kmerCnt;

   goto errCleanUp_fun13_sec05_sub04;

   /*****************************************************\
   * Fun13 Sec05 Sub04:
   *   - clean up return varible for errors
   \*****************************************************/

   errCleanUp_fun13_sec05_sub04:;

   freeHeapAry_kmerCnt(
      retHeapAryST,
      *numSeqUI
   );

   retHeapAryST = 0;

   goto cleanUp_fun13_sec05_sub05;

   /*****************************************************\
   * Fun13 Sec05 Sub05:
   *   - clean up variables always cleaned up
   *     (error/no error)
   \*****************************************************/

   cleanUp_fun13_sec05_sub05:;

   freeStack_seqST(&seqStackST);

   if(
         faFILE
      && faFILE != stdin
      && faFILE != stdout
   ) fclose(faFILE);

   faFILE = 0;

   return retHeapAryST;
} /*faToKmerCnt_kmerCnt*/

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
