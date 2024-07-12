/*########################################################
# Name: samEntryStruct
#  - Holds structer to hold a sam file entry. This also
#    includes the functions needed to support this
#    structure.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - header guards and definitions
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
'  o fun06: makeSamEntry
'    - Makes an heap allocated samEntry structure
'  o .h fun07: qhistToMed_samEntry
'    - Gets the median q-score for an histogram of
'       q-scores in a samStruct
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
'  o fun13: p_samEntryAsFastq
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
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Header gaurds and definitions
\-------------------------------------------------------*/

#ifndef SAMENTRYSTRUCT_H
#define SAMENTRYSTRUCT_H

#define def_ajdQ_samEntry 33 /*offest to get q-score of 0*/
#define def_maxQ_samEntry 94 /*highest possible q-score*/

/*-------------------------------------------------------\
| ST01: samEntry
|  - Holds a single samfile entry
\-------------------------------------------------------*/
typedef struct samEntry
{ /*samEntry*/
    /*Buffers for storing strings*/
    char qryIdStr[128];/*Holds the query id/name*/
    unsigned char lenQryIdUC;  /*Length of query id/name*/

    char refIdStr[128];/*Holds the reference id/name*/
    unsigned char lenRefIdUC;  /*Length of ref id/name*/

    char *cigTypeStr;  /*Holds the cigar type entry*/
    int *cigValAryI;   /*Holds the cigar number*/
    unsigned int lenCigUI;     /*Length of cigar entry*/
    unsigned int lenCigBuffUI;/*# bytes malloc to cigStr*/

    char rNextStr[128];/*Holds the rNext entry*/
    unsigned char lenRNextUC;  /*Length of rNextStr*/

    /*The sequence and q-scores entries have lengths
    ` stored in readLenUI
    */
    char *seqStr;      /*Holds the sequence entry*/
    unsigned int lenSeqBuffUI;/*bytes malloc seqStr/qStr*/

    char *qStr;        /*Holds the q-score entry*/
    unsigned int lenQBuffUI; /*number q-score entries*/

    char *extraStr;    /*Extra entries (columns 12-end)*/
    unsigned int lenExtraUI;   /*Length of extraStr*/
    unsigned int lenExtraBuffUI;/*# bytes malloc to qStr*/

    /*Flags/single numeric values in the sam entry*/
    unsigned char mapqUC;      /*Holds mapping quality*/
    unsigned short flagUS;     /*Holds the flag*/
    int pNextI;        /*PNEXT postion of mate/next read*/
    int tLenI;         /*TLEN observed template length*/

    /*Stats*/

    /*The q-scores; includes the softmasked regions*/
    float medianQF;    /*holds median read q-score*/
    float meanQF;      /*holds mean read q-score*/
 
    unsigned int refStartUI;  /*First reference base*/
    unsigned int refEndUI;    /*Last reference base*/
    unsigned int readLenUI;   /*Holds read length*/
    unsigned int alnReadLenUI;/*Number ref bases aligned*/

    unsigned int numMatchUI;/*Holds number of matches*/
    unsigned int numSnpUI;  /*Holds number of mismatches*/
    unsigned int numInsUI;  /*Holds number of insertions*/
    unsigned int numDelUI;  /*number of deletions*/
    unsigned int numMaskUI; /*number soft masked bases*/

    /*These variables are used in finding the q-scores*/
    unsigned int qHistUI[def_maxQ_samEntry];
    unsigned long sumQUL;             /*Total for mean Q*/
}samEntry;

/*-------------------------------------------------------\
| Fun01: blank_samEntry
| Use:
|  - Sets all values to 0, or for c-strings to '\0'
| Input:
|  - samSTPtr:
|    - Pointer to samEntry structure to blank
| Output:
|  - Modifies:
|    o Sets every variable to 0
|    o the first character in every c-string to be '\0'
\-------------------------------------------------------*/
#define \
blank_samEntry( \
   samSTPtr \
){\
    unsigned int iIter = 0;\
    \
    (samSTPtr)->qryIdStr[0] = '\0';/*query id/name*/\
    (samSTPtr)->lenQryIdUC = 0;    /*Length of query id*/\
    \
    (samSTPtr)->refIdStr[0] = '\0';/*reference id/name*/\
    (samSTPtr)->lenRefIdUC = 0;  /*Length; reference id*/\
    \
    (samSTPtr)->lenCigUI = 0;     /*Length; cigar entry*/\
    \
    for( \
       iIter=0; \
       iIter < (samSTPtr)->lenCigBuffUI; \
       ++iIter \
    ){ /*Loop: Clear all cigar entries*/\
       (samSTPtr)->cigTypeStr[iIter] = '\0';\
       (samSTPtr)->cigValAryI[iIter] = 0;\
    } /*Loop: Clear all cigar entries*/\
    \
    (samSTPtr)->rNextStr[0] = '\0'; /*rNext entry*/\
    (samSTPtr)->lenRNextUC = 0;     /*Length of rNext*/\
    \
    if((samSTPtr)->seqStr)\
       (samSTPtr)->seqStr[0] = '\0';\
    \
    if((samSTPtr)->qStr)\
       (samSTPtr)->qStr[0] = '\0';\
    \
    if((samSTPtr)->extraStr)\
       (samSTPtr)->extraStr[0] = '\0';\
    \
    (samSTPtr)->lenExtraUI = 0;\
    \
    /*Flags/single numeric values in the sam entry*/\
    (samSTPtr)->mapqUC = 0;      /*mapping quality*/\
    (samSTPtr)->flagUS = 0;      /*flag in samEntryCStr*/\
    (samSTPtr)->pNextI = 0;      /*PNEXT postion*/\
    (samSTPtr)->tLenI = 0;       /*TLEN*/\
    \
    /*Stats*/\
    \
    /*The q-scores excluded the softmasked regions*/\
    (samSTPtr)->medianQF = 0;    /*median read q-score*/\
    (samSTPtr)->meanQF = 0;      /*mean read q-score*/\
    \
    (samSTPtr)->refStartUI = 0;  /*First reference base*/\
    (samSTPtr)->refEndUI = 0;    /*Last reference base*/\
    (samSTPtr)->readLenUI = 0;   /*Holds read length*/\
    (samSTPtr)->alnReadLenUI = 0;/*Aligned length*/\
    \
    (samSTPtr)->numMatchUI = 0;  /*number of matches*/\
    (samSTPtr)->numSnpUI = 0;    /*number of mismatches*/\
    (samSTPtr)->numInsUI = 0;    /*number of insertions*/\
    (samSTPtr)->numDelUI = 0;    /*number of deletions*/\
    (samSTPtr)->numMaskUI = 0;   /*number soft masked*/\
    \
    /*These variables are used in finding the q-scores*/\
    for(iIter = 0; iIter < def_maxQ_samEntry; ++iIter)\
       (samSTPtr)->qHistUI[iIter] = 0;\
    \
    (samSTPtr)->sumQUL = 0;\
} /*blank_samEntry*/

/*-------------------------------------------------------\
| Fun02: init_samEntry
|  - Initializes a samEntry structure for use. This 
|    function should only ever be called once per
|    structure or after freeStack_samEntry has been used.
|    Use blank_samEntry to reset the values in a samEntry
|    struct after initialization.
| Input:
|  - samSTPtr:
|    o pointer to samEntry struct to initialize
| Output:
|  - Modifies:
|    o all variables in samEntry to be 0 or null
|  - Returns:
|    o 0 for success
|    o 1 for memory error
\-------------------------------------------------------*/
#define \
init_samEntry( \
   samSTPtr \
){ \
    (samSTPtr)->seqStr = 0; \
    (samSTPtr)->lenSeqBuffUI = 0; \
    \
    (samSTPtr)->qStr = 0; \
    (samSTPtr)->lenQBuffUI = 0; \
    \
    (samSTPtr)->cigTypeStr = 0; \
    (samSTPtr)->lenCigBuffUI = 0; \
    (samSTPtr)->cigValAryI = 0; \
    \
    (samSTPtr)->extraStr = 0; \
    (samSTPtr)->lenExtraBuffUI = 0; \
    \
    blank_samEntry(samSTPtr); \
} /*init_samEntry*/

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
);

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
);

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
void
freeHeap_samEntry(
   struct samEntry **samSTPtr
);

/*-------------------------------------------------------\
| Fun06: makeSamEntry
|  - Makes an heap allocated samEntry structure
| Input:
| Output:
|  - Returns:
|    o Pointer to a samEntry structure
|    o 0 if had an memory error
\-------------------------------------------------------*/
struct samEntry * makeSamEntry();

/*-------------------------------------------------------\
| Fun07: qhistToMed_samEntry
|   - Gets the median q-score for an histogram of q-scores
|     in a samStruct
| Input:
|   - samSTPTr:
|     o Pointer to samEntry struct with a q-score
|       histogram to find the median q-score for
| Output:
|   - Modifies:
|     o samSTPtr->medianQF to have the median q-score
\-------------------------------------------------------*/
#define \
qhistToMed_samEntry( \
   samSTPtr \
){\
    unsigned int numBasesUI = 0;\
    unsigned int midPointUL = (samSTPtr)->readLenUI >> 1;\
    unsigned int uiQ = 0;\
    \
    for(uiQ = 0; uiQ < def_maxQ_samEntry; uiQ++)\
    { /*Loop: through q-score histogram for midpoint*/\
        numBasesUI += (samSTPtr)->qHistUI[uiQ];\
        \
        if(numBasesUI >= midPointUL)\
        { /*if found the midpoint, then find the median*/\
            if(numBasesUI > midPointUL || numBasesUI & 1)\
                (samSTPtr)->medianQF = (float) uiQ;\
            \
            else\
            { /*Else: even & 2 differnt Q's at midpoint*/\
                numBasesUI = uiQ;\
                ++uiQ;\
                \
                while((samSTPtr)->qHistUI[uiQ]==0) ++uiQ;\
                \
                (samSTPtr)->medianQF =\
                    (numBasesUI + uiQ) / ((float) 2);\
            } /*Else: even & 2 differnt Q's at midpoint*/\
            \
            break;\
        } /*if found the midpoint, then find the median*/\
    } /*Loop: through q-score histogram for midpoint*/\
} /*qhistToMed_samEntry*/

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
);

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
);

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
|    o samFILE to be on the next line
|    o buffStr to hold a sam file line (resized if needed)
|    o lenBuffUL to hold the resized length of buffStr
|  - Returns:
|    o 0 for success
|    o 1 for EOF (End Of File)
|    o 64 for memory errors
\-------------------------------------------------------*/
char get_samLine(
   struct samEntry *samSTPtr,
   char **buffStr,
   unsigned long *lenBuffUL,
   void *samVoidFILE
);

/*-------------------------------------------------------\
| Fun11: findRefPos_samEntry
|   - Find an reference coordinate in an sequence in
|     an sam entry structure
| Input:
|   - samSTPtr:
|     o Pointer to samEntry structu with sequence to find
|       the reference position in
|   - siCig:
|     o Index of the cigar entry I am on
|   - cigBaseOnSI:
|     o Number of bases left in the cigar entry
|   - targPosSI:
|     o Reference position I wish to find
|   - refPosSI:
|     o Position I am currently at in the reference
|       sequence
|   - seqPosSI:
|     o Position I am currently at in the sequence in
|       samSTPtr
| Output:
|   - Modifies:
|     o siCig to point to the next open cigar entry
|       - will be > samSTPtr->lenCigUI when the sequence
|         does not end at at targPosSI
|     o cigBaseOnSI to have the number of bases remianing
|       in the current siCig entry
|     o refPosSI to be the index to the position in targSI
|     o seqPosSI to be the index of the sequence base at
|       refPosSI (check cigar to see if deletion)
|   - Returns:
|     o Index of last cigar entry I incurmented
\-------------------------------------------------------*/
#define \
findRefPos_samEntry(\
   samSTPtr,    /*Sam entry structure with cigar*/\
   siCig,       /*Current cigar entry on*/\
   cigBaseOnSI, /*Number of base left in cigar entry*/\
   targPosSI,   /*Position looking for*/\
   refPosSI,    /*Current reference position*/\
   seqPosSI     /*Current sequence position*/\
)({/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ` Fun11 TOC: findRefPos_samEntry\
   `   - Find an reference coordinate in an sequence in\
   `     an sam entry structure\
   `   o fun11 sec01:\
   `     - Start loop and check insertions/soft masking\
   `   o fun11 sec02:\
   `     - Move position in deletion cases\
   `   o fun11 sec03:\
   `     - Move position for snp/match cases\
   `   o fun11 sec04:\
   `     - Move to the next cigar entry\
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/\
   \
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:\
   ^   - Start loop and check insertions/soft masking\
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
   \
   int lastCigMacSI = siCig;\
   \
   /*Check if I did a full cigar entry move*/\
   while(refPosSI < targPosSI)\
   { /*Loop: till I am on the target base*/\
      lastCigMacSI = (siCig);\
      \
      switch((samSTPtr)->cigTypeStr[(siCig)])\
      { /*Switch: check what the next entry is*/\
         case 'S':\
         case 'I':\
         /*Case: Softmasking or insertions*/\
            (seqPosSI) += (cigBaseOnSI);\
            ++(siCig);\
            (cigBaseOnSI) = 0;\
            break;\
         /*Case: Softmasking or insertions*/\
         \
         /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
         ^ Fun11 Sec02:\
         ^   - Move position in deletion cases\
         \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
         \
         case 'D':\
         /*Case: Deletion*/\
            (refPosSI) += (cigBaseOnSI);\
            \
            /*This is the case most of the time, so\
            `   branchless will be slower than if\
            */\
            if((refPosSI) <= (targPosSI))\
            { /*If: I have not found target position*/\
               ++(siCig);\
               (cigBaseOnSI) = 0;\
            } /*If: I have not found target position*/\
            \
            else \
            { /*Else: I overshot the target*/\
               (cigBaseOnSI) =\
                  (int) (\
                       (refPosSI)\
                     - (targPosSI)\
                  ); /*Find how many bases overshot by*/\
               \
               /*Make corrections for overshooting*/\
               (refPosSI) -= (cigBaseOnSI);\
            } /*Else: I overshot the target*/\
            \
            break;\
         /*Case: Deletion*/\
         \
         /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
         ^ Fun11 Sec03:\
         ^   - Move position for snp/match cases\
         \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
         \
         case 'M':\
         case 'X':\
         case '=':\
         /*Case: match (M or =) or snp (M or X)*/\
            (refPosSI) += (cigBaseOnSI);\
            (seqPosSI) += (cigBaseOnSI);\
            \
            /*This is the case most of the time, so\
            `   branchless will be slower than if\
            */\
            if((refPosSI) <= (targPosSI))\
            { /*If: I have not found target position*/\
               (cigBaseOnSI) = 0;\
               ++(siCig);\
            } /*If: I have not found target position*/\
            \
            else \
            { /*Else: I overshot the target*/\
               (cigBaseOnSI) =\
                  (int) (\
                       (refPosSI)\
                     - (targPosSI)\
                  ); /*Find how many bases overshot by*/\
               \
               /*Make corrections for overshooting*/\
               (refPosSI) -= (cigBaseOnSI);\
               (seqPosSI) -= (cigBaseOnSI);\
            } /*Else: I overshot the target*/\
            \
            break;\
         /*Case: match (M or =) or snp (M or X)*/\
         \
         default: \
         /*Case: hard mask of some kind*/ \
            ++(siCig);\
            (cigBaseOnSI) = 0;\
         /*Case: hard mask of some kind*/ \
      } /*Switch: check what the next entry is*/\
      \
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun11 Sec03:\
      ^   - Move to the next cigar entry\
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\
      \
      if((siCig) >= (samSTPtr)->lenCigUI)\
         break; /*End of the sequence*/\
      \
      /*This case will be true most of the time, unless\
      `   the start has already been found\
      */\
      if((cigBaseOnSI) == 0)\
         (cigBaseOnSI) = (samSTPtr)->cigValAryI[(siCig)];\
   } /*Loop: till I am on the target base*/\
   \
   lastCigMacSI;\
}) /*findRefPos_samEntry*/

/*-------------------------------------------------------\
| Fun12: p_samEntry
| Use:
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
char p_samEntry(
   struct samEntry *samSTPtr,
   char **buffStr,
   unsigned long *lenBuffUL,
   char pNoNewLineBl,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun13: p_samEntryAsFastq
| Use:
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
void p_samEntryAsFastq(
   struct samEntry *samSTPtr,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun14: pfa_samEntry
| Use:
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
);

/*-------------------------------------------------------\
| Fun15: pstats_samEntry
| Use:
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
);

#endif

/*-------------------------------------------------------\
| Note01:
|   - Notes about the sam file format from the sam file
|     pdf
\-------------------------------------------------------*/

/*
Sam file table for first 11 columns (all sam files have)
| Col | Field |  Type  |        Brief description              |
|:---:|:-----:|:------:|:-------------------------------------:|
|  1  | QNAME | String |       Query template NAME             |
|  2  | FLAG  |  Int   |          bitwise FLAG                 |
|  3  | RNAME | String |     Reference sequence NAME           |
|  4  |  POS  |  Int   |  1-based leftmost mapping POSition    |
|  5  | MAPQ  |  Int   |          MAPping Quality              |
|  6  | CIGAR | String |            CIGAR string               |
|  7  | RNEXT | String | Reference name of the mate/next read  |
|  8  | PNEXT |  Int   |   Position of the mate/next read      |
|  9  | TLEN  |  Int   |      observed Template LENgth         |
| 10  |  SEQ  | String |          segment SEQuence             |
| 11  | QUAL  | String | ASCII of Phred-scaled base Quality+33 |
*/

/*eqx cigar entry
`   matches are '='
`   mimsatches are 'X'
'   non-eqx merges matches and mismatches into 'M'
`   deletions are 'D'
`   insertions are 'I'
`   soft masks are 'S'
`   Numbers come before entry (=, X, D, I, or S) & show
`     number of times the entry is repeated
`   Everything else is a hard mask (was removed) or
'      removed the bases from the sequence
`   EX: 10S5=1X701 (10 soft masked bases, 5 matches,
`       1 SNP, 701 matches)
*/

/* Sam file flag values tables
| Bit  | FLAG  |                        Description                                 |
|:----:|:-----:|:------------------------------------------------------------------:|
| 1    | 0x1   | template having multiple segments in sequencing                    |
| 2    | 0x2   | each segment properly aligned according to the aligner             |
| 4    | 0x4   | segment unmapped                                                   |
| 8    | 0x8   | next segment in the template unmapped                              |
| 16   | 0x10  | SEQ being reverse complemented                                     |
| 32   | 0x20  | SEQ of the next segment in the template being reverse complemented |
| 64   | 0x40  | the first segment in the template                                  |
| 128  | 0x80  | the last segment in the template                                   |
| 256  | 0x100 | secondary alignment                                                |
| 512  | 0x200 | not passing filters, such as platform/vendor quality controls      |
| 1024 | 0x400 | PCR or optical duplicate                                           |
| 2048 | 0x800 | supplementary alignment                                            |
*/

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
