/*#######################################################\
\#######################################################*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samClust.h"

/*.h files only*/
#include "../generalLib/dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .h  #include "../generalLib/ulCp.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef PLAN9
void
#else
int
#endif
main(
){
   schar *samFileStr =
      (schar *)
      "/home/reason/files/testing/diflu/map-cal2009-filt.sam";

   uint numReadsUI = 0;
   FILE *samFILE = 0;

   samFILE = 
      fopen(
         (char *) samFileStr,
         "r"
      );
       
   numReadsUI = getNumReads_samClust(samFILE);
   printf("reads: %u\n", numReadsUI);
   exit(0); 
}
