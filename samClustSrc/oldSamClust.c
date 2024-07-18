/*#######################################################\
# Name: samClust
#   - holds functions and structures to support read
#     clustering
#   - for know i am keeping this simple and only doing
#     snp and large deletion events (no large insertions)
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o .h st01: event_samClust
'     - is an event (base/ins/del) in a sam clust sequence
'   o .h st02: samClust
'     - holds the event sequences i am clustering
'   o .h fun01: blank_event_samClust
'     - sets all variables to 0 in a event_samClust struct
'   o .h Fun02: init_event_samClust
'     - initializes (sets severything to null) an
'        event_samClust structure
'   o fun03: freeStack_event_samClust
'     - frees all variables in an event_samClust structure 
'   o fun04: freeHeap_event_samClust
'     - frees an event_samClust structure
'   o fun05: freeHeapList_event_samClust
'     - frees an list of heap allocated event structures
'   o fun06: mk_event_samClust
'     - make an heap allocated event_samClust structure
'   o .h fun07: blank_samClust
'     - blanks a samClust struct
'   o .h fun08: init_samClust
'     - initializes a samClust struct (null/0)
'   o fun09: freeStack_samClust
'     - frees variables in an samClust structure
'   o fun10: freeHeap_samClust
'     - frees a samClust structure
'   o fun11: setup_samClust
'     - allocates memory for a samClust structure
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

#include "samClust.h"

#include "../generalLib/samEntry.h

/*These have no .c files*/
#include "dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .h  #include "../generalLib/base10str.h"
!   o .h  #include "../generalLib/ulCp.h"
!   o .h  #include "../generalLib/numToStr.h"
!   o .h  #include "../generalLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun03: freeStack_event_samClust
|   - frees all variables in an event_samClust structure 
| Input:
|   - eventSTPtr:
|     o pointer to an event_samClust struct with variables
|       to free
| Output:
|   - Frees:
|     o nothing; only variables are list related
|   - Modifies:
|     o blanks evenSTPtr with blank_event_samClust (fun01)
\-------------------------------------------------------*/
void
freeStack_event_samClust(
   struct event_samClust *eventSTPtr
){
   if(eventSTPtr)
   { /*If: a structure was input*/
      blank_event_samClust(eventSTPtr);
   } /*If: a structure was input*/
} /*freeStack_event_samClust*/

/*-------------------------------------------------------\
| Fun04: freeHeap_event_samClust
|   - frees an event_samClust structure
| Input:
|   - eventSTPtr:
|     o pointer to an event_samClust struct to free
| Output:
|   - Frees:
|     o frees eventSTPtr
| Warning:
|   - this does not handle linked lists
\-------------------------------------------------------*/
void
freeHeap_event_samClust(
   struct event_samClust *eventSTPtr
){
   freeStack_event_samClust(eventSTPtr);
   free(eventSTPtr);
} /*freeHeap_event_samClust*/

/*-------------------------------------------------------\
| Fun05: freeHeapList_event_samClust
|   - frees an list of heap allocated event structures
| Input:
|   - eventSTPtr:
|     o pointer to an event_samClust struct list to free
|   - freeSeqBl:
|     o 1: free sequence events
|     o 0: only free alternative events
| Output:
|   - Frees:
|     o frees all structures in eventSTPtr list
| Warning:
|   - this will generate an infinite loop for circular
|     lists
\-------------------------------------------------------*/
void
freeHeapList_event_samClust(
   struct event_samClust *eventSTPtr,
   schar freeSeqBl
){
   if(eventSTPtr)
   { /*If: have an structure to free*/
      if(
            eventSTPtr->nextEvent
         && freeSeqBl
      ){ /*If: there is another event*/
         freeHeapList_event_samClust(
            eventSTPtr->nextEvent,
            1
         );

         eventSTPtr->nextEvent = 0; /*no longer present*/
      } /*If: there is another event*/

      if(eventSTPtr->altEvent)
      { /*If: there is an alternative event*/
         freeHeapList_event_samClust(
            eventSTPtr->altEvent,
            freeSeqBl
         );

         eventSTPtr->altEvent = 0; /*no longer present*/
      } /*If: there is an alternative event*/

      freeStack_event_samClust(eventSTPtr);
      free(eventSTPtr);
   } /*If: have an structure to free*/
} /*freeHeapList_event_samClust*/

/*-------------------------------------------------------\
| Fun06: mk_event_samClust
|   - make an heap allocated event_samClust structure
| Input:
| Output:
|   - Returns:
|     o pointer to new initialized event_samClust struct
|     o 0 for memory errors
\-------------------------------------------------------*/
struct event_samClust
mk_event_samClust(
){
   struct event_samClust *retHeapST = 0;

   retHeapST = malloc(sizeof(struct event_samClust));

   if(retHeapST)
      init_event_samClust(retHeapST);
      
   return retHeapST;
} /*mk_event_samClust*/

/*-------------------------------------------------------\
| Fun09: freeStack_samClust
|   - frees variables in an samClust structure
| Input:
|   - samClustSTPtr:
|     o pointer to an event_samClust struct with variables
|       to free
| Output:
|   - Frees:
|     o frees numEventsAryUI and eventAryST in samClustPtr
|   - Modifies:
|     o sets lenRefUI in samClustPtr to 0
\-------------------------------------------------------*/
void
freeStack_samClust(
   struct samClust *samClustSTPtr
){
   uint uiPos = 0;

   if(samClustSTPtr)
   { /*If: have a samClust structure*/
      free(samClustSTPtr->numEventsAryUI);
      samClustSTPtr->numEventsAryUI = 0;

      if(samClustSTPtr->eventAryST)
      { /*If: have an event matrix to free*/
         for(
            uiPos = 0;
            uiPos < samClustSTPtr->lenRefUI;
            ++uiPos
         ){ /*Loop: free the event list*/
            freeHeapList_event_samClust(
               samClustSTPtr->eventAryST->altEvent
            ); /*free the event linked list*/

            freeStack_event_samClust(
               samClustSTPtr->eventAryST
            ); /*free the array event*/
         } /*Loop: free the event list*/

         free(samClustSTPtr->eventAryST);
         samClustSTPtr->eventAryST = 0;
      } /*If: have an event matrix to free*/
   } /*If: have a samClust structure*/

   blank_samClust(samClustSTPtr);
} /*freeStack_samClust*/

/*-------------------------------------------------------\
| Fun10: freeHeap_samClust
|   - frees a samClust structure
| Input:
|   - samClustSTPtr:
|     o pointer to an event_samClust struct to free
| Output:
|   - Frees:
|     o frees samClustSTPtr and its internal variables
\-------------------------------------------------------*/
void
freeHeap_samClust(
   struct samClust *samClustSTPtr
){
   freeStack_samClust(samClustSTPtr);
   free(samClustSTPtr);
} /*freeHeap_samClust*/

/*-------------------------------------------------------\
| Fun11: setup_samClust
|   - allocates memory for a samClust structure
| Input:
|   - samClustSTPtr:
|     o pointer to an event_samClust struct to setup
|   - refLenUI:
|     o length of the reference sequence this is for
| Output:
|   - Modifies:
|     o lenRefUI in samClustSTPtr to have the reference
|       length
|     o numEventsAryUI in samClustSTPtr to be the length
|       of reference
|     o evetnAryST in samClustSTPtr to be the length
|       of reference
|   - Returns:
|     o 0 for no errors
|     o def_memErr_samClust for memory errors
\-------------------------------------------------------*/
signed char
setup_samClust(
   struct samClust *samClustSTPtr,
   unsigned int refLenUI
){
   schar errSC = 0;
   uint uiPos = 0;

   samClustSTPtr->lenRefUI = refLenUI;

   samClustSTPtr->numEventsAryUI =
      calloc(refLenUI, sizeof(uint));

   if(! samClustSTPtr->numEventsAryUI)
      goto memErr_fun11;

   samClustSTPtr->eventAryST =
      malloc(refLenUI * sizeof(struct event_samClust));

   if(! samClustSTPtr->eventAryST)
      goto memErr_fun11;

   for(
      uiPos = 0;
      uiPos < refLenUI;
      ++uiPos
   )init_event_samClust(samClustSTPtr->eventAryST[uiPos]);

   errSC = 0;
   goto ret_fun11;

   memErr_fun11:;
   errSC = def_memErr_samClust;
   goto ret_fun11;

   ret_fun11:;

   return errSC
} /*freeHeap_samClust*/

/*-------------------------------------------------------\
| Fun12: addRead_samClust
|   - adds a read to a samClust structure
| Input:
|   - samClustSTPtr:
|     o pointer to an event_samClust struct to add read to
|   - samSTPtr:
|     o pointer to samEntr structure with read to add
| Output:
|   - Modifies:
|     o samClustSTPtr numEventsAryUI and eventAryST to
|       have the new read.
|   - Returns:
|     o 0 for no errors
|     o def_memErr_samClust for memory errors
\-------------------------------------------------------*/
signed char
addRead_samClust(
   struct samClust *samClustSTPtr, /*add read to*/
   struct samEntry *samSTPtr       /*has read to add*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun12 TOC:
   '   - adds a read to a samClust structure
   '   o fun12 sec01:
   '     - variable declarations
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint uiCig = 0;   /*cigar entry on*/
   uint uiEvent = 0; /*event on the cigar*/
   uint uiNt = 0;    /*nucleotide on*/

   struct event_samClust *eventSTPtr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec02:
   ^   - scan for snps and dels 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   
} /*addRead_samClust*/

