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
'     - guards, forward declerations, defined variables
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
|   - guards, forward declerations, defined variables
\-------------------------------------------------------*/

#ifndef SAM_FILE_READ_CLUSTER_H
#define SAM_FILE_READ_CLUSTER_H

typedef struct samEntry samEntry;

#define def_memErr_samClust 2
#define def_fileErr_samClust 4

/*event types DO NOT CHANGE*/
#define def_noType_samClust 0
#define def_a_samClust 1
#define def_t_samClust 2
#define def_g_samClust 4
#define def_c_samClust 8
#define def_indel_samClust 16

/*-------------------------------------------------------\
| ST01: event_samClust
|   - is an event (base/ins/del) in a sam clust sequence
\-------------------------------------------------------*/
typedef
struct
event_samClust
{
   unsigned int posUI;              /*position in ref*/
   unsigned int numReadsUI;         /*# reads supporting*/
   signed char typeSC;              /*mutation type*/
   struct event_samClust *nextEvent;/*next event*/
   struct event_samClust *altEvent; /*alt events at pos*/
      /*other alternative events at this position*/
}event_samClust;

/*-------------------------------------------------------\
| ST02: samClust
|   - holds the event sequences i am clustering
\-------------------------------------------------------*/
typedef
struct
samClust
{
   uint lenRefUI = 0;    /*length of ref this is for*/
   uint *numEventsAryUI; /*number of events per position*/
   struct event_samClust *eventAryST;
      /*this will be the length of the reference. It is
      `   here for quick lookups
      */
}samClust;

/*-------------------------------------------------------\
| Fun01: blank_event_samClust
|   - sets all variables to 0 in an event_samClust struct
| Input:
|   - eventSTPtr:
|     o pointer to an event_samClust structure to blank
| Output:
|   - Modifies:
|     o posUI in eventSTptr to be 0
|     o typeSC in eventSTptr to be def_noType_samClust
\-------------------------------------------------------*/
#define \
blank_event_samClust( \
   eventSTPtr
){ \
   (eventSTPtr)->posUI = 0; \
   (eventSTPtr)->typeSC = 0; \
   (eventSTPtr)->numReadsUI = 0; \
} /*blank_event_samClust*/

/*-------------------------------------------------------\
| Fun02: init_event_samClust
|   - initializes (sets severything to null) an
|      event_samClust structure
| Input:
|   - eventSTPtr:
|     o pointer to an event_samClust struc to initialize
| Output:
|   - Modifies:
|     o All variables in eventSTPtr to by null/0/defaults
\-------------------------------------------------------*/

#define \
blank_event_samClust( \
   eventSTPtr
){ \
   (eventSTPtr)->nextEvent = 0; \
   (eventSTPtr)->altEvent = 0; \
   blank_event_samClust((eventSTPtr)); \
} /*blank_event_samClust*/

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
);

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
);

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
);

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
);
/*-------------------------------------------------------\
| Fun07: blank_samClust
|   - blanks a samClust struct (non-pointer variables set
|     to 0)
| Input:
|   - samClustSTPtr:
|     o pointer to a samClust structer to blank
| Output:
|   - Modifies:
|     o lenRefUI in samClustSTPtr to be 0
\-------------------------------------------------------*/
#define \
blank_samClust( \
   samClustSTPtr \
){ \
   (samClustSTPtr)->lenRefUI = 0; \
} /*blank_samClust*/

/*-------------------------------------------------------\
| Fun08: init_samClust
|   - initializes a samClust struct (null/0)
| Input:
|   - samClustSTPtr:
|     o pointer to a samClust structer to initialize
| Output:
|   - Modifies:
|     o all pointers and variables to by null
\-------------------------------------------------------*/
#define \
init_samClust( \
   samClustSTPtr \
){ \
   (samClustSTPtr)->numEventsAryUI = 0; \
   (samClustSTPtr)->eventAryST = 0; \
   blank((samClustSTPtr)); \
} /*init_samClust*/

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
);

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
);

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

#endif

