CC=cc
PREFIX=/usr/local/bin

CFLAGS=\
   -static \
   -O3 \
   -Wall \
   -Wextra \
   -std=c89

MACCFLAGS=\
   -O3

DEBUGCFLAGS=\
   -static \
   -O0 \
   -Wall \
   -Wextra \
   -std=c89 \
   -g \
   -ggdb

primFindSOURCE=\
   kmerFind.c \
   memwater/seqST.c \
   memwater/alnSetST.c \
   memwater/memwater.c \
   inputPrimFind.c \
   primFind.c

SOURCE=\
   kmerFind.c \
   memwater/seqST.c \
   memwater/alnSetST.c \
   memwater/memwater.c \
   fluST.c \
   inputGetDIIds.c \
   getDIIds.c

all:
	$(CC) $(CFLAGS) $(SOURCE) -o getDIIds

mac:
	$(CC) $(MACCFLAGS) $(SOURCE) -o getDIIds

check:
	$(CC) $(DEBUGCFLAGS) $(SOURCE) -o getDIIds

buildPrim:
	$(CC) $(CFLAGS) $(primFindSOURCE) -o primFind

buildPrimMac:
	$(CC) $(CFLAGS) $(primFindSOURCE) -o primFind

checkPrimFind:
	$(CC) $(DEBUGCFLAGS) $(primFindSource) -o primFind

install:
	mv getDIIds $(PREFIX)
	chmod a+x $(PREFIX)/getDIIds

primFindInstall:
	mv primFind $(PREFIX)
	chmod a+x $(PREFIX)/primFind
