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

SOURCE=\
   ../generalLib/seqST.c \
   ../generalLib/samEntry.c \
   ../generalLib/kmerCnt.c \
   ../generalAln/alnSet.c \
   ../generalAln/dirMatrix.c \
   ../water/water.c \
   ../getDICoordsSrc/diCoords.c \
   diScan.c \
   findDIFrag.c

NAME=findDIFrag

all:
	$(CC) $(CFLAGS) $(SOURCE) -o $(NAME)

mac:
	$(CC) $(MACCFLAGS) $(SOURCE) -o $(NAME)

check:
	$(CC) $(DEBUGCFLAGS) $(SOURCE) -o $(NAME)

bug:
	$(CC) $(DEBUGCFLAGS) $(SOURCE) -o $(NAME)
	gdb -x bug-findDIFrag.txt $(NAME)

install:
	mv $(NAME) $(PREFIX)
	chmod a+x $(PREFIX)/$(NAME)
