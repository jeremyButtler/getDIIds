CC=cc
PREFIX=/usr/local/bin

CFLAGS=\
   -O3 \
   -std=c89 \
   -static \
   -Wall \
   -Wextra

MACCFLAGS=\
   -O3

BUGCFLAGS=\
   -O0 \
   -std=c89 \
   -static \
   -Wall \
   -Wextra \
   -ggdb \
   -g

SOURCE=\
   ../generalAln/alnSet.c \
   ../generalLib/seqST.c \
   ../generalLib/samEntry.c \
   ../generalLib/kmerCnt.c \
   ../generalAln/dirMatrix.c \
   memwater.c \
   alnMemwater.c

NAME=alnMemwater

all:
	$(CC) $(CFLAGS) $(SOURCE) -o $(NAME)

mac:
	$(CC) $(MACCFLAGS) $(SOURCE) -o $(NAME)

check:
	$(CC) $(BUGCFLAGS) $(SOURCE) -o $(NAME)

bug:
	$(CC) $(BUGCFLAGS) $(SOURCE) -o $(NAME)
	gdb -x bug-alnMemwater.txt $(NAME)

install:
	mv alnwater $(PREFIX)
	chmod $(PREFIX)/alnwater
