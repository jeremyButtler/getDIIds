CC=cc
PREFIX=/usr/local/bin
NAME=samClust

CFLAGS=\
   -O3 \
   -static \
   -std=c89 \
   -Wall \
   -Wextra

MACCFLAGS=\
   -O3

BUGCFLAGS=\
   -O0 \
   -static \
   -std=c89 \
   -Wall \
   -Wextra \
   -g \
   -ggdb

SOURCE=\
   samClust.c \
   mainClust.c

all:
	$(CC) -o $(NAME) $(CFLAGS) $(SOURCE)

mac:
	$(CC) -o $(NAME) $(MACCFLAGS) $(SOURCE)

check:
	$(CC) -o $(NAME) $(BUGCFLAGS) $(SOURCE)

bug:
	$(CC) -o $(NAME) $(BUGCFLAGS) $(SOURCE)
	gdb -x bug-$(NAME).txt $(NAME)

install:
	cp $(NAME) $(PREFIX)
	chmod a+x $(PREFIX)/$(NAME)
