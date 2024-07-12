CC=cc
PREFIX=/usr/local/bin

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
   ../generalLib/samEntry.c \
   diCoords.c \
   getDICoords.c

NAME=getDICoords

all:
	$(CC) $(CFLAGS) $(SOURCE) -o $(NAME)

mac:
	$(CC) $(MACCFLAGS) $(SOURCE) -o $(NAME)

check:
	$(CC) $(BUGCFLAGS) $(SOURCE) -o $(NAME)

bug:
	$(CC) $(BUGCFLAGS) $(SOURCE) -o $(NAME)
	gdb -x bug-diCoords.txt getDICoords

install:
	mv $(NAME) $(PREFIX)
	chmod a+x $(PREFIX)/$(NAME)

updatePath:
	printf "export PATH=%s:%s" "$$PATH" "$(PREFIX)" >> "$$HOME/.bashrc";

ppath:
	(printf "%s" "$$PATH" | tr ':' '\n' | grep "$(PREFIX)") || printf "install location (%s) not in path\n" "$(PREFIX)"
