CC=cc
PREFIX=/usr/local/bin

CFLAGS=\
   -O3

STATICCFLAGS=\
   -O3 \
   -static \
   -std=c89 \
   -Wall \
   -Wextra

core:
	make CC=$(CC) MACCFLAGS=$(CFLAGS) mac -C getDIIdsSrc
	make CC=$(CC) MACCFLAGS=$(CFLAGS) mac -C getDICoordsSrc
	make CC=$(CC) MACCFLAGS=$(CFLAGS) mac -C findDIFragSrc

side:
	make CC=$(CC) MACCFLAGS=$(CFLAGS) mac -C memwater
	make CC=$(CC) MACCFLAGS=$(CFLAGS) mac -C water

staticCore:
	make CC=$(CC) CFLAGS=$(STATICCFLAGS) -C getDIIdsSrc
	make CC=$(CC) CFLAGS=$(STATICCFLAGS) -C getDICoordsSrc
	make CC=$(CC) CFLAGS=$(STATICCFLAGS) -C findDIFragSrc

staticSide:
	make CC=$(CC) CFLAGS=$(STATICCFLAGS) mac -C memwater
	make CC=$(CC) CFLAGS=$(STATICCFLAGS) mac -C water

all: core side

static: staticCore staticSide

# update path if needed
path:
	(printf "%s" "$$PATH" | tr ':' '\n' | grep "$(PREFIX)") || (printf "export PATH=%s:%s" "$$PATH" "$(PREFIX)" >> "$$HOME/.bashrc")

# this is so the path is updated last
install: coreInstall path

coreInstall:
	mv getDIIdsSrc/getDIIds $(PREFIX)
	mv getDICoordsSrc/getDICoords $(PREFIX)
	mv findDIFragSrc/findDIFrag $(PREFIX)

	chmod a+x $(PREFIX)/getDIIdsSrc/getDIIds
	chmod a+x $(PREFIX)/getDICoordsSrc/getDICoords
	chmod a+x $(PREFIX)/findDIFragSrc/findDIFrag

# this is so the path is updated last
allInstall: allInstallSub path

allInstallSub:
	mv getDIIdsSrc/getDIIds $(PREFIX)
	mv getDICoordsSrc/getDICoords $(PREFIX)
	mv findDIFragSrc/findDIFrag $(PREFIX)
	mv memwater/alnMemwater $(PREFIX)
	mv water/alnwater $(PREFIX)

	chmod a+x $(PREFIX)/getDIIds
	chmod a+x $(PREFIX)/getDICoords
	chmod a+x $(PREFIX)/findDIFrag
	chmod a+x $(PREFIX)/alnMemwater
	chmod a+x $(PREFIX)/alnwater
