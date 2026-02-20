# GNU Automake config


# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
if WITH_UMPKFCVM_LARGEDATA_DIR
	echo "USING WITH_UMPKFCVM_LARGEDATA_DIR in uwpkfcvm"
else
	if [ ! -d "./uwpkfcvm" ]; then \
		./make_data_files.py; \
	fi
endif


install:
if WITH_UMPKFCVM_LARGEDATA_DIR
	if [ ! -d "${prefix}/data/uwpkfcvm" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		$(LN_S) -f $(UMPKFCVM_LARGEDATA_DIR) ${prefix}/data/uwpkfcvm; \
		echo "USING WITH_UMPKFCVM_LARGEDATA_DIR in uwpkfcvm"; \
	fi
	cp config ${prefix}/data
else
	if [ ! -d "${prefix}/data/uwpkfcvm" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		mv uwpkfcvm ${prefix}/data; \
		$(LN_S) -f ${prefix}/data/uwpkfcvm .; \
	fi
	cp config ${prefix}/data
endif

clean:
	echo "clean"

