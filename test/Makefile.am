# Autoconf/automake file

bin_PROGRAMS = test_uwpkfcvm

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${PROJ_INCLUDES} -I${UCVM_SRC_PATH}/src/ucvm -I../src
AM_LDFLAGS = ${LDFLAGS} ${PROJ_LDFLAGS}  -L../src -luwpkfcvm -lm

objects = test_api.o
TARGETS = $(bin_PROGRAMS)

.PHONY = run_unit

all: $(bin_PROGRAMS)

############################################
# Executables
############################################

install: all
	mkdir -p ${prefix}/test
	cp test_uwpkfcvm ${prefix}/test

test_uwpkfcvm: $(objects)
	$(CC) -o $@ $^ $(AM_LDFLAGS) 

run_unit : test_uwpkfcvm
	./test_uwpkfcvm

clean: 
	rm -rf test_uwpkfcvm test_api.o
