# Where the executable files will be copied
destination := $(HOME)/bin

#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
opt :=  -O0 -fno-inline-functions      # For debugging

# Flags to determine the warning messages issued by the compiler
warn := \
 -Wall \
 -Wcast-align \
 -Wcast-qual \
 -Wmissing-declarations \
 -Wmissing-prototypes \
 -Wnested-externs \
 -Wpointer-arith \
 -Wstrict-prototypes \
 -Wno-unused-parameter \
 -Wno-unused-function \
 -Wshadow \
 -Wundef \
 -Wwrite-strings

CFLAGS := -g -std=gnu99 $(warn) $(opt)
lib := -L/usr/local/lib -lgsl -lgslcblas -lpthread -lm

.c.o:
	$(CC) $(CFLAGS) $(incl) -c -o ${@F}  $<

# test jobqueue.c
XJOBQUEUE := xjobqueue.o misc.o 
xjobqueue : $(XJOBQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XJOBQUEUE) $(lib)

xjobqueue.o : jobqueue.c jobqueue.h
	$(CC) $(CFLAGS) -c -DTEST -o $@ jobqueue.c

# Make dependencies file
depend : *.c *.h
	echo '#Automatically generated dependency info' > depend
	$(CC) -MM $(incl) *.c >> depend

clean :
	rm -f *.a *.o *~ 

include depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: clean

