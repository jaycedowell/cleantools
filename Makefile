# $ID: build_beam, v1.0 2008/07/24 Exp $

CFLAGS= -c -O3 -fPIC -m64 -D_REENTRANT -I/usr/local/itt/idl/external/include/
LDFLAGS= -shared -lm
SRCS=build_beam.c alfa_clean.c lib_version.c
OBJS=build_beam.o alfa_clean.o lib_version.o

all: clean_tools.so

$(SRCS):
	$(CC) $(CFLAGS) -c $*.c

clean_tools.so: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

# Clean up
clean:
	-rm -f clean_tools.so *.o
