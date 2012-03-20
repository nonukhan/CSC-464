# Makefile is courtesy of Marc Baudoin -- babafou@babafou.eu.org
# On many systems I prefer the gcc compiler over the native compiler
#
     CC = cc # or gcc
 CFLAGS = -O
LDFLAGS = -lm
     RM = rm -f

all : fastDNAml

fastDNAml : fastDNAml.o
	$(CC) $(CFLAGS) -o fastDNAml fastDNAml.o $(LDFLAGS)

fastDNAml.o : fastDNAml.c fastDNAml.h
	$(CC) $(CFLAGS) -c fastDNAml.c

clean :
	$(RM) fastDNAml.o
