###########################################################
# The makefile for the MPEG program.
# There should be no changes for most UNIX compilers.
###########################################################

MFLAGS = -O3 -pg -g
DEFS = system.h globals.h marker.h
BASELINE =  mpeg.o codec.o huffman.o io.o chendct.o lexer.o marker.o me.o mem.o stat.o stream.o transform.o y4m_input.o vidinput.o

.c.o:
	$(CC) $(MFLAGS) -I. -c $*.c 

.c.ln:
	lint -c $*.c 

all: mpeg

clean:
	rm *.o mpeg

mpeg: $(BASELINE)
	$(CC) $(MFLAGS) $(BASELINE) -lm -o mpeg

mpeg.o: mpeg.c $(DEFS) vidinput.h ogg/os_types.h
codec.o: codec.c $(DEFS)
marker.o: marker.c $(DEFS) marker.h
huffman.o: huffman.c $(DEFS) huffman.h
io.o: io.c $(DEFS) vidinput.h ogg/os_types.h
chendct.o: chendct.c $(DEFS)
lexer.o:lexer.c
mem.o: mem.c 
me.o: me.c 
stat.o: stat.c 
stream.o: stream.c $(DEFS)
transform.o: transform.c $(DEFS) dct.h
vidinput.o: vidinput.c vidinput.h ogg/os_types.h
y4m_input.o: y4m_input.c vidinput.h ogg/os_types.h	

lcheck: mpeg.ln codec.ln marker.ln huffman.ln io.ln chendct.ln lexer.ln mem.ln me.ln stat.ln stream.ln transform.ln
	lint  mpeg.ln codec.ln marker.ln huffman.ln io.ln chendct.ln lexer.ln mem.ln me.ln stat.ln stream.ln transform.ln


#
# Should you wish to modify the interpreter 
# or the Huffman tables, use these commands.
#
#
#ctables.h: tables.h htable
#	rm -f ctables.h
#	htable <tables.h >ctables.h
#
# Make sure that you move the first comment from the start to 
# within the braces %{ /* */ %} in htable.lex, lexer.l,
# otherwise it will give you a error (definitions too long). 
# Caution: the htable name may conflict with /usr/etc/htable
# on some machines.
#
#htable:htable.lex
#	lex htable.lex
#	mv lex.yy.c htable.c
#	cc htable.c -o htable -ll
#
#
#  Caution: Sometimes -ll is required.
#
#lexer.c:lexer.l
#	lex lexer.l
#	mv lex.yy.c lexer.c
#
