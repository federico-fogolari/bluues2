CC     = cc 
CFLAGS = -g -O3 -Wno-unused-result 
LFLAGS = -lm 

install: bluues2 
	chmod 755 bluues2; cp bluues2 /usr/local/bin/

bluues2: bluues2.c 
	$(CC) bluues2.c -o $@ $(CFLAGS) $(LFLAGS)  

