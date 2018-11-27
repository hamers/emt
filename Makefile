all: emtfunctions.c triplefunctions.c
	$(CC) -fPIC -shared -O3 -o emtlibrary.so -lm emtfunctions.c triplefunctions.c
clean: 
	$(RM) emtlibrary.so
	
	
