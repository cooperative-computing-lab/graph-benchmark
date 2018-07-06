CC=mpic++
CFLAGS=-03 -fopenmp.

LIBS=-lm

DEPS = csrSpMV.h distribution.h
OBJ = csrSpMV.o distribution.o 


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

spmv: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 