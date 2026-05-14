PGMNAME = gfa
#CC = icc
CC = gcc -lm
OBJ = $(PGMNAME).o \
cdna.o     findIR.o \
findAPR.o  findMR.o    print_gff_file.o  rcdna.o \
findDR.o   findSTR.o   is_subset.o  print_tsv_file.o \
findGQ.o   findZDNA.o  print_usage.c     read_mult_fasta.o
# SIMD primitives in simd_match.h are static-inline so they fold into the
# IR/MR/DR inner loops. Default build enables SSSE3 (universal on x86_64
# since 2006) which gives 16-byte vector compare + reverse. Override with
# CFLAGS="-O2 -march=native" (or `-mavx2 -mssse3`) for the 32-byte AVX2
# path on modern CPUs.
CFLAGS = -O2 -mssse3
#CFLAGS = -O2 -ansi
LFLAGS = 
$(PGMNAME): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LFLAGS)

clean:
	rm *.o
