#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "gfa.h"
/**********************************************************
 *
 *      Procedure to read the dna file and put into array
 *
 *
 **********************************************************/
/*
 * Defined in main module file to avoid stack overflow problem
 * can easily be 300Mb
 */
extern char *dna;

//Scan the file once: count records and compute the longest record's sequence
//length (bases only, comments/whitespace skipped). The returned max is used to
//size the heap dna[] buffer to the actual need rather than a 300 Mbp hard cap.
int get_fasta_count_ex(FILE *dna_file, int *max_seq_len_out);

int get_fasta_count(FILE *dna_file) {
	int max_seq_len = 0;
	return get_fasta_count_ex(dna_file, &max_seq_len);
}

int get_fasta_count_ex(FILE *dna_file, int *max_seq_len_out) {
	enum { BUFSZ = 1 << 16 };
	static unsigned char buf[BUFSZ];
	size_t n;
	int fasta_count = 0;
	int cur_seq_len = 0;
	int max_seq_len = 0;
	int in_header = 0;
	int first_byte_seen = 0;

	while ((n = fread(buf, 1, BUFSZ, dna_file)) > 0) {
		size_t i;
		if (!first_byte_seen) {
			if (buf[0] != '>') {
				fprintf(stderr, " base read =%c.\n", buf[0]);
				fprintf(stderr,
						" FATAL ERROR: Sequence file NOT in FASTA format\n");
				exit(88);
			}
			first_byte_seen = 1;
		}
		for (i = 0; i < n; i++) {
			unsigned char c = buf[i];
			if (c == '>') {
				if (fasta_count > 0 && cur_seq_len > max_seq_len) {
					max_seq_len = cur_seq_len;
				}
				cur_seq_len = 0;
				fasta_count++;
				in_header = 1;
			} else if (in_header) {
				if (c == '\n') in_header = 0;
			} else if (isalpha(c)) {
				cur_seq_len++;
			}
		}
	}
	if (cur_seq_len > max_seq_len) max_seq_len = cur_seq_len;

	if (!first_byte_seen) {
		fprintf(stderr, " End of Sequence file!\n");
		return (0);
	}
	*max_seq_len_out = max_seq_len;
	return (fasta_count);
}

int read_mult_fasta(FILE *dna_file, int fasta, char fasta_title[]) {
	register int i;
	int base, fasta_len;
	char line[MAX_LINE + 1];
	BOOLEAN start;
	int fasta_count = 0;

	/* Declare Procedure */
	//void nulls(char line[], int n);

	fasta_len = MAX_FASTA_SIZE;
	start = TRUE;
	i = 0;
	//fprintf(stderr, "\n fasta number =%d \n", fasta);
	if (fasta > 1) {
		//advance to correct fasta section
		while (fasta_count <= (fasta-1)) {
			base = getc(dna_file);
			if (base == '>') {
				fasta_count++;
			}
		}
		ungetc(base, dna_file);
	}

	while ((base = getc(dna_file)) != EOF) {
		switch (base) {
			case '>': /* if (start) get header information
			 else break and return */
			if (start) {
				if (fgets(line, MAX_LINE, dna_file) != NULL) {
					/* crack some part of header and make title */
					memset((char *) fasta_title, '\0', fasta_len);
					//nulls(fasta_title, fasta_len);
					if (strlen(line) < fasta_len) fasta_len = strlen(line) - 1;
					strncpy(fasta_title, &line[0], fasta_len);
				}
				start = FALSE;
			}
			else {
				ungetc(base, dna_file);
				fprintf(stderr, " INFO: Program read %d bases \n", i);
				return (i);
			}
				break;
			default: /* look for DNA sequence */
			if (isalpha(base)) {
				dna[i] = tolower(base);
				i++;
			}
				break;
		}
	}
	//fprintf(stderr, " INFO: Program read %d bases \n", i);
	return (i);
}
