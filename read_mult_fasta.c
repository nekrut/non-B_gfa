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

/* Lookup table: maps ASCII bytes to (lowercase letter) for A-Za-z, 0 otherwise.
 * Lets the body-read pass replace per-byte isalpha+tolower with a single
 * load from this table. Built lazily on first use. */
static unsigned char base_lut[256];
static int base_lut_ready = 0;

static void build_base_lut(void) {
	int c;
	for (c = 0; c < 256; c++) base_lut[c] = 0;
	for (c = 'A'; c <= 'Z'; c++) base_lut[c] = (unsigned char)(c | 0x20);
	for (c = 'a'; c <= 'z'; c++) base_lut[c] = (unsigned char) c;
	base_lut_ready = 1;
}

/* Buffered FASTA reader. Replaces a getc-per-byte tokenizer that called
 * isalpha+tolower on every byte. Reads 64 KB at a time, uses memchr to
 * locate '>' and '\n', and a 256-byte lookup table to filter+downcase
 * bases in bulk. Semantics preserved:
 *   - fasta_title gets the first MAX_FASTA_SIZE chars of the header line
 *     (between '>' and the trailing newline), null-padded.
 *   - dna[] is filled with lowercase A-Za-z bytes for the chosen record;
 *     other characters (digits, whitespace, punctuation) are skipped.
 *   - returns base count for the record.
 */
int read_mult_fasta(FILE *dna_file, int fasta, char fasta_title[]) {
	enum { BUFSZ = 1 << 16 };
	static unsigned char buf[BUFSZ];
	size_t n = 0, pos = 0;
	int i = 0;                /* base index into dna[] */
	int records_seen = 0;     /* number of '>' chars consumed so far */
	int target = fasta;       /* 1-based record we want */
	int in_header = 0;        /* currently reading header bytes */
	int header_written = 0;   /* bytes already copied into fasta_title */
	int header_cap = MAX_FASTA_SIZE;

	if (!base_lut_ready) build_base_lut();
	memset((char *) fasta_title, '\0', header_cap);

	for (;;) {
		if (pos >= n) {
			n = fread(buf, 1, BUFSZ, dna_file);
			pos = 0;
			if (n == 0) break;        /* EOF */
		}

		if (in_header) {
			/* Copy header bytes up to the next '\n'. */
			unsigned char *nl = memchr(buf + pos, '\n', n - pos);
			size_t avail = nl ? (size_t)(nl - (buf + pos)) : (n - pos);
			if (records_seen == target && header_written < header_cap) {
				size_t take = avail;
				if (header_written + (int) take > header_cap)
					take = header_cap - header_written;
				memcpy(fasta_title + header_written, buf + pos, take);
				header_written += (int) take;
			}
			pos += avail;
			if (nl) {
				pos += 1;             /* consume the newline */
				in_header = 0;
			}
			continue;
		}

		/* Sequence body: scan for the next '>' which marks end of record. */
		unsigned char *gt = memchr(buf + pos, '>', n - pos);
		size_t body_end = gt ? (size_t)(gt - buf) : n;

		if (records_seen == target) {
			/* This is the record we want; filter+downcase its bytes into dna[]. */
			size_t k;
			for (k = pos; k < body_end; k++) {
				unsigned char b = base_lut[buf[k]];
				if (b) dna[i++] = b;
			}
		}

		pos = body_end;
		if (gt) {
			pos += 1;                 /* consume '>' */
			records_seen++;
			if (records_seen > target) {
				/* We've finished the target record. */
				return i;
			}
			in_header = 1;
		}
	}

	/* EOF reached while reading the target record (or before reaching it). */
	if (records_seen < target) {
		fprintf(stderr, " WARNING: requested FASTA record %d but file ended after %d\n",
				target, records_seen);
	}
	return i;
}
