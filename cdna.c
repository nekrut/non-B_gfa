#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

#if defined(__SSSE3__)
#  include <tmmintrin.h>
#endif

/* Build the forward complement strand into dna3[].
 *
 * Mapping (lowercase): a<->t, c<->g, n<->n. The original code used a
 * per-byte switch; here we use a SIMD shuffle-table lookup keyed on the
 * low nibble of each ASCII byte, which uniquely identifies the five
 * letters of interest (a=1, c=3, t=4, g=7, n=14). _mm_shuffle_epi8
 * applies the lookup 16 bytes per instruction. Other letters map to
 * zero, matching the calloc'd dna3 buffer the original code left for
 * non-acgtn bytes.
 */
void cdna(int ndna) {

	int i = 0;

#if defined(__SSSE3__)
	const __m128i LUT = _mm_setr_epi8(
			0,           /*  0 */
			(char) 't',  /*  1 'a'->'t' */
			0,           /*  2 */
			(char) 'g',  /*  3 'c'->'g' */
			(char) 'a',  /*  4 't'->'a' */
			0,           /*  5 */
			0,           /*  6 */
			(char) 'c',  /*  7 'g'->'c' */
			0,           /*  8 */
			0,           /*  9 */
			0,           /* 10 */
			0,           /* 11 */
			0,           /* 12 */
			0,           /* 13 */
			(char) 'n',  /* 14 'n'->'n' */
			0            /* 15 */);
	const __m128i NIBBLE = _mm_set1_epi8(0x0F);

	for (; i + 16 <= ndna; i += 16) {
		__m128i in = _mm_loadu_si128((const __m128i *) (dna + i));
		__m128i idx = _mm_and_si128(in, NIBBLE);
		__m128i out = _mm_shuffle_epi8(LUT, idx);
		/* The low-nibble LUT collides on non-acgtn letters that share
		 * a nibble with one of a/c/g/t (e.g. 's' has nibble 3 like
		 * 'c', so would be miscomplemented to 'g'). Build an acgtn
		 * equality mask and AND it in so any other lowercase byte
		 * stays zero -- matching the original switch's behaviour of
		 * leaving the dna3 slot untouched (calloc'd zero) for non-
		 * acgtn input. */
		__m128i eqa = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 'a'));
		__m128i eqc = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 'c'));
		__m128i eqg = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 'g'));
		__m128i eqt = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 't'));
		__m128i eqn = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 'n'));
		__m128i valid = _mm_or_si128(
				_mm_or_si128(_mm_or_si128(eqa, eqc),
						_mm_or_si128(eqg, eqt)),
				eqn);
		out = _mm_and_si128(out, valid);
		_mm_storeu_si128((__m128i *) (dna3 + i), out);
	}
#endif
	for (; i < ndna; i++) {
		switch (dna[i]) {
			case 'a': dna3[i] = 't'; break;
			case 'c': dna3[i] = 'g'; break;
			case 'g': dna3[i] = 'c'; break;
			case 't': dna3[i] = 'a'; break;
			case 'n': dna3[i] = 'n'; break;
		}
	}
	fprintf(stderr, "Complement Finished:%d\n", i);
	return;
} /* END */
