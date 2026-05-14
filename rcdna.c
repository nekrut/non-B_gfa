#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

#if defined(__SSSE3__)
#  include <tmmintrin.h>
#endif

/* Build the reverse complement strand into dna2[].
 *
 * dna2[ndna-i-1] = complement(dna[i]). The SIMD path loads 16 bytes
 * from dna[i..i+15], applies the same nibble-indexed complement LUT
 * used in cdna(), then reverses the 16 bytes in-register before
 * storing them at dna2[ndna-i-16..ndna-i-1]. Reversal uses one more
 * shuffle with a [15..0] mask.
 */
void rcdna(int ndna) {

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
	const __m128i REV = _mm_setr_epi8(
			15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

	for (; i + 16 <= ndna; i += 16) {
		__m128i in = _mm_loadu_si128((const __m128i *) (dna + i));
		__m128i idx = _mm_and_si128(in, NIBBLE);
		__m128i comp = _mm_shuffle_epi8(LUT, idx);
		/* Mask out non-acgtn bytes that share a low nibble with one
		 * of a/c/g/t (e.g. 's'/'w'/'d'/'q'); see cdna.c for the same
		 * fixup. The original scalar switch left dna2 slots untouched
		 * (calloc'd zero) for these. */
		__m128i eqa = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 'a'));
		__m128i eqc = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 'c'));
		__m128i eqg = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 'g'));
		__m128i eqt = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 't'));
		__m128i eqn = _mm_cmpeq_epi8(in, _mm_set1_epi8((char) 'n'));
		__m128i valid = _mm_or_si128(
				_mm_or_si128(_mm_or_si128(eqa, eqc),
						_mm_or_si128(eqg, eqt)),
				eqn);
		comp = _mm_and_si128(comp, valid);
		__m128i rev = _mm_shuffle_epi8(comp, REV);
		_mm_storeu_si128((__m128i *) (dna2 + ndna - i - 16), rev);
	}
#endif
	for (; i < ndna; i++) {
		int k = ndna - i - 1;
		switch (dna[i]) {
			case 'a': dna2[k] = 't'; break;
			case 'c': dna2[k] = 'g'; break;
			case 'g': dna2[k] = 'c'; break;
			case 't': dna2[k] = 'a'; break;
			case 'n': dna2[k] = 'n'; break;
		}
	}
	fprintf(stderr, "Reverse Complement Finished:%d\n", i);
	return;
} /* END */
