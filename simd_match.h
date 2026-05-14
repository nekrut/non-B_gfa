#ifndef SIMD_MATCH_H_
#define SIMD_MATCH_H_

/*
 * SIMD-accelerated extension primitives used by findIR/findMR/findDR.
 *
 * Both helpers return the length k of the longest prefix that satisfies the
 * loop's match-and-bounds conditions. The original byte-at-a-time `while`
 * loops compute the same k; these helpers process 32 bytes per AVX2
 * iteration, 16 bytes per SSSE3 iteration, with a scalar tail.
 *
 * Defined as `static inline` so the find-loops can inline the helper bodies
 * and avoid an indirect-call overhead per (strti, sp) pair. AVX2 / SSSE3 /
 * SSE2 paths are gated on compile-time __AVX2__ / __SSSE3__ / __SSE2__
 * macros (set by gcc when -mavx2 / -mssse3 / -msse2 is on the command line).
 *
 * Default Makefile builds with `-msse2 -mssse3` (universal on x86_64 since
 * 2006). Pass CFLAGS="-O2 -mavx2 -mssse3 -msse2" or `-march=native` to
 * enable the AVX2 path.
 *
 * Bounds policy: callers must pass `max_len` no larger than what the original
 * loop's index-bounds checks would have permitted. The DNA buffers in gfa.c
 * are calloc'd with 64 bytes of trailing zero slack so 32-byte loads near the
 * high end of the buffer remain in the allocation; trailing zero bytes do
 * not equal any base in {a,c,g,t,n} and so they don't change k.
 */

#include <stdint.h>

#if defined(__AVX2__)
#  include <immintrin.h>
#endif
#if defined(__SSSE3__)
#  include <tmmintrin.h>
#endif
#if defined(__SSE2__)
#  include <emmintrin.h>
#endif

/* Forward-forward match used by findDR: a and b both walk forward.
 * Stops at the first index k where:
 *   - a[k] != b[k], OR
 *   - a[k] == 'n', OR
 *   - k == max_len. */
static inline int forward_match_n_on_a(const unsigned char *a,
                                       const unsigned char *b,
                                       int max_len) {
	int k = 0;
#if defined(__AVX2__)
	{
		const __m256i nv = _mm256_set1_epi8('n');
		const __m256i ones = _mm256_set1_epi8(-1);
		while (k + 32 <= max_len) {
			__m256i va = _mm256_loadu_si256((const __m256i *)(a + k));
			__m256i vb = _mm256_loadu_si256((const __m256i *)(b + k));
			__m256i eq = _mm256_cmpeq_epi8(va, vb);
			__m256i is_n = _mm256_cmpeq_epi8(va, nv);
			__m256i stop = _mm256_or_si256(
					_mm256_xor_si256(eq, ones), is_n);
			uint32_t mask = (uint32_t) _mm256_movemask_epi8(stop);
			if (mask) return k + __builtin_ctz(mask);
			k += 32;
		}
	}
#endif
#if defined(__SSE2__)
	{
		const __m128i nv = _mm_set1_epi8('n');
		const __m128i ones = _mm_set1_epi8(-1);
		while (k + 16 <= max_len) {
			__m128i va = _mm_loadu_si128((const __m128i *)(a + k));
			__m128i vb = _mm_loadu_si128((const __m128i *)(b + k));
			__m128i eq = _mm_cmpeq_epi8(va, vb);
			__m128i is_n = _mm_cmpeq_epi8(va, nv);
			__m128i stop = _mm_or_si128(_mm_xor_si128(eq, ones), is_n);
			uint32_t mask = (uint32_t) _mm_movemask_epi8(stop);
			if (mask) return k + __builtin_ctz(mask);
			k += 16;
		}
	}
#endif
	while (k < max_len && a[k] == b[k] && a[k] != (unsigned char) 'n') k++;
	return k;
}

/* Reverse-forward match used by findIR/findMR: left walks DOWN from left_end,
 * right walks UP from `right`. The match at step k compares left_end[-k]
 * against right[k]. Stops at the first k where:
 *   - left_end[-k] != right[k], OR
 *   - right[k] == 'n', OR
 *   - k == max_len. */
static inline int reverse_forward_match_n_on_right(const unsigned char *left_end,
                                                   const unsigned char *right,
                                                   int max_len) {
	int k = 0;
#if defined(__AVX2__)
	{
		/* Per-lane byte reverse + lane swap = full 32-byte reverse. */
		const __m256i nv = _mm256_set1_epi8('n');
		const __m256i ones = _mm256_set1_epi8(-1);
		const __m256i lane_rev = _mm256_setr_epi8(
				15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
				15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
		while (k + 32 <= max_len) {
			__m256i vl = _mm256_loadu_si256(
					(const __m256i *)(left_end - 31 - k));
			__m256i rev_per_lane = _mm256_shuffle_epi8(vl, lane_rev);
			__m256i rev = _mm256_permute2x128_si256(rev_per_lane, rev_per_lane, 0x01);
			__m256i vr = _mm256_loadu_si256((const __m256i *)(right + k));
			__m256i eq = _mm256_cmpeq_epi8(rev, vr);
			__m256i is_n = _mm256_cmpeq_epi8(vr, nv);
			__m256i stop = _mm256_or_si256(
					_mm256_xor_si256(eq, ones), is_n);
			uint32_t mask = (uint32_t) _mm256_movemask_epi8(stop);
			if (mask) return k + __builtin_ctz(mask);
			k += 32;
		}
	}
#endif
#if defined(__SSSE3__)
	{
		const __m128i nv = _mm_set1_epi8('n');
		const __m128i ones = _mm_set1_epi8(-1);
		const __m128i rev_mask = _mm_setr_epi8(
				15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
		while (k + 16 <= max_len) {
			__m128i vl = _mm_loadu_si128(
					(const __m128i *)(left_end - 15 - k));
			__m128i rev = _mm_shuffle_epi8(vl, rev_mask);
			__m128i vr = _mm_loadu_si128((const __m128i *)(right + k));
			__m128i eq = _mm_cmpeq_epi8(rev, vr);
			__m128i is_n = _mm_cmpeq_epi8(vr, nv);
			__m128i stop = _mm_or_si128(_mm_xor_si128(eq, ones), is_n);
			uint32_t mask = (uint32_t) _mm_movemask_epi8(stop);
			if (mask) return k + __builtin_ctz(mask);
			k += 16;
		}
	}
#endif
	while (k < max_len && left_end[-k] == right[k]
			&& right[k] != (unsigned char) 'n')
		k++;
	return k;
}

/* Forward-forward match without the 'n' stop. Used by findSTR where the
 * original strncmp-based loop matches 'n' == 'n' as equal (rather than
 * stopping). Returns the length of the longest prefix where a[k] == b[k],
 * stopping at the first mismatch or k == max_len.
 */
static inline int forward_match(const unsigned char *a,
                                const unsigned char *b,
                                int max_len) {
	int k = 0;
#if defined(__AVX2__)
	{
		while (k + 32 <= max_len) {
			__m256i va = _mm256_loadu_si256((const __m256i *)(a + k));
			__m256i vb = _mm256_loadu_si256((const __m256i *)(b + k));
			__m256i eq = _mm256_cmpeq_epi8(va, vb);
			uint32_t mask = ~(uint32_t) _mm256_movemask_epi8(eq);
			if (mask) return k + __builtin_ctz(mask);
			k += 32;
		}
	}
#endif
#if defined(__SSE2__)
	{
		while (k + 16 <= max_len) {
			__m128i va = _mm_loadu_si128((const __m128i *)(a + k));
			__m128i vb = _mm_loadu_si128((const __m128i *)(b + k));
			__m128i eq = _mm_cmpeq_epi8(va, vb);
			uint32_t mask = (~(uint32_t) _mm_movemask_epi8(eq)) & 0xFFFFu;
			if (mask) return k + __builtin_ctz(mask);
			k += 16;
		}
	}
#endif
	while (k < max_len && a[k] == b[k]) k++;
	return k;
}

/* Broadcast-match: compare a single byte `a` against up to 16 bytes at `b`,
 * returning a 16-bit bitmask where bit k is set iff b[k] == a AND b[k] != 'n'
 * AND k < n. n must be in [0, 16]. b must be a valid 16-byte load (the DNA
 * buffer's 64-byte trailing-zero slack makes that safe near the high end of
 * the allocation).
 *
 * Used by findDR to flatten the sp byte-compare loop: spMax+1 is typically
 * ~11, so one SIMD load+cmpeq+movemask replaces an 11-iteration scalar loop
 * of byte loads and conditional branches, and the caller iterates only the
 * sp positions that pass the first-byte check.
 */
static inline uint32_t broadcast_match_mask_16(unsigned char a,
                                               const unsigned char *b,
                                               int n) {
	if (n <= 0) return 0;
	if (a == (unsigned char) 'n') return 0;
#if defined(__SSE2__)
	{
		const __m128i av = _mm_set1_epi8((char) a);
		__m128i vb = _mm_loadu_si128((const __m128i *) b);
		__m128i eq = _mm_cmpeq_epi8(vb, av);
		uint32_t mask = (uint32_t) _mm_movemask_epi8(eq);
		if (n < 16) mask &= ((uint32_t) 1 << n) - 1;
		return mask;
	}
#else
	{
		uint32_t mask = 0;
		int k;
		int lim = n < 16 ? n : 16;
		for (k = 0; k < lim; k++) {
			if (b[k] == a) mask |= (uint32_t) 1 << k;
		}
		return mask;
	}
#endif
}

#endif /* SIMD_MATCH_H_ */
