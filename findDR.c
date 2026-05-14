#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "gfa.h"
#include "simd_match.h"
#include <time.h>

/******************************
 *  Direct Repeat Finder   ****
 *****************************
 *
 * Per-strti structure (post-restructuring):
 *
 *   1. Precompute LCP[D] = longest common prefix of dna[strti..] and
 *      dna[strti+D..], for every distance D in [mindir, maxdir+dspacer].
 *      The first-byte filter is run in 16-D-wide SIMD chunks so the
 *      ~75% of D values whose first byte mismatches dna[strti] cost
 *      only a single SIMD load+cmpeq+movemask, no extension call.
 *
 *   2. Early exit: if max(LCP[D]) < mindir, no DR can possibly exist
 *      at this strti (the size of any DR equals the LCP at the chosen
 *      D, and size must be at least mindir). Skip the rest. On uniform
 *      random DNA this fires for ~99.97% of strti.
 *
 *   3. Iterate (size, sp) in the original priority order (size DESC,
 *      then sp ASC) and use LCP[size+sp] as an O(1) replacement for
 *      the per-pair extension call. First hit wins; for sp==0 the BTR
 *      tail is computed with a single SIMD extension from the matched
 *      position.
 *
 * The original (size, sp, extension) triple loop is preserved exactly
 * for ordering and tie-break semantics; extension results are just
 * served from a precomputed table instead of being recomputed.
 */
int findDR(int mindir, int maxdir, int dspacer, int total_bases) {

	time_t rawtime;
	time(&rawtime);

	int strti = 0;
	int ndx = 0;
	int size = 0;
	int lasti;
	int sizeMin = 0;
	int spMin;
	int spMax;
	int totlen;
	int end = 0;

	/*******************************************
	 * Start looking for direct repeats ********
	 *******************************************/
	lasti = total_bases - (mindir * 2); //last i value that needs to be examined

	const int MIN_D = mindir;
	const int MAX_D_GLOBAL = maxdir + dspacer;
	const int LCP_LEN = MAX_D_GLOBAL - MIN_D + 1;
	int *LCP = (int *) calloc((size_t) LCP_LEN, sizeof(int));
	if (!LCP) {
		fprintf(stderr, "FATAL: findDR could not allocate LCP table\n");
		exit(23);
	}

	for (strti = 0; strti <= lasti; strti++) {
		while (dna[strti] == 'n') {
			strti++;
		}
		if (strti >= lasti) {
			break;
		}

		unsigned char b0 = (unsigned char) dna[strti];
		if (b0 == (unsigned char) 'n') continue;

		/* Cap MAX_D for this strti to stay inside the sequence. */
		int MAX_D = MAX_D_GLOBAL;
		if (MAX_D > total_bases - strti - 1)
			MAX_D = total_bases - strti - 1;
		if (MAX_D < MIN_D) continue;

		int n_D = MAX_D - MIN_D + 1;

		/* Phase 1: clear and populate LCP[D] for D in [MIN_D..MAX_D]. */
		memset(LCP, 0, (size_t) n_D * sizeof(int));
		int max_lcp = 0;

		for (int D_base = MIN_D; D_base <= MAX_D; D_base += 16) {
			int chunk_n = MAX_D - D_base + 1;
			if (chunk_n > 16) chunk_n = 16;
			uint32_t mask = broadcast_match_mask_16(
					b0,
					(const unsigned char *) &dna[strti + D_base],
					chunk_n);
			while (mask) {
				int off = __builtin_ctz(mask);
				mask &= mask - 1;
				int D = D_base + off;
				int max_lcp_for_D = total_bases - strti - D;
				if (max_lcp_for_D > maxdir) max_lcp_for_D = maxdir;
				if (max_lcp_for_D < mindir) continue;
				int lcp = forward_match_n_on_a(
						(const unsigned char *) &dna[strti],
						(const unsigned char *) &dna[strti + D],
						max_lcp_for_D);
				LCP[D - MIN_D] = lcp;
				if (lcp > max_lcp) max_lcp = lcp;
			}
		}

		/* Phase 2: early exit if no D has a long-enough LCP. */
		if (max_lcp < mindir) continue;

		/* Phase 3: iterate (size, sp) in original priority order, using
		 * LCP[D] as an O(1) extension lookup. First hit wins per strti. */
		int found = 0;
		for (size = maxdir; size >= mindir; size--) {
			if (((size * 2) + dspacer) <= (end - strti)) {
				/* Original used `size = sizeMin; continue;` here, which
				 * exits the size loop on the next for-step. break is
				 * equivalent and clearer. */
				break;
			}
			spMin = max(0, ((end - strti) - (size * 2)) + 2);
			spMax = min(dspacer, lasti - strti);
			if (spMax < spMin) continue;

			for (int sp = spMin; sp <= spMax; sp++) {
				int D = size + sp;
				int idx = D - MIN_D;
				if (idx < 0 || idx >= n_D) continue;
				if (LCP[idx] < size) continue;

				/* DR found at (size, sp). For sp==0 (Big Tandem Repeat)
				 * extend past `size` with a SIMD forward match from the
				 * already-matched stretch. */
				totlen = size;
				int j = strti + size + sp + size;
				if (sp == 0) {
					int btr_max = total_bases - j;
					int i_btr = strti + size;
					if (btr_max > 0) {
						int btr_k = forward_match_n_on_a(
								(const unsigned char *) &dna[i_btr],
								(const unsigned char *) &dna[j],
								btr_max);
						totlen += btr_k;
						j += btr_k;
					}
				}

				drep[ndx].start = strti + 1;
				drep[ndx].len = size;
				drep[ndx].loop = sp;
				drep[ndx].num = totlen / drep[ndx].len;
				drep[ndx].end = j;
				drep[ndx].sub = totlen % drep[ndx].len;
				drep[ndx].strand = 0;
				ndx++;
				end = j - 1;
				found = 1;
				break;
			}
			if (found) break;
		}
		(void) sizeMin;
	}

	free(LCP);
	return (ndx);
}/* END of direct*/
