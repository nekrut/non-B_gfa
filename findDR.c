#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"
#include "simd_match.h"
#include <time.h>

/******************************
 *  Direct Repeat Finder   ****
 *****************************/
int findDR(int mindir, int maxdir, int dspacer, int total_bases) {

	time_t rawtime;
	time(&rawtime);

	register int j, i, k, sp, end;
	int strti = 0;
	int ndx = 0;
	int size = 0;
	int lasti;
	int sizeMin = 0;
	int spMin;
	int spMax;
	int totlen;
	i = j = k = sp = size = end = 0;

	/*******************************************
	 * Start looking for direct repeats ********
	 *******************************************
	 /*/
	lasti = total_bases - (mindir * 2); //last i value that needs to be examined

	for (strti = 0; strti <= lasti; strti++) {
		while (dna[strti] == 'n') {//skip n's
			strti++;
		}
		if (strti >= lasti) {
			break;
		}

		for (size = maxdir; size >= mindir; size--) {
			if (((size * 2) + dspacer) <= (end - strti)) {// sizes are not big enough to escape prev.
				sp = dspacer;
				size = sizeMin;
				continue;
				//break;
			}
			//only examine spacers that give large enough repeats to escape previous
			spMin = max(0,((end-strti)-(size*2))+2);
			//watch for end of sequence
			spMax = min(dspacer,lasti-strti);
			if (spMax < spMin) continue;

			/* Broadcast-match: compute the bitmask of sp values in
			 * [spMin..spMax] where dna[strti+size+sp] == dna[strti]. One
			 * SIMD load+cmpeq replaces the spMax-spMin+1 byte-compare
			 * iterations the scalar inner loop would do; we then iterate
			 * only the set bits. Requires spMax-spMin < 16. */
			unsigned char b0 = (unsigned char) dna[strti];
			int n_sp = spMax - spMin + 1;
			uint32_t mask;
			if (n_sp <= 16) {
				mask = broadcast_match_mask_16(
						b0,
						(const unsigned char *) &dna[strti + size + spMin],
						n_sp);
			} else {
				mask = ~(uint32_t) 0;  // very rare: spacer range > 16
			}

			while (mask) {
				int sp_offset = __builtin_ctz(mask);
				mask &= mask - 1;
				sp = spMin + sp_offset;
				int max_len = size;
				int rhs_cap = total_bases - strti - size - sp;
				if (rhs_cap < max_len) max_len = rhs_cap;
				if (max_len < 0) max_len = 0;
				k = forward_match_n_on_a(
						(const unsigned char *) &dna[strti],
						(const unsigned char *) &dna[strti + size + sp],
						max_len);
				i = strti + k;
				j = strti + size + sp + k;
				if (k == size) {//DR found!
					totlen = k;
					if (sp == 0) {
						while (dna[i] == dna[j]) {//expand Big Tandem Repeat (BTR)
							totlen++;
							j++;
							i++;
						}
					}
					//Set new DR
					drep[ndx].start = strti + 1;
					drep[ndx].len = size;
					drep[ndx].loop = sp;
					drep[ndx].num = totlen / drep[ndx].len; //repeats
					drep[ndx].end = j;
					drep[ndx].sub = (totlen % drep[ndx].len); //remainder
					drep[ndx].strand = 0;

					ndx++;
					end = j - 1;
					sp = dspacer;
					size = sizeMin;
					/* Exit the broadcast-mask loop: the original scalar code
					 * relied on setting size=sizeMin to break the outer
					 * size-loop on the NEXT for-step. With the mask loop
					 * iterating bits, we have to break out explicitly --
					 * otherwise the next mask bit triggers an extension at
					 * size=0 and the divisor in `totlen / size` is zero. */
					break;
				}
			}
		}
	}
	return (ndx);
}/* END of direct*/
