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
			for (sp = spMin; sp <= spMax; sp++) {
				/* Inline byte-0 fast-path: most (size, sp) pairs on real DNA
				 * fail at the first byte. The static-inline SIMD helper is
				 * cheap, but its m128-constant setup still costs cycles. A
				 * scalar byte compare here lets us skip the SIMD entry for
				 * the ~3/4 of positions where dna[strti] doesn't match. */
				unsigned char b0 = (unsigned char) dna[strti];
				if (b0 == (unsigned char) 'n'
						|| (unsigned char) dna[strti + size + sp] != b0) {
					k = 0;
					i = strti;
					j = strti + size + sp;
				} else {
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
				}
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
				}
			}
		}
	}
	return (ndx);
}/* END of direct*/
