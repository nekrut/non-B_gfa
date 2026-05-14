#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "gfa.h"
#include "simd_match.h"

/*******************************************************************
 *  findIR eXplorer (crux)                                      *
 *  Program to locate possible cruciforms in nucleic acid sequence *
 *******************************************************************/

void delIRep(int nreps, int toRemove) {//shifts stack down, but does not reset ndx
	int i = 0;
	for (i = toRemove; i < nreps; i++) {
		irep[i] = irep[i + 1];
	}
	//fprintf(stderr, " removed start %d stop %d  \n", mrep[toRemove].start,mrep[toRemove].end);

}

int findIR(int mincrf, int cspacer, int cut, int shortSpacer, int total_bases) {

	register int i, j, k, sp;
	int strti = 0;
	int cBack = 0;//counter to look at previous elements for overlap etc.
	int maxcBack = 10;//
	int ndx = 0;
	int tmpStart = 0;
	int tmpStop = 0;
	int maxSP = 0;
	i = j = k = sp = 0;
	BOOLEAN rightShifted = FALSE;
	BOOLEAN leftShifted = FALSE;

	/*******************************************
	 * Start looking for inverted repeats*******
	 *******************************************
	 *
	 * Per-strti structure: precompute LCP[sp] = longest reverse-forward
	 * match starting at (dna[strti], dna3[strti+sp+1]) for each sp in
	 * [0, cspacer], using broadcast_match_mask_16 to filter the ~75% of
	 * sp values whose first byte mismatches dna[strti]. If max(LCP) <
	 * mincrf, no IR is possible at this strti and the sp loop is
	 * skipped. Otherwise the original sp loop runs unchanged but reads
	 * k from the table instead of recomputing it. */
	int LCP_CAP = cspacer + 1;
	int *LCP = (int *) calloc((size_t) LCP_CAP, sizeof(int));
	if (!LCP) {
		fprintf(stderr, "FATAL: findIR could not allocate LCP table\n");
		exit(23);
	}

	for (strti = mincrf; strti <= (total_bases - mincrf); strti++) {
		maxSP = min(cspacer,(total_bases-(strti+mincrf)));
		if (maxSP < 0) continue;

		unsigned char b0 = (unsigned char) dna[strti];
		if (b0 == (unsigned char) 'n') continue;

		memset(LCP, 0, (size_t) (maxSP + 1) * sizeof(int));
		int max_lcp = 0;

		for (int chunk = 0; chunk <= maxSP; chunk += 16) {
			int chunk_n = (maxSP + 1) - chunk;
			if (chunk_n > 16) chunk_n = 16;
			uint32_t mask = broadcast_match_mask_16(
					b0,
					(const unsigned char *) &dna3[strti + 1 + chunk],
					chunk_n);
			while (mask) {
				int off = __builtin_ctz(mask);
				mask &= mask - 1;
				int sp_i = chunk + off;
				int max_len = strti + 1;
				int rhs_cap = total_bases - strti - sp_i - 1;
				if (rhs_cap < max_len) max_len = rhs_cap;
				if (max_len < mincrf) continue;
				int lcp = reverse_forward_match_n_on_right(
						(const unsigned char *) &dna[strti],
						(const unsigned char *) &dna3[strti + sp_i + 1],
						max_len);
				LCP[sp_i] = lcp;
				if (lcp > max_lcp) max_lcp = lcp;
			}
		}

		if (max_lcp < mincrf) continue;

		for (sp = 0; sp <= maxSP; sp++) {
			k = LCP[sp];
			i = strti - k;
			j = strti + sp + 1 + k;
			if (k >= mincrf) {
				if ((k <= cut) && (sp > shortSpacer)) {//check for short IR spacers
					continue;
				}
				tmpStart = ((strti - k) + 2);// in ncbi coordinates (+1 to array cords)
				tmpStop = (strti + k + sp + 1);// in ncbi coordinates (+1 to array cords)
				if ((ndx == 0)) {//first one, can't compare current to prev if prev doesn't exist
					rightShifted = FALSE;
					leftShifted = FALSE;
					irep[ndx].start = tmpStart;
					irep[ndx].sub = tmpStop;//min loop boundary set to end by default
					irep[ndx].len = k;
					irep[ndx].loop = sp;
					irep[ndx].num = 1;
					irep[ndx].end = tmpStop;
					irep[ndx].strand = 0;
					ndx++;
				}
				else {//Not first one

					//check for immediate inclusions
					//old within new, new larger looped
					while ((irep[ndx - 1].end <= tmpStop)
							&& (irep[ndx - 1].start >= tmpStart) && irep[ndx
							- 1].len < k && ((ndx - 1) >= 0)) {
						//old within new, new better
						ndx--;//replace previous
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					//new within old, new better
					while ((irep[ndx - 1].end >= tmpStop)
							&& (irep[ndx - 1].start <= tmpStart) && irep[ndx
							- 1].len < k && ((ndx - 1) >= 0)) {
						ndx--;//replace previous
						rightShifted = FALSE;
						leftShifted = FALSE;

					}
					//old within new, old better
					if ((irep[ndx - 1].end <= tmpStop) && (irep[ndx - 1].start
							>= tmpStart) && irep[ndx - 1].len > k) {
						//don't add new
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					//new within old, old better
					if ((irep[ndx - 1].end >= tmpStop) && (irep[ndx - 1].start
							<= tmpStart) && irep[ndx - 1].len > k) {
						//don't add new
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					else if ((tmpStop == irep[ndx - 1].end) && (k == irep[ndx
							- 1].len) && (!rightShifted)) {
						leftShifted = TRUE;//to check for alternate shifting
						rightShifted = FALSE;
						ndx--;
						irep[ndx].num = irep[ndx].num + 1;//add one to Permutation count
						irep[ndx].sub = tmpStart;//adjust minimum loop boundary
						ndx++;
					}
					else if ((tmpStart == irep[ndx - 1].start) && (k
							== irep[ndx - 1].len) && (!leftShifted)) {
						rightShifted = TRUE;//to check for alternate shifting
						leftShifted = FALSE;

						//need check as rightShfited grows that it doesn't swallow previous
						if ((irep[ndx - 2].end <= tmpStop)
								&& (irep[ndx - 2].start >= tmpStart)
								&& irep[ndx - 2].len < k) {
							//old within new, new better
							//replace old with current shifting
							irep[ndx - 2] = irep[ndx - 1];
							ndx--;//replace previous
						}

						ndx--;
						irep[ndx].num = (irep[ndx].num + 1);//add one to Permutation count
						irep[ndx].end = tmpStop;//adjust end
						//mrep[ndx].sub = -12;
						ndx++;
					}
					else {//neither shifted, add new repeat
						rightShifted = FALSE;
						leftShifted = FALSE;
						irep[ndx].start = tmpStart;
						irep[ndx].sub = tmpStop;//min loop boundary set to end by default
						irep[ndx].len = k;
						irep[ndx].loop = sp;
						irep[ndx].num = 1;
						irep[ndx].end = tmpStop;
						irep[ndx].strand = 0;
						ndx++;

						for (cBack = 1; cBack <= maxcBack; cBack++) {
							while ((((irep[ndx - (1 + cBack)].end >= irep[ndx
									- 1].end) && (irep[ndx - (1 + cBack)].start
									<= irep[ndx - 1].start)) || ((irep[ndx - (1
									+ cBack)].end <= irep[ndx - 1].end)
									&& (irep[ndx - (1 + cBack)].start
											>= irep[ndx - 1].start)))
									&& ((cBack + 1) <= ndx)) {
								//maximize stem length, then minimize loop length
								if ((irep[ndx - (1 + cBack)].len == irep[ndx
										- 1].len)) {//if stems are equal, keep shortest loop
									//if previous loop is larger, delete it
									if ((irep[ndx - (1 + cBack)].loop
											> irep[ndx - 1].loop)) {
										delIRep((ndx - 1), (ndx - (1 + cBack)));
									}
								}
								//otherwise keep longest stem

								else if ((irep[ndx - (1 + cBack)].len
										< irep[ndx - 1].len)) {//previous stem is shorter = delete it
									delIRep((ndx - 1), (ndx - (1 + cBack)));
								}
								//reset everything and remove end ndx
								rightShifted = FALSE;
								leftShifted = FALSE;
								--ndx;
								cBack = 1;
							}
						}//while
					}//for cBack
				}//else first
			}//if ndx > maxcBack + 2
		}//if maxcBack>0
	}//if k>minmir
	free(LCP);
	return (ndx);
}/* END of findIR*/
