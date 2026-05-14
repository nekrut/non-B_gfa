#ifndef GFA_H_
#define GFA_H_
#define MAX_REPS 2500000
//#define MAX_REPSIZE 1500
//MAX_DNA used to be a hard 300 Mbp cap on a static dna[] buffer; dna/dna2/dna3
//are now heap-allocated to the longest record's actual length and the cap is
//gone.
//#define MAX_DNA 300000000
//#define MAX_BINS 50000
#define MAXCOL 101
#define MAX_FASTA_SIZE 80
#define MAX_LINE 256

#undef FALSE
#undef TRUE
typedef enum {
	FALSE, TRUE
} BOOLEAN;

/*   Structure definition for REP   **************************
 *    typedef struct REP
 *        {
 *        int start;Starting position of repeat
 *		  int loop; Size of Loop (also kv score for Z-DNA)
 *        int pos;  Position of "i", internal use only
 *        int len;  Length of the repeat pattern
 *        int num;  Number of times pattern is repeated (permutations for MR)
 *        int end;  End position of the motif
 *        int ID;   Unique ID for each repeat for the SQL database
 *        char typ; Type of repeat: tandem, mirror, cruciform
 *        int sub;  holds type for str, Islands for G4, MinLoop
 *        	for IR and MR, and remainder for DR
 *        int special; tag for if sequence is Slippped, Cruciform, Triplex etc.
 *
 *
 *        }REP;
 **************************************************************/

typedef struct REP {
	int start;
	int loop;
//	int pos;
	int len;
	int num;
	int end;
	int sub;
	int strand; //0 for plus, 1 for minus
	int special;//1 for yes, 0 for no
} REP;



/******************************
 *        A-Tract        *****
 *****************************/
typedef struct A_Tract {
	int strt;
	short len;
} A_Tract;

/******************************
 *  potential Bent DNA      *****
 *****************************/
typedef struct potential_Bent_DNA {
	double a_center;
	int strt;
	int end;
} potential_Bent_DNA;

/******************************
 *        G-Island        *****
 *****************************/
typedef struct G_Island {
	int strt;
	int len;
	//int gap;
} G_Island;

/******************************
 *  potential G-quad      *****
 *****************************/
//typedef struct potential_G_Quads {
//	int island;
//	int islands;
//	int runSize;
//	int pruns1;
//	int pruns2;
//	int pruns3;
//} potential_G_Quads;

//these are global so findGQ can process forward and reverse strand in one pass
extern int nGisls;
extern int nCisls;

//Heap-allocated by main() at startup. dna/dna2/dna3 are sized to the longest
//record in the input; REP and island buffers are sized to MAX_REPS.
extern char *dna;
extern char *dna2; //reverse complement DNA
extern char *dna3; //complement DNA

extern G_Island *gisle;
extern G_Island *rcgisle;

extern potential_Bent_DNA *pAPRs;

extern REP *irep;
extern REP *mrep;
extern REP *drep;
extern REP *grep;
extern REP *zrep;
extern REP *srep;
extern REP *arep;

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

#endif

#endif /* GLOBVARS_H_ */

