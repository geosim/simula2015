// riga commentata

//#ifndef EXE_VERS
//   #define MyPrintf acutPrintf
//#else
//   #define MyPrintf printf
//#endif

#define EXTRA	20	/* NS - 01/03/06: abbassato per recuerare memoria */

#define MAX_LEN_HANDLE		16 + 1

#define CARICO				0
#define PORTATA				1
#define CARICO_E_PORTATA	2
#define CARICO_O_PORTATA	3
#define MIN_SCABR_PRV		40
#define PG					(float) 3.14159		/* PI GRECO */
#define DENSA				(float) 1.293		/* DENSITA` ASSOLUTA */

#define ACQUA	1		/* valore tipo fluido */
#define GASB    2       /* valore tipo fluido */
#define GASM    3       /* valore tipo fluido */

#define FALSO	0
#define VERO	1

//#define GS_CAN  2

#include <math.h>

/*
**		PROTOTIPI
*/
int     simula(void);
int     finp(void);
int     finpserb(void);
int     finph24(void);
int     finp_con(FILE *fp);
int     carico_np(int stato_lato);
int     snola(int stato_lato);
int     riam(void);
int     riam24(void);
int     distrhn(void);
int     fout(void);
int     fout_con(void);
int     control_prv(void);
int     contr_serb(void);
int     calcolo_D_prv(void);
int     mem_array(void);
int     ciclo_maglie(void);
int     veribb(void);
long    ext(long *dati, long *punt, long riga, long col);
int     w_ext(long *dati, long *punt, long riga, long col, long val);
int     error(int line, int int1, int int2, double float1, double float2);
int     error24(int line, int int1, double float1, double float2);
long    cod_prv(long punt);
int     cod_vv(long punt);
int     res(long punt, float *y, float *w);
int     resh2o(long punt, float *y, float *w);
int     resgas(long punt, float *y, float *w);
double  vcol(double a, double b, double c, double d, float *e);
float   interp(float port,int pomp);
int     inter(long *ier, long maglia);
int     nnnew(long i, long j);
int     sono_serb(int nodo);
int     ripristino_del(void);
int     carico_serb_start(void);
TCHAR   *td_cad2(int passo);
int     bumbum(void);
long    new_nodo(long i);
int     ripesco_nodo(void);
void    orologio(void);
float   dtocomm(float diamcalc);
float   diamott(long i);
int     dfromv(void);
float   diam_gas(long punt, float q);
int     continuo(double car, double port);
int     chiudo_prv(void);
int     print_label(void);
long    bu(long npart);
long    codice_prv_taglio(void);
int     new_nodo(int i);
int     riallineo(void);
void	VariazioneDiametro(long reg, int DimAum, int *FineRegolazione);
int		RicalcoloDiametriRegolatori(float *MyQSum);
void	CalcolaQsum(float *QSum);
int		ChkRegolatori(void);

#ifndef EXE_VERS
	// Versione per ARX
	int	simula(void);
#else
	// Versione per EXE
	int	simula(int argc, char* argv[]);
#endif

/*  protezione */
#ifdef DOS
	int protez(void);
#endif
