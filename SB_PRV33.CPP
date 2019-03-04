/******************************************************************************/
/*Name: HY_PRV33.CPP
                                
	Module description: simulatore reti acqua e gas (bassa e media pressione)
         
	Author: 
                                      
	(c) Copyright 1998-2010 by IREN Acqua Gas S.P.A

   21/2/2005 Roby - Aggiunto codice GEOsim in lettura e scrittura lati
                    (Vedi vettore IDs)
   01/09/2008 Stefano - modifiche al codice delle PRV per instabilità

/******************************************************************************/

#include "stdafx.h"

#ifndef EXE_VERS
   #include <gs_def.h>
#endif

#include "hydro33.h"


#ifndef EXE_VERS
   #define MyPrintf acutPrintf
#else
   #define MyPrintf wprintf
#endif


/*************************************************************************/
/*  GLOBAL VARIABLES                                                     */
/*************************************************************************/

bool  Aborted = false;
float    Err_q, Err_h, *Qe, *Carico, *L, *D, *Dorig, *Q, *Qd, *A, *B, *QuotaZ,
      	 *Scabr, Cmq, Dens, Cstp, Pprec, Ppapp, *Ypprec, *Carico_serb,
	     *Pq[11], *Ph[11],
	     *H_sfioro24, *Volume24, *Base24, *Q_riemp24, *Z_serb24, *Tot_erog_serb24,
	     *Carico_prv, *Delta_prv, *D_prv,
	     Min_vel_prv = (float) 0.5, Min_q_prv = (float) 0.0005,
	     *H_monte, *H_valle, *Diamcomm;
float	Veloc,		/*  */
		*QMaxReg,	/* Portata massima fissata per un regolatore */
		Max_toll_D;	/*  */

// StatoRegolatore:  0 regolatore a posto
//                   1 regolatore nel quale sono in corso la diminuzione dei diametri dei lati afferenti e la portata è superiore al max consentito
// StatoRegolatore viene posto = 1 quando si diminuisce il diam di un lato che afferisce ad un regolatore. Da quel momento in poi
// quel lato può diminuire di diam e non aumentare (elimina instabilità)
int      Contr_q_h, Fluido, Ictrp, Cambio_maglie=0, Primo_giro,
	     Contr_prv=0, Cambio_prv, Chiuse_prv=0, *Nl,
        *StatoRegolatore;

long     *Nnn, *Ni, *Nf, *Np, *La, Nodi=0, Lati=0, *Coll_Serb_start,
	     Maglie=0, Lati_f=0, Nc_princ, Nc, *No, Num_pomp=0, Num_serb=1,
	     *P_nola, *D_nola, *Ivet, **D_k, Tot_num_serb=0, Nodo_partenza, 
	     *Serb24, *Quota_sfioro24, Nserb=0, Num_prv=0, *Lato_prv, Maglie_max=0;
long    *IDs; // Vettore per gli identificativi GEOsim dei lati
long     Tot_mem=0, Ftell_nodi=0, Ftell_lati=0, Ftell_cons=0;
long     Last, Last2, Last3; // Last delimita le valvole, Last2 delimita le pompe, Last3 delimita le PRV
int      Ora, Ora_iniz=0, Ora_fin=0, Ora_errore, Ver24;
int      Lancetta=-1;
int      Diametri, Calcolo_diametri=FALSO, Num_loop_diam=1;
int      Connessioni=FALSO, Debug_on=FALSO, Rete, Ciclo, f_Vincoli=1;
double   Perd_unit;
TCHAR    *Ftemp;
char     *Latodel, *Nododel, *Nodo_start, *Dtocalc;
char     **Handle;
// QMaxRegTot serve nel gas. Memorizza la portata totale regolata da regolatori con limitazione 
// di portata (se =0 vuol dire che non esistono regolatori con limitazione di portata)
float	 QMaxRegTot=0;
// GlobalIter memorizza il contatore iter di veribb. Serve globale per utilizzarlo in ResGas
long     GlobalIter=0;
// in chiudo_prv() controllo:
// se GlobalIter <= MaxGlobalIter
// allora Contr_prv = 1;              /* non stabile */
// in resgas() controllo;
// se *(Nl+punt) == -1 && GlobalIter > MaxGlobalIter allora
// calcolo le prv come tali altrimenti le valuto come lati normali
#define  MaxGlobalIter -1
// PrvAcceleratoreConvergenza smorza la convergenza delle prv
#define  PrvAcceleratoreConvergenza 0.75

FILE     *Filedatiinp,*Fileoutlati,*Fileoutnodi,*File_errore,*File_esecuzione;
FILE     *FileDEBUG, *FileDEBUG_CSV;

#ifndef EXE_VERS

// roby test
void Stampa_valori_prv(long contatore_iter);

///////////////////////////////////////////////////////////////////////////////
// inizio - Versione per ARX
int simula(void)
{
	struct resbuf *p;
	int    argc;
	TCHAR  *argv;

  	p = ads_getargs();
	argc = p->resval.rint;
	p = p->rbnext;
	argv = p->resval.rstring;

	if (argc == 3) Debug_on = VERO;

	Ftemp  =  argv;
	Primo_giro = 1;
	Ciclo = 0;

	// ads_textscr();

	int ret=GS_GOOD;

   // 14/05/98 PAOLA DEPRATI - "do ... while" fittizio affinchè venga chiuso il file <File_errore>
   do
   {
      for (Ora=Ora_iniz; Ora<=Ora_fin; Ora++) 
      {
         if (acedUsrBrk() == 1) { Aborted = true; return RTERROR; } // per uscita forzata dall'utente

         Tot_mem = 0;
		   Ora_errore = 0;  /* serve per controllare la scrittura in errore.out */
		   Cambio_prv = 0;  /* azzero flag controllo PRV */
		   if(Cambio_maglie) Rete++;
		   else Rete = 1;

		   Ciclo++;

		   if((Ora == Ora_iniz) && !(Cambio_maglie)) {
			   if((ret=finp())==GS_BAD) break;              /* lettura dati iniziali e allocazione memoria */
			   carico_np('0');                              /* carico vettore Np */
			   if((ret=snola('0'))==GS_BAD) break;          /* definisco vettore Nola */
			   if(!bumbum()) {                              /* controllo topologia */
				   if(Connessioni) {                         /* esco per controllo connessioni */
					   fout_con();
					   break;
				   }
		         // -- inizio sezione inserita da Stefano (01/08/2006) per calcolare portata e carico in reti non magliate
               if(!Cambio_prv) {
			         if((ret=distrhn())==GS_BAD) break;
		         }
   			   if((ret=riam())==GS_BAD) break;              /* assegno portate di prima approssimazione nei rami */
		         // -- fine sezione inserita da Stefano (01/08/2006) per calcolare portata e carico in reti non magliate
				   if((ret=fout())==GS_BAD) break;
				   continue;
			   }
			   if((ret=snola('1'))==GS_BAD) break;          /* CONTROLLARE !! */
			   if((ret=riam())==GS_BAD) break;              /* assegno portate di prima approssimazione nei rami */
			   if((ret=ciclo_maglie())==GS_BAD) break;      /* interprete di rete */
		   }
		   else {
			   finph24();
			   carico_np('0');                              /* carico vettore Np */
			   if((ret=snola('0'))==GS_BAD) break  ;        /* ridefinisco Nola */
			   if(!bumbum()) {                              /* controllo topologia */
				   if((ret=fout())==GS_BAD) break;
				   continue;
			   }
			   if((ret=snola('1'))==GS_BAD) break  ;        /* CONTROLLARE !! */
			   if(Cambio_maglie) {
				   if((ret=riam())==GS_BAD) break  ;         /* assegno portate di prima approssimazione nei rami */
				   if((ret=ciclo_maglie())==GS_BAD) break;   /* interprete di rete */
			   }
			   else
				   riam24();
		   }

		   if(Maglie > 0) 
            if (veribb() == FALSO) break;

		   if(!Cambio_prv) {
			   if((ret=distrhn())==GS_BAD) break;
		   }
		   if((ret=fout())==GS_BAD) break;

		   MyPrintf(_T("\n"));
	   }
   }
   while(0);

	// chiudo il file di errore se e` stato aperto
	if(File_errore != NULL) fclose(File_errore);
	if(Debug_on) fclose(FileDEBUG);
	if(Debug_on) fclose(FileDEBUG_CSV);

   if (Filedatiinp) fclose(Filedatiinp);
   if (Fileoutlati) fclose(Fileoutlati);
   if (Fileoutnodi) fclose(Fileoutnodi);

   if (Aborted)
   	MyPrintf(td_cad2(60));   /* "\nVerifica abortita\n" */
   else
	   MyPrintf(td_cad2(2));   /* "\nVerifica conclusa\n" */

	return ret;
}
// fine - Versione per ARX
///////////////////////////////////////////////////////////////////////////////

#else

///////////////////////////////////////////////////////////////////////////////
// inizio - Versione per EXE
int _tmain(int argc, _TCHAR* argv[])
{
/*
** Parametri di avvio:
** primo = valore fittizio: 3=debug
** secondo= path file dati di input e output
*/
   if (argc == 3) Debug_on = VERO;
   if (argc < 2) return 0;

	int ret=VERO;

	// 02/08/96 non serve Ftemp = (char *) calloc(150, sizeof(char));
	Ftemp  =  argv[2];
   TCHAR   ftempinp[150];
	//Ftemp = (char *) calloc(150, sizeof(char));
	//strcpy(Ftemp, "c:/programmi/topobase3/bin/");
	Primo_giro = 1;
	Ciclo = 0;

	File_esecuzione = _wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("esecuzione.txt")), _T("w+"))

	// 14/05/98 PAOLA DEPRATI - "do ... while" fittizio affinchè venga chiuso il file <File_errore>
	do
	{
		for (Ora=Ora_iniz; Ora<=Ora_fin; Ora++) {
			Tot_mem = 0;
			Ora_errore = 0;  // serve per controllare la scrittura in errore.out 
			Cambio_prv = 0;  // azzero flag controllo PRV 
			if(Cambio_maglie) Rete++;
			else Rete = 1;

			Ciclo++;

			if((Ora == Ora_iniz) && !(Cambio_maglie)) {
				if((ret=finp())==FALSO) break;				// lettura dati iniziali e allocazione memoria
				carico_np('0');								// carico vettore Np
				if((ret=snola('0'))==FALSO) break;			// definisco vettore Nola
				if(!bumbum()) {								// controllo topologia 
					if(Connessioni) {						// esco per controllo connessioni 
						fout_con();
						break;
					}
					if((ret=fout())==FALSO) break;
					continue;
				}
				if((ret=snola('1'))==FALSO) break;          // CONTROLLARE !!
				if((ret=riam())==FALSO) break;              // assegno portate di prima approssimazione nei rami
				if((ret=ciclo_maglie())==FALSO) break;      // interprete di rete
			}
			else {
				finph24();
				carico_np('0');								// carico vettore Np
				if((ret=snola('0'))==FALSO) break  ;        // ridefinisco Nola
				if(!bumbum()) {								// controllo topologia
					if((ret=fout())==FALSO) break;
					continue;
				}
				if((ret=snola('1'))==FALSO) break  ;        // CONTROLLARE !!
				if(Cambio_maglie) {
					if((ret=riam())==FALSO) break  ;		// assegno portate di prima approssimazione nei rami
					if((ret=ciclo_maglie())==FALSO) break;	// interprete di rete
				}
				else
					riam24();
			}

		   if(Maglie > 0) 
            if (veribb() == FALSO) break;

			if(!Cambio_prv) {
				if((ret=distrhn())==FALSO) break;
			}
			if((ret=fout())==FALSO) break;

			MyPrintf(_T("\n"));
		}
	}
	while(0);

	// chiudo il file di errore se e` stato aperto
   if (File_errore) fclose(File_errore);
   if (File_esecuzione) fclose(File_esecuzione);
   if (Debug_on && FileDEBUG) fclose(FileDEBUG);
   if (Debug_on && FileDEBUG_CSV) fclose(FileDEBUG_CSV);
	
	MyPrintf(td_cad2(2));   // "\nVerifica conclusa\n"
	fwprintf(File_esecuzione, td_cad2(2));

	return ret;
}
// fine - Versione per EXE
///////////////////////////////////////////////////////////////////////////////
#endif


/* =================================
** FINPSERB - gestisce l'input dei serbatoi (h24)
*/
int finpserb()
{
	long kkk;

	fscanf(Filedatiinp,"%d",&Ora_fin);         /* ora inizio e fine simul, numero serbatoi */
	kkk=1;
	while(kkk<=Nserb) {
		fscanf(Filedatiinp,"%d,%f,%f,%f,%d", &Serb24[kkk], &H_sfioro24[kkk], &Volume24[kkk],
											&Z_serb24[kkk], &Quota_sfioro24[kkk]);
		if (*(H_sfioro24+kkk)<=0){
			if(error(18, *(Serb24+kkk), 0, 0, 0)==FALSO) return FALSO;
			return FALSO;
		}
		if (*(Volume24+kkk)<=0){
			if(error(19, *(Serb24+kkk), 0, 0, 0)==FALSO) return FALSO;
			return FALSO;
		}
		if (*(Z_serb24+kkk)<0){
			if(error(20, *(Serb24+kkk), 0, 0, 0)==FALSO) return FALSO;
			return FALSO;
		}
		*(Base24 + kkk) = *(Volume24 + kkk) / *(H_sfioro24 + kkk);
		kkk++;
	}

	MyPrintf(td_cad2(0));  /* "\n\n   GEOsim PRESSURE NETWORK simulator h24 v.1\n\n" */
   #ifdef EXE_VERS
	   fwprintf(File_esecuzione, td_cad2(0));
   #endif

	return VERO;
}
/* ============================ 
**
** Calcola il valore di Qsum
** QSum = Sommatoria delle portate erogate da ogni regolatore
** QSum = -1 ==>  Errore
**
*/
void CalcolaQsum(float *QSum)
{	   
	long	li;

	for(li=1; li<=Num_serb; li++) *(QSum+li)=0;	// azzeramento QSum
	
	li = 1;
	while (li<=Lati) {	// calcolo della QSum per i lati che ha il regolatore
		if((*(Ni+li) <= Num_serb) && (*(Q+li) < 0)) {
			/* Ni è Serb e Q < 0 */
			*(QSum + *(Ni+li))  = *(QSum + *(Ni+li)) + fabs(*(Q + li));
		}

		if((*(Ni+li) <= Num_serb) && (*(Q+li) > 0)) {
			/* Ni è Serb e Q > 0 */
			*(QSum + *(Ni+li))  = *(QSum + *(Ni+li)) - fabs(*(Q + li));
		}

		if((*(Nf+li) <= Num_serb) && (*(Q+li) < 0)) {	
			/* Nf è Serb e Q < 0 */
			*(QSum + *(Nf+li))  = *(QSum + *(Nf+li)) - fabs(*(Q + li) - *(Qd + li));
		}

		if((*(Nf+li) <= Num_serb) && (*(Q+li) > 0)) {	
			/* Nf è Serb e Q > 0 */
			*(QSum + *(Nf+li))  = *(QSum + *(Nf+li)) + fabs(*(Q + li) - *(Qd + li));
		}

		li++;
	}	
	return;
}
/* ============================
**
** Subroutine VariazioneDiametro:
** Modifico il diametro dei lati che convergono nel regolatore che supera QMaxReg imposta
** DiminuisciAumenta =0 diminuisce diam
** DiminuisciAumenta =1 aumenta diam
** 
*/
void VariazioneDiametro(long regolatore, int DiminuisciAumenta, int *FineRegolazione)
{
	long	li;
	float Range=0.0005f; // % con cui varia il diametro (precedente +- precedente * Range)
   Range=0.05f;
	
	li = 1;
	
	while (li<=Lati) {
		if(*(Ni+li) == regolatore || *(Nf+li) == regolatore){
			if (DiminuisciAumenta == 0) {
				*(D+li) = *(D+li) - (*(D+li) * Range);
				if (*(D+li) > (*(Dorig+li)))
					*(D+li) = *(Dorig+li);
				else
					*FineRegolazione = 0; // Ho variato il diametro, quindi continuo
			}
			else {
				*(D+li) = *(D+li) + (*(D+li) * Range);
				if (*(D+li) > (*(Dorig+li)))
					*(D+li) = *(Dorig+li);
				//else
				//	*FineRegolazione = 0; // Ho variato il diametro, quindi continuo
			}
		}
		li++;
	}
	return;
}

/* ============================
**
** Subroutine RicalcoloDiametriRegolatori:
** Ricalcola i diametri dei lati che convergono nel regolatore che supera QMaxReg imposta
** 1 = Non ho soddisfatto le condizioni dei regolatori... faccio altri giri...
** 0 = OK tutti i regolatori a posto.
**
*/
int RicalcoloDiametriRegolatori(float *MyQSum)
{
	long	li=1;	// nodo
	int		FineRegolazione; // 1 non ho variato niente; 0 devo continuare perchè ho variato diam
	int		Aumento=1;		// Aumento diametro
	int		Diminuisco=0;	// Diminuisco il diametro
	double	Range=0.049;		// precisione di controllo portata
	
	FineRegolazione = 1;

	while (li<=Num_serb) { // per ogni regolatore che supera il valore imposto
		if(*(Nododel+li) == '0' &&	// nodo che è connesso alla rete
			*(MyQSum+li) < (-1)* *(QMaxReg + li) && 
			(*(QMaxReg + li)>0) &&
			(
			(*(MyQSum+li) <= (-1)* (*(QMaxReg + li) + Range)) || 
			(*(MyQSum+li) >= (-1)* (*(QMaxReg + li) - Range))
			)
			)
		{
			// Il regolatore li supera QMaxreg e diminuisco la sezione dei lati che convergono
			// nel regolatore 
			VariazioneDiametro(li, Diminuisco, &FineRegolazione);	// Vario la sezione...
         *(StatoRegolatore + li) = 1;
			//FineRegolazione = 0;
			}

		if(*(Nododel+li) == '0' &&	// nodo che è connesso alla rete
         *(StatoRegolatore + li) < 1 &&
			*(MyQSum+li) > (-1)* *(QMaxReg + li) && 
			(*(QMaxReg + li)>0) &&
			(
			(*(MyQSum+li) <= (-1)* (*(QMaxReg + li) + Range)) || 
			(*(MyQSum+li) >= (-1)* (*(QMaxReg + li) - Range))
			)
			)
		{
			// Il regolatore li ha la portata inferiore a QMaxreg e aumento la sezione dei lati che convergono
			// nel regolatore 
			VariazioneDiametro(li, Aumento, &FineRegolazione);	// Vario la sezione...
			//FineRegolazione = 0;
			}

		li++;
	}

	return(FineRegolazione);	
}
/* ============================
**
**  Subroutine ChkRegolatori:
**	Controlla se la portata di tutti i regolatori (GAS) e valuta se spalmare eventuali
**	eccedenze sui regolatori successivi che non superano i loro valore massimo.
**	consentito
**
*/
int ChkRegolatori()
{	   
	/* long	li;	 Nodo che supera il QMaxreg */
	float	*MyQSum;
	int FineRegolaz=0; // Non ho finito... continuo iterazione su variazione diametri...

	/* Alloco lo spazio per MyQSum = Sommatoria delle portate erogate da ogni regolatore */
	MyQSum = (float *) calloc((Num_serb+1) , sizeof(float));
	if(MyQSum == NULL) return(-1);	// Errore... fine

	/* Calcolo la sommatoria delle portate (MyQSum) erogate da ogni regolatore */
	CalcolaQsum(MyQSum);

	FineRegolaz = RicalcoloDiametriRegolatori(MyQSum);

	free (MyQSum);

	return (FineRegolaz);	/* FineRegolaz= 1 => Tutto OK (soddisfo le condiz del regolatore)*/
}

/* ====================================
** FINP: leggo i dati di input
*/
int finp()
{
	int		j, i3, rc;
	long	li, lj, li1, li2, lik, ljj3, temp_num_prv=0, id;
	float	f1, f2, f3, f4, qtot= (float) 0.0, deltah_prv;
	TCHAR	ftempinp[150];
	char	c1;
	FILE	*fp;

   int   ret=VERO;
   
   // 14/05/98 PAOLA DEPRATI - "do ... while" fittizio affinchè venga chiuso il file <Filedatiinp>
   do
   {
      if((Filedatiinp=_wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("data.inp")),_T("r")))==NULL) 
      {
		   MyPrintf(_T("\nfmgr error data.inp\n"));
		   ret=FALSO;
         break;
	   }

	   fscanf(Filedatiinp,"%d,%d",&Fluido,&Ver24); /* Debug_on */

	   if(Debug_on) FileDEBUG = _wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("debug.txt")), _T("w+"));
	   if(Debug_on) FileDEBUG_CSV = _wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("debug_test.csv")), _T("w+"));

	   if(Ver24 == 3) {  /* controllo connessioni */
		   Ver24 = FALSO;
		   Connessioni = VERO;
		   finp_con(Filedatiinp);
         break;
      }

      fscanf(Filedatiinp,"%ld,%ld,%ld,%ld,%ld,%d,%ld", &Nodi,&Lati,&Num_serb,&Num_pomp,&Num_prv,&Contr_q_h, &Last);
	   fscanf(Filedatiinp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &Err_q,&Err_h,&Min_vel_prv,&Min_q_prv,
						   &Cmq,&Cstp,&Pprec,&Ppapp,&Veloc,&Dens); /*Perd_unit, Max_toll_D*/
	   Perd_unit = (float) 4.6;
	   Max_toll_D = (float) 0.0;
	   if(Fluido == GASM) Perd_unit = pow((double) Perd_unit, (double) 2.0);

						   // Last valvole
	   Last2 = 2 * Last;	// pompe
	   Last3 = 3 * Last;	// Prv

	   if(Fluido == ACQUA) {
		   Err_q = Err_q / (float) 1000.0;          /* da l/s in mc/s */
		   Min_q_prv = Min_q_prv / (float) 1000.0;  /* da l/s in mc/s */
	   }

	   if(Ver24 == 2) {  /* calcolo diametri */
		   Ver24 = FALSO;
		   Calcolo_diametri = VERO;
		   fscanf(Filedatiinp,"%d", &Num_loop_diam);
	   }

	   Lati_f = Num_serb - 1;
	   Maglie = Lati - Nodi + 1 + Lati_f;
      /*
      Stefano 08/08/2011
      Il calcolo delle maglie è errato. I dati che vengono in questo momento processati non afferiscono ad una rete
      connessa. I dati possono far riferimento a più reti. Il calcolo del numero delle Maglie può quindi non essere corretto.
      La funzione riallineo() chiamata in seguito ricalcola il numero delle Maglie prendendo in considerazione solo
      i nodi e lati connessi. Per questo motivo forzo il numero delle Maglie a 1 per non variare la struttura del programma
      ben sapendo che poi, in riallineo(), il numero corretto di Maglie viene ricalcolato
      */
      Maglie = 1;
	   Nserb = Num_serb - Num_pomp;  /* serbatoi con volume */

	   #ifdef DOS      /* RAM DISPONIBILE - SPECIFICA DEL DOS */
         /*if(Ora == Ora_iniz) per HIGHC
           MyPrintf(td_cad2(3), coreleft());*/ /* "\nMemoria disponibile : %ld bytes" */
         #ifdef EXE_VERS
		   //fwprintf(File_esecuzione, td_cad2(3), coreleft());
         #endif
	   #endif

	   rc = mem_array();   /* DEFINIZIONE DINAMICA DEGLI ARRAY */
	   if (rc)
      {
		   /*  "\nMEMORIA INSUFFICIENTE interruzione al passo %d, allocati %ld bytes" */
		   MyPrintf(td_cad2(4), rc, Tot_mem);
			#ifdef EXE_VERS
            fwprintf(File_esecuzione, td_cad2(4), rc, Tot_mem);
         #endif
		   ret=FALSO;
         break;
	   }

	   if(Ver24) finpserb();
	   else 
		{
			MyPrintf(td_cad2(5));  /* "\n\n   GEOsim PRESSURE NETWORK simulator h1 v.5\n\n" */
			#ifdef EXE_VERS
            fwprintf(File_esecuzione, td_cad2(5));		
         #endif
		}

	   Nc = 1;
	   Nodo_partenza = 1;
	   Nc_princ = 1;

	   /* teoricamente Z(i) per i serbatoi non servirebbe: il programma lo legge
          in cons.inp. */
	   li=1;
	   while (li<=Nodi) {	/* NS - Gestione regolatori */
			if((Fluido == ACQUA) || (li > Num_serb)) {
				fscanf(Filedatiinp,"%ld,%f", &li1, &f1);
			}
			else {	/* solo per GASB, GASM: leggo in QMaxReg portata max fissata x ogni regolatore */
				fscanf(Filedatiinp,"%ld, %f, %f", &li1, &f1, &f4);
				*(QMaxReg+li) = f4;		/* portata imposta al regolatore, definiamo -1 = regol. senza vicoli */
				/*if (f4 == -1.0)
					f_Vincoli = 0;	// ho un regolatore senza limitazioni
				else
					QMaxRegTot = QMaxRegTot + f4;*/
			}
			*(Nnn+li) = li1;
			*(QuotaZ+li) = f1;		
			*(Qe+li) = (float) 0.0;	/* serbatoi senza consumo */
			li++;
	   }

	   li=1;
	   while (li<=Lati)
      {
		   fscanf(Filedatiinp,"%ld,%ld,%ld,%f,%f,%c,%d,%f",
		                       &id,&li1,&li2,&f1,&f2,&c1,&i3,&f3);

		   if(li1 == li2) {
			   /* Numerazione lati errata - nodo iniziale %d ripetuto anche come finale */
			   if(error(2, li1, 0, 0, 0)==FALSO) ret=FALSO;
			   ret=FALSO;
            break;
		   }
		   lik=1;
		   while (lik<=Nodi) {
			   if(*(Nnn+lik) == li1) *(Ni+li)=lik;
			   if(*(Nnn+lik) == li2) *(Nf+li)=lik;
			   lik++;
		   }

		   if(*(Ni+li) == 0) {
			   /* "\nErrore: il lato %d %d punta al nodo %.0f che non esiste." */
			   if(error(22, li1, li2, (float) li1, 0)==FALSO) ret=FALSO;
			   ret=FALSO;
            break;
		   }
		   if(*(Nf+li) == 0) {
			   /* "\nErrore: il lato %d %d punta al nodo %.0f che non esiste." */
			   if(error(22, li1, li2, (float) li2, 0)==FALSO) ret=FALSO;
			   ret=FALSO;
            break;
		   }

		   if(*(Ni+li) > (*(Nf+li))) {
			   lik = *(Ni+li);
			   *(Ni+li) = *(Nf+li);
			   *(Nf+li) = lik;
		   }

		   *(IDs + li) = id; // memorizzo ID GEOsim del lato
		   *(Nl  + li) = i3;
		   *(L   + li) = f1;
		   *(D   + li) = f2;
		   *(Dorig + li) = f2;	/* Carico i diametri originali */

		   if(*(D+li) <= 0.0) {
			   /* "\nErrore: il lato %d %d ha il diametro = 0." */
			   if(error(23, li1, li2, 0, 0)==FALSO) ret=FALSO;
			   ret=FALSO;
            break;
		   }

		   if(Calcolo_diametri == VERO)
			   *(Dtocalc+li) = c1;
		   if(Fluido == ACQUA) {
			   *(Scabr+li) = f3;
			   *(D+li) = *(D+li) * (float) 0.001;         /* da mm in m */
		   }
		   else
			   *(L+li) = *(L+li) * (float) 0.001;         /* da m in km */

		   switch (*(Nl+li)) {
		   case -1:  /* valvole PRV */
			   ++temp_num_prv;
            *(Lato_prv+temp_num_prv) = li;
            if(*(Nnn+ *(Nf+li)) >= Last3) {
				   deltah_prv = *(QuotaZ+ *(Nf+li));
				   *(QuotaZ+ *(Nf+li)) = *(QuotaZ+ *(Ni+li));
            }
            else {
				   deltah_prv = *(QuotaZ+ *(Ni+li));
				   *(QuotaZ+ *(Ni+li)) = *(QuotaZ+ *(Nf+li));
            }
				/* GAS_PRV: in deltah_prv trovo la press di regolaz + la quota allora tolgo la
				quota e la reinserisco in parti di gas */
				if(Fluido == GASB)
					deltah_prv = (deltah_prv - *(QuotaZ+ *(Nf+li))) + *(QuotaZ+ *(Nf+li))*DENSA*(Dens - (float) 1.0);
				else if(Fluido == GASM)
					deltah_prv = (deltah_prv - *(QuotaZ+ *(Nf+li))) + *(QuotaZ+ *(Nf+li))*DENSA*(Dens - (float) 1.0) / (float) 10000.0;

            *(Carico_prv+temp_num_prv) = deltah_prv;
            *(Delta_prv+temp_num_prv) = (float) 0.0;
		      if(Fluido == ACQUA)
   			   *(L+li) = (float) 1.0;
		      else
			      *(L+li) = (float) 0.001;         /* da m in km */
            *(D+li) = (float) 0.0005;
			   if(Fluido == ACQUA) *(Scabr+li) = (float) 80.0;
            break;
         case -2:  /* valvole ON/OFF */
		      if(Fluido == ACQUA)
   			   *(L+li) = (float) 1.0;
		      else
			      *(L+li) = (float) 0.001;         /* da m in km */
            *(D+li) = (float) 1.0e10;
			   if(Fluido == ACQUA) *(Scabr+li) = (float) 120.0;
            break;
		   }
		   li++;
	   }

      if (ret==FALSO)
         break;

      li=2;
	   while(li<=Num_serb) {
		   lj = Lati+li-1;
		   *(Ni+lj) = 1;
		   *(Nf+lj) = li;
		   *(D+lj)  = (float) 1.0e10;
		   if(Fluido == ACQUA)
   			*(L+lj) = (float) 1.0;
		   else
			   *(L+lj) = (float) 0.001;         /* da m in km */
		   if(Fluido == ACQUA)  *(Nl+lj) = 0;
		   li++;
	   }

	   fclose(Filedatiinp); Filedatiinp = NULL;

		/* Lettura consumi */
	   if((Filedatiinp=_wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("cons.inp")), _T("r")))==NULL) {
		   MyPrintf(_T("\nfmgr error cons.inp\n"));
			ret=FALSO;
         break;
	   }

	   li=1;
	   while (li<=Num_serb) {
		   if(Ver24) {
			   fscanf(Filedatiinp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
                                &f1,&f2,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3);
			   *(Carico+li) = f1;
			   if(li<=Nserb) {
				   *(Q_riemp24+li) = f2;
				   /* calcolo volumi iniziali serbatoi (solo se h24) */
				   *(Volume24+li) = (*(Base24+li) * (*(Carico+li) - *(Z_serb24+li)));
			   }
		   }
		   else {
			   fscanf(Filedatiinp,"%f", &f1);
			   *(Carico+li) = f1;
		   }
		   /* nel caso del gas arriva nel carico la pressione di regolazione
		    a cui va AGGIUNTA la quota (in parti di gas) */
		   if(Fluido == GASB)
			   *(Carico+li) = *(Carico+li)+*(QuotaZ+li)*DENSA*(Dens - (float) 1.0);
		   else if(Fluido == GASM)
			   *(Carico+li) = *(Carico+li)+*(QuotaZ+li)*DENSA*(Dens - (float) 1.0) / (float) 10000.0;
		   li++;
	   }

	   li=1;
	   while (li<=Num_serb) {
		   *(Carico_serb+li) = *(Carico+li);
		   if(li>1) {
			   lj = Lati+li-1;
			   if(Fluido == GASM) /* bug eliminato 27/12/95 */
				   *(Scabr+lj) = (*(Carico+Nc)-*(Carico+li)) * (*(Carico+Nc)+*(Carico+li));
			   else
				   *(Scabr+lj) = (*(Carico+Nc)-*(Carico+li));
			   *(Q+lj) = (float) 0.0;
		   }
		   li++;
	   }

	   li=Num_serb+1;
	   while (li<=Nodi) {
		   if(Ver24)
			   fscanf(Filedatiinp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
							    &f1,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3);
		   else
			   fscanf(Filedatiinp,"%f", &f1);
		   *(Qe+li) = f1 * Cmq;
		   if(Fluido == ACQUA) *(Qe+li) = *(Qe+li) / (float) 1000.0;  /* da l/s in mc/s */
		   qtot = qtot + (*(Qe+li));
		   li++;
	   }

	   li=1;
	   while(li<=Lati) {
		   if(Ver24)
		    fscanf(Filedatiinp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
							   &f1,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3,&f3);
		   else
		    fscanf(Filedatiinp,"%f", &f1);

		   switch (*(Nl+li)) {
           case -1:  /* valvole PRV */
               break;
		   case -2:  /* valvole ON/OFF (VV) */
               *(Qd+li) = f1; /* -1 chiusa 0 aperta */
               /* azzero La ovvero elimino lati da costruzione maglie */
               *(La+li) = (f1 == -1) ? -10000 : 0;
               break;
           default:
               if(*(Nl+li)>20) {  /* pompa */
				   *(Qd+li) = f1; /* -1 non attiva 0 attiva */
               }
               else {
				   *(Qd+li) = f1 * Cmq;
				   if(Fluido == ACQUA) *(Qd+li) = *(Qd+li) * (float) 0.001;  /* da l/s in mc/s */
				   if(*(Qd+li) < 0.0) *(Qd+li) = (float) 0.0;
				   qtot = qtot + (*(Qd+li));
			   }
		   }
		   *(Q+li) = (float) 0.0;
		   li++;
	   }

	   Ftell_cons = ftell(Filedatiinp); /* serve per h24 */

	   ripristino_del();

	   /* "\nDati input:\nNodi: %d   Lati: %d"  */
	   MyPrintf(td_cad2(7), Nodi, Lati);
		#ifdef EXE_VERS
         fwprintf(File_esecuzione, td_cad2(7), Nodi, Lati);
      #endif

	   if(Fluido == ACQUA) {
		   /* "\nStazioni di pompaggio: %3d\nValvole PRV: %3d" */
		   MyPrintf(td_cad2(8),Num_pomp,Num_prv);
		   /* "\nErogazione totale: %9.2f l/s\n"  */
		   MyPrintf(td_cad2(10),qtot*1000.0);
         #ifdef EXE_VERS
            fwprintf(File_esecuzione, td_cad2(8), Num_pomp, Num_prv);
            fwprintf(File_esecuzione, td_cad2(10), qtot*1000.0);
         #endif
	   }
	   else {
		   /* "\nErogazione totale: %9.2f mc/h\n" */
		   MyPrintf(td_cad2(11),qtot);
         #ifdef EXE_VERS
            fwprintf(File_esecuzione, td_cad2(11), qtot);
         #endif
	   }

	   if(Ver24)
      {
		   MyPrintf(td_cad2(12),Ora); /* "\nSimulazione ora: %3d\n" */
         #ifdef EXE_VERS
   			fwprintf(File_esecuzione, td_cad2(12), Ora);
         #endif
	   }

	   if(Num_serb == 0) {
		   /* Non esistono serbatoi/regolatori o punti di immissione in rete */
		   if(error(3, 0, 0, 0, 0)==FALSO) return FALSO;
		   return FALSO;
	   }

	   calcolo_D_prv();  /* assegno diametri a prv */

	   *(Qe+1) = (qtot > 0.0) ? (-qtot) : (float) 0.0;

	   /*
	   ** lettura curve caratteristiche stazioni pompaggio (solo h2o)
	   */

	   if(Num_pomp > 0) {
 		   if((fp=_wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("stp.inp")), _T("r")))==NULL) {
			   MyPrintf(_T("\nfmgr error stp.inp\n"));
			   return FALSO;
		   }
		   li=1;
		   while(li<=Num_pomp) {
			   ljj3 = 0;
			   lj=1;
			   while(lj<=Lati) {
				   if(*(Nl+lj) == (20+li)) {
					   ljj3 = lj;
					   break;
				   }
				   lj++;
			   }
			   for(j=1; j<=10; j++) {
				   fscanf(fp,"%f,%f",&f1, &f2);
				   Pq[j][li] = f1 * ((*(Nnn+*(Nf+ljj3)) > Last2) ? (float) -0.001 : (float) 0.001);
				   Ph[j][li] = f2;
			   }
			   *(Ypprec+li) = Ph[10][li] + (float) fabs(Ph[1][li] - Ph[10][li]) / (float) 2.0;
			   li++;
		   }
		   fclose(fp);
	   }

		/* 
		** NS - Se almeno uno dei regolatori non ha vincoli (cioè per convenzione è posto a -1) 
		** allora non faccio controlli altrimenti verifico che la portata massima dei regolatori
		** non superi il consumo totale della rete (qtot)
		*/
	   /*if (f_Vincoli){
		   if (qtot > QMaxRegTot) {
				MyPrintf("\n*** Fine **** \nLimitazione imposta ai regolatori: %9.2f", QMaxRegTot);
				MyPrintf(" minore del consumo: %9.2f \n*************", qtot);
				if (Filedatiinp) fclose(Filedatiinp);
				return FALSO;
			}
		}*/

	   /*   primo controllo dati immessi  */

	   li=1;
	   while(li<=Lati+Lati_f) {
		   ljj3 = 0;
		   lj=1;
		   while(lj<=Lati+Lati_f) {
			   if(*(Ni+lj) == *(Ni+li) && *(Ni+lj) != 0 && *(Ni+li) != 0)
				   if(*(Nf+lj) == *(Nf+li) && *(Nf+lj) != 0 && *(Nf+lj) != 0) {
					   ljj3++;
					   if(ljj3 > 1) {
						   /* Asola, vedi lati con nodo iniziale %d e nodo finale %d */
						   if(error(6, *(Nnn+(*(Ni+lj))), *(Nnn+(*(Nf+lj))), 0, 0)==FALSO) return FALSO;
						   if(li>Lati || lj>Lati)
							   /* L'asola e' composta dal lato fra due serbatoi\ne il corrispondente lato fittizio */
							   if(error(25, 0, 0, 0, 0)==FALSO) return FALSO;
						   return FALSO;
					   }
			   }
			   lj++;
		   }
		   li++;
	   }

   }  
   while(0);

   if (Filedatiinp) { fclose(Filedatiinp); Filedatiinp = NULL; }

   return ret;
}

/* =====================
** FINPH24 - gestisce l'input da database nel caso di h24
**           ripristina l'input anche per PRV chiusa
*/

int finph24()
{
	long    li, lj;
	float   cons[24], qtot = (float) 0.0;

	if(Ver24) {
		Nc = 0;     /* riazzero serbatoio principale (inserito per chiusura serb) */
		Lati_f = 0;     /* riazzero numero lati fittizi */
		fseek(Filedatiinp,Ftell_cons,SEEK_SET);

		li=1;
		while (li<=Nodi) {
			if(!Cambio_maglie) {
				fscanf(Filedatiinp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
					                &cons[1],&cons[2],&cons[3],&cons[4],&cons[5],&cons[6],&cons[7],&cons[8],&cons[9],&cons[10],&cons[11],
						            &cons[12],&cons[13],&cons[14],&cons[15],&cons[16],&cons[17],&cons[18],&cons[19],&cons[20],&cons[21],&cons[22],&cons[23]);
				*(Qe+li) = cons[Ora] * Cmq;
				if(Fluido == ACQUA) *(Qe+li) = *(Qe+li) / (float) 1000.0;  /* da l/s in mc/s */
			}

			if(li<=Nserb) {
				if(*(Serb24+li) > 0) { /* cioe` serbatoio aperto (se c_serb < 0 => serb chiuso) */
					if (!Cambio_maglie) *(Q_riemp24+li) = *(Qe+li)/Cmq;
					*(Qe+li) = (float) 0.0;  /* consumo serb = 0.0 */
					*(Carico+li) = *(Volume24+li) / *(Base24+li) + *(Z_serb24+li);
					*(QuotaZ+li) = *(Carico+li);
					if(Nc != 0) {
						Lati_f++;
						*(Ni+(Lati+Lati_f)) = Nc;
						*(Nf+(Lati+Lati_f)) = li;
					}
					else
						Nc = li;
				}
				else      /* cioe' serb chiuso */
					*(Qe+li) = (float) 0.0;     /* azzero erogazione preced. (consumo serb = 0.0) */
			}

			if((*(Nnn+li)>Last2) && (*(Nnn+li)<Last3)) { /* cioe` pompa */
			if(Nc != 0) {
				Lati_f++;
				*(Ni+(Lati+Lati_f)) = Nc;
				*(Nf+(Lati+Lati_f)) = li;
			}
			else
				Nc = li;
			}
			if(li > Num_serb) *(Carico+li) = (float) 0.0; /* nodi normali */
			qtot = qtot + (*(Qe+li));
			li++;
		}
	}
	else {
		li=1;
		while (li<=Num_serb) {
			*(Carico+li) = *(Carico_serb+li);
			li++;
		}
		li=Num_serb+1;
		while (li<=Nodi) {
			*(Carico+li) = (float) 0.0;
			qtot = qtot + (*(Qe+li));
			li++;
		}
	}

	Maglie = Lati - Chiuse_prv - Nodi + 1 + Lati_f;

	li=1;
	while(li<=Lati) {
		if(!Cambio_maglie) {
			fscanf(Filedatiinp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
								&cons[1],&cons[2],&cons[3],&cons[4],&cons[5],&cons[6],&cons[7],&cons[8],&cons[9],&cons[10],&cons[11],
								&cons[12],&cons[13],&cons[14],&cons[15],&cons[16],&cons[17],&cons[18],&cons[19],&cons[20],&cons[21],&cons[22],&cons[23]);

			switch (*(Nl+li)) {
            case -1:  /* valvole PRV */
				break;
            case -2:  /* valvole ON/OFF (VV) */
				*(Qd+li) = cons[Ora]; /* -1 chiusa 0 aperta */
				/* azzero La ovvero elimino lati da costruzione maglie */
				*(La+li) = (cons[Ora] == -1) ? -10000 : 0;
				break;
            default:
				if(*(Nl+li)>20) {  /* pompa */
					*(Qd+li) = cons[Ora]; /* -1 non attiva 0 attiva */
					if(cons[Ora] == -1) {
						lj=1;
						while(lj<=Tot_num_serb) {
							if(*(Coll_Serb_start+lj) == *(Nf+li)) {
								*(Coll_Serb_start+lj) = -1 * *(Coll_Serb_start+lj);
								break;
							}
							lj++;
						}
					}
				}
				else {
					*(Qd+li) = cons[Ora] * Cmq;
					if(Fluido == ACQUA) *(Qd+li) = *(Qd+li) * (float) 0.001;  /* da l/s in mc/s */
					if(*(Qd+li) < 0.0) *(Qd+li) = (float) 0.0;
					qtot = qtot + (*(Qd+li));
				}
			}
		}
		*(Q+li) = (float) 0.0;
		li++;
	}

	if(!Cambio_maglie) 
   {
		MyPrintf(td_cad2(12),Ora); /* "\nSimulazione ora: %3d\n" */

      #ifdef EXE_VERS
         fwprintf(File_esecuzione, td_cad2(12), Ora);
      #endif

		if(Fluido == ACQUA)
      {
			MyPrintf(td_cad2(10),qtot*1000.0); /* "\nErogazione totale: %9.2f l/s\n"  */
         #ifdef EXE_VERS
            fwprintf(File_esecuzione, td_cad2(10), qtot*1000.0);
         #endif
		}
		else 
      {
			MyPrintf(td_cad2(11),qtot); /* "\nErogazione totale: %9.2f mc/h\n" */
         #ifdef EXE_VERS
   			fwprintf(File_esecuzione, td_cad2(11), qtot);
         #endif
		}
	}

	*(Qe+Nc) = (-qtot);  /* assegno erogazione a serbatoio principale */

	li=Lati+1;
	while(li<=Lati+Lati_f) {
		*(Scabr+li) = (*(Carico+*(Ni+li))-*(Carico+*(Nf+li))); /* sbagliato se GASM */
		*(Q+li) = (float) 0.0;
		li++;
	}

	li=1;
	while(li<=Num_prv) {
		*(D+labs(*(Lato_prv+li))) = *(D_prv+li);
		*(Delta_prv+li) = (float) 0.0;
		li++;
	}

	/* azzero La ovvero elimino lati da costruzione maglie */
	li=1;
	while(li<=Lati+Num_serb-1) {
		*(La+li) = ((cod_prv(li) < 0)||(cod_vv(li) == 1)) ? -10000 : 0;
		li++;
	}

	return VERO;
}

/* ===========================
**  FOUT.C - Gestisce l'output a database
**           Scrittura su file dei risultati finali
*/
int fout()
{
//	int     i;
	long    li, ljl, lj;
	float   v = (float) 0, qi, qf, qp, p, qie, qfe, pp, Dens_e;
	TCHAR   ftempinp[150];

   if (Aborted) return FALSO;

	if(Lati_f != 0) {		/* RIDISTRIBUISCO PORTATE LATI FITTIZI */
		li=Lati+1;
		while(li<=Lati+Lati_f) {
			ljl = *(Nf+li);
			*(Qe+Nc) = *(Qe+Nc) + (*(Q+li));
			*(Qe+ljl) = *(Qe+ljl) - (*(Q+li));
			li++;
		}
	}

	//GAS_PRV: if(Fluido == ACQUA) {
		control_prv();      /* controllo valv PRV con velocita` minima */
		if(Cambio_maglie) {
			ripristino_del();
			return VERO; /* esco subito se trovo PRV da chiudere */
		}
	if(Fluido == ACQUA) { //GAS_PRV: spostato qua
		li=1;
		while(li<=Nodi) {  /* stefano 05/01/96 */
			/* if(*(Carico + i) == 0 && *(Nnn+i) < Last && *(Nododel+i) == '0') stefano 10/01/96*/
			if(*(Carico + li) < 0 && *(Nnn+li) < Last && *(Nododel+li) == '0') {
				/* Il nodo %d ha carico < 0; \nATTENZIONE, i risultati della simulazione possono non essere attendibili! */
				if(error(24, *(Nnn+li), 0, 0, 0)==FALSO) return FALSO;
			}
			li++;
		}
	}

	if(Ora == Ora_iniz && !Fileoutlati) {
		if((Fileoutlati=_wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("lati.out")), _T("w+")))==NULL) {
			MyPrintf(_T("\nfmgr error lati.out\n"));
			return FALSO;      /* FINE */
		}
		if((Fileoutnodi=_wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("nodi.out")), _T("w+")))==NULL) {
         MyPrintf(_T("\nfmgr error nodi.out\n"));
         return FALSO;      /* FINE */
		}
	}

	if(Fluido == ACQUA) {
		li=1;
		while(li<=Lati) {
			if(*(Latodel+li) == '1' || *(Nl+li) == -1) { li++; continue; }
			if((*(Nl+li)<20)&&(*(Nl+li)>=0)) {  /* SCRIVO SE IL LATO NON E' UN PEZZO SPECIALE */
				qi = *(Q+li) * (float) 1000.0;
				qp = (*(Nl+li)>=0) ? (*(Qd+li) * (float) 1000.0) : (float) 0.0;
				qf = qi - qp;
				p = (float) fabs(*(Carico+(*(Ni+li))) - (*(Carico+(*(Nf+li)))));
				if(p < (-99.0)) p = (float) (-99.0);
				v = (float) (fabs(qi) > fabs(qf)) ? *(Q+li) : ((qi / (float) 1000.0)-(qp / (float) 1000.0));
				v = (float) 4.0 * (float) fabs(v) / (PG * (float) pow((double) *(D+li), (double) 2.0));

            // Scrivo i dati dei lati acqua
            // Codice lato, nodo iniziale, nodo finale, portata iniziale, portata finale, velocità, pressione
				fprintf(Fileoutlati, "%5d,%5d,%5d,%10.3f,%10.3f,%6.3f,%7.3f,%6.1f\n",
				      					*(IDs + li), *(Nnn+(*(Ni+li))), *(Nnn+(*(Nf+li))), qi, qf, v, p, (*(D+li) * 1000.0));
			}
			else if(*(Nl+li)>=20) {  /* scrivo se il lato e` una pompa */
				/* La stazione di pompaggio %d funziona in dissipazione */
				if(((*(Q+li)>0)&&(*(Nnn+*(Ni+li))<Last2)) || ((*(Q+li)<0)&&(*(Nnn+*(Nf+li))<Last2)))
					if(error(11, ((*(Nnn+*(Ni+li))>Last2) ? *(Nnn+*(Nf+li)) : *(Nnn+(*(Ni+li)))), 0, 0, 0)==FALSO) return FALSO;
			else if((*(Q+li) > Pq[10][(*(Nl+li)-20)]) || (*(Q+li) < Pq[1][(*(Nl+li)-20)]))
				/* La stazione di pompaggio %d non puo` erogare %f */
//				if(error(12,((*(Nnn+*(Ni+i))>Last2) ? *(Nnn+*(Nf+i)) : *(Nnn+(*(Ni+i)))),0,(*(Q+i)*1000.0),0)==FALSO) return FALSO;
				if(error(12,((*(Nnn+*(Ni+li))>Last2) ? *(Nnn+*(Nf+li)) : *(Nnn+(*(Ni+li)))),0,(*(Q+li)*1000.0),0)==FALSO) return FALSO;
			}
			++li;
		}
	}
	else {
		li=1;
		while (li<=Lati) {
			if(*(Latodel+li) == '1' || *(Nl+li) == -1) { li++; continue;}
			qi = *(Q+li);
			qp = (*(Nl+li)>=0) ? *(Qd+li) : (float) 0.0;
			qf = qi - qp;
			p = (float) fabs(*(Carico+(*(Ni+li))) - (*(Carico+(*(Nf+li)))));
			if(p < (-99.0)) p = (float) (-99.0);

			/***** Anna: calcolo le portate effettive - inizio *****/
			if (Fluido==GASM) {
				/* pp = (fabs(*(Carico + *(Ni+i))) + fabs(*(Carico + *(Nf+i))))/2.0; */
				/* pp =  fabs(*(Carico + *(Nf+i))); */

				pp = (float) ( fabs(*(Carico + *(Nf+li))) < fabs(*(Carico + *(Ni+li))) ) ?
					 (float) fabs(*(Carico + *(Nf+li))):
					 (float) fabs(*(Carico + *(Ni+li)));

				Dens_e = Dens * (float) 1.293 * pp / (float) 1.034 * (float) 0.9479;  /* 0.554 */
				if(Dens_e == 0) Dens_e = (float) 0.00001; /* stefano 05/01/96 */

				qie = qi*Dens / Dens_e / (float) 0.8625;  /* 0.554 */ /* 0.8625 per Azzolin */
				qfe = qf*Dens / Dens_e / (float) 0.8625;  /* 0.554 */
			}
			else {
				qie = qi;
				qfe = qf;
			}
			v = (fabs(qie) > fabs(qfe)) ? qie / (float) 3600.0 : qfe / (float) 3600.0;
			v = (float) 4.0 * (float) fabs(v) / (float) (PG * pow((double) (*(D+li) / (float) 1000.0), (double) 2.0));
			if(v > 999) v = (float) 999;

         // Scrivo i dati dei lati gas
         // Codice lato, nodo iniziale, nodo finale, portata iniziale, portata finale, velocità, pressione
			fprintf(Fileoutlati, "%5d,%5d,%5d,%9.1f,%9.1f,%6.3f,%7.3f,%6.1f\n",
					               *(IDs + li), *(Nnn+(*(Ni+li))), *(Nnn+(*(Nf+li))), qi, qf, v, p, *(D+li));
			/***** Anna: calcolo le portate effettive - fine *****/
			++li;
		}
	}

	if(Fluido == ACQUA) {
		/*
		**  assegno l'erogazione del nodo di aspirazione delle pompe al nodo
		**  di mandata e nel carico del nodo mandata metto la prevalenza
		*/
		li=1;
		while(li<=Lati) {
			if(*(Latodel+li) == '1') { li++; continue;}
			if(*(Nl+li) > 20) {
				if(*(Nnn+*(Ni+li))<Last2)
					*(Qe+(*(Ni+li))) = *(Qe+(*(Nf+li)));
				else
					*(Qe+(*(Nf+li))) = *(Qe+(*(Ni+li)));
			}
			++li;
		}

		li=1;
		while(li<=Nodi) {
			if((*(Nododel+li) == '1') || ((*(Nnn+li)>Last) && (*(Nnn+li)<Last3))) { li++; continue; }
			qp = (float) 1000.0 * (*(Qe+li));
			p = *(Carico+li) - (*(QuotaZ+li));
			if(p < (float) (-3.0)) p = (float) (-3.0);
			fprintf(Fileoutnodi, "%5d,%8.2f,%6.2f,%6.2f\n", *(Nnn+li), qp, (*(Carico+li)<(-3)) ? -3 : *(Carico+li), p);
			++li;
		}
	}
	else {
		li=1;
		while(li<=Nodi) {
			if((*(Nododel+li) == '1') || ((*(Nnn+li)>Last) && (*(Nnn+li)<Last3))) { li++; continue; }
			qp = *(Qe+li);
			if(Fluido == GASB)
				p = *(Carico+li) - (*(QuotaZ+li) * DENSA * (Dens - (float) 1.0));
			else
				p = *(Carico+li) - (*(QuotaZ+li) * DENSA * (Dens - (float) 1.0) / (float) 10000.0);
			if(p < (-3.0)) p = (float) (-3.0);
			fprintf(Fileoutnodi, "%5d,%9.1f,%9.4f,%9.4f\n", *(Nnn+li), qp, (*(Carico+li)<(-3)) ? -3 : *(Carico+li), p);
			++li;
		}
	}

	if(Calcolo_diametri && (Num_loop_diam>0)) dfromv(); /* determino diametri */

	/* memorizzo l'erogazione dei serbatoi */
	if(Ver24) {
		li=1;
		while(li<=Nserb) {
			if(*(Nododel+li) == '0') *(Tot_erog_serb24+li) = *(Tot_erog_serb24+li) + *(Qe+li);
			++li;
		}
	}

   /*
   ///////////////////// roby test inizio
	li=1;
   FILE *miof;
   miof=fopen("C:\\GS\\GS8\\TEMP\\TEST1.TXT", "w");
   fprintf(miof, "--- inizio ciclo Tot_num_serb\n");
	while(li<=Tot_num_serb) {
      fprintf(miof, "%d, %c\n", *(Nnn+abs(*(Coll_Serb_start+li))),
                             *(Nododel+abs(*(Coll_Serb_start+li))));
		li++;
	}
	li=1;
   fprintf(miof, "\n--- inizio ciclo Nodi\n");
   while(li<=Nodi) {
      fprintf(miof, "%d, %c\n", *(Nnn+li),
		                        *(Nododel+li));
		li++;
	}
   fclose(miof);
   ///////////////////// roby test fine
   */

	li=1;
	while(li<=Tot_num_serb) {
		if(*(Nododel+abs(*(Coll_Serb_start+li))) == '0' && *(Coll_Serb_start+li) > 0)
			*(Coll_Serb_start+li) = -1 * *(Coll_Serb_start+li);
		li++;
	}

   /* controllo se ho verificato tutte le reti */
	li=1;
	while(li<=Tot_num_serb) {
		if(*(Nododel+abs(*(Coll_Serb_start+li))) == '0' && *(Coll_Serb_start+li) > 0)
			*(Coll_Serb_start+li) = -1 * *(Coll_Serb_start+li);
		li++;
	}

	li=1;
	while(li<=Tot_num_serb) {
		if(*(Coll_Serb_start+li) > 0) { /* mancano reti da verificare */
			Nodo_partenza = *(Coll_Serb_start+li);
			ripristino_del();
			Cambio_maglie = 1;          /* cambio topologia */
			--Ora;                      /* ritorno indietro di un'ora */
			return VERO;
		}
		li++;
	}

	/* controllo serbatoi per h24 */
	if(Ver24) contr_serb();   /* controllo serbatoi da chiudere */

	li=1;
	while(li<=Num_prv) {
		if(*(Lato_prv+li) < 0) Chiuse_prv--; /* ripristino lato */
		*(Lato_prv+li) = labs(*(Lato_prv+li)); /* apro */
		li++;
	}

	li=1;
	while(li<=Tot_num_serb) {
		*(Coll_Serb_start+li) = abs(*(Coll_Serb_start+li));
		li++;
	}

	li=1;
	while(li<=Lati) {
		if(*(Nl+li)>20 && *(Qd+li) == -1.0)
      { // roby 21/2/05
			lj=1;
			while(lj<=Tot_num_serb) {
				if(*(Coll_Serb_start+lj) == *(Nf+li)) {
					*(Coll_Serb_start+lj) = -1 * *(Coll_Serb_start+lj);
					break;
				}
				lj++;
			}
      } // roby 21/2/05
		li++;
	}

	Nodo_partenza = Nc_princ; /* ripristino nodo partenza originale */

	/* ripristino nodi e lati eliminati per bumbum */
	ripristino_del();

	if(Calcolo_diametri)
		if(Num_loop_diam > 0) {
			Num_loop_diam--;
			Cambio_maglie = 1;          /* ricalcolo diametri */
			--Ora;                      /* ritorno indietro di un'ora */
		}

	if(Num_pomp > 0) {
		li=1;
		while(li<=Num_pomp) {
			*(Ypprec+li) = Ph[10][li] + (float) fabs(Ph[1][li] - Ph[10][li]) / (float) 2.0;
			li++;
		}
	}

	if(Cambio_maglie && Ver24) { /* possibile attivazione in contr_serb || 28/02/96 stefano - aggiunto Ver24 */
		fseek(Fileoutnodi,Ftell_nodi,SEEK_SET);
		fseek(Fileoutlati,Ftell_lati,SEEK_SET);
		return VERO;
	}

	if(Ver24) {
		fprintf(Fileoutnodi, "%d\n",-1*(Ora+1)); /* segno ripetizione ora */
		fprintf(Fileoutlati, "%d\n",-1*(Ora+1)); /* segno ripetizione ora */
		Ftell_nodi = ftell(Fileoutnodi);
		Ftell_lati = ftell(Fileoutlati);
	}

	if(Ora == Ora_fin) {
		fclose(Fileoutlati); Fileoutlati = NULL;
		fclose(Fileoutnodi); Fileoutnodi = NULL;
	}

	return VERO;
}

int ripristino_del()
{
	long i;

	i=1;
	while(i<=Nodi) {
		*(Nododel+i) = '1';
		i++;
	}
	i=1;
	while(i<=Lati+Num_serb-1) {
		*(Latodel+i) = '1';
		i++;
	}

	return VERO;
}

/* =========================
**  CONTROL_PRV
**  se Cambio_prv = 0 allora riattiva tutte le PRV
*/

int control_prv()
{
	Cambio_maglie = 0;

	if(Cambio_prv) {
		Chiuse_prv++;            /* tolgo un lato */
		Cambio_maglie = 1;       /* cambio topologia */
		--Ora;                   /* ritorno indietro di un'ora */
								 /* "\nRipeto simulazione ora %d." */
                                 /* error24(5, (Ora+1), 0.0, 0.0); */
   }

   return VERO;
}

/* ========================
**  CONTR_SERB
*/
int contr_serb()
{
	long   i;
	float vol_temp = (float) 0.0;

	Cambio_maglie = 0;  /* flag per indicare che se 1 bisogna ricalcolare topol maglie */

	i=1;
	while(i<=Nserb) {
		if(*(Serb24+i) < 0) {   /* serbatoio chiuso */
			/* controllo se lasciarlo ancora chiuso o aprirlo */
			if(*(Carico+i) < (*(Z_serb24+i)+*(H_sfioro24+i))) {
				*(Volume24+i) = (*(Base24+i) * *(H_sfioro24+i)); /* necessario ?? */
				Cambio_maglie = 1;
				*(Serb24+i) = *(Nnn+i);
				/* "\nApro il serbatoio %d.\nPiezometrica inferiore al pelo libero (%.1f/%.1f)." */
				if(error24(3, *(Nnn+i), *(Carico+i), (*(Z_serb24+i)+*(H_sfioro24+i)))==FALSO) return FALSO;
			}
		}
		else {
			/* volume invasato all'inizio simulazione */
			vol_temp = *(Volume24+i) + (*(Q_riemp24+i) + *(Tot_erog_serb24+i)) * (float) 3600.0;
  
			/* controllo su pelo libero > h sfioro */
  
			if(vol_temp > (*(Base24+i) * *(H_sfioro24+i))) {
				/* "\nIl serbatoio %d sfiora %.1f mc."  */
				if(error24(2, *(Nnn+i), vol_temp - (*(Base24+i)* *(H_sfioro24+i)), 0.0)==FALSO) return FALSO;
				/* controllo il segno dell'erogazione (scarto il caso di serbatoio
				pieno perche' entra piu' di quanto esce) e che la quota parte di
				volume sfiorato si maggiore di quanto impostato in Quota_sfioro */
				if((*(Tot_erog_serb24+i) > 0) && (vol_temp - (*(Base24+i)* *(H_sfioro24+i))) >
					((float)*(Quota_sfioro24+i) * 0.01 * (*(Base24+i)* *(H_sfioro24+i)))) {
					*(Serb24+i) = (-1) * *(Serb24+i);  /* chiudo */
					Cambio_maglie = 1;
					/* "\nChiudo il serbatoio %d." */
					if(error24(4, *(Nnn+i), 0.0, 0.0)==FALSO) return FALSO;
				}
			}
		}
		i++;
	}

	if(Cambio_maglie) {
		--Ora; /* ritorno indietro di un'ora */
		/* "\nRipeto simulazione ora %d." */
		if(error24(5, (Ora+1), 0.0, 0.0)==FALSO) return FALSO;
	}
	else {
		i=1;
		while(i<=Nserb) {
			if(*(Serb24+i) > 0) {
				/* volume invasato all'inizio simulazione */
				*(Volume24+i) = *(Volume24+i) + (*(Q_riemp24+i) + *(Tot_erog_serb24+i)) * (float) 3600.0;

				/* controllo su volume < 0 oppure pelo libero > h sfioro */

				if(*(Volume24+i)<0.0) {
					/* riempo con la max capacita` */
					*(Volume24+i) = (*(Base24+i) * *(H_sfioro24+i));
					/* "\nIl serbatoio %d e` vuoto.\nReinizializzo con %.1f mc." */
					if(error24(1, *(Nnn+i), *(Volume24+i), 0.0)==FALSO) return FALSO;
				}
				else if(*(Volume24+i) > (*(Base24+i) * *(H_sfioro24+i)))
				*(Volume24+i) = (*(Base24+i)* *(H_sfioro24+i));
			}
			i++;
		}
	}

	/* azzero erogazione totale serbatoi */
	i=1;
	while(i<=Nserb) {
		*(Tot_erog_serb24+i) = (float) 0.0;
		i++;
	}

	return VERO;
}

/* =========================== 
**  Subroutine CARICO_NP
*/
int carico_np(int stato_lato)
{
	int   i;

	i=1;
	while(i<=Nodi) {
		Np[i] = 0;
		i++;
	}

	i=1;
	while(i<=Lati+Lati_f) {
		/* salto le PRV chiuse */
		if(((cod_prv(i) >= 0)||(cod_vv(i) == 0)) && (*(Latodel+i) != stato_lato)) {
			/* cioe` PRV o VV aperta e lato normale */
			*(Np+(*(Ni+i))) = *(Np+(*(Ni+i))) + 1;
			*(Np+(*(Nf+i))) = *(Np+(*(Nf+i))) + 1;
		}
		i++;
	}

	if(Primo_giro && (!Connessioni)) carico_serb_start();

	return VERO;
}

/* ============================
**  Subroutine SNOLA
*/
int snola(int stato_lato)
{
	int		nodo_reale;
	long	i, j, tot_ind=0, elemento=0;

	/* errore calloc */
	if(stato_lato == '0') {
		tot_ind = (Lati+Lati_f)*2;
		Tot_mem = Tot_mem+((tot_ind+1)*sizeof(int) + (Nodi+2)*sizeof(int));
		if (D_nola == NULL && P_nola == NULL)
		{
			D_nola = (long *) calloc((tot_ind+1) , sizeof(long));
			P_nola = (long *) calloc((Nodi+2) , sizeof(long));
		}
		else
		{
			D_nola = (long *) realloc(D_nola, (tot_ind+1)*sizeof(long));
		}
		if (D_nola == NULL || P_nola == NULL) {
			MyPrintf(_T("\nMEMORIA INSUFFICIENTE in nola"));
			return FALSO;
		}
	}

	/* parte precedente eliminata
	if(D_nola == NULL) {
		tot_ind = (Lati+Lati_f)*2;
		Tot_mem = Tot_mem+((tot_ind+1)*sizeof(int) + (Nodi+2)*sizeof(int));
		D_nola = (int *) calloc((tot_ind+1) , sizeof(int));
		P_nola = (int *) calloc((Nodi+2) , sizeof(int));
		if (D_nola == NULL || P_nola == NULL) {
			MyPrintf(_T("\nMEMORIA INSUFFICIENTE in nola"));
			return FALSO;
		}
	}
	*/

	*(P_nola + 1) = 1;              /* IL PRIMO PUNT. PUNTA AL PRIMO ELEMENTO */

	/* controllare ciclo  (stefano) */
	i=1;
	while(i<=Nodi) {         /* CARICAMENTO VETTORE NODI-LATI NOLA */
		nodo_reale = 0;
		j=1;
		while(j<=Lati+Lati_f) {
			/* salto le PRV o VV chiuse e i lati eliminati temporaneamente */
			if(((*(Ni+j) == i) || (*(Nf+j)==i)) &&
				 (cod_prv(j) >= 0) &&
				 (cod_vv(j) == 0) &&
				 (*(Latodel+j) != stato_lato)) {
				++elemento;
				nodo_reale = 1;
				*(D_nola+elemento) = j;
				*(P_nola + i + 1) = elemento + 1;
			}
			j++;
		}
		if(!nodo_reale) *(P_nola + i + 1) = elemento + 1;
		i++;
	}

	*D_nola = elemento;          /* ULTIMO ELEMENTO DEI DATI DI NOLA */
	*P_nola = Nodi + 1;          /* ULTIMO ELEMENTO DEI PUNTATORI A NOLA */
	*(P_nola + Nodi + 1) = (*D_nola+1); /*NEXT VAL. DA CARICARE IN D_NOLA */

	return VERO;
}

/* ================================ 
**  Subroutine DISTRHN
*/
int distrhn(void)
{
	int     ks;
	long    i, ier, nodi_reali=0, prv;
	float   y=(float) 0.0, w=(float) 0.0, portata;

	i=1;
	while(i<=Nodi) {
		if(*(Nododel+i) == '1')
			Np[i] = 1;
		else {
			Np[i] = 0;
			nodi_reali++;
		}
		i++;
	}

	*(Np+Nc) = 1;
	ier = 1;

	while(ier < nodi_reali) {      /* CAUSA LOOP SE RETI SEPARATE */
		i=1;
		while(i<=Lati+Lati_f) {
			/* controllo PRV chiuse e lati temporaneamente eliminati */
			ks = ((cod_prv(i) >= 0) && /*||*/ (cod_vv(i) == 0) && /*||*/
				(*(Latodel+i) != '1')) ? (*(Np+(*(Ni+i))) + *(Np+(*(Nf+i)))) : 0;
			if(ks == 1) {
				// GAS_PRV: prima di gestione gas_prv era: if(Fluido == ACQUA && *(Nl+i) == -1 && Maglie == 0) {
				if(*(Nl+i) == -1 && Maglie == 0) {
					ier++;
					portata = *(Q+i);
					prv = cod_prv(i);
					if(*(Np+(*(Ni+i))) == 0) {
						if((*(Nnn+ *(Ni+i))>=Last3) && (portata > 0.0))
							if(*(Carico+(*(Nf+i))) <= *(Carico_prv+prv))
								*(Carico+*(Ni+i)) = *(Carico+*(Nf+i));
							else {
								*(Lato_prv+prv) = (-1) * *(Lato_prv+prv); /* chiudo */
								Cambio_prv = 1;  /* attivo flag => ripeto simulazione */
								return VERO;
							}
							else if((*(Nnn+ *(Nf+i))>=Last3) && (portata < 0.0))
								if(*(Carico+(*(Nf+i))) <= *(Carico_prv+prv))
									*(Carico+*(Ni+i)) = *(Carico+*(Nf+i));
								else {
									*(Lato_prv+prv) = (-1) * *(Lato_prv+prv); /* chiudo */
									Cambio_prv = 1;  /* attivo flag => ripeto simulazione */
									return VERO;
								}
							else if(*(Nnn+ *(Ni+i)) >= Last3)
								*(Carico+*(Ni+i)) = *(Carico_prv+prv);
							else if(*(Nnn+ *(Nf+i)) >= Last3) {
								/* Errore non codificato */
								if(error(7, 0, 0, 0, 0)==FALSO) return FALSO;
								return FALSO;
							}
						*(Np+(*(Ni+i))) = 1;
					}
					else {
						if((*(Nnn+ *(Nf+i))>=Last3) && (portata < 0.0))
							if(*(Carico+(*(Ni+i))) <= *(Carico_prv+prv))
								*(Carico+*(Nf+i)) = *(Carico+*(Ni+i));
							else {
								*(Lato_prv+prv) = (-1) * *(Lato_prv+prv); /* chiudo */
								Cambio_prv = 1;  /* attivo flag => ripeto simulazione */
								return VERO;
							}
						else if((*(Nnn+ *(Ni+i))>=Last3) && (portata > 0.0))
							if(*(Carico+(*(Ni+i))) <= *(Carico_prv+prv))
								*(Carico+*(Nf+i)) = *(Carico+*(Ni+i));
							else {
								*(Lato_prv+prv) = (-1) * *(Lato_prv+prv); /* chiudo */
								Cambio_prv = 1;  /* attivo flag => ripeto simulazione */
								return VERO;
							}
						else if(*(Nnn+ *(Nf+i)) >= Last3)
							*(Carico+*(Nf+i)) = *(Carico_prv+prv);
						else if(*(Nnn+ *(Ni+i)) >= Last3) {
							/* Errore non codificato */
							if(error(7, 0, 0, 0, 0)==FALSO) return FALSO;
							return FALSO;
						}
						*(Np+(*(Nf+i))) = 1;
					}
					continue;
				}
				if(Fluido == ACQUA && *(Nl+i) > 20) {	/* Riguarda specificatamente le pompe */
					/* lascio nel nodo fittizio che simula il serbatoio il carico
						corrispondente alla quota di prelievo della pompa.
						Questo per evitare errori nella sub res (stefano 12-03-92)*/
					if((*(Nnn+(*(Ni+i))) > Last2) && (*(Nnn+(*(Ni+i))) < Last3)) {
						*(Carico+(*(Nf+i))) = *(QuotaZ+(*(Nf+i)))+interp(*(Q+i),*(Nl+i)-20);
						*(Carico+(*(Ni+i))) = *(QuotaZ+(*(Ni+i)));
					}
					else {
						*(Carico+(*(Ni+i))) = *(QuotaZ+(*(Ni+i)))+interp(*(Q+i),*(Nl+i)-20);
						*(Carico+(*(Nf+i))) = *(QuotaZ+(*(Nf+i)));
					}
					if(*(Np+(*(Ni+i))) == 0) {
						ier++;
						*(Np+(*(Ni+i))) = 1;
					}
					if(*(Np+(*(Nf+i))) == 0) {
						ier++;
						*(Np+(*(Nf+i))) = 1;
					}
				}
				else {
					res(i, &y, &w);
					ier++;
					if(Fluido != GASM) {
						if(*(Np+(*(Ni+i))) == 0) {
							*(Carico+(*(Ni+i))) = *(Carico+(*(Nf+i))) + (y);
							*(Np+(*(Ni+i))) = 1;
						}
						else {
							*(Carico+(*(Nf+i))) = *(Carico+(*(Ni+i))) - (y);
							*(Np+(*(Nf+i))) = 1;
						}
					}
					else {
						if(*(Np+(*(Ni+i))) == 0) {
							/*
                     if((pow((double)*(Carico+(*(Nf+i))),(double) 2)+(double) y) > 0.0)
								*(Carico+(*(Ni+i))) = sqrt((float) fabs(pow((double)*(Carico+(*(Nf+i))),(double) 2)+(double) y));
							else
								*(Carico+(*(Ni+i))) = (float) -1.0; // 0.0; eliminato - stefano 10/01/96
                     */
                     /* STEFANO 05/05/2014 inizio */
							if(y > 0.0)
								*(Carico+(*(Ni+i))) = sqrt((float) fabs(pow((double)*(Carico+(*(Nf+i))),(double) 2)+(double) y));
							else  {
                        if ((pow((double)*(Carico+(*(Nf+i))),(double) 2)-(double) fabs(y)) < 0) {
								   if(error(27, *(Nnn+(*(Ni+i))), *(Nnn+(*(Nf+i))), 0, 0)==FALSO) return FALSO;
								   return FALSO;
                        }
                        else
								   *(Carico+(*(Ni+i))) = sqrt((float) fabs(pow((double)*(Carico+(*(Nf+i))),(double) 2)-(double) fabs(y)));
                     }
                     /* STEFANO 05/05/2014 fine */
							*(Np+(*(Ni+i))) = 1;
						}
						else {
							/*
                     if((pow((double)*(Carico+(*(Ni+i))),(double) 2)-(double) y) > 0.0)
								*(Carico+(*(Nf+i))) = sqrt((float) fabs(pow((double)*(Carico+(*(Ni+i))),(double) 2)-(double) y));
							else
								*(Carico+(*(Nf+i))) = (float) -1.0; // 0.0; eliminato - stefano 10/01/96
                     */
                     /* STEFANO 05/05/2014 inizio */
							if(y > 0.0) {
                        if ((pow((double)*(Carico+(*(Ni+i))),(double) 2)-(double) y) < 0) {
								   if(error(27, *(Nnn+(*(Ni+i))), *(Nnn+(*(Nf+i))), 0, 0)==FALSO) return FALSO;
								   return FALSO;
                        }
                        else
								   *(Carico+(*(Nf+i))) = sqrt((float) fabs(pow((double)*(Carico+(*(Ni+i))),(double) 2)-(double) y));
                     }
							else
							   *(Carico+(*(Nf+i))) = sqrt((float) fabs(pow((double)*(Carico+(*(Ni+i))),(double) 2)+(double) fabs(y)));
                     /* STEFANO 05/05/2014 fine */
							*(Np+(*(Nf+i))) = 1;
						}
					}
				}
			}
			i++;
		}
		if(ier == 1) {
			/* Errore non codificato */
			if(error(7, 0, 0, 0, 0)==FALSO) return FALSO;
			return FALSO;
		}
	}

	return VERO;
}

/* ===============================
**  Subroutine RIAM
*/
int riam(void)
{
	long	nodo1, nodo2, ier=1, j, iaa, ih, jj;
	float	qtrasp;

	while(ier > 0) {
		ier = 0;
		nodo1=1;
		while(nodo1<=Nodi) {
			if(*(Np+nodo1) != 1 || *(Nododel+nodo1) == '1') {
				nodo1++;
				continue;
			}
			qtrasp = (-1) * *(Qe+nodo1);
			ier = 1;
			iaa = 0;
			j=1;
			while(j<=(*(P_nola+nodo1+1)-*(P_nola+nodo1))) {  /* stefano 28-12-92 */
				ih = ext(D_nola, P_nola, nodo1, j);
				nodo2 = *(Nf+abs(ih));
				if(nodo2 == nodo1) nodo2 = *(Ni+abs(ih));
				if(ih != 0) {      /* CORREZIONE PER LE Q DISTRIBUITE */
					if(ih < 0) {
						ih = (-ih);
						if(nodo2 < nodo1) qtrasp = qtrasp + *(Q+ih) - ((*(Nl+ih)>=0) ? *(Qd+ih) : (float) 0.0);
						else qtrasp = qtrasp - *(Q+ih);
					}
					else {
						if(*(Np+nodo2) != 0) {
							if(iaa > 0) {
								/* dato che Np+nodo1 = 1 non si puo' passare per due nodo2 diversi
								   e quindi avere gia' memorizzato un iaa */
								/* Controllare il lato %d (%d %.0f) */
								if(error(10, iaa, *(Nnn+nodo1), (float) *(Nnn+nodo2), 0)==FALSO) return FALSO;
								return FALSO;
							}
							iaa = ih;
							jj = nodo2;
						}
						else {
							if(nodo2 < nodo1) qtrasp = qtrasp + *(Q+ih) - ((*(Nl+ih)>=0) ? *(Qd+ih) : (float) 0.0);
							else qtrasp = qtrasp - *(Q+ih);
						}
					}          /* fine else ih < 0  */
				}            /* fine if ih != 0   */
				j++;
			}              /* fine for */
			ih = iaa;
			nodo2 = jj;
			if(nodo2 < nodo1) *(Q+ih) = (-qtrasp) + ((*(Nl+ih)>=0) ? *(Qd+ih) : (float) 0.0);
			if(nodo1 < nodo2) *(Q+ih) = qtrasp;
			(*(Np+nodo1))--;
			(*(Np+nodo2))--;
			nodo1++;
		}                /* fine for i<=Nodi */
	}                  /* fine while ier > 0 */
	return VERO;
}

int nnnew(long i, long j)  /* ===================== */
{
	long nlt=1, lat=1, val=0;

	while(lat != 0) {
		lat = abs(ext(D_nola, P_nola, i, nlt));
		if(lat == 0) {
			val = 0;
			break;
		}

		if((*(Ni+lat) == j) || (*(Nf+lat) == j)) {
			val = lat;  // ext(D_nola, P_nola, i, nlt); perchè non tornare l'abs se lo utilizzo come abs?
			break;
		}
		nlt++;
   }

   return(val);
}

/* ============================= 
**  Subroutine RIAM24
*/
int riam24(void)
{
	long	j, i, nodoi, nodof, lato_d_min;
	float	dmin;

	if(riam()==FALSO) return FALSO;       /* assegno portate nei rami */

	j=1;
	while(j<=Maglie) {
		dmin = (float) 1.0e10;              /* taglio lato con D minimo */

		i=1;
		while(i<=*(D_k[j])) {
			nodoi = abs(*(D_k[j]+i));
			if(*(D+nodoi) < dmin) {
				dmin = *(D+nodoi);
				lato_d_min = nodoi;         /* memorizzo il lato con diam minimo */
			}
			i++; 
		}
		nodoi = *(Ni+lato_d_min);
		nodof = *(Nf+lato_d_min);
		*(Np+nodoi) = *(Np+nodoi) - 1;
		*(Np+nodof) = *(Np+nodof) - 1;
		if(riam()==FALSO) return FALSO;	/* assegno portate nella maglia aperta */
		j++;
	}
	return VERO;
}

/* ================================= 
**  Subroutine CICLO_MAGLIE
*/
int ciclo_maglie(void)
{
	long	jj, i, ierl, npart=0, npart_temp, io, jo, rc;

	if(Maglie == 0) {      /* SE NON ESISTONO MAGLIE CALCOLO I CARICHI ED ESCO */
		MyPrintf(td_cad2(1));  /* "\nNessuna maglia" */
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(1));
      #endif
		return VERO;
	}

	npart = codice_prv_taglio();  /* cerco PRV */
	if(npart == 0) {
		i=1;
		while(i<=Nodi) {
			if(*(Nododel+i) == '0') {
   				npart = i;
   				break;
   			}
			i++;
		}
	}
	MyPrintf(td_cad2(6));  /* "\nMaglia:\n" */
   #ifdef EXE_VERS
	   fwprintf(File_esecuzione, td_cad2(6));
   #endif

	jj=1;
	while(jj<=Maglie)
   {
		ierl = bu(npart);
		if(inter(&ierl, jj)==FALSO) return FALSO;

      #ifndef EXE_VERS
         if (acedUsrBrk() == 1) { Aborted = true; return FALSO; }
      #endif

		npart = *(Ni+ierl);  /* nodo partenza per bu */
		*(La+ierl) = (int) (-10000);

		MyPrintf(_T("%3d    \r"),jj);     /* STAMPA DEL NUM. MAGLIA IN ELABORAZIONE */

		io = *(Ni+ierl);
		jo = *(Nf+ierl);

		if(io >= jo) {  /* non dovrebbe mai succedere, tantomeno = */
			MyPrintf(_T("\n\nINTERNAL ERROR %d %d \n\n"),io, jo);   /* da eliminare */
			io = *(Nf+ierl);
			jo = *(Ni+ierl);
			if((io == 0) || (jo == 0)) return FALSO;
		}

		i = 1;
		while(i != 0) {
   			rc = abs(ext(D_nola, P_nola, io, i));   // stefano 01/08/96 generava errore se rc < 0
   			if(rc == 0) break;
   			if( (*(Ni+rc) == jo) || (*(Nf+rc) == jo)) {
   				rc =  w_ext(D_nola, P_nola, io, i, -rc);
   				if(rc) {
   					/* Errore non codificato - contattare AMGA */
   					if(error(17, 1, 0, 0, 0)==FALSO) return FALSO;
   					return FALSO;
   				}
   			}
   			++i;
		}

		i = 1;
		while(i != 0) {
			rc = abs(ext(D_nola, P_nola, jo, i));   // stefano 01/08/96 generava errore se rc < 0
   			if(rc == 0) break;
   			if( (*(Ni+rc) == io) || (*(Nf+rc) == io)) {
				rc =  w_ext(D_nola, P_nola, jo, i, -rc);
				if(rc) {
   					/* Errore non codificato - contattare AMGA */
					if(error(17, 2, 0, 0, 0)==FALSO) return FALSO;
					return FALSO;
		   		}
   			}
   			++i;
		}
		/* elimino il lato io jo */
		*(Np+io) = *(Np+io) - 1;
		*(Np+jo) = *(Np+jo) - 1;

		if(riam()==FALSO) return FALSO;

		npart_temp = npart;
		npart = codice_prv_taglio();  /* cerco fra PRV */
		if(npart == 0) npart = npart_temp;
		jj++;
	}
	if(Ora == Ora_iniz && !(Cambio_maglie)) {
		/* "\n\nMemoria utilizzata  : %ld bytes" */
		MyPrintf(td_cad2(40),Tot_mem);
      #ifdef EXE_VERS
         fwprintf(File_esecuzione, td_cad2(40), Tot_mem);
      #endif
	}

	return VERO;
}

long bu(long npart)    /* =============================================== */
{
	long i, in, iter_bu, ierl;

	i=1;
	while(i<=Nodi) *(No+i++) = 0;
	i=1;
	while(i<=Lati+Lati_f) {
	   if(*(La+i) != -10000) *(La+i)=0;
	   i++;
	}

	*(No+npart) = 1000;
	ierl = 0;
	in = 1000;
	iter_bu = 1;

	while(1) {
		i=1;
		while(i<=Lati+Lati_f) {
			if(*(La+i) == -10000 || *(Latodel+i) == '1') { i++; continue; }
			/* test se lato gia' esaminato - esegue se no */
			if((*(No+(*(Ni+i))) == in) && (*(No+(*(Nf+i))) != in - 1)) {

				/* esegue e punta il lato e il nodo testando
					se si chiude la maglia */

				if(*(No+(*(Nf+i))) == 0) {                 /* NON SI CHIUDE */
					*(La+i) = iter_bu;
					*(No+(*(Nf+i))) = in + 1;
				}
				else {       /* si chiude la maglia */
					if(*(No+(*(Nf+i))) == in) *(La+i) = (-1); /* NUMERO LATI MAGLIA DISPARI */
					else *(La+i) = iter_bu;               /* NUMERO LATI MAGLIA PARI */

					ierl = i;   /* puntatore lato di chiusura e successivo taglio */
					return(ierl);  /* trovata maglia - puntatore lato */
				}
			}

			/* ciclo identico al precedente ma con invertiti ni e nf */
			if((*(No+(*(Nf+i))) == in) && (*(No+(*(Ni+i))) != in - 1)) {
				if(*(No+(*(Ni+i))) == 0) {
					*(La+i) = iter_bu;
					*(No+(*(Ni+i))) = in + 1;
				}
				else {
					if(*(No+(*(Ni+i))) == in) *(La+i) = (-1);
					else *(La+i) = iter_bu;
					ierl = i;
					return(ierl);  /* trovata maglia - puntatore lato */
				}
			}
			i++;
		}

		/* iterazione successiva - maglia non ancora chiusa */
		in++;
		iter_bu++;
	}

   return(ierl);
}

/* =================================
**  Subroutine INTER
*/
int inter(long *lato, long maglia)
{
	long    i, j, nodp, nlma, iniz1, iniz2, iniz3, iniz4, n1, n2,
			nn, nn1, ko, nnma, nodoi, nodof, lato_d_min, latom, ks;
	float   dmin;

	nodp = (*(No+(*(Ni+*lato))) < *(No+(*(Nf+*lato)))) ? *(No+(*(Nf+*lato))) : *(No+(*(Ni+*lato)));
	nlma = (nodp - 1000) * 2;
	iniz3 = nodp - 1000;
	iniz4 = iniz3 + 1;

	if(*(La+*lato) == -1) iniz3 = iniz4;
	if(*(La+*lato) == -1) nlma++;
	nnma = nlma+1;
	iniz2 = nnma / 2 + 1;
	iniz1 = iniz2 - 1;
	if(nnma / 2 * 2 != nnma) iniz1 = iniz2;
	if(*(La+*lato) == -1) {
		n1 = *(Ni+*lato);
		n2 = *(Nf+*lato);
		--iniz3;
		++iniz4;
		nn = *(No+(*(Ni+*lato)));
	}
	else {
		nn = 1000 + *(La+*lato);
		if(*(No+(*(Ni+*lato))) == nn) n1 = *(Ni+*lato);
		else n1 = *(Nf+*lato);
		n2 = n1;
	}
	*(Ivet+iniz1) = n1;             /* PRIMI NODI  TROVATI */
	*(Ivet+iniz2) = n2;
	--iniz1;
	++iniz2;

	do    /* CICLO SUL LATO SINISTRO */
	{
		nn1 = nn - 1000;
		j = 1;
		while(j != 0) {
			i = abs(ext(D_nola, P_nola, n1, j)); // 02/08/96 inserito abs()
			if(i==0)break;
			if(*(La+i) == nn1) {
				if(n1 == *(Ni+i)) {
					n1 = *(Nf+i);
					goto lab2000;
				}
				if(n1 == *(Nf+i)) {
					n1 = *(Ni+i);
					goto lab2000;
				}
			}
			++j;
		}

		lab2000:           /* LATO PUNTATO */
		*(La+i) = (-2);    /* se i = 0 e' lo stesso */

		j = 1;
		while(j != 0) {
			i = abs(ext(D_nola, P_nola, n2, j)); // 02/08/96 inserito abs()
			if(i == 0) break;
			if(*(La+i) == nn1) {
				if(n2 == *(Ni+i)) {
					n2 = *(Nf+i);
					break;
				}
				if(n2 == *(Nf+i)) {
					n2 = *(Ni+i);
					break;
				}
			}
			++j;
		}
		*(La+i) = (-2);    /* se i = 0 e' lo stesso */

		/*
		** AGGIORNAMENTO INDICI E INSERIMENTO NODI TROVATI IN Ivet
		*/

		--nn;
		--iniz3;
		++iniz4;
		*(Ivet+iniz1) = n1;
		*(Ivet+iniz2) = n2;
		--iniz1;
		++iniz2;

		/* test se la maglia si richiude prima del previsto (presenza di rami) */
		/*
		**    test di errore - non ha trovato la maglia o nessun percorso
		*/
		if((n1 != n2) && (nn == 1000)) {
			/* Errore nella costruzione delle maglie; lato %d nodi %d e %.0f */
			if(error(14, i, *(Nnn+n1), (float) *(Nnn+n2), 0)==FALSO) return FALSO;
			return FALSO;
		}
	}
	while(n1 != n2);    /* ?????? stefano  */
	nnma = iniz2 - iniz1 - 1;
	nlma = nnma - 1;

	if(D_k[maglia] == NULL)	{
		D_k[maglia] = (long *) calloc(nlma+1, sizeof(long));
		Tot_mem = (Tot_mem + (nlma+1 * sizeof(long)));
	}
	else {
		D_k[maglia] = (long *) realloc(D_k[maglia], (nlma+1)*sizeof(long));
		Tot_mem = (Tot_mem + (nlma+1 * sizeof(long)));
	}

	if(D_k[maglia] == NULL) {
		MyPrintf(_T("\nMEMORIA INSUFFICIENTE nella gestione delle maglie\n"));
		return FALSO;
	}

	*(D_k[maglia]) = nlma;            /* lati nella maglia */

	dmin = (float) 1.0e5;        /* TAGLIO LATO DELLA MAGLIA TROVATA (ierl) */
	ko=1;
	while (ko <= nnma-1) {
		nodoi = *(Ivet+(ko+iniz1));
		nodof = *(Ivet+(ko+1+iniz1));

		latom = nnnew(nodoi, nodof);
		if(*(D+latom) < dmin) {
			dmin = *(D+latom);
			lato_d_min = latom;  /* memorizzo il lato con diam minimo */
		}
		ks = 1;
		if (nodoi > nodof) ks = -1;

		*(D_k[maglia]+ko) = ks*latom;      /* *(k[ij]+jj) = ks; */
		ko++;
	}
	*lato = lato_d_min;

	return VERO;   /* ritorno lato diametro minimo */
}

/* ================================
**  Subroutine VERIBB
*/
int veribb(void)
{
	int	esco=0;
	long	i, ii, ks, n_lato, iter;
	float	cper, y, w, ba, bm ,xm;
	float	SumQ=0;
   //long  MagliaMaxBi;

	cper = (float) 0.75; /* acceleratore di convergenza (.7 - .8) */
	iter = (int) 0;
   GlobalIter=iter;

	print_label();

	if(Fluido != ACQUA)
	   for(i=1; i<=Num_serb; i++) *(StatoRegolatore+i)=0;	// azzeramento StatoRegolatore

   /*
   /////////////////////// test inizio
   FILE *MyF;
   MyF=fopen("C:\\GS\\GS8\\TEMP\\lati_di_maglia_750.TXT","w");

	i=1;
	while(i<=*(D_k[750]))
   {
		fprintf(MyF, "%ld\n", *(IDs + abs(*(D_k[750]+i)))); // ID del lato

		i++; 
	}
   fclose(MyF);

   MyF=fopen("C:\\GS\\GS8\\TEMP\\TEST.TXT","a");

	fprintf(MyF, "-------------------------------------\n");
	long MyLi=1;
	while(MyLi<=Nodi)
   {
		fprintf(MyF, "%5d;%c\n", *(Nnn+MyLi), *(Nododel+MyLi));
		++MyLi;
	}

   MyLi=1;
	while(MyLi<=Lati+Lati_f)
   {
		fprintf(MyF, "%ld;%5d;%5d;%c\n", 
              *(IDs + MyLi), *(Nnn+*(Ni+(MyLi))), *(Nnn+*(Nf+(MyLi))), *(Latodel+MyLi));
		++MyLi;
	}

   fclose(MyF);
   /////////////////////// test fine
   */

	do
   {
      #ifndef EXE_VERS
   		if (acedUsrBrk() == 1) { Aborted = true; return FALSO; }
      #endif

		bm = (float) 0.0;
		xm = (float) 0.0;
		Ictrp = 0;  // azzero controllo pompe in curva - stefano 30-01-92

		// Scrittura sul file debug_test.csv per analisi dei dati di convergenza nei test...
		if(Debug_on) {		
			long li, reg;
			reg=1;
			while (reg<=Num_serb) {
				li=1;
				while (li<=Lati) {
					if (((*(Ni+li)==reg) || (*(Nf+li)==reg)) && (*(QMaxReg + reg)>0))
						fprintf(FileDEBUG_CSV, "%ld, %ld, %ld, %f , %f | ", reg, *(Nnn + *(Ni+li)), *(Nnn + *(Nf+li)), *(Q+li), *(D+li));
					li++;
				}
				reg++;
			}
			fprintf(FileDEBUG_CSV, "\n");
		}

		//-------------------------------------------------
      // test
      //FILE *fp;

      //fp = fopen("c://gs//gs8//temp//maglie.txt", "a");

      //float *ValoriY;
      //ValoriY = (float *) calloc((Lati + Lati_f + 1), sizeof(float));

		i=1;
		while(i<=Maglie)
      {
			*(B+i) = (float) 0.0;
			*(A+i) = (float) 0.0;

			ii=1;
			while(ii<=*D_k[i]) {        // stefano 24-12-92
				n_lato = abs(*(D_k[i]+ii));
				ks = (int) n_lato / (*(D_k[i]+ii));
				res(n_lato, &y, &w);
            
            //ValoriY[n_lato] = y; // test

				*(B+i) = *(B+i) + (float) ks * (y);
				*(A+i) = *(A+i) + (w);
				ii++;
			}

			ba = (float) fabs(*(B+i));   // memorizzo il massimo sbilanciamento
			if(ba > bm)
         {
            //MagliaMaxBi = i; // test
            bm = ba;
         }

         // test
         //fprintf(fp, "%f;", *(B+i));

			if((*(A+i) == 0) || (*(B+i) == 0))
				*(B+i) = (float) 0.0001;
			else
				*(B+i) = -(*(B+i)) / *(A+i);

			// b(i) e' diventata una portata (mc/s)

			if(Fluido == ACQUA) {
				if(*(B+i) < -1) *(B+i) = (float) (-1);
				if(*(B+i) > 1) *(B+i) = (float) 1;
			}
			else {
				if(*(B+i) < -50000.0) *(B+i) = (float) (-50000.0); // gasb
				if(*(B+i) > 50000.0) *(B+i) = (float) 50000.0;   // gasb 
			}

			ba = (float) fabs(*(B+i));
			if(ba > xm) xm = ba;
			*(B+i) = *(B+i) * cper;

			// applica la portata correttiva *(B+i) ad ogni lato della maglia
         ii=1;
			while(ii<=*D_k[i]) {
				n_lato = abs(*(D_k[i]+ii));
				ks = (int) n_lato / (*(D_k[i]+ii));
				if(ks < 0) *(Q+n_lato) = *(Q+n_lato) - *(B+i);
				else if(ks > 0) *(Q+n_lato) = *(Q+n_lato) + *(B+i);
				ii++;
			}
			i++; 
		}

      // test
      //fprintf(fp, "\n");
      //fclose(fp);

      /*
      // test inizio
      FILE *fp;

      fp = fopen("c://gs//gs8//temp//magliaMax.txt", "a");
		ii=1;
		while(ii<=*D_k[MagliaMaxBi])
      {
			n_lato = abs(*(D_k[MagliaMaxBi]+ii));
        
         fprintf(fp, "%d;%d;%d;%d;%d;%f;%f;%d\n", iter, MagliaMaxBi, n_lato, 
                                            *(Nnn+*(Ni+n_lato)), *(Nnn+*(Nf+n_lato)),
                                            ValoriY[n_lato], ValoriY[n_lato] / *(L+n_lato), *(IDs+n_lato));

			ii++;
		}
      fclose(fp);
      free(ValoriY);
      // test fine
      */

		iter++;
      GlobalIter=iter;

		if(Num_prv > 0)  /* calcolo la perdita di carico per le PRV */
			esco = chiudo_prv();

		switch(Contr_q_h) {
		case CARICO:
			if(Fluido == ACQUA)
				MyPrintf(_T(" %4d %12.5f %3d %3d    \r"), iter, bm, Ictrp, Contr_prv);
			else
				MyPrintf(_T(" %4d %12.5f %3d   \r"), iter, bm, Contr_prv);
			break;
		case PORTATA:
			if(Fluido == ACQUA)
				MyPrintf(_T(" %4d %12.5f %3d %3d    \r"), iter, xm*1000.0, Ictrp, Contr_prv);
			else
				MyPrintf(_T(" %4d %12.5f %3d   \r"), iter, xm, Contr_prv);
			break;
		default:
			if(Fluido == ACQUA)
				MyPrintf(_T(" %4d %12.5f %12.5f %3d %3d    \r"), iter, bm, xm*1000.0, Ictrp, Contr_prv);
			else
				MyPrintf(_T(" %4d %12.5f %12.5f %3d   \r"), iter, bm, xm, Contr_prv);
			break;
		}

      #ifndef EXE_VERS
         if (acedUsrBrk() == 1) return FALSO; // per uscita forzata dall'utente
      #endif

      // test Stefano 21/9/06
      // Stampo in un file valori di PRV
      //Stampa_valori_prv(iter);

		if(esco) return VERO;

	}
	while (continuo(bm,xm));

	return VERO;
}

/* =============================
**  Subroutine CHIUDO_PRV
*/
int chiudo_prv()
{
	long	prv, prv_da_chiudere2=0, prv_da_chiudere1=0, continuo=FALSO;
	float	prec_delta_prv, car_monte, car_valle,
			velocita, min_velocita=(float) 1e8,
			portata_ass, min_portata_ass=(float) 1e8, portata;
	float pp, Dens_e;
   // 09092008 int   VersoQ;

	Contr_prv = 0;

	/* controllo se in questa rete esistono PRV da verificare */
	prv=1;
	while(prv<=Num_prv) {
		if((*(Lato_prv+prv) > 0) && (*(Latodel+labs(*(Lato_prv+prv))) == '0'))
			continuo = VERO;
		prv++;
	}
	if(!continuo) return 0;

	/* memorizzo carichi iterazione precedente */
	prv=1;
	while(prv<=Num_prv) {
		if((*(Lato_prv+prv) > 0) && (*(Latodel+*(Lato_prv+prv)) == '0'))
			if(*(Nnn + *(Ni+*(Lato_prv+prv))) > *(Nnn + *(Nf+*(Lato_prv+prv)))) {
				*(H_monte+prv) = *(Carico+ *(Nf+ *(Lato_prv+prv)));
				*(H_valle+prv) = *(Carico+ *(Ni+ *(Lato_prv+prv)));
			}
			else {
				*(H_monte+prv) = *(Carico+ *(Ni+ *(Lato_prv+prv)));
				*(H_valle+prv) = *(Carico+ *(Nf+ *(Lato_prv+prv)));
			}
		prv++;
	}

	if(distrhn()==FALSO) return FALSO;     /* calcolo carichi */

	prv=1;
	while(prv<=Num_prv) {
		if((*(Lato_prv+prv) > 0) && (*(Latodel+*(Lato_prv+prv)) == '0')) {
  			if(*(Nnn + *(Ni+*(Lato_prv+prv))) > *(Nnn + *(Nf+*(Lato_prv+prv)))) {
				car_monte = *(Carico+ *(Nf+ *(Lato_prv+prv)));
				car_valle = *(Carico+ *(Ni+ *(Lato_prv+prv)));
            // 09092008 // se la Q è positiva vuol dire che la valvola non può regolare (Q che scorre nel verso opposto)
            // 09092008 if (*(Q+ *(Lato_prv+prv)) > 0)
            // 09092008    VersoQ = -1;
            // 09092008 else
            // 09092008    VersoQ = 1;
			}
			else {
				car_monte = *(Carico+ *(Ni+ *(Lato_prv+prv)));
				car_valle = *(Carico+ *(Nf+ *(Lato_prv+prv)));
            // 09092008 // se la Q è negativa vuol dire che la valvola non può regolare (Q che scorre nel verso opposto)
            // 09092008 if (*(Q+ *(Lato_prv+prv)) < 0)
            // 09092008    VersoQ = -1;
            // 09092008 else
            // 09092008    VersoQ = 1;
			}

			prec_delta_prv = *(Delta_prv+prv); // perdita di carico della PRV del ciclo precedente

			// if((car_monte <= *(Carico_prv+prv)) && (car_valle <= *(Carico_prv+prv))) /* carico in rete <nominale valv */
         // carico monte: carico del nodoiniziale che identifica la valvola PRV
         // carico valle: carico del nodofinale che identifica il nodo fittizio inserito.
         // Il nodo fittizio si trova sempre fra la PRV e il nodo che è stato selezionato (NODO_REGO) per indicare la direzione di regolazione

         // se il carico di monte è minore del carico da controllare è inutile che regoli quindi imposto la perdita di carico a zero
			if(car_monte <= *(Carico_prv+prv)) /* carico in rete < nominale valv */
				*(Delta_prv+prv) = (float) 0.0;        /* valvola tutta aperta */
            // 09092008 controllo anche il verso della Q. Se opposto al verso di funzionamento della valvola applico una perdita di carico
            // 09092008 if (VersoQ > 0)
				// 09092008    *(Delta_prv+prv) = (float) 0.0;        /* valvola tutta aperta */
            // 09092008 else
				// 09092008    *(Delta_prv+prv) = (float) prec_delta_prv + PrvAcceleratoreConvergenza * 
            // 09092008                         (fabs(*(Carico_prv+prv) - car_monte) - prec_delta_prv);
         // se il carico di valle è maggiore del carico di monte è inutile che regoli dato che la PRV viene attraversata da una Q in senso inverso e quindi imposto la perdita di carico a zero
         // 08092008 forse. Ripristino la situazione precedente (inserendo lo smorzatore PrvAcceleratoreConvergenza) altrimenti le PRV attraversate da Q opposta vengono chiuse per prime
			else if(car_valle >= car_monte)   // ??????
				   //*(Delta_prv+prv) = (float) 0.0;        /* valvola tutta aperta */
					// 08092008 *(Delta_prv+prv) = (float) fabs(*(Carico_prv+prv) - car_valle);
					*(Delta_prv+prv) = (float) prec_delta_prv + (float) PrvAcceleratoreConvergenza * 
                                  (fabs(*(Carico_prv+prv) - car_valle) - prec_delta_prv);
				else
					// stefano 22082008 *(Delta_prv+prv) = (float) fabs(*(Carico_prv+prv) - car_monte);
               // calcolo la perdita di carico come quella del ciclo precedente + il 75% di quella che sarebbe necessaria
					*(Delta_prv+prv) = (float) prec_delta_prv + (float) PrvAcceleratoreConvergenza * 
                                  (fabs(*(Carico_prv+prv) - car_monte) - prec_delta_prv);

			/* controllo se tutto il sistema e` stabile */
			/*if((fabs(prec_delta_prv - *(Delta_prv+prv)) > 0.2) &&*/
			if( ((fabs(*(H_monte+prv)-car_monte) > 0.2) && 
              (fabs(*(H_valle+prv)-car_valle) > 0.2) 
             ) &&
             *(D+ *(Lato_prv+prv)) > 0.001
             || GlobalIter <= MaxGlobalIter
           )
				Contr_prv = 1;              /* non stabile */
		}
		prv++;
	}

	/* se tutte le PRV sono in funzione allora controllo
       la velocita` o la portata per chiudere.
       DEVO SCEGLIERE QUELLA PIU` PROSSIMA ALLA CHIUSURA =>
       min vel or min q or verso q inverso */

	if(Contr_prv == 0) {
		prv=1;
		while(prv<=Num_prv) {
			if((*(Lato_prv+prv) > 0) && (*(Latodel+*(Lato_prv+prv)) == '0')) {  // se la prv è attiva
				portata = *(Q+*(Lato_prv+prv));
				if(((*(Nnn+ *(Ni+ *(Lato_prv+prv)))>=Last3) && (portata > 0.0)) ||
					((*(Nnn+ *(Nf+ *(Lato_prv+prv)))>=Last3) &&
					(portata < 0.0))) {  // se la portata è opposta alla drirezione di controllo della prv
					if(((*(Nnn+ *(Ni+ *(Lato_prv+prv)))>=Last3) &&
						(fabs(*(H_valle+prv) - *(Carico+*(Ni+*(Lato_prv+prv)))) < 1.0)
						) ||
						((*(Nnn+ *(Nf+ *(Lato_prv+prv)))>=Last3) &&
						(fabs(*(H_valle+prv) - *(Carico+*(Nf+*(Lato_prv+prv)))) < 1.0)
						))
						prv_da_chiudere1 = prv;
				}
				if(Fluido == ACQUA)
               velocita = (float) 4.0 * (float) fabs(*(Q+*(Lato_prv+prv))) / (PG *
						      (float) pow((double) *(D+*(Lato_prv+prv)), (double) 2.0));
            else {
			      if (Fluido==GASM) {
				      pp = (float) ( fabs(*(Carico + *(Nf+*(Lato_prv+prv)))) < fabs(*(Carico + *(Ni+*(Lato_prv+prv)))) ) ?
					      (float) fabs(*(Carico + *(Nf+*(Lato_prv+prv)))):
					      (float) fabs(*(Carico + *(Ni+*(Lato_prv+prv))));

				      Dens_e = Dens * (float) 1.293 * pp / (float) 1.034 * (float) 0.9479;  /* 0.554 */
				      if(Dens_e == 0) Dens_e = (float) 0.00001; /* stefano 05/01/96 */
				      velocita = *(Q+*(Lato_prv+prv)) * Dens / Dens_e / (float) 0.8625;  /* 0.554 */ /* 0.8625 per Azzolin */
			      }
			      else {
				      velocita = *(Q+*(Lato_prv+prv));
			      }
			      velocita = velocita / (float) 3600.0;
			      velocita = (float) 4.0 * (float) fabs(velocita) / (float) (PG * pow((double) (*(D+*(Lato_prv+prv)) / (float) 1000.0), (double) 2.0));
            }
	         portata_ass = (float) fabs(portata);
		       /* controllo velocita` e valore assoluto portata */
             /* stefano 08092008 tolgo il caso di Q=0 cioè non chiudo se Q=0 (&& portata_ass > 0) */
			    if(((velocita < Min_vel_prv) ||  /* !!!!!!!!!!!! */
				   (portata_ass < Min_q_prv))
               && portata_ass > 0) {
					if((velocita < min_velocita) || (portata_ass < min_portata_ass)) {
						prv_da_chiudere2 = prv;
						min_velocita = velocita;
						min_portata_ass = portata_ass;
					}
				}
			}
			prv++;
		}
	}

	/* do piu` importanza alla portata in verso opposto */
	if(prv_da_chiudere1 > 0) {
		/* "Chiudo PRV %d\r" */
		if(Ver24) {
			if(error24(6,
					(*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere1)))>=Last3) ?
					*(Nnn+ *(Nf+ *(Lato_prv+prv_da_chiudere1))) :
					*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere1))),
					0.0, 0.0)==FALSO) return FALSO;
		}
		else {
			MyPrintf(td_cad2(47), // "                                                   Chiudo PRV %d\r"
				(*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere1)))>=Last3) ?
				*(Nnn+ *(Nf+ *(Lato_prv+prv_da_chiudere1))) :
				*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere1))));
         #ifdef EXE_VERS
			   fwprintf(File_esecuzione, td_cad2(47),
      		         (*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere1)))>=Last3) ?
				         *(Nnn+ *(Nf+ *(Lato_prv+prv_da_chiudere1))) :
				         *(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere1))));
         #endif
		}


		*(Lato_prv+prv_da_chiudere1) = (-1) * *(Lato_prv+prv_da_chiudere1); /* chiudo */
		Cambio_prv = 1;  /* attivo flag => ripeto simulazione */
		return 1;
	}

	if(prv_da_chiudere2 > 0) {
		if(Ver24) {
			if(error24(6,
					(*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere2)))>=Last3) ?
					*(Nnn+ *(Nf+ *(Lato_prv+prv_da_chiudere2))) :
					*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere2))),
					0.0, 0.0)==FALSO) return FALSO;
		}
		else {
			MyPrintf(td_cad2(47), // "                                                   Chiudo PRV %d\r"
				   (*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere2)))>=Last3) ?
				   *(Nnn+ *(Nf+ *(Lato_prv+prv_da_chiudere2))) :
				   *(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere2))));
         #ifdef EXE_VERS
			   fwprintf(File_esecuzione, td_cad2(47),
			            (*(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere2)))>=Last3) ?
				         *(Nnn+ *(Nf+ *(Lato_prv+prv_da_chiudere2))) :
				         *(Nnn+ *(Ni+ *(Lato_prv+prv_da_chiudere2))));
         #endif
		}

		*(Lato_prv+prv_da_chiudere2) = (-1) *
									 *(Lato_prv+prv_da_chiudere2); /* chiudo */
		Cambio_prv = 1;  /* attivo flag => ripeto simulazione */
		return 1;
	}

	return 0;
}

/* =============================
**  Subroutine CONTINUO
*/
int continuo(double car, double port)
{
	int esco=0;	/* Tutto OK, fine */

	if((Contr_q_h == CARICO) && (car > Err_h))
		esco = 1;
	else if((Contr_q_h == PORTATA) && (port > Err_q))
		esco = 1;
	else if((Contr_q_h == CARICO_E_PORTATA) && ((port > Err_q) || (car > Err_h)))
		esco = 1;
	else if((Contr_q_h == CARICO_O_PORTATA) && (port > Err_q) && (car > Err_h))
		esco = 1;
	/* GAS_PRV: prima...: if((Fluido == ACQUA) && (Contr_prv == 1)) */
	if(Contr_prv == 1) esco = 1;
	if((Fluido == ACQUA) && (Ictrp == 1) && (Cstp == 1)) esco = 1;

	if(Fluido != ACQUA && QMaxRegTot > 0 && esco == 0)
      // controllo se i regolatori con limitatore di portata sono tutti sotto soglia
      // il controllo lo faccio solo a rete bilanciata (esco = 0)
		if (ChkRegolatori() != 1)
		//if (ChkRegolatori() == 1)
			//esco=0;	 // così esco dal ciclo in veribb
		//else
			esco=1;	 // continuo a calcolare

	return esco;
}


/* ================================= 
**  Subroutine PRINT_LABEL
**  stampa intestazione di veribb
*/
int print_label()
{
	/*if(Cambio_maglie) {
	  MyPrintf(_T("\n"));
	  return VERO;
	}*/

	switch(Contr_q_h) {
	case CARICO:
		/* Controllo sul carico (%f m mm kg/cmq)*/
		MyPrintf(td_cad2(42), Err_h);
		/* "\n\n Passo     Carico" */
		MyPrintf(td_cad2(41));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(42), Err_h);
   		fwprintf(File_esecuzione, td_cad2(41));
      #endif
		break;
	case PORTATA:
		/* Controllo sulla portata (%f l/s mc/h)*/
		MyPrintf(td_cad2(44), (Fluido == ACQUA) ? Err_q*1000.0 : Err_q);
		/* "\n\n Passo     Portata" */
		MyPrintf(td_cad2(43));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(44), (Fluido == ACQUA) ? Err_q*1000.0 : Err_q);
   		fwprintf(File_esecuzione, td_cad2(43));
      #endif
		break;
	case CARICO_E_PORTATA:
		/* Controllo su carico (%f m mm kg/cmq) and portata (%f l/s mc/h)*/
		MyPrintf(td_cad2(46),Err_h, _T("and"), (Fluido == ACQUA) ? Err_q*1000.0 : Err_q);
		/* "\n\n Passo     Carico    Portata" */
		MyPrintf(td_cad2(45));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(46), Err_h, _T("and"), (Fluido == ACQUA) ? Err_q*1000.0 : Err_q);
	      fwprintf(File_esecuzione, td_cad2(45));
      #endif
		break;
	case CARICO_O_PORTATA:
		/* Controllo su carico (%f m mm kg/cmq) or portata (%f l/s mc/h)*/
		MyPrintf(td_cad2(46),Err_h, _T("or"), (Fluido == ACQUA) ? Err_q*1000.0 : Err_q);
		/* "\n\n Passo     Carico    Portata" */
		MyPrintf(td_cad2(45));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(46), Err_h, _T("or"), (Fluido == ACQUA) ? Err_q*1000.0 : Err_q);
	   	fwprintf(File_esecuzione, td_cad2(45));
      #endif
		break;
	}
	if(Fluido == ACQUA)
   {
		MyPrintf(_T("  P  PRV\n\n"));  // switch controllo pompe e prv
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, _T("  P  PRV\n\n"));
      #endif
	}
	else
   {
		MyPrintf(_T("  PRV\n\n"));	// GAS_PRV
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, _T("  PRV\n\n"));
      #endif
	}

	return VERO;
}

/* ================================================================== */
/*
**  MEM_ARRAY : Allocazione dinamica della memoria per gli array usati
**  Il codice di ritorno indica in che passo e` stata interrotta
**  l'allocazione della memoria oppure 0 per SUCCESS
**
**  mem : Indica i bytes di memoria allocata.
*/
int mem_array(void)
{
	long    i;
	FILE    *fp;
	TCHAR   ftempinp[150];
	float   dd;

	if(Connessioni) {
		Handle = (char **) calloc((Lati+1+EXTRA), sizeof(char *));
		i=1;
		while(i<(Lati+1+EXTRA)) {
			Handle[i] = (char *) calloc(MAX_LEN_HANDLE, sizeof(char));
			i++;
		}
		if(Handle == NULL) return(1);
	}
	else {
		A = (float *) calloc((Maglie+1+EXTRA) , sizeof(float));
		if(A == NULL) return(2);
		Tot_mem = Tot_mem + ((Maglie+1+EXTRA) * sizeof(float));

		B = (float *) calloc((Maglie+1+EXTRA) , sizeof(float));
		if(B == NULL) return(3);
		Tot_mem = (Tot_mem + ((Maglie+1+EXTRA) * sizeof(float)));

		/* aggiunta per errore calloc */
		D_k = (long **) calloc((Maglie+1+EXTRA) , sizeof(long *));
		if(D_k == NULL) return(3);
		Tot_mem = Tot_mem + ((Maglie+1+EXTRA) * sizeof(long *));
		Maglie_max = Maglie;

		Qe = (float *) calloc((Nodi+1+EXTRA) , sizeof(float));
		if(Qe == NULL) return(4);
		Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(float)));

		Carico = (float *) calloc((Nodi+1+EXTRA) , sizeof(float));
		if(Carico == NULL) return(5);
		Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(float)));

		QuotaZ = (float *) calloc((Nodi+1+EXTRA) , sizeof(float));
		if(QuotaZ == NULL) return(6);
		Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(float)));

		No = (long *) calloc((Nodi+1+EXTRA) , sizeof(long));
		if(No == NULL) return(8);
		Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(long)));

		Ivet = (long *) calloc((Nodi+1+EXTRA) , sizeof(long));
		if(Ivet == NULL) return(9);
		Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(long)));

		Carico_serb = (float *) calloc((Num_serb+EXTRA), sizeof(float));
		if(Carico_serb == NULL) return(13);
		Tot_mem = (Tot_mem + ((Num_serb+EXTRA) * sizeof(float)));

		if(Num_prv > 0) {
			Lato_prv = (long *) calloc((Num_prv+1+EXTRA) , sizeof(long));
			if(Lato_prv == NULL) return(14);
			Tot_mem = (Tot_mem + ((Num_prv+1+EXTRA) * sizeof(long)));

			Carico_prv = (float *) calloc((Num_prv+1+EXTRA) , sizeof(float));
			if(Carico_prv == NULL) return(15);
			Tot_mem = (Tot_mem + ((Num_prv+1+EXTRA) * sizeof(float)));

			D_prv = (float *) calloc((Num_prv+1+EXTRA) , sizeof(float));
			if(D_prv == NULL) return(15);
			Tot_mem = (Tot_mem + ((Num_prv+1+EXTRA) * sizeof(float)));

			H_monte = (float *) calloc((Num_prv+1+EXTRA) , sizeof(float));
			if(H_monte == NULL) return(16);
			Tot_mem = (Tot_mem + ((Num_prv+1+EXTRA) * sizeof(float)));

			H_valle = (float *) calloc((Num_prv+1+EXTRA) , sizeof(float));
			if(H_valle == NULL) return(17);
			Tot_mem = (Tot_mem + ((Num_prv+1+EXTRA) * sizeof(float)));

			Delta_prv = (float *) calloc((Num_prv+1+EXTRA) , sizeof(float));
			if(Delta_prv == NULL) return(18);
			Tot_mem = (Tot_mem + ((Num_prv+1+EXTRA) * sizeof(float)));
		}

      // Alloco vettore per memorizzare gli identificatori GEOsim dei lati
		IDs = (long *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(long));
		if(IDs == NULL) return(20);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(long)));

		L = (float *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(float));
		if(L == NULL) return(19);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(float)));

		D = (float *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(float));
		if(D == NULL) return(20);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(float)));

		/* NS - 01/03/06: Salviamo i diametri originali dei lati */
		Dorig = (float *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(float));
		if(Dorig == NULL) return(21);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(float)));

		Q = (float *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(float));
		if(Q == NULL) return(22);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(float)));

		Qd = (float *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(float));
		if(Qd == NULL) return(23);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(float)));
		
		/* NS - 6/12/05 */
		QMaxReg = (float *) calloc((Num_serb+1) , sizeof(float));
		if(QMaxReg == NULL) return(24);
		Tot_mem = (Tot_mem + ((Num_serb+1) * sizeof(float)));

		/* stefano 11/03/2010 */
   	if(Fluido != ACQUA) {
         StatoRegolatore = (int *) calloc((Num_serb+1) , sizeof(int));
		   if(StatoRegolatore == NULL) return(31);
		   Tot_mem = (Tot_mem + ((Num_serb+1) * sizeof(int)));
      }

		if(Calcolo_diametri) {
			Dtocalc = (char *) calloc((Lati+Lati_f+1+EXTRA), sizeof(char));
			if(Dtocalc == NULL) return(25);
			Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(char)));

			if((fp=_wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("diam.inp")), _T("r")))!=NULL) {
				fscanf(fp,"%d",&Diametri);
				Diamcomm = (float *) calloc((Diametri+1+EXTRA) , sizeof(float));
				if(Diamcomm == NULL) return(25);
				Tot_mem = (Tot_mem + ((Diametri+1+EXTRA) * sizeof(float)));

				i=1;
				while(i<=Diametri) {
				   fscanf(fp,"%f",&dd);
				   *(Diamcomm+i) = dd;
				   i++;
				}

				fclose(fp);

				if(Fluido == ACQUA) {
					i=1;
					while(i<=Diametri) {
						*(Diamcomm+i) = (*(Diamcomm+i)) * (float) 0.001;     /* diam in m */
						i++;
					}
				}
			}
			else
				Diametri = 0; /* calcola diametri reali */
		}

		Scabr = (float *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(float));
		if(Scabr == NULL) return(21);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(float)));

		La = (long *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(long));
		if(La == NULL) return(29);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(long)));

		Nl = (int *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(long));
		if(Nl == NULL) return(30);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(long)));

		if((Fluido == ACQUA)&&(Num_pomp > 0)) {

			Ypprec = (float *) calloc((Num_pomp+1+EXTRA) , sizeof(float));
			if(Ypprec == NULL) return(31);
			Tot_mem = (Tot_mem + ((Num_pomp+1+EXTRA) * sizeof(float)));

			i=0;
			while(i<11) {
				Pq[i] = (float *) calloc((Num_pomp+1) , sizeof(float));
				if(Pq == NULL) return(32);
				Tot_mem = (Tot_mem + ((Num_pomp+1) * sizeof(float)));

				Ph[i] = (float *) calloc((Num_pomp+1), sizeof(float));
				if(Ph == NULL) return(33);
				Tot_mem = (Tot_mem + ((Num_pomp+1) * sizeof(float)));
				i++;
			}
		}
	}

	Nnn = (long *) calloc((Nodi+1+EXTRA) , sizeof(long));
	if (Nnn == NULL) return(7);
	Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(long)));

	/* modificato con nodi per variazione precedente (ciclo_maglie) */
	Np = (long *) calloc((Nodi+1+EXTRA) , sizeof(long));
	if(Np == NULL) return(10);
	Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(long)));

	Nododel = (char *) calloc((Nodi+1+EXTRA), sizeof(char));
	if(Nododel == NULL) return(11);
	Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(char)));

	Nodo_start = (char *) calloc((Nodi+1+EXTRA), sizeof(char));
	if(Nodo_start == NULL) return(12);
	Tot_mem = (Tot_mem + ((Nodi+1+EXTRA) * sizeof(char)));

	Latodel = (char *) calloc((Lati+Lati_f+1+EXTRA), sizeof(char));
	if(Latodel == NULL) return(24);
	Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(char)));

	#ifdef DOS
		Ni = (long *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(long));
		if(Ni == NULL) return(26);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(long)));
	#else
		Ni = (long *) calloc((Lati+Lati_f+1000) , sizeof(long));
		if(Ni == NULL) return(27);
		Tot_mem = (Tot_mem + ((Lati+Lati_f+1000) * sizeof(long)));
	#endif

	Nf = (long *) calloc((Lati+Lati_f+1+EXTRA) , sizeof(long));
	if(Nf == NULL) return(28);
	Tot_mem = (Tot_mem + ((Lati+Lati_f+1+EXTRA) * sizeof(long)));

	if(Ver24) {
		Serb24 = (long *) calloc((Nserb+1) , sizeof(long));          /* numero serb */
		H_sfioro24 = (float *) calloc((Nserb+1) , sizeof(float));  /* quota sfioro (da Z_serb) */
		Tot_erog_serb24 = (float *) calloc((Nserb+1) , sizeof(float));  /* erogazione totale serbatoio */
		Volume24 = (float *) calloc((Nserb+1) , sizeof(float));    /* volume serb */
		Q_riemp24 = (float *) calloc((Nserb+1) , sizeof(float));   /* portata riempimento (o svuotam) serb */
		Z_serb24 = (float *) calloc((Nserb+1) , sizeof(float));    /* quota slm fondo serb */
		Base24 = (float *) calloc((Nserb+1) , sizeof(float));      /* area base serb */
		Quota_sfioro24 = (long *) calloc((Nserb+1) , sizeof(long));  /* percentuale per attivare chiusura serb (0 non attivo) */
	}

	return (0);             /* TUTTO OK */
}

int error(int line, int int1, int int2, double float1, double float2)
{
	TCHAR ftempinp[150];

	if(File_errore == NULL)
		if((File_errore=_wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("errore.out")), _T("w+")))==NULL) 
      {
			MyPrintf(_T("\nfmgr error errore.out\n"));
			return FALSO;    /* FINE */
		}

	if(line > 0)
   {
		fprintf(File_errore, "\nError N. %2d",line);

		/* "\nErrore N. %d\n>>> controllare file errore.out <<<\n\n" */
		MyPrintf(td_cad2(27),line);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(27), line);
      #endif
	}

	switch (line) {
	case -1:
		MyPrintf(td_cad2(54),int1, float1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(54), int1, float1);
      #endif
		fwprintf(File_errore, td_cad2(54), int1, float1);
		break;
	case 2:
		MyPrintf(td_cad2(16),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(16), int1);
      #endif
		fwprintf(File_errore, td_cad2(16),int1);
		break;
	case 3:
		MyPrintf(td_cad2(17));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(17));
      #endif
		fwprintf(File_errore, td_cad2(17));
		break;
	case 6:
		MyPrintf(td_cad2(19),int1,int2);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(19), int1, int2);
      #endif
		fwprintf(File_errore, td_cad2(19),int1,int2);
		break;
	case 7:
		MyPrintf(td_cad2(22));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(22));
      #endif
		fwprintf(File_errore, td_cad2(22));
		break;
	case 10:
		MyPrintf(td_cad2(23),int1,int2,float1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(23), int1, int2, float1);
      #endif
		fwprintf(File_errore, td_cad2(23),int1,int2,float1);
		break;
	case 11:
		MyPrintf(td_cad2(20),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(20), int1);
      #endif
		fwprintf(File_errore, td_cad2(20),int1);
		break;
	case 12:
		MyPrintf(td_cad2(21),int1,float1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(21), int1, float1);
      #endif
		fwprintf(File_errore, td_cad2(21),int1,float1);
		break;
	case 14:
		MyPrintf(td_cad2(25),int1,int2,float1);
      #ifdef EXE_VERS
		   fwprintf(File_esecuzione, td_cad2(25), int1, int2, float1);
      #endif
		fwprintf(File_errore, td_cad2(25),int1,int2,float1);
		break;
	case 16:
		MyPrintf(td_cad2(26));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(26));
      #endif
		fwprintf(File_errore, td_cad2(26));
		break;
	case 17:
		MyPrintf(td_cad2(22));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(22));
      #endif
		fwprintf(File_errore, td_cad2(22));
		break;
	case 18:
		MyPrintf(td_cad2(28),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(28), int1);
      #endif
		fwprintf(File_errore, td_cad2(28),int1);
		break;
	case 19:
		MyPrintf(td_cad2(29),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(29), int1);
      #endif
		fwprintf(File_errore, td_cad2(29),int1);
		break;
	case 20:
		MyPrintf(td_cad2(30),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(30), int1);
      #endif
		fwprintf(File_errore, td_cad2(30),int1);
		break;
	case 21:
		MyPrintf(td_cad2(48),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(48), int1);
      #endif
		fwprintf(File_errore, td_cad2(48),int1);
		break;
	case 22:
		MyPrintf(td_cad2(52),int1,int2,float1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(52), int1, int2, float1);
      #endif
		fwprintf(File_errore, td_cad2(52),int1,int2,float1);
		break;
	case 23:
		MyPrintf(td_cad2(53),int1,int2);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(53), int1, int2);
      #endif
		fwprintf(File_errore, td_cad2(53),int1,int2);
		break;
	case 24:
		/* Il nodo %d ha carico = 0; \nATTENZIONE i risultati della simulazione possono non essere attendibili! */
		MyPrintf(td_cad2(55),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(55), int1);
      #endif
		fwprintf(File_errore, td_cad2(55),int1);
		break;
	case 25:
		/* L'asola e' composta dal lato fra due serbatoi\ne il corrispondente lato fittizio */
		MyPrintf(td_cad2(56));
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(56));
      #endif
		fwprintf(File_errore, td_cad2(56));
		break;
	case 26:
		//\nLimitazione imposta ai regolatori: %9.2f minore del consumo: %9.2f \n", qtot);
		MyPrintf(td_cad2(58), float1, float2);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(58), float1, float2);
      #endif
		fwprintf(File_errore, td_cad2(58), float1, float2);
		break;
      /* STEFANO 05/05/2014 inizio */
	case 27:
		//\nPerdita di carico troppo elevata nel lato %d - %d. Legge resistenza non applicabile.\n"
		MyPrintf(td_cad2(61), int1, int2);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(58), float1, float2);
      #endif
		fwprintf(File_errore, td_cad2(61), int1, int2);
		break;
      /* STEFANO 05/05/2014 fine */
	default:
		MyPrintf(_T(" "));
		break;
	}

	return VERO;
}

/* ================================================================== */
/*
**  ERROR24 - Gestione dell'errore per h24
**  line            : linea del file messaggi errori.dat
**  int1            : parametro 2 INTERO
**  float1          : parametro 3 FLOATING
*/
int error24(int line, int int1, double float1, double float2)
{
	TCHAR ftempinp[150];

	if(File_errore == NULL)
	if((File_errore=_wfopen(wcscat(wcscpy(ftempinp,Ftemp), _T("errore.out")), _T("w+")))==NULL) 
   {
		MyPrintf(_T("\nfmgr error errore.out\n"));
		return FALSO;    /* FINE */
	}

	if(Ora_errore == 0) {     /* Scrivo solo 1 volta */
		fwprintf(File_errore,td_cad2(12),Ora); /* "\nSimulazione ora: %3d\n" */
		Ora_errore = 1;
	}

	switch (line) {
	case 1:
		MyPrintf(td_cad2(13),int1,float1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(13), int1, float1);
      #endif
		fwprintf(File_errore, td_cad2(13),int1,float1);
		break;
	case 2:
		MyPrintf(td_cad2(15),int1,float1);
      #ifdef EXE_VERS
		   fwprintf(File_esecuzione, td_cad2(15), int1, float1);
      #endif
		fwprintf(File_errore, td_cad2(15),int1,float1);
		break;
	case 3:
		MyPrintf(td_cad2(31),int1,float1,float2);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(31), int1, float1, float2);
      #endif
		fwprintf(File_errore, td_cad2(31),int1,float1,float2);
		break;
	case 4:
		MyPrintf(td_cad2(32),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(32), int1);
      #endif
		fwprintf(File_errore, td_cad2(32),int1);
		break;
	case 5:
		MyPrintf(td_cad2(33),int1);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(33), int1);
      #endif
		fwprintf(File_errore, td_cad2(33),int1);
		break;
	case 6:
		MyPrintf(td_cad2(50),int1); // "                                                   Chiudo PRV %d\r"
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(50), int1);
      #endif
		fwprintf(File_errore, _T("\n\n"));
		fwprintf(File_errore, td_cad2(50),int1);
		break;
	default:
		MyPrintf(_T(" "));
		break;
	}

	return VERO;
}

/* ================================================================== */
/*
**
**  Subroutine BUMBUM
**
*/
int bumbum()
{
	long  nodo_segno, i;

	/* non considero lati fittizi */
	i=Lati+1;
	while(i<=Lati+Lati_f) {
		(*(Np+*(Ni+i)))--;
		(*(Np+*(Nf+i)))--;
		i++;
	}

	/* azzero vettore Nodo_start */
	i=1;
	while(i<=Nodi) {
		*(Nodo_start+i) = '0';
		i++;
	}
	/* in nodo_segno il numero di nodo che devo iniziare a seguire */
	nodo_segno = Nodo_partenza;

	MyPrintf(td_cad2(49));  /*"\nControllo topologia  |" */
   #ifdef EXE_VERS
	   fwprintf(File_esecuzione, td_cad2(49));
   #endif

	while(nodo_segno != 0) {
		orologio();
		while(nodo_segno != 0) {
			if(*(Np+nodo_segno) < 2)
				/* segno il nodo per indicarlo attivo */
				*(Nododel+nodo_segno) = '0';

			/* diminuisco le partizioni del nodo */
			(*(Np+nodo_segno))--;
			/* 0 = nodo non piu` da seguire, 1 si */
			*(Nodo_start+nodo_segno) = (*(Np+nodo_segno) == 0) ? '0' : '1';
			/* cerco nodo collegato a nodo_segno */
			nodo_segno = new_nodo(nodo_segno);
		}
		/* cerco nodo lasciato marcato */
		nodo_segno = ripesco_nodo();
	}

	MyPrintf(td_cad2(51));  /*"\b\bcompletato." */
   #ifdef EXE_VERS
   	fwprintf(File_esecuzione, td_cad2(51));
   #endif

	return riallineo();
}

/* ===================================================== */
/*
**  funzione new_nodo
**  restituisce nodo in o fi di lato che e` diventato ramo e quindi
**  viene eliminato da nola.
**  A differenza della versione di nnnew da verifica, il parametro
**  passato e` uno solo (numero nodo) e restituisce un nodo invece
**  che un lato
*/

long new_nodo(long i)
{
	long nlt = 1, lat = 1, val = 0;

	while(lat != 0) {
		lat = ext(D_nola, P_nola, i, nlt);
		if(lat == 0) {
			val = 0;
			break;
		}

		if((lat > 0) && (*(Latodel+lat) == '1') && (cod_prv(lat) >= 0) && (cod_vv(lat) == 0)) {
			if((*(Ni+lat) == i) && (*(Np+(*(Nf+lat))) > 0)) {
				*(Latodel+lat) = '0';   /* lato concatenato */
				val = *(Nf+lat);
				if(sono_serb(val)) {
					*(Nododel+val) = '0';
					val = 0;
				}
				break;
			}
			if((*(Nf+lat) == i) && (*(Np+(*(Ni+lat))) > 0)) {
				*(Latodel+lat) = '0';   /* lato concatenato */
				val = *(Ni+lat);
				if(sono_serb(val)) {
					*(Nododel+val) = '0';
					val = 0;
				}
				break;
			}
		}
		nlt++;
	}
	return(val);
}

/* =============================================================== */
int ripesco_nodo()
{
	long i, val = 0;

	i=1;
	while(i<=Nodi) {
		if((*(Nodo_start+i) != '0') && (!sono_serb(i))) {
			val = i;
			break;
		}
		i++;
	}

	return(val);
}

int riallineo()
{
	long	nodi_reali=0, lati_reali=0, num_prv_reali=0, nstp_reali=0,
			temp_Lati_f=0, i, j, jj3;
	float	qtot=(float) 0.0;

	if(Connessioni) return (0);

	*(Qe+Nc) = (float) 0.0;
	Nc = 0;         /* riazzero serbatoio principale  */
	f_Vincoli=1;	// controllo se esistono regolatori con portata controllata
	QMaxRegTot=0;	// portata totale imposta ai regolatori della rete connessa

	/* azzero erogazioni nodi pompe */
	if(Num_pomp > 0) {
		i=1;
		while(i<=Lati) {
			if(*(Nl+i) > 20) {
				*(Qe+(*(Ni+i))) = (float) 0.0;
				*(Qe+(*(Nf+i))) = (float) 0.0;
			}
			i++;
		}
	}

	i=1;
	while(i<=Nodi) {
		if(*(Nododel+i) == '0') {
			if(sono_serb(i)) {
				*(Qe+i) = (float) 0.0;
				if(Nc != 0) {
					temp_Lati_f++;
					*(Ni+(Lati+temp_Lati_f)) = Nc;
					*(Nf+(Lati+temp_Lati_f)) = i;
					/* *(D+(Lati+temp_Lati_f))  = (float) 1.0e10;
					if(Fluido == ACQUA)  *(Nl+(Lati+temp_Lati_f)) = 0; */
				}
				else Nc = i;

				/* portata imposta al regolatore, definiamo -1 = regol. senza vicoli */
				if (*(QMaxReg+i) == -1.0)
					f_Vincoli = 0;	/* ho un regolatore senza limitazioni */
				else
					QMaxRegTot = QMaxRegTot + *(QMaxReg+i);
			}
			qtot = qtot + (*(Qe+i));
			nodi_reali++;
		}
		i++;
	}

	/* 
	** NS - Se almeno uno dei regolatori non ha vincoli (cioè per convenzione è posto a -1) 
	** allora non faccio controlli altrimenti verifico che la portata massima dei regolatori
	** non superi il consumo totale della rete (qtot)
	*/
	if (f_Vincoli && Fluido != ACQUA){
		if (qtot > QMaxRegTot) {
			//\nLimitazione imposta ai regolatori minore del consumo: %9.2f \n
			if(error(26, 0, 0, QMaxRegTot, qtot)==FALSO) return FALSO;
			// exit(EXIT_FAILURE); // stefano 23/03/2006 genera errore in ACAD: chiude ACAD
         return FALSO;
		}
	}

	if(Nc == 0) {
		/* Non esistono serbatoi/regolatori o punti di immissione in rete */
		if(error(3, 0, 0, 0, 0)==FALSO) return FALSO;
		// exit(EXIT_FAILURE); // stefano 23/03/2006 genera errore in ACAD: chiude ACAD
      return FALSO;
	}

	/* Nodo_partenza = Nc;  controllare */
	i=1;
	while(i<=Lati) {
		if(*(Latodel+i) == '0') {
			qtot = qtot + ((*(Nl+i)>=0) ? *(Qd+i) : (float) 0.0);
			lati_reali++;
			if(*(Nl+i) == -1) ++num_prv_reali;
			if(*(Nl+i) > 20) ++nstp_reali;
		}
		i++;
	}

	/* assegno erogazione a serbatoio principale */
	*(Qe+Nc) = (qtot > 0.0) ? (-qtot) : (float) 0.0;

	Lati_f = temp_Lati_f;
	// 05/02/2007 commentata da Stefano Maglie = lati_reali - nodi_reali + 1 + Lati_f;
	// la riga sotto consente di superare un problema: in caso di reti non connesse
	// il numero delle maglie viene calcolato non corretto. A volte può essere < 0
	// Per non indagare a fondo e rivedere la logica si decide di inserire un valore
	// che consente di superare il problema a scapito dell'utilizzo di RAM e rimandare a
	// tempi migliori la soluzione.
	//PAOLETTA		dopo la modifica errore nella funzione
	//				long bu(long npart)
	//Maglie = nodi_reali + 1 + Lati_f;		//lati_reali - nodi_reali + 1 + Lati_f;
	Maglie = lati_reali - nodi_reali + 1 + Lati_f;

	/* parte aggiunta perche' Maglie puo' aumentare */
	if(Maglie > Maglie_max)
	{
		Maglie_max = Maglie;
      // roby 22/12/2011 il realloc alloca nuovi elementi senza porli = 0
      D_k = (long **) _recalloc(D_k, (Maglie+1+EXTRA), sizeof(long *)); // roby 
		if(D_k == NULL) exit(1);
		Tot_mem = Tot_mem + ((Maglie+1+EXTRA) * sizeof(long *));
      
		// 02/08/96 inizio - aggiunta perchè anche A e B vengono modificati
		A = (float *) realloc(A, (Maglie+1+EXTRA) * sizeof(float));
		if(A == NULL) exit(1);
		Tot_mem = Tot_mem + ((Maglie+1+EXTRA) * sizeof(float));

		B = (float *) realloc(B, (Maglie+1+EXTRA) * sizeof(float));
		if(B == NULL) exit(1);
		Tot_mem = (Tot_mem + ((Maglie+1+EXTRA) * sizeof(float)));
		// 02/08/96 fine - aggiunta perchè anche A e B vengono modificati
	}

	i=Lati+1;
	while(i<=Lati+Lati_f) {
		if(Fluido == GASM) /* bug eliminato 27/12/95 */
			*(Scabr+i) = (*(Carico+*(Ni+i))-*(Carico+*(Nf+i))) * (*(Carico+*(Ni+i))+*(Carico+*(Nf+i)));
		else
			*(Scabr+i) = (*(Carico+*(Ni+i))-*(Carico+*(Nf+i)));
		*(Q+i) = (float) 0.0;
		*(Latodel+i) = '0';
		i++;
	}

	/* "\nRete connessa:\n
		Nodi: %d   Lati: %d   Maglie: %d   Maglie fittizie: %3d"  */
	MyPrintf(td_cad2(9), nodi_reali, lati_reali, Maglie,Lati_f);
   #ifdef EXE_VERS
   	fwprintf(File_esecuzione, td_cad2(9), nodi_reali, lati_reali, Maglie, Lati_f);
   #endif

	if(Fluido == ACQUA)
   {
		/* "\nStazioni di pompaggio: %3d\nValvole PRV: %3d" */
		MyPrintf(td_cad2(8),nstp_reali,num_prv_reali);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(8), nstp_reali, num_prv_reali);
      #endif
	}
	else
   {
		/* GAS_PRV: "\nValvole PRV: %3d" */
		MyPrintf(td_cad2(59), num_prv_reali);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(59), num_prv_reali);
      #endif
	}

	if(Debug_on) {
		fwprintf(FileDEBUG, _T("nodi connessi:\n"));
		i=1;
		while(i<=Nodi) {
			if(*(Nododel+i) == '0') fprintf(FileDEBUG,"%d\n",*(Nnn+i));
			i++;
		}

		fprintf(FileDEBUG,"lati connessi:\n");
		i=1;
		while(i<=Lati+Num_serb-1) {
			if(*(Latodel+i) == '0')
				fprintf(FileDEBUG,"%d,%d\n", *(Nnn+*(Ni+i)),*(Nnn+*(Nf+i)));
			i++;
		}
	}

	carico_np('1');  /* ricarico il vettore np */

	/* controllo dati immessi nuova rete */

	i=1;
	while(i<=Lati+Lati_f) {
		jj3 = 0;
		j=1;
		while(j<=Lati+Lati_f) {
			if(*(Ni+j) == *(Ni+i) && *(Ni+j) != 0 && *(Ni+i) != 0)
				if(*(Nf+j) == *(Nf+i) && *(Nf+j) != 0 && *(Nf+j) != 0) {
					jj3++;
					if(jj3 > 1) {
						/* Asola, vedi lati con nodo iniziale %d e nodo finale %d */
						if(error(6, *(Nnn+(*(Ni+j))), *(Nnn+(*(Nf+j))), 0, 0)==FALSO) return FALSO;
						if(i>Lati || j>Lati)
							/* L'asola e' composta dal lato fra due serbatoi\ne il corrispondente lato fittizio */
							if(error(25, 0, 0, 0, 0)==FALSO) return FALSO;
						exit(1);
					}
				}
			j++;
		}
		i++;
	}

	return ((nodi_reali < 2) || (lati_reali < 1)) ? 0 : 1;
}

int sono_serb(int nodo)
{
   return (((!Ver24 && (nodo <= Num_serb)) ||
		   (Ver24 && (nodo <= Num_serb) &&
                     ((*(Serb24+nodo) > 0) || (*(Nnn+nodo) > Last2))))
          ? 1 : 0);
}

/* ==================================================== */
int carico_serb_start()
{
	long i,j;

	i=1;
	while(i<=Num_serb) {
		Tot_num_serb = Tot_num_serb + *(Np+i);
		i++;
	}

	/*Tot_num_serb = Tot_num_serb - (Lati_f*2);*/    /* partizioni totali serb */

	Coll_Serb_start = (long *) calloc((Tot_num_serb+EXTRA), sizeof(long));
	if(Coll_Serb_start == NULL) {
		MyPrintf(_T("\nMEMORIA INSUFFICIENTE in carico_serb_start"));
		exit(1);
	}

	Tot_num_serb = 0;

	i=1;
	while(i<=Num_serb) {
		j=1;
		while(j<=Lati) {
			if(*(Ni+j) == i) {
				Tot_num_serb++;
				*(Coll_Serb_start+Tot_num_serb) = *(Nf+j);
			}
			if(*(Nf+j) == i) {
				Tot_num_serb++;
				*(Coll_Serb_start+Tot_num_serb) = *(Ni+j);
			}
			j++;
		}
		i++;
	}

	Primo_giro = 0;

	i=1;
	while(i<=Lati) {
		if(*(Nl+i)>20 && *(Qd+i) == -1.0)
			j=1;
			while(j<=Tot_num_serb) {
				if(*(Coll_Serb_start+j) == *(Nf+i)) {
					*(Coll_Serb_start+j) = -1 * *(Coll_Serb_start+j);
					break;
				}
				j++;
			}
		i++;
	}

	j=1;
	while(j<=Tot_num_serb) {
		if(*(Coll_Serb_start+j) > 0) {
			Nodo_partenza = *(Coll_Serb_start+j);
			break;
		}
		j++;
	}

	if(Nodo_partenza == 0) {
		/* Non esistono serbatoi/regolatori o punti di immissione in rete */
		if(error(3, 0, 0, 0, 0)==FALSO) return FALSO;
		exit(1);
	}

	return VERO;
}

/* ================================================================== */
/*
**  EXT:
*/
long ext(long *dati, long *punt, long riga, long col)
{
	long   val;

	if ( riga<=0 || riga>=*punt ||
		 col<=0  || *(punt+riga)>*dati ||
		 col > (*(punt+riga+1) - *(punt+riga)) )
		val = (long) NULL;                           /* ERRORE */
	else
		val = *(dati + (*(punt + riga)) + col - 1);   /* TUTTO OK */

	return(val);
}

/* ================================================================== */
/*
**
**  W_EXT:
**  Scrive un valore VAL nel vettore DATI il cui puntatore e` PUNT
**  nella posizione riga, col (corrispondente alla posizione
**  dell'elemento nella matrice che ha originato DATI).
**  In caso di errore ritorna 1.
**
*/
int w_ext(long *dati, long *punt, long riga, long col, long val)
{
	long   rc;

	if ( riga<=0 || riga>=*punt ||
		 col<=0  || *(punt+riga)>*dati ||
		 col > (*(punt+riga+1) - *(punt+riga)) )
		rc = (int) 1;                                       /* ERRORE */
	else {
		*(dati + (*(punt + riga)) + col - 1) = val;   /* SCRIVO IL VALORE */
		rc = (int) 0;
	}

	return(rc);
}

/* ==================================================== */
void orologio()
{
	Lancetta++;
	switch (Lancetta) {
	case 0:
		MyPrintf(_T("\b|"));
		break;
	case 1:
		MyPrintf(_T("\b/"));
		break;
	case 2:
		MyPrintf(_T("\b-"));
		break;
	case 3:
		MyPrintf(_T("\b\\"));
		Lancetta = -1;
		break;
	}
	return;
}

long cod_prv(long punt)
{
	long   rc=0, i;

	if(!Connessioni) { /* 15-02-96 - stefano - bug per valvole chiuse/aperte */
		i=1;
		while(i<=Num_prv) {
			if(labs(*(Lato_prv+i)) == punt) {
				rc = labs(*(Lato_prv+i)) / *(Lato_prv+i) * i;
				break;
			}
			i++;
		}
	}

	return(rc);
}

int cod_vv(long punt)
{
	if(!Connessioni) { /* 15-02-96 - stefano - bug per valvole chiuse/aperte */
		if((*(Nl+punt) == -2 ||  /* valvola VV */
			*(Nl+punt) > 20) &&   /* pompa */
			*(Qd+punt) == -1.0) return 1;  /* -1 chiusa 0 aperta */
	}

	return 0;
}

long codice_prv_taglio()
{
	long quale=0, i;

   i=1;   // stefano - 12/10/98 - bug: spostato controllo se il lato prv è > 0 sopra (altrimenti cerca in np un codice negativo)
	while(i<=Num_prv) {
		if((*(Lato_prv+i) > 0) &&
         (*(Np+ *(Ni+ *(Lato_prv+i))) > 0) &&
			(*(Np+ *(Nf+ *(Lato_prv+i))) > 0) &&
			(*(Latodel+ *(Lato_prv+i)) == '0')) {
			quale = *(Ni+ *(Lato_prv+i));
			break;
		}
		i++;
	}

	return quale;
}

int calcolo_D_prv()
{
	long i, j;

	i=1;
	while(i<=Num_prv) {
      if(*(Lato_prv+i) > 0) {   // stefano 12/10/98 - inserito per bug - non deve cercare nei vettori per *(Lato_prv+i) < 0
		   j=1;
		   while(j<=Lati) {
			   if((*(Ni+*(Lato_prv+i)) < Last3) &&
				   (*(Ni+*(Lato_prv+i)) == *(Ni+j)) ||
				   (*(Ni+*(Lato_prv+i)) == *(Nf+j))) {
				   *(D_prv+i) = *(D+j);
				   *(D+*(Lato_prv+i)) = *(D+j);
				   break;
			   }
			   else if((*(Nf+*(Lato_prv+i)) < Last3) &&
				   (*(Nf+*(Lato_prv+i)) == *(Ni+j)) ||
				   (*(Nf+*(Lato_prv+i)) == *(Nf+j))) {
				   *(D_prv+i) = *(D+j);
				   *(D+*(Lato_prv+i)) = *(D+j);
				   break;
			   }
			   j++;
		   }
		}
		i++;
	}

   return VERO;
}

/* ================================================================== */
/*
** DFROMV  trova il diametro ottimale fissata la velocita`
*/

int dfromv()
{
	long i;
	float prec_D, curr_toll, max_toll_D = (float) 0.0;

	if(Debug_on) fprintf(FileDEBUG,"\nCalcolo diam\nlato  prec  att  max_toll\n");
	i=1;
	while(i<=Lati) {
		prec_D = *(D+i);
		if((*(Dtocalc+i)=='1')&&(*(Latodel+i)=='0')) *(D+i) = dtocomm(diamott(i));
		if(*(D+i) <= 0.0) {
			if(Fluido == ACQUA) *(D+i) = (float) 0.01;
			else *(D+i) = (float) 100.0;
		}
		curr_toll = (float) fabs(*(D+i)-prec_D);
		if(curr_toll > max_toll_D) max_toll_D = curr_toll;
		if(Debug_on) fprintf(FileDEBUG,"%4d %6.1f %6.1f %6.1f\n",i,prec_D,*(D+i),max_toll_D);
		i++;
	}
	if(Fluido == ACQUA) max_toll_D = max_toll_D * (float) 1000.0;
	/* Iterazione %3d: massimo scostamento diametri => %6.1f mm */
	if(error(-1, Ciclo, 0, max_toll_D, 0)==FALSO) return FALSO;

	/* controllo diametri con la tolleranza indicata */
	if(max_toll_D <= Max_toll_D) Num_loop_diam = 1;

	return VERO;
}

/* =========================================================== */
/*
** DIAMOTT  trova il diametro ottimale fissata la velocita`
*/

float diamott(long i)
{
	float diam, qi, qf;

	qi = *(Q+i);
	qf = qi - *(Qd+i);
	if (fabs(qi) < fabs(qf)) qi = qi - qf;

	if(Veloc > 0.0) {
		if(Fluido == ACQUA) {
			/* scrivo se il lato non e' un pezzo speciale */
			if((*(Nl+i)<20)&&(*(Nl+i)>0))
				/* calcolo diametro ottimale data la velocita` */
				diam = (float) sqrt((float) fabs(((float) 4 * qi)/(PG * Veloc)));
			else
				diam = *(D+i);
		}
		else // GAS
			/* GAS_PRV: scrivo se il lato non e' un pezzo speciale */
			if(*(Nl+i)>=0)
				/* calcolo diametro ottimale data la velocita` */
				diam = (float) 1000 * (float) sqrt(fabs(((float) 4 * qi / (float) 3600.0)/(PG * Veloc)));
			else
				diam = *(D+i);
	}
   else {
      /*if(Fluido == ACQUA)
         diam = diam_acqua(i, qi);
      else*/
         diam = diam_gas(i, qi);
   }

   return ((float) diam);
}

/* =========================================================== */
/*
** DTOCOMM  trova il diametro commerciale dato un diametro qualsiasi
*/

float dtocomm(float diamcalc)

{
	long  i;
	float diam = diamcalc;

	if (Diametri > 0) {
		diam = (float) 0.0;
		if(diamcalc < *(Diamcomm+1))
			diam = *(Diamcomm+1);
		else {
			i=1;
			while(i<=(Diametri-1)) {
				if(diamcalc < *(Diamcomm+(i+1))) {
				   diam = (fabs((double) (diamcalc-*(Diamcomm+i))) <
						   fabs((double) (*(Diamcomm+(i+1))-diamcalc)))
						  ? *(Diamcomm+i) : *(Diamcomm+(i+1));
				   break;
				}
				i++;
			}
		}
		if(diam == 0) diam = *(Diamcomm+Diametri);
	}

	if(diam == 0.0) diam = diamcalc;   /* diametro non trovato = diametro calcolato */
 
	return diam;
}

/* -------------  modificare dopo --------------------*/
/* =========================================================== */
/*
** finp_con
*/
int finp_con(FILE *fp)
{
	long    i, in, i1, i2, ik, nodo_start, rc;
	char    c1[10];

	Nodo_partenza = 1;
	Num_serb = 1;   /* fittizio per ripristino_del() */

	fscanf(fp,"%d",&nodo_start);
	fscanf(fp,"%d",&Nodi);
	fscanf(fp,"%d",&Lati);
	Lati_f = 0;   /* inserito 21/12/95 */

	rc = mem_array();   /* DEFINIZIONE DINAMICA DEGLI ARRAY */
	if(rc) 
   {
		/*  "\nMEMORIA INSUFFICIENTE interruzione al passo %d, allocati %ld bytes" */
		MyPrintf(td_cad2(4), rc, Tot_mem);
      #ifdef EXE_VERS
   		fwprintf(File_esecuzione, td_cad2(4), rc, Tot_mem);
      #endif
		exit(1);
	}

	i=1;
	while(i<=Nodi) {
		fscanf(fp,"%d", &in);
		*(Nnn+i) = in;
		if(*(Nnn+i) == nodo_start) Nodo_partenza = i;
		i++;
	}

	i=1;
	while(i<=Lati) {
		fscanf(fp,"%d,%d,%s",&i1,&i2,c1);
		ik=1;
		while(ik<=Nodi) {
			if(*(Nnn+ik) == i1) *(Ni+i)=ik;
			if(*(Nnn+ik) == i2) *(Nf+i)=ik;
			ik++;
		}
		strcpy(*(Handle+i),c1);
		i++;
	}

	/*MyPrintf(td_cad2(57));*/  /* "terminato.\n" */
   #ifdef EXE_VERS
	   //fwprintf(File_esecuzione, td_cad2(57));
   #endif
	/* fclose(fp);    eliminato 21/12/95 perche' gia' chiamato sopra */

	ripristino_del();

	return VERO;
}

/* ================================================================== */
/*
**  fout_con
*/
int fout_con()
{
	FILE   *fp; 
	long   i;
	TCHAR  ftempinp[150];

   if (Aborted) return FALSO;

	if((fp=_wfopen(wcscat(wcscpy(ftempinp,Ftemp),_T("hydro2.u28")),_T("w")))==NULL)
   {
		MyPrintf(_T("\nfmgr out error hydro2.u28\n"));
		exit(1);
	}

	i=1;
	while(i<=Lati) {
		if(*(Latodel+i) == '0') fprintf(fp,"%s\n", *(Handle+i));
		i++;
	}

	fclose(fp);

	return VERO;
}

/* ============================
**  Subroutine RES
*/
int res(long punt, float *y, float *w)
{
	switch(Fluido) {
	case ACQUA:
		resh2o(punt,y,w);
		break;
	case GASB:
		resgas(punt,y,w);
		break;
	case GASM:
		resgas(punt,y,w);
		break;
	}

	return VERO;
}

int resh2o(long punt, float *y, float *w)
{
	long    ipunp, kk, ko, ihx, prv, iter_res;
	float   fff, yyy, rrr, errcol=(float) 1, qq, vo, ee, sl, v, vl, z, a, r;
	// register float f[8], g[8];
	float f[8], g[8];

	f[1] = (float) 1.31e-6;         /* Colebrook */
	f[3] = (float) 2.0;             /* Strickler */
	g[3] = (float) (16.0/3.0);      /* Strickler */
	f[4] = (float) 1.85;            /* Hazen - Williams */
	g[4] = (float) 4.87;            /* Hazen - Williams */
	f[5] = (float) 1.785714;        /* Bazin */
	g[5] = (float) 4.785714;        /* Bazin */
	f[6] = (float) 2.0;             /* Darcy */
	g[6] = (float) 5.0;             /* Darcy */
	f[7] = (float) 1.83;            /* Orsi */
	g[7] = (float) 4.87;            /* Orsi */

	if(*(Nl+punt) > 20) {           /* STAZIONE POMPAGGIO */
		ipunp = *(Nl+punt) - 20;              /* PUNTATORE POMPA */
		*y = interp(*(Q+punt),ipunp);
		*w = (float) 0.0;
		fff = *y;

		/* cerco stabilita` */

		*y = (*(Ypprec+ipunp) + Ppapp * ((*y) - *(Ypprec+ipunp)));
		fff = (float) fabs(*y - *(Ypprec+ipunp));
		*(Ypprec+ipunp) = (*y);
		kk = 1;
		qq = Pq[1][ipunp];
		if(Ph[1][ipunp] <= (*y)) goto lab12;
		kk=2;
		while(kk<=10) {
			if(Ph[kk][ipunp] <= *y) {
				yyy = ((*y) - Ph[kk-1][ipunp]) / (Ph[kk][ipunp] - Ph[kk-1][ipunp]);
				rrr = Pq[kk][ipunp]-Pq[kk-1][ipunp];
				qq = Pq[kk - 1][ipunp] + rrr * yyy;
				goto lab12;
			}
			kk++;
		}
		qq = Pq[10][ipunp];

		/*
		** controllo che la portata nel lato e quella della pompa (calcolata
		** sulla curva) non differiscano piu' di Pprec(isione)
		*/
		lab12:

		/*  if(fabs(*(Q+punt) - qq) > Pprec) Ictrp = 1;*/
		/*  if(fabs(fff - *y) > 0.1) Ictrp = 1;*/
		/*  if(fff > 1.0) Ictrp = 1;*/
		if(fff > Pprec) Ictrp = 1;  /* 06/05/96 stefano per Aosta */

		*y = (Pq[1][ipunp] < 0) ? *y : (-1 * *y); /* nodo > Last2 come nf */

		return VERO;
	}          /* fine if stazione pompaggio */

	ko = *(Nl+punt);
	ihx = ko + 1;

	switch(ihx) {
    case -1:                                    /* valvola VV */
		*y = (float) 0.0;
		*w = (float) 0.0;
		break;
    case 0:                                     /* valvola PRV */
        prv = cod_prv(punt);
		/* calcolo w con Strickler */
        if(*(Q+punt) == 0) {
			*y = (float) 0;
			*w = (float) 0;
        }
        else {
			*y = ((int) (*(Q+punt) / fabs(*(Q+punt)))) * *(Delta_prv+prv);
			/* *(L+punt) = 1; lunghezza PRV = 1m */
			r = (float) 10.29359 / (float) pow((double) *(Scabr+punt), (double) 2);
         // se la perdita di carico è diversa da zero calcolo il diametro
         // altrimenti lascio il precedente
			if(*y != 0.0)
			//if(*y == 0.0)
			//	*(D+punt) = *(D_prv+prv);
			//else
				*(D+punt) = (float) pow((double) fabs(r / *y *
			                           *(Q+punt) * fabs(*(Q+punt))), ((double) 3.0 / (double) 16.0));
			   r = (float) 10.29359 / (float) pow((double) *(Scabr+punt), (double)2) / (float) pow((double) *(D+punt), (double)g[3]);
			   *w = f[3] * r * (float) pow((double) fabs(*(Q+punt)),(double)(f[3]-(float) 1.0));
        }
        break;
    case 1:
		*y = *(Scabr+punt);
		*w = (float) 0.0;
		break;
    case 2:
    case 3:
		vo = (float) 4.0 * *(Q+punt) / PG / (float) pow((double) *(D+punt), 2);
		ee = *(Scabr+punt) * (float) 0.001 / *(D+punt);
		sl = (float) 6523.955148 / *(D+punt) * (f[*(Nl+punt)] / *(D+punt)) * (f[*(Nl+punt)] / *(D+punt));
		vl = (float) vcol(sl, ee, *(D+punt), f[*(Nl+punt)], &z);
		*y = sl;
		a = (float) 0.0;
		iter_res = (int) 0;
		while(errcol > 0.001) {
			*y = *y + a;
			if(fabs(*y) > sl)
				v=(float) vcol(*y, ee, *(D+punt), f[*(Nl+punt)], w);
			else {
				*w = vl / sl;
				v = (*w) * (*y);
			}
			a = (vo - v) / (*w);
			errcol = (float) fabs(vo - v);
			iter_res++;
			if(iter_res > 30) {
				/* Piu` di trenta iterazioni */
				if(error(16, 0, 0, 0, 0)==FALSO) return FALSO;
				exit(1);                   /* FINE */
			}
		}
		*w = vl / sl;
		if(fabs(*y) > sl) v = (float) vcol(*y, ee, *(D+punt), f[*(Nl+punt)], w);
		*y = *y * *(L+punt);
		*w = *(L+punt) * (float) 4.0 / PG / (float) (pow((double) *(D+punt), 2)) / (*w);
		break;
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
		v = *(Scabr+punt);
		if(ko == 3) v = (float) 10.29359 / (float) pow((double) *(Scabr+punt), 2);
		else if(ko == 4) v = (float) 1./((float) .0941 * (float)pow((double) *(Scabr+punt), 1.85));
		else if(ko == 5) v = (float) .000857 * (float) (pow((double) (1 + 2 * *(Scabr+punt) / pow((double) *(D+punt), 0.5)), 2));
		else if(ko == 6) v = (float) .000648 * (float) (pow((double) (1 + 2 * *(Scabr+punt) / pow((double) *(D+punt), 0.5)), 2));
		else if(ko == 7) v = v * (float) 0.0009862;
		r = v * *(L+punt) / (float) pow((double) *(D+punt), (double)g[ko]);
		if(*(Qd+punt) > 0.0) {
			if(ko <= 2) {
				if(error(15, 0, 0, 0, 0)==FALSO) return FALSO;
				exit(1);
			}
			*w = r / *(Qd+punt) * (*(Q+punt) * (float) pow((double) fabs(*(Q+punt)), (f[ko] - 1.0)) - (*(Q+punt)-*(Qd+punt))* (float) pow((double)
			                      (float) fabs(*(Q+punt)-*(Qd+punt)), (f[ko] - 1.0)));
			*y = r / *(Qd+punt) / (f[ko] + (float)1) * (float) (pow((double) fabs(*(Q+punt)), (f[ko] + 1.0)) - pow((double) fabs(*(Q+punt)-*(Qd+punt)),
			                                            (double)(f[ko] + (float) 1.0)));
		}
		else {
			*w = f[*(Nl+punt)] * r * (float) pow((double) fabs(*(Q+punt)), (f[*(Nl+punt)] - 1.0));
			*y = r * *(Q+punt) * (float) pow((double) fabs(*(Q+punt)), (f[*(Nl+punt)] - 1.0));
		}
		break;
	default:
		*y = *(Scabr+punt);
		*w = (float) 0.0;
		break;
	}

	return VERO;
}

/*
** VCOL - Ritorna anche il parametro z modificato!!!
*/
double vcol(double y, double e, double d, double v, float *z)
{
   float a=(float) 0, sr=(float) 0;

   sr = (float) (sqrt(19.62 * d) / sqrt(fabs(y)));
   a = (float) (-2.0) * (float) log10 (e / (float) 3.715 + 2.51 / d / (sr * fabs(y)) * v);
   *z = (float) 0.5 * sr * (a + (float) 0.868589 / (float) ( 1 + e * d / v * (sr * fabs(y)) / (float) 9.32465));

   return(y * a * sr);
}

int resgas(long punt, float *y, float *w)
{
	float   r,costante;
	long    prv;

	if(Fluido == GASB) costante = (float) 2.32e8;
	else costante = (float) 48600.0;

   /* if(*(Scabr+punt)==0.0) { bug 21/12/95 */
	if(punt <= Lati) { 		/* ovvero non lato fittizio */
// GAS_PRV: inizio ----------------
		if (*(Nl+punt) == -1 && GlobalIter > MaxGlobalIter){ /* GAS_PRV */ /* Stefano 08092008 */
			prv = cod_prv(punt);
			if(*(Q+punt) == 0) {
				*y = (float) 0;
				*w = (float) 0;
			}
			else {
				*y = ((int) (*(Q+punt) / fabs(*(Q+punt)))) * *(Delta_prv+prv);
				/* *(L+punt) = 1; lunghezza PRV = 1m */
            // *(L+punt) inserita perchè 0.001 m (nel gas le lunghezze sono in km)
   			r = costante * Dens * *(L+punt) * *(Q+punt) * (float)pow((double) fabs(*(Q+punt)),(double) 0.82);
            // se la perdita di carico è diversa da zero calcolo il diametro
            // altrimenti lascio il precedente
				if(*y != 0.0)
				//if(*y == 0.0)
				//	*(D+punt) = *(D_prv+prv);
				//else
 					*(D+punt) = (float) pow((double) (r / *y), ((double) 1.0 / (double) 4.82));

            // *(L+punt) inserita perchè 0.001 m (nel gas le lunghezze sono in km)
            r = costante * Dens * *(L+punt) / (float) pow((double) *(D+punt),(double)4.82);
   			*w = (float) 1.82 * r * (float) pow((double) fabs(*(Q+punt)),(double) 0.82);
			}
        }
		else
		{
// GAS_PRV: fine --------------
			r = costante * Dens * *(L+punt) / (float) pow((double) *(D+punt),(double)4.82);
   		if (*(Qd+punt) > 0.0) {
				*w = r / *(Qd+punt) * (*(Q+punt) * (float) pow((double) fabs(*(Q+punt)),(double) 0.82)-(*(Q+punt) - *(Qd+punt)) * (float) pow(fabs((double)
					 (*(Q+punt) - *(Qd+punt))),(double) 0.82));
				*y = r / *(Qd+punt) / (float) 2.82 * (float) (pow((double) fabs(*(Q+punt)),(double) 2.82)-pow((double) fabs(*(Q+punt)-*(Qd+punt)),(double) 2.82));
			}
			else {
   				*w = (float) 1.82 * r * (float) pow((double) fabs(*(Q+punt)),(double) 0.82);
   				*y = r * *(Q+punt) * (float) pow((double) fabs(*(Q+punt)),(double) 0.82);
			}
		}
	}
	else {
   		*w = (float) 0.0;
   		*y = *(Scabr+punt);
	}

   // inizio test
   //FILE *fp;
   //if (*(IDs+punt)== 3414)
   //   if ((fp = fopen("c://gs//gs8//temp//3414.txt", "a")) != NULL)
   //   {
	  //    fprintf(fp, "%ld;%f;%f\n", GlobalIter, *(Q+punt), *y);
   //      fclose(fp); 
   //   }

	return VERO;
}

float interp(float port,int pomp)
/* dato il valore della portata e il numero identificativo della
** pompa restituisce la prevalenza
*/
{
	int     kk;
	float   yyy, rrr;

	if(((Pq[1][pomp] < 0) && (port >= Pq[1][pomp])) ||
		((Pq[1][pomp] > 0) && (port <= Pq[1][pomp])) ||
		((Pq[1][pomp] == 0) && (port == Pq[1][pomp])))  /* stefano 10/01/96 */
		return (Ph[1][pomp]);

	for(kk=2; kk<=10; kk++) {
		//if(((Pq[1][pomp] < 0) && (port >= Pq[kk][pomp])) ||
		//((Pq[1][pomp] > 0) && (port <= Pq[kk][pomp]))) { /* stefano 10/01/96 */
		if(fabs(port) <= fabs(Pq[kk][pomp])) {
			yyy = (fabs(port) - fabs(Pq[kk-1][pomp]))/fabs(Pq[kk][pomp] - Pq[kk-1][pomp]);
			rrr = Ph[kk - 1][pomp] - Ph[kk][pomp];
			return (Ph[kk - 1][pomp] - rrr * yyy);
		}
	}

	return (Ph[10][pomp]);
}

float diam_gas(long punt, float q)
{
	float   costante, diam;

	if(Fluido == GASB) costante = (float) 2.32e8;
	else costante = (float) 48600.0;

	/* if(*(Scabr+punt)==0.0) bug 21/12/95 */
	if(punt > Lati)
		diam = (float) pow((double) (costante * Dens * 1.0 / (float) Perd_unit *
		                   pow((double) fabs(q),(double) 1.82) ),(double) (1.0/4.82));

	return diam;
}

TCHAR *td_cad2(int passo)
{
	switch (passo) {
	case 0:
		return _T("\n\n     GEOsim PRESSURE NETWORK simulator h24  v.1\n\n");
	case 1:
		return _T("\nNessuna maglia");
	case 2:
		return _T("\nVerifica conclusa\n");
	case 3:
		return _T("\nMemoria disponibile : %ld bytes");
	case 4:
		return _T("\nMEMORIA INSUFFICIENTE interruzione al passo %d, allocati %ld bytes");
	case 5:
		//return "\n\n     GEOsim PRESSURE NETWORK simulator h1  v.7  (21/02/05)\n\n";
		//return "\n\n     GEOsim PRESSURE NETWORK simulator h1  v.7.2  (31/07/08)\n\n";
		return _T("\n\n     GEOsim PRESSURE NETWORK simulator h1  v.7.3  (11/05/2009)\n\n");
	case 6:
      return _T("\nMaglia:\n");
	case 7:
		return _T("\nDati input:\nNodi: %d   Lati: %d");
	case 8:
		return _T("\nStazioni di pompaggio: %3d\nValvole PRV: %3d");
	case 9:
		return _T("\nRete connessa:\nNodi: %d   Lati: %d   Maglie: %d   Maglie fittizie: %3d");
	case 10:
		return _T("\nErogazione totale: %9.2f l/s\n");
	case 11:
		return _T("\nErogazione totale: %9.2f mc/h\n");
	case 12:
		return _T("\nSimulazione ora: %3d");
	case 13:
		return _T("\nIl serbatoio %d e` vuoto.\nReinizializzo con %.1f mc.");
	case 15:
		return _T("\nIl serbatoio %d sfiora %.1f mc.");
	case 16:
		return _T("\nNumerazione lati errata - nodo iniziale %d ripetuto anche come finale.");
	case 17:
		return _T("\nNon esistono serbatoi/regolatori o punti di immissione in rete.");
	case 19:
		return _T("\nAsola, vedi lati con nodo iniziale %d e nodo finale %d.");
	case 20:
		return _T("\nLa stazione di pompaggio %d funziona in dissipazione.");
	case 21:
		return _T("\nLa stazione di pompaggio %d non puo` erogare %f l/s.");
	case 22:
		return _T("\nErrore non codificato - contattare AMGA");
	case 23:
		return _T("\nControllare il lato %d (%d %.0f).");
	case 25:
		return _T("\nErrore nella costruzione delle maglie; lato %d nodi %d e %.0f");
	case 26:
		return _T("\nPiu` di trenta iterazioni");
	case 27:
		return _T("\nErrore N. %d\n>>> controllare file errore.out <<<\n\n");
	case 28:
		return _T("\nIl serbatoio %d ha un'altezza di sfioro <=0.");
	case 29:
		return _T("\nIl serbatoio %d ha un volume <=0.");
	case 30:
		return _T("\nLa quota minima del serbatoio %d deve essere >=0.");
	case 31:
		return _T("\nApro il serbatoio %d.\nPiezometrica inferiore al pelo libero (%.1f/%.1f).");
	case 32:
		return _T("\nChiudo il serbatoio %d.");
	case 33:
		return _T("\nRipeto simulazione ora %d.");
	case 34:
		return _T("\nCalcolo diametri concluso.");
	case 35:
		return _T("\nCalcolo diametri :\nVelocita` del fluido: %4.2f");
	case 36:
		return _T("\n%d lati su %d con diametri da calcolare\n");
	case 37:
		return _T("\nnodo = %d");
	case 38:
		return _T("\nErrore N. %2d");
	case 39:
		return _T("\nMEMORIA INSUFFICIENTE");
	case 40:
		return _T("\nMemoria utilizzata  : %ld bytes");
	case 41:
		return _T("\n\n Passo     Carico");
	case 42:
		if(Fluido == ACQUA)
			return _T("\n\nControllo sul carico (%f m)");
		else if(Fluido == GASB)
			return _T("\n\nControllo sul carico (%f mm)");
		else if(Fluido == GASM)
			return _T("\n\nControllo sul carico (%f kg/cmq)");
	case 43:
		return _T("\n\n Passo     Portata");
	case 44:
		if(Fluido == ACQUA)
			return _T("\n\nControllo sulla portata (%f l/s)");
		else
			return _T("\n\nControllo sulla portata (%f mc/h)");
	case 45:
		return _T("\n\n Passo     Carico       Portata ");
	case 46:
		if(Fluido == ACQUA)
			return _T("\n\nControllo su carico (%f m) %s portata (%f l/s)");
		else if(Fluido == GASB)
			return _T("\n\nControllo su carico (%f mm) %s portata (%f mc/h)");
		else if(Fluido == GASM)
			return _T("\n\nControllo su carico (%f kg/cmq) %s portata (%f mc/h)");
	case 47:
		return _T("                                                   Chiudo PRV %d\r");
	case 48:
		return _T("\n\nFatal Error %d: PRV %d");
	case 49:
		return _T("\nControllo topologia  |");
	case 50:
		return _T("                                                   Chiudo PRV %d\r");
	case 51:
		return _T("\b\bcompletato.");
	case 52:
		return _T("\nErrore: il lato %d %d punta al nodo %.0f che non esiste.\n");
	case 53:
		return _T("\nErrore: il lato %d %d ha il diametro = 0.\n");
	case 54:
		return _T("\nIterazione %3d: massimo scostamento diametri => %6.1f mm\n");
	case 55:
		return _T(" Il nodo %d ha carico < 0; \nATTENZIONE i risultati della simulazione possono non essere attendibili!");
	case 56:
		return _T("L'asola e' composta dal lato fra due serbatoi\ne il corrispondente lato fittizio");
	case 57:
		return _T("terminato.\n");
	case 58:
		return _T("\nLimitazione imposta ai regolatori: %9.2f minore del consumo: %9.2f \n");
	case 59:
		return _T("\nValvole PRV: %3d");
	case 60:
		return _T("\nVerifica abortita\n");
   /* STEFANO 05/05/2014 inizio */
   case 61:
      return _T("Perdita di carico troppo elevata nel lato %d - %d. Legge resistenza non applicabile.\n");
   /* STEFANO 05/05/2014 fine */
	default:
		return _T(" ");
	}
}


//////////////////////////////
// funzione di debug per PRV test
void Stampa_valori_prv(long contatore_iter)
{
	FILE    *fp, *fp1;
   long	prv = 1;
   long Nodo_iniz, Nodo_final;
   double car_monte=-1, car_valle=-1, Portata, Diametro;
   double Prev_car_monte, Prev_car_valle, deltaprv;
   char connesso;
   long OnOff;

   if ((fp = fopen("c://gs//gs8//temp//TestPrv.txt", "a")) == NULL) return;
   if ((fp1 = fopen("c://gs//gs8//temp//TestPrv1.txt", "a")) == NULL) return;
   //fprintf(fp, "Contatore = %ld\n", contatore_iter);
   //fprintf(fp, "Nodo_iniz, Nodo_final,car_monte, car_valle, Portata, Diametro\n");

	while(prv<=Num_prv)
   {
		if((*(Lato_prv+prv) > 0) && (*(Latodel+*(Lato_prv+prv)) == '0')) {
         Nodo_iniz = *(Nnn + *(Ni+*(Lato_prv+prv)));
         Nodo_final = *(Nnn + *(Nf+*(Lato_prv+prv)));

			if(*(Nnn + *(Ni+*(Lato_prv+prv))) > *(Nnn + *(Nf+*(Lato_prv+prv)))) {
				car_monte = *(Carico+ *(Nf+ *(Lato_prv+prv)));
				car_valle = *(Carico+ *(Ni+ *(Lato_prv+prv)));
			}
			else {
				car_monte = *(Carico+ *(Ni+ *(Lato_prv+prv)));
				car_valle = *(Carico+ *(Nf+ *(Lato_prv+prv)));
			}

		   Portata  = *(Q+ *(Lato_prv+prv));
		   Diametro = *(D+ *(Lato_prv+prv));

         Prev_car_monte = *(H_monte+prv);
         Prev_car_valle = *(H_valle+prv);
         deltaprv       = *(Delta_prv+prv);
         connesso       = *(Latodel+*(Lato_prv+prv));
         OnOff          = *(Lato_prv+prv);

		   fprintf(fp, "%ld;%f;", 
                                                         Nodo_iniz, 
                                                         Diametro);
		   //fprintf(fp, "%ld;%ld;%ld;%f;%f;%f;%f;%f;%f;%f;%c;%ld", contatore_iter,
     //                                                    Nodo_iniz, Nodo_final,
     //                                                    car_monte, car_valle,
     //                                                    Portata, Diametro,
     //                                                    Prev_car_monte, Prev_car_valle,
     //                                                    deltaprv, connesso, 
     //                                                    OnOff);
		   fprintf(fp1, "%ld;%ld;%ld;%f;%f;%f;%f;%f;%f;%f;%c;%ld\n", contatore_iter,
                                                         Nodo_iniz, Nodo_final,
                                                         car_monte, car_valle,
                                                         Portata, Diametro,
                                                         Prev_car_monte, Prev_car_valle,
                                                         deltaprv, connesso,
                                                         OnOff);
      }
      prv++;
   }
   fprintf(fp, "\n");
   fclose(fp); 
   fclose(fp1); 
}



/* 03/01/96 - inserimento modifiche per velocita' gas media */
/* 05/01/96 - correzione bug velocita' (Dens_e = 0) => inserito controllo carichi = 0 in fout() *
/* 10/01/96 - correzione in interp per portate == 0 */
/* 10/01/96 - correzione in fout per nodi con carico < 0 e non == 0 */
/* 10/01/96 - inserito carico = -1 e non 0.0 in distrhn */
/* 29/01/96 - inserito controllo per asola fra due serbatoi (errore nel prg se si inserisce un nodo fra due serb e non se si disegnano due lati) */
/* 15-02-96 - stefano - bug per valvole chiuse/aperte */
/* 28-02-96 - stefano - bug per reti separate (mi sembrava di averlo gia' eliminato if(Cambio_maglie && Ver24) { possibile attivazione in contr_serb || 28/02/96 stefano - aggiunto Ver24 */
/* 07/05/96 - stefano - reinserito Pprec nel controllo delle pompe  => if(fff > Pprec) Ictrp = 1;  06/05/96 stefano per Aosta */
/* 01/08/96 - stefano - inserito abs() in ciclo_maglie (non sono sicuro)*/
/* 21/11/97 - inizio trasformazione da int a long */
