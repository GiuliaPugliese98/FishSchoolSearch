/**************************************************************************************
 *
 * CdL Magistrale in Ingegneria Informatica
 * Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
 *
 * Progetto dell'algoritmo Fish School Search 221 231 a
 * in linguaggio assembly x86-32 + SSE
 *
 * Fabrizio Angiulli, aprile 2019
 *
 **************************************************************************************/

/*
 *
 * Software necessario per l'esecuzione:
 *
 *    NASM (www.nasm.us)
 *    GCC (gcc.gnu.org)
 *
 * entrambi sono disponibili come pacchetti software
 * installabili mediante il packaging tool del sistema
 * operativo; per esempio, su Ubuntu, mediante i comandi:
 *
 *    sudo apt-get install nasm
 *    sudo apt-get install gcc
 *
 * potrebbe essere necessario installare le seguenti librerie:
 *
 *    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
 *    sudo apt-get install libc6-dev-i386
 *
 * Per generare il file eseguibile:
 *
 * nasm -f elf32 fss32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o fss32.o fss32c.c -o fss32c -lm && ./fss32c $pars
 *
 * oppure
 *
 * ./runfss32
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define type float
#define MATRIX type *
#define VECTOR type *

#define EPSILON 0.00000001

typedef struct
{
	MATRIX x;	  // posizione dei pesci
	VECTOR xh;	  // punto associato al minimo di f, soluzione del problema
	VECTOR c;	  // coefficienti della funzione
	VECTOR r;	  // numeri casuali
	int np;		  // numero di pesci, quadrato del parametro np
	int d;		  // numero di dimensioni del data set
	int iter;	  // numero di iterazioni
	type stepind;     // parametro stepind
	type stepvol;     // parametro stepvol
	type wscale;      // parametro wscale
	int display;
	int silent;
} params;

typedef struct
{
	VECTOR w;
	VECTOR deltaf;
	VECTOR deltax;
	VECTOR baricentro;
	int rand;
	type fishesWeight;
	type stepindIni; // parametro stepind
	type stepvolIni; // parametro stepvol
} var;

int np;
int d;
int iter;
int effettuato;
int allineamento = 4;
int allineamentoPerfetto = 0;
/*
 *
 *	Le funzioni sono state scritte assumendo che le matrici siano memorizzate
 * 	mediante un array (float*), in modo da occupare un unico blocco
 * 	di memoria, ma a scelta del candidato possono essere
 * 	memorizzate mediante array di array (float**).
 *
 * 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
 * 	matrici per righe (row-major order) o per colonne (column major-order).
 *
 * 	L'assunzione corrente � che le matrici siano in row-major order.
 *
 */

void *get_block(int size, int elements)
{
	return _mm_malloc(elements * size, 16);
}

void free_block(void *p)
{
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols)
{
	return (MATRIX)get_block(sizeof(type), rows * cols);
}

void dealloc_matrix(MATRIX mat)
{
	free_block(mat);
}

/*
 *
 * 	load_data
 * 	=========
 *
 *	Legge da file una matrice di N righe
 * 	e M colonne e la memorizza in un array lineare in row-major order
 *
 * 	Codifica del file:
 * 	primi 4 byte: numero di righe (N) --> numero intero
 * 	successivi 4 byte: numero di colonne (M) --> numero intero
 * 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
 *
 *****************************************************************************
 *	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
 * 	della matrice.
 *****************************************************************************
 *
 */
MATRIX load_data(char *filename, int *n, int *k)
{
	FILE *fp;
	int rows, cols, status, i;

	fp = fopen(filename, "rb");

	if (fp == NULL)
	{
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	MATRIX data = alloc_matrix(rows, cols);
	status = fread(data, sizeof(type), rows * cols, fp);
	fclose(fp);

	*n = rows;
	*k = cols;

	return data;
}

/*
 * 	save_data
 * 	=========
 *
 *	Salva su file un array lineare in row-major order
 *	come matrice di N righe e M colonne
 *
 * 	Codifica del file:
 * 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
 * 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
 * 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
 */
void save_data(char *filename, void *X, int n, int k)
{
	FILE *fp;
	int i;
	fp = fopen(filename, "wb");
	if (X != NULL)
	{
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++)
		{
			fwrite(X, sizeof(type), k, fp);
			// printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type) * k;
		}
	}
	else
	{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

// PROCEDURE ASSEMBLY

//extern void prova(params *input);
extern void addVettori(VECTOR v1, VECTOR v2, VECTOR ris, int dim);
extern void subVettori(VECTOR v1, VECTOR v2, VECTOR ris, int dim);
extern type distEuclidea(VECTOR v1, VECTOR v2, int dim);
extern VECTOR copyVector(VECTOR v, int inizio, int dim);
extern type pesoTot(VECTOR v, int dim);
extern void prodVet_x_Scalare(VECTOR v1, type s, VECTOR ris, int dim);
extern void prodVet_x_ScalareUn(VECTOR v1, type s, VECTOR ris, int dim);
extern type prodScalare(VECTOR v1, VECTOR v2, int dim);

///////////

// copyVector
VECTOR copyVector(VECTOR v, int inizio, int dim)
{
	// Controlla se l'indice di inizio è un multiplo dell'allineamento.
	if (inizio % allineamento == 0)
		return v + inizio;

	// Alloca un blocco di memoria di dimensione `dim` elementi, di tipo `type`.
	VECTOR ret = get_block(sizeof(type), dim);

	// Definisce il numero di elementi da copiare in un ciclo unrollato.
	int unrolling = 16;

	// Indice del ciclo for.
	int i = 0;

	// Ciclo unrollato che copia `unrolling` elementi alla volta.
	for (i = 0; i + unrolling < dim; i += unrolling)
	{
		ret[i] = v[i + inizio];
		ret[i + 1] = v[i + inizio + 1];
		ret[i + 2] = v[i + inizio + 2];
		ret[i + 3] = v[i + inizio + 3];
		ret[i + 4] = v[i + inizio + 4];
		ret[i + 5] = v[i + inizio + 5];
		ret[i + 6] = v[i + inizio + 6];
		ret[i + 7] = v[i + inizio + 7];
		ret[i + 8] = v[i + inizio + 8];
		ret[i + 9] = v[i + inizio + 9];
		ret[i + 10] = v[i + inizio + 10];
		ret[i + 11] = v[i + inizio + 11];
		ret[i + 12] = v[i + inizio + 12];
		ret[i + 13] = v[i + inizio + 13];
		ret[i + 14] = v[i + inizio + 14];
		ret[i + 15] = v[i + inizio + 15];
	}

	// Ciclo che copia gli elementi rimanenti.
	for (; i < dim; i++)
		ret[i] = v[i + inizio];

	// Restituisce il puntatore al vettore copiato.
	return ret;
} // copyVector

/*
void addVettori(VECTOR v1,VECTOR v2, VECTOR ris, int dim){
	for(int i =0; i<dim; i++)
		ris[i]=v1[i]+v2[i];
}

/*
type distEuclidea(VECTOR v1, VECTOR v2, int dim){
	type v=0;
	for(int i=0;i<dim;i++){
		v+= ((v2[i]-v1[i])*(v2[i]-v1[i]));
	}
	return (type)sqrt(v);
}


type pesoTot(VECTOR v, int dim){
	type tmp=0;
	for(int i=0; i<dim;i++){
		tmp+=v[i];
	}
	return tmp;
}


void prodVet_x_Scalare(VECTOR v1, type s, VECTOR ris, int dim){
	for(int i =0; i<dim; i++){
		ris[i]=v1[i]*s;
	}
}


type prodScalare(VECTOR v1, VECTOR v2,int dim){
	type ris=0.0;
	for(int i =0; i<dim; i++){
		ris+=v1[i]*v2[i];
	}
	return ris;
}


void subVettori(VECTOR v1,VECTOR v2, VECTOR ris, int dim){
	for(int i =0; i<dim; i++){
		ris[i]=v1[i]-v2[i];
	}
}
*/

// funzioneObiettivo
// f(x) = e^x + x^2 − c ◦ x
type funzioneObiettivo(VECTOR x, params *input, int dim)
{
	// Calcola il prodotto scalare tra il vettore x e se stesso.
	type x2 = prodScalare(x, x, dim);

	// Calcola l'esponenziale del prodotto scalare.
	type ex2 = (type)exp(x2);

	// Calcola il prodotto scalare tra il vettore x e il vettore c.
	type cx = prodScalare(input->c, x, dim);

	// Restituisce la funzione obiettivo.
	return ex2 + x2 - cx;
} // funzioneObiettivo

// funzioneMatrix
type funzioneMatrix(MATRIX matrice, params *input, int inizio, int dim)
{
	// Copia dalla matrice un vettore di dimensione dim a partire dall'indice inizio sulla variabile vettore
	VECTOR vettore = copyVector(matrice, inizio, dim);

	// Calcola il valore della funzione obiettivo su questo vettore
	type ret = funzioneObiettivo(vettore, input, dim);

	// Se la riga non è allineata, libera la memoria allocata per il vettore
	if ((inizio % allineamento) != 0)
		free_block(vettore);

	// Restituisce il valore calcolato della funzione obiettivo
	return ret;
} // funzioneMatrix

type getRand(params* input,int it, int pesce, int dim,int individuale)
{
    
	return input->r[it*pesce*(d+1)+dim+individuale];
}

// addMatriceVettore
void addMatriceVettore(MATRIX matrix, VECTOR vector, int riga, int d)
{
	if ((riga * d % allineamento) == 0)
	{
		// Somma vector alla riga riga di matrix con allineamento
		addVettori(matrix + riga * d, vector, matrix + riga * d, d);
	}
	else
	{
		int unroll = 4;
		int i = 0;

		for (i = 0; i + unroll < d; i += unroll)
		{
			// Somma vector alla riga riga di matrix senza allineamento
			matrix[riga * d + i] += vector[i];
			matrix[riga * d + i + 1] += vector[i + 1];
			matrix[riga * d + i + 2] += vector[i + 2];
			matrix[riga * d + i + 3] += vector[i + 3];
		}

		for (; i < d; i++)
		{
			// Ripeti il calcolo per le componenti rimanenti
			matrix[riga * d + i] += vector[i];
		}
	}
} // addMatriceVettore

// addMatriceVettoreBroad
void addMatriceVettoreBroad(MATRIX matrix, VECTOR vector, int d)
{
	int unroll = 4;
	int pesce = 0;
	for (pesce = 0; pesce + unroll < np; pesce += unroll)
	{
		// Somma vector alla riga pesce di matrix
		addMatriceVettore(matrix, vector, pesce, d);
		// Somma vector alla riga pesce + 1 di matrix
		addMatriceVettore(matrix, vector, pesce + 1, d);
		// Somma vector alla riga pesce + 2 di matrix
		addMatriceVettore(matrix, vector, pesce + 2, d);
		// Somma vector alla riga pesce + 3 di matrix
		addMatriceVettore(matrix, vector, pesce + 3, d);
	}
	for (; pesce < np; pesce++)
	{ 
		// Somma vector alla riga pesce di matrix
		addMatriceVettore(matrix, vector, pesce, d);
	}
}// addMatriceVettoreBroad

// prodTrasMatVet
void prodTrasMatVet(MATRIX matrix, VECTOR vector, VECTOR ris, int righe, int dim)
{
	// Alloca memoria per un vettore parz
	VECTOR parz = get_block(sizeof(type), dim);
	int unroll = 4;

	// Inizializza il vettore ris a zero
	for (int i = 0; i < dim; i++)
	{
		ris[i] = 0;
	}

	int pesce = 0;
	// Verifica se è possibile ottimizzare il calcolo con allineamento perfetto
	if (allineamentoPerfetto)
	{
		// Loop ottimizzato per allineamento perfetto
		for (pesce = 0; pesce + unroll < righe; pesce += 4)
		{
			// Calcola il prodotto tra matrix e vector[pesce] e memorizza il risultato in parz
			// matrix + pesce * dim sposta il puntatore matrix di pesce * dim elementi nella memoria
			prodVet_x_Scalare(matrix + pesce * dim, vector[pesce], parz, dim);
			// Somma parz a ris
			addVettori(ris, parz, ris, dim);
			//Come sopra per pesce +1, +2 e +3
			prodVet_x_Scalare(matrix + (pesce + 1) * dim, vector[pesce + 1], parz, dim);
			addVettori(ris, parz, ris, dim);
			prodVet_x_Scalare(matrix + (pesce + 2) * dim, vector[pesce + 2], parz, dim);
			addVettori(ris, parz, ris, dim);
			prodVet_x_Scalare(matrix + (pesce + 3) * dim, vector[pesce + 3], parz, dim);
			addVettori(ris, parz, ris, dim);
		}
	}
	else
	{
		// Loop per allineamento non perfetto
		for (pesce = 0; pesce + unroll < righe; pesce += 4)
		{

			if ((pesce * dim % allineamento) == 0)
			{
				// Calcola il prodotto con allineamento per matrix[pesce] e memorizza il risultato in parz
				prodVet_x_Scalare(matrix + pesce * dim, vector[pesce], parz, dim);
				// Somma parz a ris
				addVettori(ris, parz, ris, dim);
				//Come sopra per pesce +1, +2 e +3
				prodVet_x_ScalareUn(matrix + (pesce + 1) * dim, vector[pesce + 1], parz, dim);
				addVettori(ris, parz, ris, dim);
				prodVet_x_ScalareUn(matrix + (pesce + 2) * dim, vector[pesce + 2], parz, dim);
				addVettori(ris, parz, ris, dim);
				prodVet_x_ScalareUn(matrix + (pesce + 3) * dim, vector[pesce + 3], parz, dim);
				addVettori(ris, parz, ris, dim);
			}
			else
			{
				// Calcola il prodotto senza allineamento per matrix[pesce] e memorizza il risultato in parz
				prodVet_x_ScalareUn(matrix + pesce * dim, vector[pesce], parz, dim);
				// Somma parz a ris
				addVettori(ris, parz, ris, dim);
				//Come sopra per pesce +1, +2 e +3
				prodVet_x_ScalareUn(matrix + (pesce + 1) * dim, vector[pesce + 1], parz, dim);
				addVettori(ris, parz, ris, dim);
				prodVet_x_ScalareUn(matrix + (pesce + 2) * dim, vector[pesce + 2], parz, dim);
				addVettori(ris, parz, ris, dim);
				prodVet_x_ScalareUn(matrix + (pesce + 3) * dim, vector[pesce + 3], parz, dim);
				addVettori(ris, parz, ris, dim);
			}
		}
	}
	for (; pesce < righe; pesce++)
	{
		if ((pesce * dim % allineamento) == 0)
		{
			// Calcola il prodotto con allineamento per matrix[pesce] e memorizza il risultato in parz
			prodVet_x_Scalare(matrix + pesce * dim, vector[pesce], parz, dim);
		}
		else
		{
			// Calcola il prodotto senza allineamento per matrix[pesce] e memorizza il risultato in parz
			prodVet_x_ScalareUn(matrix + pesce * dim, vector[pesce], parz, dim);
		}

		// Somma parz a ris
		addVettori(ris, parz, ris, dim);
	}

	// Libera la memoria allocata per parz
	free_block(parz);
}// prodTrasMatVet

// minimoVettore
type minimoVettore(VECTOR vector, int dim)
{

	// Dichiara una variabile index per memorizzare l'indice dell'elemento minimo del vettore.
	int index = 0;

	// Inizializza la variabile min con il valore del primo elemento del vettore.
	type min = vector[0];

	// Cicla su tutti gli elementi del vettore.
	for (int i = 0; i < dim; i++)
	{

		// Se l'elemento corrente del vettore è minore del valore minimo attuale,
		// allora aggiorna l'indice dell'elemento minimo e il valore minimo.
		if (vector[i] < min)
		{
			index = i;
			min = (vector[i]);
		}
	}

	// Restituisce il valore minimo del vettore.
	return vector[index];
} // minimoVettore

// zeroRowMatrix
void zeroRowMatrix(MATRIX matrix, int riga, int dim)
{

	// Azzera tutti gli elementi della riga di matrice specificata da riga.
	for (int i = 0; i < dim; i++)
	{
		matrix[riga * dim + i] = 0;
	} // for
} // zeroRowMatrix

// replaceMatrixRowVector
void replaceMatrixRowVector(MATRIX matrix, VECTOR vector, int riga, int dim)
{

	// Sostituisce la riga di matrice specificata da riga con il vettore vector.
	for (int i = 0; i < dim; i++)
	{
		matrix[riga * dim + i] = vector[i];
	} // for
} // replaceMatrixRowVector

// generaPosizioneCasuale
void generaPosizioneCasuale(VECTOR pos, MATRIX posPesci, int pesce, int dim, params *input, var *vars, int it)
{
	for (int i = 0; i < dim; i++)
	{
		// Calcola la nuova posizione del pesce.
		// yi (j) = xi (j) + rand(−1, 1) · step ind
		pos[i] = posPesci[pesce * dim + i] + (getRand(input, it, pesce, i, 0) * 2 - 1) * input->stepind;
	} // for
}// generaPosizioneCasuale

///////////

// movimentoIndividuale
void movimentoIndividuale(params *input, var *vars, int pesce, int it)
{
	// Dichiara una variabile newPosition di tipo VECTOR e alloca in memoria d elementi di tipo type.
	VECTOR newPosition = get_block(sizeof(type), d);

	// Genera una posizione casuale per il pesce.
	generaPosizioneCasuale(newPosition, input->x, pesce, d, input, vars, it);

	// Calcola la variazione della funzione deltaf tra la nuova posizione e la posizione corrente.
	// ∆fi = f (yi) − f (xi)
	type deltaf = funzioneObiettivo(newPosition, input, d) - funzioneMatrix(input->x, input, pesce * d, d);

	// Se la variazione della funzione deltaf è negativa, il pesce si sposta alla nuova posizione.
	if (deltaf < 0)
	{

		// Segna che il pesce si è mosso.
		effettuato = 1;

		// Aggiorna la variazione della funzione deltaf del pesce.
		vars->deltaf[pesce] = deltaf;

		// Copia la posizione corrente del pesce in una variabile temporanea currentFishPosition.
		VECTOR currentFishPosition = copyVector(input->x, pesce * d, d);

		// Calcola lo spostamento del pesce.
		VECTOR deltaCurrentFishPosition = get_block(sizeof(type), d);
		// ∆xi = yi − xi
		subVettori(newPosition, currentFishPosition, deltaCurrentFishPosition, d);

		// Aggiorna la matrice deltax con lo spostamento del pesce.
		replaceMatrixRowVector(vars->deltax, deltaCurrentFishPosition, pesce, d);

		// Aggiorna la matrice x con la nuova posizione del pesce.
		replaceMatrixRowVector(input->x, newPosition, pesce, d);

		// Se l'indice di inizio della posizione del pesce non è un multiplo dell'allineamento, libera la memoria allocata per la variabile temporanea currentFishPosition.
		if ((pesce * d % allineamento) != 0)
			free_block(currentFishPosition);

		// Libera la memoria allocata per la variabile deltaCurrentFishPosition.
		free_block(deltaCurrentFishPosition);
	} // if

	// Altrimenti, il pesce rimane nella sua posizione corrente e la sua variazione della funzione deltaf viene azzerata.
	else
	{
		// ∆xi = yi − xi = 0
		zeroRowMatrix(vars->deltax, pesce, d);
		// ∆fi = f (yi) − f (xi) = 0
		vars->deltaf[pesce] = 0;
	} // else

	// Libera la memoria allocata per la variabile newPosition.
	free_block(newPosition);

} // movimentoIndividuale

// alimentazione
void alimentazione(params *input, var *vars)
{

	// Calcola il peso totale dei pesci.
	vars->fishesWeight = pesoTot(vars->w, np);

	// Se è stato effettuato almeno uno spostamento, allora si alimenta la popolazione.
	if (effettuato)
	{

		// Trova il valore massimo della variazione della funzione.
		// Ovvero il minimo cioè migliore valore della f dell'algoritmo FSS.
		type max = -minimoVettore(vars->deltaf, np);

		// Se il valore massimo della variazione della funzione è maggiore di una certa soglia, allora si alimentano i pesci.
		if (max > EPSILON)
		{

			// Alloca un vettore temporaneo ris.
			VECTOR ris = get_block(sizeof(type), np);

			// Calcola il prodotto scalare tra un vettore e un numero scalare, calcola quindi la normalizzazione della variazione della funzione.
			//(∆fi/maxj(∆fj))
			prodVet_x_Scalare(vars->deltaf, (type)1.0 / max, ris, np);

			// Aggiorna i pesi dei pesci.
			// Wi = Wi + (∆fi/maxj(∆fj))
			subVettori(vars->w, ris, vars->w, np);

			// Libera la memoria allocata per il vettore ris.
			free_block(ris);
		} // if
	}	  // if
} // alimentazione

// movimentoIstintivo
// i pesci che hanno incontrato un maggiore miglioramento attirano i pesci nella propria posizione
void movimentoIstintivo(params *input, var *vars)
{

	// Verifica se almeno un pesce si è spostato
	if (effettuato)
	{

		// Alloca memoria per il vettore I di dimensione d
		VECTOR I = get_block(sizeof(type), d);

		// Alloca memoria per un vettore num di dimensione d
		VECTOR num = get_block(sizeof(type), d);

		// Calcola il prodotto tra il vettore deltax, deltaf e mette il risultato in num
		// sommatoria per i che va da 1 a np di: ∆x · ∆fi --> numeratore
		prodTrasMatVet(vars->deltax, vars->deltaf, num, np, d);

		// Calcola il reciproco del peso totale di deltaf e memorizzalo in denom
		// 1 / sommatoria per i che va da 1 a np di: ∆fi --> denominatore
		type denom = (type)1.0 / pesoTot(vars->deltaf, np);

		// Verifica se denom è maggiore di EPSILON o minore di -EPSILON
		if (denom > EPSILON || denom < -EPSILON)
		{

			// Moltiplica num per denom e memorizza il risultato in I
			// I = num*denom
			prodVet_x_Scalare(num, denom, I, d);

			// Aggiungi il vettore I alla matrice input->x
			// xi = xi + I
			addMatriceVettoreBroad(input->x, I, d);
		}

		// Libera la memoria allocata per I
		free_block(I);

		// Libera la memoria allocata per num
		free_block(num);
	}
} // mov istintivo

// baricentro
void baricentro(params *input, var *vars)
{
	// Calcola il reciproco del peso totale di w e lo memorizza in denom
	// 1 / sommatoria per i che va da 1 a np di: Wi
	type denom = (type)1.0 / pesoTot(vars->w, np);

	// Calcola il prodotto tra la matrice x e il vettore w e memorizza il risultato in baricentro
	// sommatoria per i che va da 1 a np di: xi · Wi
	prodTrasMatVet(input->x, vars->w, vars->baricentro, np, d);

	// Moltiplica baricentro per denom e memorizza il risultato in baricentro
	// baricentro = baricento * denom --> B = (sommatoria per i che va da 1 a np di: xi · Wi) / (sommatoria per i che va da 1 a np di: Wi)
	prodVet_x_Scalare(vars->baricentro, denom, vars->baricentro, d);
} // baricentro

// movimentoVolitivo
void movimentoVolitivo(params *input, var *vars, int it)
{
	// Inizializza una variabile segno a 1
	int segno = 1;

	// Verifica se il peso totale di w è maggiore del valore precedentemente assunto, salvato in fishesWeight.
	if (pesoTot(vars->w, np) > vars->fishesWeight)
	{

		// Se la condizione è vera, cambia il segno a -1
		segno = -1;
		// printf(pesoAumentato.);
	}

	// Alloca memoria per un vettore xiMinusB di dimensione d
	MATRIX xiMinusBMatrix = alloc_matrix(np, d);

	// Alloca memoria per un vettore movVolitivo di dimensione d
	MATRIX movVolitivoMatrix = alloc_matrix(np, d);
	
	// Loop che itera su ciascun pesce
	#pragma omp parallel for 
	for (int pesce = 0; pesce < np; pesce++)
	{
		VECTOR xiMinusB = copyVector(xiMinusBMatrix, pesce, d);
		VECTOR movVolitivo = copyVector(movVolitivoMatrix, pesce, d);

		// Genera un numero casuale
		type rnd = getRand(input, it, pesce, d, 1);

		// Copia il vettore x_i dall'input
		VECTOR x_i = copyVector(input->x, pesce * d, d);

		// Calcola la distanza euclidea tra x_i e baricentro
		// dist(xi , B)
		type distanza = distEuclidea(x_i, vars->baricentro, d);

		// Calcola (segno · stepvol · rand(0, 1)) / dist(xi , B)
		type scalare = (input->stepvol * rnd * (type)segno) / distanza;

		// Sottrai baricentro da x_i e memorizza il risultato in xiMinusB
		// xi − B
		subVettori(x_i, vars->baricentro, xiMinusB, d);

		// Calcola il prodotto tra xiMinusB e scalare e memorizza il risultato in movVolitivo
		//((segno · stepvol · rand(0, 1)) · (xi − B))/ dist(xi , B)
		prodVet_x_Scalare(xiMinusB, scalare, movVolitivo, d);

		// Aggiungi movVolitivo alla matrice input->x per il pesce corrente
		// xi + movVolitivo
		addMatriceVettore(input->x, movVolitivo, pesce, d);
			
		// Se la riga corrente non è allineata, libera la memoria allocata per x_i
		if ((pesce * d % allineamento) != 0)
		{
			free_block(x_i);
		}	
	} // for

	// Libera la memoria allocata per xiMinusB
	dealloc_matrix(xiMinusBMatrix);

	// Libera la memoria allocata per movVolitivo
	dealloc_matrix(movVolitivoMatrix);
} // movimentoVolitivo

// minimo
void minimo(params *input)
{
	// Inizializza il valore minimo con il valore della funzione obiettivo per il primo vettore contenuto nella matrice
	type valore_minimo = funzioneMatrix(input->x, input, 0, d);
	// Inizializza l'indice associato al valore minimo
	int index = 0;

	// Loop per iterare attraverso tutti i vettori
	for (int i = 0; i < np; i++)
	{
		// Calcola il valore della funzione obiettivo per il vettore corrente
		type valore_tmp = funzioneMatrix(input->x, input, i * d, d);

		// Verifica se il valore corrente è minore del valore minimo
		if (valore_tmp < valore_minimo)
		{
			// Se sì, aggiorna il valore minimo e l'indice associato
			valore_minimo = valore_tmp;
			index = i;
		}
	}

	// Copia il vettore con il valore minimo nell'array xh
	for (int i = 0; i < d; i++)
	{
		input->xh[i] = input->x[index * d + i];
	}
} // minimo

// init
void init(params *input, var *vars)
{
	// Assegna il valore di d dalla struttura input a una variabile globale d
	d = input->d;

	// Assegna il valore di np dalla struttura input a una variabile globale np
	np = input->np;

	// Verifica se d è un multiplo di allineamento e imposta allineamentoPerfetto di conseguenza
	if (d % allineamento == 0)
	{
		allineamentoPerfetto = 1;
	}

	// Assegna il valore di iter dalla struttura input a una variabile globale iter
	iter = input->iter;

	// Alloca memoria per il vettore w con dimensione np
	vars->w = get_block(sizeof(type), np);

	// Alloca memoria per la matrice deltax con dimensioni np x d
	vars->deltax = alloc_matrix(np, d);

	// Alloca memoria per il vettore deltaf con dimensione np
	vars->deltaf = get_block(sizeof(type), np);

	// Assegna il valore di stepind dalla struttura input a stepindIni
	vars->stepindIni = input->stepind;

	// Assegna il valore di stepvol dalla struttura input a stepvolIni
	vars->stepvolIni = input->stepvol;

	// Alloca memoria per il vettore baricentro con dimensione d
	vars->baricentro = get_block(sizeof(type), input->d);

	// Alloca memoria per il vettore xh con dimensione d
	input->xh = get_block(sizeof(type), input->d);

	// Inizializza la variabile rand di vars a 0
	vars->rand = 0;

	// Inizializza la variabile fishesWeight di vars a 0
	vars->fishesWeight = 0;

	// Inizializza il vettore dei pesi w con wscale/2 e calcola il peso totale fishesWeight
	for (int i = 0; i < input->np; i++)
	{
		vars->w[i] = input->wscale / 2;
		vars->fishesWeight += vars->w[i];
	}

	// Inizializza i vettori baricentro e xh a 0
	for (int i = 0; i < input->d; i++)
	{
		vars->baricentro[i] = 0;
		input->xh[i] = 0;
	}
} // init

// aggiornaParametri
void aggiornaParametri(params *input, var *vars)
{
	// Aggiorna il valore di stepind: step ind - (step ind (initial) / Itmax )
	input->stepind = input->stepind - (vars->stepindIni / iter);

	// Aggiorna il valore di stepvol: step vol - (step vol (initial) / Itmax )
	input->stepvol = input->stepvol - (vars->stepvolIni / iter);
} // aggiornaParametri

// fss
void fss(params *input)
{
	int it = 0;
	var *vars = get_block(sizeof(var), 1);

	// Inizializza la posizione xi con d valori casuali e il peso Wi a Wscale / 2
	init(input, vars);

	// Loop principale delle iterazioni
	while (it < iter)
	{

		// Inizializza effettuato a 0, la variabile che tiene traccia se un pesce ha effettuato un movimento.
		effettuato = 0;

		// Loop per eseguire il movimento individuale per ciascun pesce
		#pragma omp parallel for
		for (int pesce = 0; pesce < np; pesce++)
		{
			// Esegue il movimento individuale per il pesce corrente
			movimentoIndividuale(input, vars, pesce, it);
		}

		// Applica l’operatore di alimentazione
		alimentazione(input, vars);
		// Esegue il movimento istintivo
		movimentoIstintivo(input, vars);
		// Calcola il baricentro
		baricentro(input, vars);
		// Esegue il movimento volitivo
		movimentoVolitivo(input, vars, it);
		// Aggiorna i parametri
		aggiornaParametri(input, vars);
		it += 1;
	}

	// Restituisce la posizione del pesce associata al minimo valore di f
	minimo(input);
} // fss

int main(int argc, char **argv)
{

	char fname[256];
	char *coefffilename = NULL;
	char *randfilename = NULL;
	char *xfilename = NULL;
	int i, j, k;
	clock_t t;
	float time;

	//
	// Imposta i valori di default dei parametri
	//

	params *input = malloc(sizeof(params));

	input->x = NULL;
	input->xh = NULL;
	input->c = NULL;
	input->r = NULL;
	input->np = 25;
	input->d = 2;
	input->iter = 350;
	input->stepind = 1;
	input->stepvol = 0.1;
	input->wscale = 10;

	input->silent = 0;
	input->display = 0;

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if (argc <= 1)
	{
		printf("%s -c <c> -r <r> -x <x> -np <np> -si <stepind> -sv <stepvol> -w <wscale> -it <itmax> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tc: il nome del file ds2 contenente i coefficienti\n");
		printf("\tr: il nome del file ds2 contenente i numeri casuali\n");
		printf("\tx: il nome del file ds2 contenente le posizioni iniziali dei pesci\n");
		printf("\tnp: il numero di pesci, default 25\n");
		printf("\tstepind: valore iniziale del parametro per il movimento individuale, default 1\n");
		printf("\tstepvol: valore iniziale del parametro per il movimento volitivo, default 0.1\n");
		printf("\twscale: valore iniziale del peso, default 10\n");
		printf("\titmax: numero di iterazioni, default 350\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc)
	{
		if (strcmp(argv[par], "-s") == 0)
		{
			input->silent = 1;
			par++;
		}
		else if (strcmp(argv[par], "-d") == 0)
		{
			input->display = 1;
			par++;
		}
		else if (strcmp(argv[par], "-c") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing coefficient file name!\n");
				exit(1);
			}
			coefffilename = argv[par];
			par++;
		}
		else if (strcmp(argv[par], "-r") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing random numbers file name!\n");
				exit(1);
			}
			randfilename = argv[par];
			par++;
		}
		else if (strcmp(argv[par], "-x") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing initial fish position file name!\n");
				exit(1);
			}
			xfilename = argv[par];
			par++;
		}
		else if (strcmp(argv[par], "-np") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing np value!\n");
				exit(1);
			}
			input->np = atoi(argv[par]);
			par++;
		}
		else if (strcmp(argv[par], "-si") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing stepind value!\n");
				exit(1);
			}
			input->stepind = atof(argv[par]);
			par++;
		}
		else if (strcmp(argv[par], "-sv") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing stepvol value!\n");
				exit(1);
			}
			input->stepvol = atof(argv[par]);
			par++;
		}
		else if (strcmp(argv[par], "-w") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing wscale value!\n");
				exit(1);
			}
			input->wscale = atof(argv[par]);
			par++;
		}
		else if (strcmp(argv[par], "-it") == 0)
		{
			par++;
			if (par >= argc)
			{
				printf("Missing iter value!\n");
				exit(1);
			}
			input->iter = atoi(argv[par]);
			par++;
		}
		else
		{
			printf("WARNING: unrecognized parameter '%s'!\n", argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//

	if (coefffilename == NULL || strlen(coefffilename) == 0)
	{
		printf("Missing coefficient file name!\n");
		exit(1);
	}

	if (randfilename == NULL || strlen(randfilename) == 0)
	{
		printf("Missing random numbers file name!\n");
		exit(1);
	}

	if (xfilename == NULL || strlen(xfilename) == 0)
	{
		printf("Missing initial fish position file name!\n");
		exit(1);
	}

	int x, y;
	input->c = load_data(coefffilename, &input->d, &y);
	input->r = load_data(randfilename, &x, &y);
	input->x = load_data(xfilename, &x, &y);

	if (input->np < 0)
	{
		printf("Invalid value of np parameter!\n");
		exit(1);
	}

	if (input->stepind < 0)
	{
		printf("Invalid value of si parameter!\n");
		exit(1);
	}

	if (input->stepvol < 0)
	{
		printf("Invalid value of sv parameter!\n");
		exit(1);
	}

	if (input->wscale < 0)
	{
		printf("Invalid value of w parameter!\n");
		exit(1);
	}

	if (input->iter < 0)
	{
		printf("Invalid value of it parameter!\n");
		exit(1);
	}

	//
	// Visualizza il valore dei parametri
	//

	if (!input->silent)
	{
		printf("Coefficient file name: '%s'\n", coefffilename);
		printf("Random numbers file name: '%s'\n", randfilename);
		printf("Initial fish position file name: '%s'\n", xfilename);
		printf("Dimensions: %d\n", input->d);
		printf("Number of fishes [np]: %d\n", input->np);
		printf("Individual step [si]: %f\n", input->stepind);
		printf("Volitive step [sv]: %f\n", input->stepvol);
		printf("Weight scale [w]: %f\n", input->wscale);
		printf("Number of iterations [it]: %d\n", input->iter);
	}

	// COMMENTARE QUESTA RIGA!
	// prova(input);
	//

	//
	// Fish School Search
	//

	t = clock();
	fss(input);
	t = clock() - t;
	time = ((float)t) / CLOCKS_PER_SEC;

	if (!input->silent)
		printf("FSS time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato di xh
	//
	sprintf(fname, "xh32_%d_%d_%d.ds2", input->d, input->np, input->iter);
	save_data(fname, input->xh, 1, input->d);
	if (input->display)
	{
		if (input->xh == NULL)
			printf("xh: NULL\n");
		else
		{
			printf("xh: [");
			for (i = 0; i < input->d - 1; i++)
				printf("%f,", input->xh[i]);
			printf("%f]\n", input->xh[i]);
		}
	}

	if (!input->silent)
		printf("\nDone.\n");

	return 0;
}
