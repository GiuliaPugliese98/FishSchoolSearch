; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf64 regression64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

section .bss			; Sezione contenente dati non inizializzati

alignb 32
stepind		resq		1
alignb 32
de		resq		1
alignb 32
risSommaEu	resq		1
alignb 32
prova1	resq		4
alignb 32
cv	resq		1

section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

global distEuclidea

distEuclidea:
    push    rbp                 ; Salva il Base Pointer
    mov     rbp, rsp            ; Il Base Pointer punta al Record di Attivazione corrente
    pushaq_distEucl             ; Salva i registri generali

    							; Inizializza i registri YMM1, RAX e RBX
    VXORPD YMM1, YMM1           ; Inizializza YMM1 a zero
    XOR RAX, RAX                ; Inizializza RAX a zero
    XOR RBX, RBX                ; Inizializza RBX a zero

cicloSommaEuclidea:
								; RDX = quanti elementi devono ancora essere processati
    SUB RDX, 16                 ; Sottrae 16 da RDX
    CMP RDX, 0                  ; Confronta RDX con zero
    JL sommaEuclidea_mezzi      ; Salta a sommaEuclidea_mezzi se RDX è negativo
    							; Carica i vettori v1 e v2 nei registri YMM0, sottrae e moltiplica
    VMOVAPD YMM0, [RDI + RAX*8]
    VSUBPD YMM0, [RSI + RAX*8]
    VMULPD YMM0, YMM0
    VADDPD YMM1, YMM0			; Somma YMM0 in YMM1
								
								; RAX = quanti elementi ho processato
								; Ripeti per i successivi 3 blocchi di 32 byte
    ADD RAX, 16                 ; Aggiunge 16 a RAX
    JMP cicloSommaEuclidea      ; Salta di nuovo a cicloSommaEuclidea

sommaEuclidea_mezzi:
    ADD RDX, 16                 ; Aggiunge 16 a RDX

cicloSommaEuclidea_mezzi:
    SUB RDX, 8                  ; Sottrae 8 da RDX
    CMP RDX, 0                  ; Confronta RDX con zero
    JL fineSommaEuclidea        ; Salta a fineSommaEuclidea se RDX è negativo
								; Carica i vettori v1 e v2 nei registri YMM0, sottrae e moltiplica
    VMOVAPD YMM0, [RDI + RAX*8]
    VSUBPD YMM0, [RSI + RAX*8]
    VMULPD YMM0, YMM0
    VADDPD YMM1, YMM0
    ADD RAX, 8                   ; Aggiunge 8 a RAX
    JMP cicloSommaEuclidea_mezzi ; Salta di nuovo a cicloSommaEuclidea_mezzi

fineSommaEuclidea:
    ADD RDX, 7                   ; Aggiunge 7 a RDX
    VXORPD YMM3, YMM3           ; Inizializza YMM3 a zero
    VXORPD YMM0, YMM0           ; Inizializza YMM0 a zero

cicloFineSommaEuclidea:
    CMP RDX, 0                  ; Confronta RDX con zero
    JL e1                        ; Salta a e1 se RDX è negativo
    ; Carica il valore di v2 in XMM2, sottrae il valore di v1, moltiplica e somma
    VMOVSD XMM2, [RSI + RAX]
    VSUBSD XMM2, [RDI + RAX]
    VMULSD XMM2, XMM2, XMM2
    VADDSD XMM3, XMM2
    SUB RDX, 1                   ; Sottrae 1 da RDX
    ADD RAX, 8                   ; Aggiunge 8 a RAX
    JMP cicloFineSommaEuclidea   ; Salta di nuovo a cicloFineSommaEuclidea

e1:
    VHADDPD YMM1, YMM1           ; Sposta dati da 128-bit inferiori a superiori in YMM1    
    VEXTRACTF128 XMM0, YMM1, 1b  ; Estrae i 128 bit superiori di YMM1 in XMM0
    VADDPD XMM0, XMM1            ; Somma i 128 bit inferiori di YMM1
    VADDPD XMM0, XMM3            ; Aggiunge il risultato precedente alla somma iniziale
    VSQRTSD XMM0, XMM0           ; Calcola la radice quadrata di XMM0
	VMOVSD [risSommaEu], XMM0    ; Salva il risultato finale in memoria

popaq_Eucl:
    							; Ripristina i registri generali
    mov rsp, rbp                ; Ripristina lo Stack Pointer
    pop rbp                     ; Ripristina il Base Pointer
    VMOVSD XMM0, [risSommaEu]
    ret                         ; Ritorna dal sottoprogramma

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global pesoTot

pesoTot:
    push        rbp                 ; Salva il Base Pointer
    mov         rbp, rsp            ; Il Base Pointer punta al Record di Attivazione corrente
    pushaq_pesoTot                  ; Salva i registri generali

    VXORPD YMM1, YMM1               ; Inizializza YMM1 a zero (somma cumulativa)
    XOR RAX, RAX                    ; Inizializza RAX a zero (indice per l'accesso ai dati)

cicloPesoTot:
    SUB RSI, 16                     ; Decrementa RSI (contatore) di 16
    CMP RSI, 0                      ; Confronta RSI con zero
    JL ciclopeso_mezzi              ; Salta a ciclopeso_mezzi se RSI è negativo
    VADDPD YMM1, [RDI + RAX*8]      ; Somma il blocco di 32 byte (4 double) a YMM1
    VADDPD YMM1, [RDI + RAX*8 + 32] ; Continua la somma con il blocco successivo
    VADDPD YMM1, [RDI + RAX*8 + 64] ; ...
    VADDPD YMM1, [RDI + RAX*8 + 96] ; ...
    ADD RAX, 16                     ; Incrementa RAX di 16
    JMP cicloPesoTot                ; Salta di nuovo a cicloPesoTot

ciclopeso_mezzi:
    ADD RSI, 16                     ; Ripristina RSI aggiungendo 16

cicloPesoTot_mezzi:
    SUB RSI, 8                      ; Decrementa RSI di 8
    CMP RSI, 0                      ; Confronta RSI con zero
    JL finePesoTot                  ; Salta a finePesoTot se RSI è negativo
    VADDPD YMM1, [RDI + RAX*8]      ; Somma il blocco di 16 byte (2 double) a YMM1
    VADDPD YMM1, [RDI + RAX*8 + 32] ; Continua la somma con il blocco successivo
    ADD RAX, 8                      ; Incrementa RAX di 8
    JMP cicloPesoTot_mezzi          ; Salta di nuovo a cicloPesoTot_mezzi

finePesoTot:
    ADD RSI, 8                      ; Ripristina RSI aggiungendo 8
    VXORPD XMM2, XMM2               ; Inizializza XMM2 a zero (somma cumulativa dei restanti elementi)

cicloFinePesoTot:
    SUB RSI, 1                      ; Decrementa RSI di 1
    JL e2                           ; Salta a e2 se RSI è negativo
    VADDSD XMM2, [RDI + RAX*8]      ; Aggiungi il singolo double a XMM2
    ADD RAX, 1                      ; Incrementa RAX di 1
    JMP cicloFinePesoTot            ; Salta di nuovo a cicloFinePesoTot

e2:
    VHADDPD YMM1, YMM1              ; Somma orizzontale degli elementi di YMM1
    VEXTRACTF128 XMM0, YMM1, 1b     ; Estrae la parte superiore di YMM1 in XMM0
    VADDPD XMM0, XMM1               ; Somma la parte inferiore di YMM1
    VADDPD XMM0, XMM2               ; Aggiungi il risultato parziale singolo in XMM2
    VMOVSD [de], XMM0               ; Memorizza il risultato finale in memoria

    popaq_pesoTot                   ; Ripristina i registri generali
    mov rsp, rbp                    ; Ripristina lo Stack Pointer
    pop rbp                         ; Ripristina il Base Pointer
    VMOVSD XMM0, [de]               ; Carica il risultato finale in XMM0
    ret                             ; Ritorna dal sotto-programma
	
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global addVettori

addVettori:
    push    rbp                 ; Salva il Base Pointer
    mov     rbp, rsp            ; Il Base Pointer punta al Record di Attivazione corrente
    pushaq_add                  ; Salva i registri generali

    ; RDI = v1
    ; RSI = v2
    ; RDX = ris
    ; RCX = dim
    
    XOR     RAX, RAX            ; i = 0

cicloaddVettori:
    SUB     RCX, 16             ; dim = dim - 4
    CMP     RCX, 0              ; dim == 0?
    JL      lavorodimezzato     ; Salta a lavorodimezzato se RCX è negativo
    VMOVAPD YMM0, [RDI + RAX*8] ; YMM0=[V1, V1, V1, V1]
    VADDPD  YMM0, [RSI + RAX*8] ; YMM0=[V1+V2, V1+V2, V1+V2, V1+V2]
        
    VMOVAPD YMM1, [RDI + RAX*8+32]
    VADDPD  YMM1, [RSI + RAX*8+32]
    
    VMOVAPD YMM2, [RDI + RAX*8+64]
    VADDPD  YMM2, [RSI + RAX*8+64]
    
    VMOVAPD YMM3, [RDI + RAX*8+96]
    VADDPD  YMM3, [RSI + RAX*8+96]
    
    VMOVAPD [RDX + RAX*8], YMM0
    VMOVAPD [RDX + RAX*8+32], YMM1
    VMOVAPD [RDX + RAX*8+64], YMM2
    VMOVAPD [RDX + RAX*8+96], YMM3
    
    ADD     RAX, 16             ; Incrementa RAX di 16
    JMP     cicloaddVettori     ; Salta di nuovo a cicloaddVettori

lavorodimezzato:
    ADD     RCX, 16             ; Ripristina RCX aggiungendo 16

cicloaddVettorimezzi:
    SUB     RCX, 8              ; Diminuisce RCX di 8
    JL      fineaddVettori      ; Salta a fineaddVettori se RCX è negativo
    VMOVAPD YMM0, [RDI + RAX*8] 
    VADDPD  YMM0, [RSI + RAX*8]
    
    VMOVAPD YMM1, [RDI + RAX*8+32]
    VADDPD  YMM1, [RSI + RAX*8+32]

    VMOVAPD [RDX + RAX*8], YMM0
    VMOVAPD [RDX + RAX*8+32], YMM1
    ADD     RAX, 8              ; Incrementa RAX di 8
    JMP     cicloaddVettorimezzi; Salta di nuovo a cicloaddVettorimezzi

fineaddVettori:
    ADD     RCX, 7              

ciclofineaddVettori:
    CMP     RCX, 0               ; dim == 0?
    JL      e3
    VMOVSD  XMM1, [RDI + RAX]   ; YMM1=[0, 0, 0, V1]
    VADDSD  XMM1, [RSI + RAX]   ; YMM1=[0, 0, 0, V1+V2]
    VMOVSD  [RDX + RAX], XMM1   ; Metto il risultato in ris
    SUB     RCX, 1              ; Diminuisce RCX di 1
    ADD     RAX, 8              ; Incrementa RAX di 8
    JMP     ciclofineaddVettori ; Salta di nuovo a ciclofineaddVettori

e3:
    popaq_add                   ; Ripristina i registri generali
    mov     rsp, rbp            ; Ripristina lo Stack Pointer
    pop     rbp                 ; Ripristina il Base Pointer
    ret                         ; Ritorna dal sotto-programma

;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
global subVettori

subVettori:
    push    rbp                 ; Salva il Base Pointer
    mov     rbp, rsp            ; Il Base Pointer punta al Record di Attivazione corrente
    pushaq_sub                  ; Salva i registri generali

    ; RDI = v1
    ; RSI = v2
    ; RDX = ris
    ; RCX = dim
    
    XOR     RAX, RAX            ; i = 0

ciclosubVettori:
    SUB     RCX, 16             ; dim = dim - 4
    CMP     RCX, 0              ; dim == 0?
    JL      lavorodimezzatosub  ; Salta a lavorodimezzatosub se RCX è negativo
    VMOVAPD YMM0, [RDI + RAX*8] ; YMM0=[V1, V1, V1, V1]
    VSUBPD  YMM0, [RSI + RAX*8] ; YMM0=[V1-V2, V1-V2, V1-V2, V1-V2]
    
    VMOVAPD YMM1, [RDI + RAX*8+32]
    VSUBPD  YMM1, [RSI + RAX*8+32]
    
    VMOVAPD YMM2, [RDI + RAX*8+64]
    VSUBPD  YMM2, [RSI + RAX*8+64]
    
    VMOVAPD YMM3, [RDI + RAX*8+96]
    VSUBPD  YMM3, [RSI + RAX*8+96]
    
    VMOVAPD [RDX + RAX*8   ], YMM0
    VMOVAPD [RDX + RAX*8+32], YMM1
    VMOVAPD [RDX + RAX*8+64], YMM2
    VMOVAPD [RDX + RAX*8+96], YMM3
    ADD     RAX, 16             ; Incrementa RAX di 16
    JMP     ciclosubVettori     ; Salta di nuovo a ciclosubVettori

lavorodimezzatosub:
    ADD     RCX, 16             ; Ripristina RCX aggiungendo 16

ciclosubVettorimezzi:
    SUB     RCX, 8              ; Diminuisce RCX di 8
    JL      finesubVettori      ; Salta a finesubVettori se RCX è negativo
    VMOVAPD YMM0, [RDI + RAX*8] 
    VSUBPD  YMM0, [RSI + RAX*8]
    
    VMOVAPD YMM1, [RDI + RAX*8+32]
    VSUBPD  YMM1, [RSI + RAX*8+32]

    VMOVAPD [RDX + RAX*8   ], YMM0
    VMOVAPD [RDX + RAX*8+32], YMM1
    ADD     RAX, 8              ; Incrementa RAX di 8
    JMP     ciclosubVettorimezzi; Salta di nuovo a ciclosubVettorimezzi

finesubVettori:
    ADD     RCX, 7              

ciclofinesubVettori:
    CMP     RCX, 0               ; dim == 0?
    JL      e4
    VMOVSD  XMM1, [RDI + RAX]   ; YMM1=[0, 0, 0, V1]
    VSUBSD  XMM1, [RSI + RAX]   ; YMM1=[0, 0, 0, V1-V2]
    VMOVSD  [RDX + RAX], XMM1   ; Metto il risultato in ris
    SUB     RCX, 1              ; Diminuisce RCX di 1
    ADD     RAX, 8              ; Incrementa RAX di 8
    JMP     ciclofinesubVettori ; Salta di nuovo a ciclofinesubVettori

e4:
    popaq_sub                   ; Ripristina i registri generali
    mov     rsp, rbp            ; Ripristina lo Stack Pointer
    pop     rbp                 ; Ripristina il Base Pointer
    ret                         ; Ritorna dal sotto-programma

;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global prodVet_x_Scalare

prodVet_x_Scalare:
    push    rbp                 			; Salva il Base Pointer
    mov     rbp, rsp            			; Il Base Pointer punta al Record di Attivazione corrente
    pushaq_prod                 			; Salva i registri generali

    ; RDI = v1
    ; XMM0 = S
    ; RSI = ris
    ; RDX = dim
    
    XOR     RAX, RAX            			; i = 0
    VBROADCASTSD YMM0, XMM0     			; Copia il valore di XMM0 in tutti gli elementi di YMM0

cicloProdVet_x_Scalare:
    SUB     RDX, 16             			; DIM -= 16
    CMP     RDX, 0              			; DIM <= 0?
    JL      prod_x_Scal_lavoromezzi 		; Se sì, salta a prod_x_Scal_lavoromezzi
    VMOVAPD YMM1, [RDI + RAX*8] 			; YMM1 = [V1, V1, V1, V1]
    VMULPD  YMM1, YMM0           			; YMM1 = [V1*S, V1*S, V1*S, V1*S]
    
    VMOVAPD YMM2, [RDI + RAX*8+32]
    VMULPD  YMM2, YMM0
    
    VMOVAPD YMM3, [RDI + RAX*8+64]
    VMULPD  YMM3, YMM0
    
    VMOVAPD YMM4, [RDI + RAX*8+96]
    VMULPD  YMM4, YMM0
    
    VMOVAPD [RSI + RAX*8], YMM1 			; Rimette a posto in memoria  
    VMOVAPD [RSI + RAX*8+32], YMM2
    VMOVAPD [RSI + RAX*8+64], YMM3
    VMOVAPD [RSI + RAX*8+96], YMM4
    ADD     RAX, 16             			; Incrementa RAX di 16
    JMP     cicloProdVet_x_Scalare 			; Salta di nuovo a cicloProdVet_x_Scalare

prod_x_Scal_lavoromezzi:
    ADD     RDX, 16             			; Ripristina RDX aggiungendo 16

cicloProdVet_x_Scalare_mezzi:
    SUB     RDX, 8              			; DIM -= 8
    CMP     RDX, 0             	 			; DIM <= 0?
    JL      fineProdVet_x_Scalare 			; Se sì, salta a fineProdVet_x_Scalare
    VMOVAPD YMM1, [RDI + RAX*8]				; YMM1 = [V1*S, V1*S, V1*S, V1*S]
    VMULPD  YMM1, YMM0
    
    VMOVAPD YMM2, [RDI + RAX*8+32]
    VMULPD  YMM2, YMM0
    
    VMOVAPD [RSI + RAX*8   ], YMM1
    VMOVAPD [RSI + RAX*8+32], YMM2
	ADD     RAX, 8              			; Incrementa RAX di 8
    JMP     cicloProdVet_x_Scalare_mezzi 	; Salta di nuovo a cicloProdVet_x_Scalare_mezzi

fineProdVet_x_Scalare:
    ADD     RDX, 7              

ciclofineProdVet_x_Scalare:
    CMP     RDX, 0              			; DIM <= 0?
    JL      e5                  			; Se sì, salta ad e5
    VMOVSD  XMM1, [RDI + RAX*8] 			; XMM1 = [0, 0, 0, V1]
    VMULSD  XMM1, XMM0          			; XMM1 = [0*S, 0*S, 0*S, V1*S]
    VMOVSD  [RSI + RAX*8], XMM1 			; Rimette a posto in memoria
    SUB     RDX, 1              			; DIM--
    ADD     RAX, 1              			; i++
    JMP     ciclofineProdVet_x_Scalare 		; Salta di nuovo a ciclofineProdVet_x_Scalare

e5:
    popaq_prod                  			; Ripristina i registri generali
    mov     rsp, rbp            			; Ripristina lo Stack Pointer
    pop     rbp                 			; Ripristina il Base Pointer
    ret                         			; Ritorna dal sotto-programma

;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global prodScalare

prodScalare:
    push        rbp                         ; salva il Base Pointer
    mov         rbp, rsp                    ; il Base Pointer punta al Record di Attivazione corrente
    pushaq_prodSc                           ; salva i registri generali
    
    ; RDI  = v1
    ; RSI  = v2
    ; RDX = dim
    
    XOR         RAX, RAX                    ; Inizializza l'indice RAX a 0
    VXORPD      YMM1, YMM1                  ; Inizializza YMM1 per la somma dei prodotti a multipli
    VXORPD      YMM2, YMM2                  ; Inizializza YMM2 per la somma singolare dei restanti 
    
cicloProdScalare:
    SUB         RDX, 16                     ; DIM-16
    CMP         RDX, 0                      ; DIM <= 0?
    JL          prodScal_lavoromezzi        ; se sì, salta a prodScal_lavoromezzi
    VMOVAPD     YMM0, [RDI + RAX*8]         ; YMM0 = [V1, V1, V1, V1]
    VMULPD      YMM0, [RSI + RAX*8]         ; YMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
    VADDPD      YMM1, YMM0                  ; Somma i prodotti parziali
    
    VMOVAPD     YMM2, [RDI + RAX*8+32]
    VMULPD      YMM2, [RSI + RAX*8+32]
    VADDPD      YMM1, YMM2
    
    VMOVAPD     YMM3, [RDI + RAX*8+64]
    VMULPD      YMM3, [RSI + RAX*8+64]
    VADDPD      YMM1, YMM3
    
    VMOVAPD     YMM4, [RDI + RAX*8+96]
    VMULPD      YMM4, [RSI + RAX*8+96]
    VADDPD      YMM1, YMM4
    
    ADD         RAX, 16                     ; Incrementa l'indice
    JMP         cicloProdScalare

prodScal_lavoromezzi:
    ADD         RDX, 16                     ; Ripristina DIM per i rimanenti 
    
cicloProdScalare_mezzi:
    SUB         RDX, 8                      ; dim = dim - 8
    CMP         RDX, 0
    JL          fineProdScalare
    VMOVAPD     YMM0, [RDI + RAX*8]         ; YMM0 = [V1, V1, V1, V1]
    VMULPD      YMM0, [RSI + RAX*8]         ; YMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
    VADDPD      YMM1, YMM0                  ; Somma i prodotti parziali
    
    VMOVAPD     YMM2, [RDI + RAX*8+32]
    VMULPD      YMM2, [RSI + RAX*8+32]
    VADDPD      YMM1, YMM2
    
    ADD         RAX, 8                      ; Incrementa l'indice
    JMP         cicloProdScalare_mezzi

fineProdScalare:
    ADD         RDX, 7                       ; Ripristina DIM per gli ultimi 
    
ciclofineProdScalare:
    CMP         RDX, 0                       ; DIM <= 0?
    JL          e6                           ; se sì, salta a e6
    VMOVSD      XMM0, [RDI + RAX]            ; XMM0 = [0, 0, 0, V1]
    VMULSD      XMM0, [RSI + RAX]            ; XMM0 = [0*V2, 0*V2, 0*V2, V1*V2]
    VADDSD      XMM2, XMM0                   ; XMM2 = [0, 0, 0, +=V1*V2]
    SUB         RDX, 1                       ; DIM--
    ADD         RAX, 8                       ; Incrementa l'indice
    JMP         ciclofineProdScalare

e6:
    VHADDPD     YMM1, YMM1                   ; Somma orizzontale di YMM1
    VEXTRACTF128 XMM3, YMM1, 1b              ; Estrae i 128 bit superiori di YMM1 in XMM3
    VADDSD      XMM1, XMM3                   ; Somma i 128 bit inferiori di YMM1
    VADDSD      XMM2, XMM1                   ; Aggiunge il risultato precedente alla somma totale
    VMOVSD      [risProdScal], XMM2          ; Memorizza il risultato finale
    
    popaq_prodScal                           ; ripristina i registri generali
    mov         rsp, rbp                     ; ripristina lo Stack Pointer
    pop         rbp                          ; ripristina il Base Pointer
    VMOVSD      XMM0, [risProdScal]
    ret                                      ; torna alla funzione C chiamante
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global prodVet_x_ScalareUn

prodVet_x_ScalareUn:
    push        rbp                         	; salva il Base Pointer
    mov         rbp, rsp                    	; il Base Pointer punta al Record di Attivazione corrente
    pushaq_prodScUn                         	; salva i registri generali
    
    ; RDI = v1
    ; XMM0 = S
    ; RSI = ris
    ; RDX = dim
    
    XOR         RAX, RAX                    	; Inizializza l'indice RAX a 0
    VBROADCASTSD YMM0, XMM0                 	; Copia il valore di XMM0 in tutti gli elementi di YMM0
    
cicloProdVet_x_ScalareUn:
    SUB         RDX, 16                     	; DIM-16
    CMP         RDX, 0                      	; DIM <= 0?
    JL          prod_x_Scal_lavoromezziUn   	; se sì, salta a prod_x_Scal_lavoromezziUn
    VMOVUPD     YMM1, [RDI + RAX*8]         	; YMM1 = [V1, V1, V1, V1]
    VMULPD      YMM1, YMM0                  	; YMM1 = [V1*S, V1*S, V1*S, V1*S]
    
    VMOVUPD     YMM2, [RDI + RAX*8+32]
    VMULPD      YMM2, YMM0
    
    VMOVUPD     YMM3, [RDI + RAX*8+64]
    VMULPD      YMM3, YMM0
    
    VMOVUPD     YMM4, [RDI + RAX*8+96]
    VMULPD      YMM4, YMM0
    
    VMOVAPD     [RSI+RAX*8], YMM1         		; Rimetto a posto in memoria  
    VMOVAPD     [RSI+RAX*8+32], YMM2
    VMOVAPD     [RSI+RAX*8+64], YMM3
    VMOVAPD     [RSI+RAX*8+96], YMM4
    ADD         RAX, 16                   		; Incrementa l'indice
    JMP         cicloProdVet_x_ScalareUn

prod_x_Scal_lavoromezziUn:
    ADD         RDX, 16                      	; Ripristina DIM per i rimanenti 
    
cicloProdVet_x_Scalare_mezziUn:
    SUB         RDX, 8                       	; dim = dim - 8
    CMP         RDX, 0
    JL          fineProdVet_x_ScalareUn
    VMOVUPD     YMM1, [RDI + RAX*8]         	; YMM1 = [V1, V1, V1, V1]
    VMULPD      YMM1, YMM0                  	; YMM1 = [V1*S, V1*S, V1*S, V1*S]
    
    VMOVUPD     YMM2, [RDI + RAX*8+32]
    VMULPD      YMM2, YMM0
    
    VMOVAPD     [RSI+RAX*8], YMM1
    VMOVAPD     [RSI+RAX*8+32], YMM2
    ADD         RAX, 8                       	; Incrementa l'indice
    JMP         cicloProdVet_x_Scalare_mezziUn
    
fineProdVet_x_ScalareUn:
    ADD         RDX, 7                        	; Ripristina DIM per gli ultimi 
    
ciclofineProdVet_x_ScalareUn:
    CMP         RDX, 0                        	; DIM <= 0?
    JL          e5                            	; se sì, fine
    VMOVSD      XMM1, [RDI + RAX*8]          	; XMM1 = [0, 0, 0, V1]
    VMULSD      XMM1, XMM0                   	; XMM1 = [0*S, 0*S, 0*S, V1*S]
    VMOVSD      [RSI+RAX*8], XMM1            	; Rimette a posto in memoria
    SUB         RDX, 1                        	; dim--
    ADD         RAX, 1                        	; i++
    JMP         ciclofineProdVet_x_ScalareUn

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------