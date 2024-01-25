;---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
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
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati

section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	stepind		resd		1
	alignb 16
	retps		resd		1
	alignb 16
	ss1		resd		4
	alignb 16
	retpt		resd		1
    alignb 16
	retde		resd		1


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
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global addVettori
	
	dim equ 20       		; Dimensione del vettore in numero di elementi
	ris equ 16       		; Offset del vettore risultato
	v2 equ 12        		; Offset del secondo vettore di input
	v1 equ 8         		; Offset del primo vettore di input

addVettori:
	
    push		ebp
	mov		ebp, esp		; Imposta il Base Pointer al Record di Attivazione corrente
	push		ebx			; Salva i registri da preservare
	push		esi
	push		edi
	
	XOR ESI,ESI         	; Inizializza un registro (ESI) utilizzato come contatore degli elementi del vettore
    MOV EDI, [EBP+dim] 		; Carica la dimensione del vettore
    MOV EBX,[EBP+v1]   		; Carica l'indirizzo del primo vettore
    MOV ECX,[EBP+v2]   		; Carica l'indirizzo del secondo vettore
    MOV EDX,[EBP+ris]  		; Carica l'indirizzo del vettore risultato

    cicloAddVettori:  							
			SUB EDI,4       				;
			CMP EDI,0
			JL fineAddVettori
			MOVAPS XMM0,[EBX + 4*ESI] 			
			ADDPS XMM0,[ECX+4*ESI]    			
			MOVAPS [EDX+4*ESI],XMM0  
			ADD ESI,4
			JMP cicloAddVettori
            
	fineAddVettori: 
        ADD EDI,3           			; Ripristina la dimensione del vettore di 3
                    
	cicloFineVettori:   				; Ciclo finale per processare gli elementi rimanenti
			CMP EDI,0        			
			JL e1            			
			
			MOVSS XMM0,[EBX+4*ESI] 		; Carica un elemento dal primo vettore in XMM0
			ADDSS XMM0,[ECX+4*ESI] 		; Aggiunge un elemento dal secondo vettore a XMM0
			MOVSS [EDX+4*ESI],XMM0 		; Memorizza il risultato nel vettore di output
			
			SUB EDI,1
			ADD ESI,1
			JMP cicloFineVettori
			
		JMP cicloFineVettori
   
	e1:
    pop	edi			; Ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; Ripristina lo Stack Pointer
	pop	ebp			; Ripristina il Base Pointer
	ret				; Torna alla funzione C chiamante

;---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global subVettori
	
	dim equ 20       	; Dimensione del vettore in numero di elementi
	ris equ 16       	; Offset del vettore risultato
	v2 equ 12        	; Offset del secondo vettore di input
	v1 equ 8         	; Offset del primo vettore di input

subVettori:
	
    push		ebp
	mov		ebp, esp	; Imposta il Base Pointer al Record di Attivazione corrente
	push		ebx		; Salva i registri da preservare
	push		esi
	push		edi

    XOR ESI,ESI         ; Inizializza un registro (ESI) utilizzato come contatore degli elementi del vettore
    MOV EDI,[EBP+dim]   ; Carica la dimensione del vettore
    MOV EBX,[EBP+v1]    ; Carica l'indirizzo del primo vettore
    MOV ECX,[EBP+v2]    ; Carica l'indirizzo del secondo vettore
    MOV EDX,[EBP+ris]   ; Carica l'indirizzo del vettore risultato
    
    										
    cicloSubVettori:  					
			SUB EDI,4
			CMP EDI,0       
			JL fineSubVettori 					
			
			MOVAPS XMM0,[EBX + 4*ESI] 		; Carica 4 elementi dal primo vettore in XMM0
			SUBPS XMM0,[ECX+4*ESI]    		; Sottrae 4 elementi dal secondo vettore a XMM0
			MOVAPS [EDX+4*ESI],XMM0  		; Memorizza i risultati nel vettore di output
			
			ADD ESI,4	 	      	           
        		JMP cicloSubVettori
            
	fineSubVettori: 
		ADD EDI,3   						
                    
	cicloFineSubVettori:                    
			CMP EDI,0						
			JL e2							
			
			MOVSS XMM0,[EBX+4*ESI]  		; Carica 4 elementi dal primo vettore in XMM0
			SUBSS XMM0,[ECX+4*ESI]			; Sottrae 4 elementi dal secondo vettore a XMM0
			MOVSS [EDX+4*ESI],XMM0			; Memorizza i risultati nel vettore di output
			
			SUB EDI,1
			ADD ESI,1
		JMP cicloFineSubVettori
   
	e2:
    	pop	edi			; Ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; Ripristina lo Stack Pointer
		pop	ebp			; Ripristina il Base Pointer
		ret				; Torna alla funzione C chiamante

;-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global prodVet_x_Scalare

dimV equ 20       ; Dimensione del vettore in numero di elementi
ris equ 16        ; Offset del vettore risultato
s equ 12          ; Offset dello scalare
v equ 8           ; Offset del vettore

prodVet_x_Scalare:
    push    ebp
    mov     ebp, esp        ; Il Base Pointer punta al Record di Attivazione corrente
    push    ebx             ; Salva i registri da preservare
    push    esi
    push    edi

    XOR     ESI, ESI         ; i=0
    MOV     EDI, [EBP+dimV]  ; EDI = dim
    MOV     EBX, [EBP+ris]   ; EBX = RIS (è l'indirizzo di UN VETTORE)
    MOVSS   XMM0, [EBP+s]    ; XMM0 = S  [0,0,0,S] (è UNO SCALARE)
    MOV     EDX, [EBP+v]     ; EDX = V1 (VETTORE DI PARTENZA)

    XORPS   XMM0, XMM0
    SHUFPS  XMM0, XMM0, 00000000b 		; XMM0 = [S, S, S, S] crea un vettore SIMD in cui tutti gli elementi sono uguali al valore dello scalare S

cicloProdVet_x_Scalare:
    SUB     EDI,4           			; DIM-4
    CMP     EDI,0             			; DIM <= 0?
    JL      prodVet_x_scalare 			; Se sì, salta a prodVet_x_scalare

    MOVAPS  XMM1, [EDX + 4*ESI]   		; XMM1 = [V1, V1, V1, V1]
    MULPS   XMM1, XMM0            		; XMM1 = [V1*S, V1*S, V1*S, V1*S]
    MOVAPS  [EBX+ 4*ESI], XMM1    		; Memorizza il risultato nel vettore di output

    ADD     ESI,4
    JMP     cicloProdVet_x_Scalare

prodVet_x_scalare:
    ADD     EDI,3
        			

ciclofineProdVet_x_Scalare:
    CMP EDI,0            			
    JL      e3                 			; Se sì, salta a e3

    XORPS XMM1, XMM1         			; XMM1 = [0,0,0,0]
    MOVSS   XMM1, [EDX+4*ESI]    		; XMM1 = [0,0,0,v1[i]]
    MULSS   XMM1, XMM0           		; XMM1 = [0*S, 0*S, 0*S, v1[i]*S]
    MOVSS   [EBX+4*ESI], XMM1    		; Memorizza il risultato nel vettore di output
    
	ADD ESI,1                  		; i++
	JMP     ciclofineProdVet_x_Scalare

e3:
    pop     edi                 		; Ripristina i registri da preservare
    pop     esi
    pop     ebx
    mov     esp, ebp            		; Ripristina lo Stack Pointer
    pop     ebp                 		; Ripristina il Base Pointer
    ret                          		; Torna alla funzione C chiamante

;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global prodScalare

dimps equ 16
v2 equ 12
v1 equ 8

prodScalare:

    push    ebp
    mov     ebp, esp        	; Il Base Pointer punta al Record di Attivazione corrente
    push    ebx             	; Salva i registri da preservare
    push    esi
    push    edi

    XOR     ESI, ESI        	; i=0
    MOV     EDI, [EBP+dimps]    ; EDI = dim
    MOV     EBX, [EBP+v1]       ; EBX = V1
    MOV     ECX, [EBP+v2]       ; ECX = V2
    							; EDX=RIS  (è un float)

    XORPS   XMM1, XMM1       	; XMM1=0

cicloProdScalare:
    SUB   EDI, 4           			; DIM-4
    CMP   EDI, 0              			; DIM <= 0?
    JL    fineProdScalare			; Se si, salta a prodScalare2
    
    MOVAPS  XMM0, [EBX + 4*ESI]  		; XMM0 = [V1, V1, V1, V1]
    MULPS   XMM0, [ECX + 4*ESI]  		; XMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
    ADDPS   XMM1, XMM0            

    ADD     ESI,4 	             		; i= i+4
    JMP     cicloProdScalare

fineProdScalare:
    ADD     EDI,3                		
    XORPS   XMM3, XMM3    			; Inizializza a zero XMM3

ciclofineProdScalare:
    CMP   EDI, 0                		; dim <= 0?
    JL    e4                    		; Se si, salta a e4
    
    MOVSS   XMM0, [EBX+4*ESI]     		; XMM0 = [0, 0, 0, v1[i]]
    MULSS   XMM0, [ECX+4*ESI]     		; XMM0 = [0, 0, 0, v1[i]*v2[i]]
    ADDSS   XMM3, XMM0            		; XMM1 = XMM1 + v1[i]*v2[i]
    
    ADD ESI,1                   		; i++
    JMP     ciclofineProdScalare

e4:
    HADDPS  XMM1, XMM1          ; Somma orizzontalmente gli elementi adiacenti in XMM1
    HADDPS  XMM1, XMM1          ; Somma di nuovo orizzontalmente gli elementi in XMM1

    ADDSS   XMM1, XMM3          ; Aggiunge il risultato accumulato XMM3 a XMM1
    MOVSS   [retps], XMM1       ; Memorizza il risultato finale in memoria

    FLD     dword [retps]       ; Carica il risultato float nello stack di coprocessore
    pop     edi                 ; Ripristina i registri
    pop     esi                 
    pop     ebx                 
    mov     esp, ebp            ; Ripristina lo Stack Pointer
    pop     ebp                 ; Ripristina il Base Pointer
    ret                         ; Torna alla funzione C chiamante
	
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global pesoTot

dimp equ 12
v equ 8

pesoTot:
    push   ebp
    mov    ebp, esp
    push   ebx
    push   esi
    push   edi

    XOR    ESI, ESI               	; i=0
    MOV    EDI, [EBP+dimp]        	; EDI = dim
    MOV    EBX, [EBP+v]           	; EBX = v
    XORPS  XMM0, XMM0             	; XMM0 = [0,0,0,0]
    XORPS  XMM1, XMM1             	; XMM1 = [0,0,0,0]

cicloPesoTot:
    SUB   EDI,4                 	; dim=dim-4
    CMP   EDI,0                  	; dim<=0?
    JL    finePesoTot        			; Se si, vado a pesoTot2
    
    ADDPS  XMM0, [EBX + ESI*4]    	; somma gli elementi del vettore EBX a XMM0
    ADD    ESI,4                	; i=i+4
    JMP    cicloPesoTot

finePesoTot:
    ADD    EDI,3

ciclofinePesoTot:
    CMP    EDI, 0                  	; dim <=0?
    JL     e5                      	; Se si, ho finito
    
    ADDSS  XMM1, [EBX + ESI*4]    	; XMM1 = [0,0,0,v]
    ADD ESI,1                     	; i++
    JMP    ciclofinePesoTot

e5:
    HADDPS XMM0, XMM0
    HADDPS XMM0, XMM0             ; somma parziale (multipli di 4)
    ADDSS  XMM0, XMM1             ; Somma delle somme parziali
    MOVSS  [retpt], XMM0

    FLD    dword [retpt]
    pop    edi
    pop    esi
    pop    ebx
    mov    esp, ebp
    pop    ebp
    ret

;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global copyAlnVector

vcopy   equ   8             ; Offset del vettore di destinazione
inizio  equ   12            ; Offset del vettore di partenza
dimcopy equ   16            ; Dimensione del vettore in numero di elementi

copyAlnVectorAsm:
    push   ebp
    mov    ebp, esp          ; Il Base Pointer punta al Record di Attivazione corrente
    push   ebx               ; Salva i registri da preservare
    push   esi
    push   edi
    
    MOV    EDX, [EBP+inizio] 	; Carica l'offset del vettore di partenza in EDX
    AND    EDX, 3             	; Calcola il resto della divisione per 4
    JZ     modulo4           	; Se resto è zero, salta a modulo4
    
    ;prints msg              	; Stampa il messaggio 'stepind:'
    ;getmem 4, [EBP+dim]      	; Non so cosa faccia getmem, la sintassi sembra errata
    
    XOR    ESI, ESI          	; i=0
    MOV    EDI, [EBP+dimcopy] 	; EDI = dim
    MOV    ECX, [EBP+inizio]  	; Carica l'offset del vettore di partenza in ECX
    SHL    ECX, 2             	; Moltiplica l'offset per 4 (dimensione di un elemento nel vettore)
    MOV    EBX, [EBP+vcopy]    	; Carica l'offset del vettore di destinazione in EBX
    ADD    EBX, ECX           	; Calcola l'indirizzo effettivo del vettore di destinazione
    		
cicloquozientecopy:
    SUB    EDI, 4              		; dim=dim-4
    CMP EDI, 0                		; Controllo se dim<=0
    JL     finecicloquozientecopy 	; Se si, salta a finecicloquozientecopy
    
	MOVUPS XMM0, [EBX+ESI]    		; Carica 4 float da [EBX+ESI] a XMM0
    MOVAPS [EAX+ESI], XMM0    		; Memorizza i 4 float in XMM0 a [EAX+ESI]
    
	ADD    ESI, 16             		; Incrementa l'indice i di 4 (spostati 4 elementi alla volta)
    JMP    cicloquozientecopy 		; Ripete il ciclo
    
finecicloquozientecopy:
    ADD    EDI, 4              ; Ripristina dim aggiungendo 4
ciclorestocopy:
    SUB    EDI, 1              ; dim=dim-1
    CMP EDI, 0                ; Controllo se dim<=0
    JL     ecopy               ; Se si, salta a ecopy
    
	MOVSS  XMM0, [ESI+EAX]     ; Carica un float da [ESI+EAX] a XMM0
    MOVSS  [ESI+EBX], XMM0     ; Memorizza il float in XMM0 a [ESI+EBX]
    
	ADD    ESI, 4              ; Incrementa l'indice i di 1
    JMP    ciclorestocopy      ; Ripete il ciclo
    
modulo4:
    ;prints msg               ; Stampa il messaggio 'stepind:'
    MOV    ECX, [EBP+inizio]  ; Carica l'offset del vettore di partenza in ECX
    SHL    ECX, 2             ; Moltiplica l'offset per 4 (dimensione di un elemento nel vettore)
    MOV    EBX, [EBP+vcopy]   ; Carica l'offset del vettore di destinazione in EBX
    ADD    EBX, ECX           ; Calcola l'indirizzo effettivo del vettore di destinazione
    MOV    EAX, EBX           ; Copia l'indirizzo del vettore di destinazione in EAX
    
ecopy:
    pop    edi                ; Ripristina i registri da preservare
    pop    esi
    pop    ebx
    mov    esp, ebp           ; Ripristina lo Stack Pointer
    pop    ebp                ; Ripristina il Base Pointer
    ret                       ; Ritorna alla funzione C chiamante

;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

global distEuclidea
	dimeucl equ 16
       	v2 equ 12
       	v1 equ 8

distEuclidea: 
	push   ebp
    	mov    ebp, esp
    	push   ebx
    	push   esi
    	push   edi

	XOR ESI,ESI
    	MOV EDI, [EBP+dimeucl]
	MOV EAX, [EBP+v1]
    	MOV EBX, [EBP+v2]
    	XORPS XMM1, XMM1  		;registro per le somme a multipli di 4
    	
cicloDistEuclidea: 
	SUB EDI,4
    	CMP EDI,0
    	JL  fineDistEuclidea

    	MOVAPS XMM0,[EAX+ESI*4]      	;XMM0=[v1, v1, v1, v1]
	SUBPS XMM0, [EBX+ESI*4]        	;XMM0=[v1-v2, v1-v2, v1-v2, v1-v2]
	MULPS XMM0, XMM0  		;XMM0=[v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2, v1-v2*v1-v2]
	ADDPS XMM1, XMM0    		;sommo le differenze dei quadrati sulle celle di XMM1

    	ADD ESI,4
    	JMP cicloDistEuclidea

fineDistEuclidea:
    	XORPS XMM2, XMM2  		;registro per le somme restanti 
	ADD EDI,3

ciclofineDistEuclidea: 
	CMP EDI, 0
	JL e7
    	
    	MOVSS XMM0, [EAX+ESI*4]  	;XMM0=[0, 0, 0, v1]
	SUBSS XMM0, [EBX+ESI*4]      	;XMM0=[0, 0, 0, (v1-v2)]
	MULSS XMM0, XMM0              	;XMM0=[0, 0, 0, (v1-v2)*(v1-v2)]
	ADDSS XMM2, XMM0              	;XMM2=[0, 0, 0, (v1-v2)*(v1-v2)]
    	ADD ESI,1
    	JMP ciclofineDistEuclidea
e7:
	HADDPS XMM1, XMM1
	HADDPS XMM1,XMM1   		;somma di tutte le differenze dei quadrati dei multipli di 4
	ADDSS XMM1, XMM0    		;somma delle differenze dei quadrati rimasti
    	ADDSS XMM1, XMM2    
   	SQRTSS XMM1, XMM1 

    	MOVSS [retde], XMM1

	FLD dword [retde]
   	pop edi 
    	pop esi
    	pop ebx
    	mov esp, ebp 
    	pop  ebp 
    	ret
;_________________        




;--------------------------

global prodVet_x_ScalareUn

	dimV equ 20
	ris equ 16
	s equ 12
	v equ 8

prodVet_x_ScalareUn:
	push		ebp
	mov		ebp, esp		; il Base Pointer punta al Record di Attivazione corrente
	push		ebx			; salva i registri da preservare
	push		esi
	push		edi

	XOR	ESI,ESI				; i=0
	MOV	EDI,[EBP+dimV]			; EDI = dim
	MOV 	EBX,[EBP+ris]			; EBX = RIS (� L'indirizzo di UN VETTORE)
	MOVSS 	XMM0,[EBP+s]			; XMM0 = S  [0,0,0,S] (� UNO SCALARE)
	MOV 	EDX,[EBP+v]			; EDX =V1 (VETTORE DI PARTENZA)
	
	;XORPS	XMM0, XMM0
	SHUFPS 	XMM0, XMM0, 00000000b	; XMM0 = [S, S, S, S]

cicloProdVet_x_ScalareUn:
	SUB 		EDI,4					; DIM-4
	CMP 		EDI,0					; DIM <= 0?
	JL 		fineProdVet_x_ScalareUn	
	
    	MOVUPS 	XMM1,[EDX + 4*ESI]			; XMM1 = [V1, V1, V1, V1]
	MULPS	XMM1, XMM0				; XMM1 = [V1*S, V1*S, V1*S, V1*S]
	MOVAPS 	[EBX+ 4*ESI], XMM1	

    	ADD 	ESI,4				
	JMP 	cicloProdVet_x_ScalareUn

fineProdVet_x_ScalareUn:
	ADD 	EDI,3					; dim=dim+4

ciclofineProdVet_x_ScalareUn:
	CMP 		EDI,0					; dim<=0?
	JL 		e3u						; se si ho finito e faccio le pop dei registri dallo stack
	;XORPS	XMM1, XMM1				; XMM1 = [0,0,0,0]
	MOVSS 	XMM1,[EDX+4*ESI]			; XMM1 = [0,0,0,v1[i]]
	MULSS 	XMM1, XMM0				; XMM1 = [0*S, 0*S, 0*S, v1[i]*S]
	MOVSS 	[EBX+4*ESI], XMM1
	ADD ESI,1						; i++
	JMP 		ciclofineProdVet_x_ScalareUn

e3u:
	pop	edi		; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp		; ripristina il Base Pointer
	ret			; torna alla funzione C chiamante
;--------------------------------------------------

