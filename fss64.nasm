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
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq_distEucl						; salva i registri generali
	
	;RDI=v1	
	;RSI=v2
	;RDX=dim
		
    VXORPD YMM1,YMM1
	XOR RAX,RAX
    XOR RBX,RBX
        
        
cicloSommaEuclidea:	
	SUB RDX,16
	CMP RDX,0
	JL sommaEuclidea_mezzi       
    VMOVAPD YMM0,[RDI + RAX*8]        
    VSUBPD YMM0,[RSI + RAX*8]
    VMULPD YMM0,YMM0
	VADDPD YMM1,YMM0
	
	VMOVAPD YMM2,[RDI + RAX*8+32]        
    VSUBPD YMM2,[RSI + RAX*8+32]
    VMULPD YMM2,YMM2
	VADDPD YMM1,YMM2
	
	VMOVAPD YMM3,[RDI + RAX*8+64]        
    VSUBPD YMM3,[RSI + RAX*8+64]
    VMULPD YMM3,YMM3
	VADDPD YMM1,YMM3
	
	VMOVAPD YMM4,[RDI + RAX*8+96]        
    VSUBPD YMM4,[RSI + RAX*8+96]
    VMULPD YMM4,YMM4
	VADDPD YMM1,YMM4
    ADD RAX,16
    JMP cicloSommaEuclidea
	
sommaEuclidea_mezzi:
	ADD RDX, 16
	
cicloSommaEuclidea_mezzi:
	SUB RDX,8
	CMP RDX,0
	JL fineSommaEuclidea
	VMOVAPD YMM0,[RDI + RAX*8]        
    VSUBPD YMM0,[RSI + RAX*8]
    VMULPD YMM0,YMM0
	VADDPD YMM1,YMM0
	
	VMOVAPD YMM2,[RDI + RAX*8+32]        
    VSUBPD YMM2,[RSI + RAX*8+32]
    VMULPD YMM2,YMM2
	VADDPD YMM1,YMM2
	ADD RAX, 8
	JMP cicloSommaEuclidea_mezzi
	
fineSommaEuclidea:
	ADD RDX,7
	VXORPD YMM3,YMM3
	VXORPD YMM0,YMM0
	
cicloFineSommaEuclidea:	
	CMP RDX,0
	JL e1
    VMOVSD XMM2,[RSI+RAX]
    VSUBSD XMM2,[RDI+RAX]
    VMULSD XMM2,XMM2,XMM2
	VADDSD XMM3,XMM2
	SUB RDX,1
	ADD RAX,8
	JMP cicloFineSommaEuclidea
	
e1:     
	VHADDPD YMM1,YMM1
    VEXTRACTF128 XMM0,YMM1,1b ;<-------
    VADDPD XMM0,XMM1
    VADDPD XMM0,XMM3
    VSQRTSD XMM0,XMM0
    VMOVSD [risSommaEu],XMM0 
       
      
	popaq_Eucl				; ripristina i registri generali
	mov		rsp, rbp	; ripristina lo Stack Pointer
	pop		rbp		; ripristina il Base Pointer
	VMOVSD XMM0,[risSommaEu] 
	ret
;_____________________________

global pesoTot 
	
	
pesoTot:
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq_pesoTot					; salva i registri generali
		
	VXORPD YMM1,YMM1
    XOR RAX,RAX
	
cicloPesoTot:	
	SUB RSI,16
	CMP RSI,0
	JL ciclopeso_mezzi      
    VADDPD YMM1,[RDI + RAX*8]

	VADDPD YMM1,[RDI + RAX*8+32]
		
	VADDPD YMM1,[RDI + RAX*8+64]
		
	VADDPD YMM1,[RDI + RAX*8+96]
		
    ADD RAX,16
    JMP cicloPesoTot
	
ciclopeso_mezzi:
	ADD RSI, 16
	
cicloPesoTot_mezzi:
	SUB RSI, 8 
	CMP RSI, 0 
	JL finePesoTot
	
	VADDPD YMM1,[RDI + RAX*8]

	VADDPD YMM1,[RDI + RAX*8+32]
	
	ADD RAX, 8
	JMP cicloPesoTot_mezzi
	
finePesoTot:
	ADD RSI,8
	VXORPD XMM2,XMM2
	
cicloFinePesoTot:	
	SUB RSI,1
	JL e2
    VADDSD XMM2,[RDI+RAX*8]	
	ADD RAX,1
	JMP cicloFinePesoTot
	
e2:     
	VHADDPD YMM1,YMM1
    VEXTRACTF128 XMM0,YMM1,1b ;<-------
    VADDPD XMM0,XMM1
    VADDPD XMM0,XMM2
    VMOVSD [de],XMM0 
       
      
	popaq_pesoTot				; ripristina i registri generali
	mov		rsp, rbp	; ripristina lo Stack Pointer
	pop		rbp		; ripristina il Base Pointer
	VMOVSD XMM0,[de] 
	ret
	
	
;---------------------------------
;cancellato copy
;--------------------
global addVettori

addVettori:
	push		rbp				
	mov			rbp, rsp			
	pushaq_add	
	
	;RDI=v1
	;RSI=v2
	;RDX=ris
	;RCX=dim
	
	XOR 			RAX,RAX			;i=0

cicloaddVettori:
	SUB 			RCX, 16			;dim=dim-4
	CMP 			RCX, 0			;dim==0?
	JL 			lavorodimezzato
	VMOVAPD 	YMM0, [RDI + RAX*8] 	;YMM0=[V1, V1, V1, V1]
	VADDPD 		YMM0, [RSI + RAX*8]	    ;YMM0=[V1+V2, V1+V2, V1+V2, V1+V2]
		
	VMOVAPD 	YMM1, [RDI + RAX*8+32]
	VADDPD 		YMM1, [RSI + RAX*8+32]	
	
	VMOVAPD 	YMM2, [RDI + RAX*8+64]
	VADDPD 		YMM2, [RSI + RAX*8+64]	
	
	VMOVAPD 	YMM3, [RDI + RAX*8+96]
	VADDPD 		YMM3, [RSI + RAX*8+96]	
	
	VMOVAPD 	[RDX + RAX*8], YMM0
	VMOVAPD 	[RDX + RAX*8+32], YMM1
	VMOVAPD 	[RDX + RAX*8+64], YMM2
	VMOVAPD 	[RDX + RAX*8+96], YMM3
	
	ADD 		RAX, 16			
	JMP 			cicloaddVettori
	
lavorodimezzato:
	ADD RCX, 16 

cicloaddVettorimezzi:
	SUB RCX, 8
	JL fineaddVettori
	VMOVAPD 	YMM0, [RDI + RAX*8] 	
	VADDPD 		YMM0, [RSI + RAX*8]	
	
	VMOVAPD 	YMM1, [RDI + RAX*8+32]
	VADDPD 		YMM1, [RSI + RAX*8+32]

	VMOVAPD 	[RDX + RAX*8], YMM0
	VMOVAPD 	[RDX + RAX*8+32], YMM1
	ADD RAX, 8
	JMP cicloaddVettorimezzi

fineaddVettori:
	ADD 		RCX, 7		
		
ciclofineaddVettori:
	CMP 			RCX, 0			;dim==0?
	JL 			e3
	VMOVSD 		XMM1, [RDI + RAX] 	;YMM1=[0, 0, 0, V1]
	VADDSD 		XMM1, [RSI + RAX]	;YMM1=[0, 0, 0, V1+V2]
	VMOVSD 		[RDX + RAX], XMM1  ;metto il risultato in ris
	SUB 		RCX, 1			;dim--
	ADD 		RAX, 8			;i++
	JMP 			ciclofineaddVettori
	
e3:
	popaq_Add					
	mov			rsp, rbp	
	pop			rbp
	ret
;---------------------------------------------------
global subVettori

subVettori:
	push		rbp				
	mov			rbp, rsp			
	pushaq_sub	
	
	;RDI=v1
	;RSI=v2
	;RDX=ris
	;RCX=dim
	
	XOR 			RAX,RAX			;i=0

ciclosubVettori:
	SUB 			RCX, 16			;dim=dim-4
	CMP 		    RCX, 0			;dim==0?
	JL 			lavorodimezzatosub
	VMOVAPD 	YMM0, [RDI + RAX*8] 	;YMM0=[V1, V1, V1, V1]
	VSUBPD 		YMM0, [RSI + RAX*8]	    ;YMM0=[V1-V2, V1-V2, V1-V2, V1-V2]
	
	VMOVAPD 	YMM1, [RDI + RAX*8+32]
	VSUBPD 		YMM1, [RSI + RAX*8+32]	
	
	VMOVAPD 	YMM2, [RDI + RAX*8+64]
	VSUBPD  	YMM2, [RSI + RAX*8+64]	
	
	VMOVAPD 	YMM3, [RDI + RAX*8+96]
	VSUBPD 		YMM3, [RSI + RAX*8+96]	
	
	VMOVAPD 	[RDX + RAX*8   ], YMM0	
	VMOVAPD 	[RDX + RAX*8+32], YMM1
	VMOVAPD 	[RDX + RAX*8+64], YMM2
	VMOVAPD 	[RDX + RAX*8+96], YMM3
	ADD 		RAX, 16			
	JMP 		ciclosubVettori
	
lavorodimezzatosub:
	ADD RCX, 16 

ciclosubVettorimezzi:
	SUB RCX, 8
	JL finesubVettori
	VMOVAPD 	YMM0, [RDI + RAX*8] 	
	VSUBPD 		YMM0, [RSI + RAX*8]	
		
	VMOVAPD 	YMM1, [RDI + RAX*8+32]
	VSUBPD 		YMM1, [RSI + RAX*8+32]

	VMOVAPD 	[RDX + RAX*8   ], YMM0
	VMOVAPD 	[RDX + RAX*8+32], YMM1
	ADD RAX, 8
	JMP ciclosubVettorimezzi

finesubVettori:
	ADD 		RCX, 7		
		
ciclofinesubVettori:
	CMP 			RCX, 0			;dim==0?
	JL 			e3
	VMOVSD 		XMM1, [RDI + RAX] 	;YMM1=[0, 0, 0, V1]
	VSUBSD 		XMM1, [RSI + RAX]	;YMM1=[0, 0, 0, V1-V2]
	VMOVSD 		[RDX + RAX], XMM1  ;metto il risultato in ris
	SUB 		RCX, 1			;dim--
	ADD 		RAX, 8			;i++
	JMP 		ciclofinesubVettori
	
e4:
	popaq_sub					
	mov			rsp, rbp	
	pop			rbp
	ret

;------------------------------------------------------------------

global prodVet_x_Scalare

prodVet_x_Scalare:
	push		rbp
	mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
	pushaq_prod							; salva i registri generali
	
	; RDI = v1
	; XMM0 = S
	; RSI = ris
	; RDX = dim
	
	XOR 		RAX,RAX
	
    VBROADCASTSD	YMM0, XMM0
		
cicloProdVet_x_Scalare:
	SUB 		RDX,16					; DIM-16
	CMP 		RDX,0					; DIM == 0?
	JL 		prod_x_Scal_lavoromezzi		; se si lavora prendendo 8 double (anzicchè 16)
	VMOVAPD YMM1, [RDI + RAX*8]			; YMM1 = [V1, V1, V1, V1]
	VMULPD	YMM1, YMM0				    ; YMM1 = [V1*S, V1*S, V1*S, V1*S]
				
	VMOVAPD YMM2, [RDI + RAX*8+32]
	VMULPD	YMM2, YMM0
	
	VMOVAPD YMM3, [RDI + RAX*8+64]
	VMULPD	YMM3, YMM0
	
	VMOVAPD YMM4, [RDI + RAX*8+96]
	VMULPD	YMM4, YMM0
	
	VMOVAPD [RSI+RAX*8   ], YMM1        ; Rimetto a posto in memoria  
	VMOVAPD [RSI+RAX*8+32], YMM2
	VMOVAPD [RSI+RAX*8+64], YMM3
	VMOVAPD [RSI+RAX*8+96], YMM4
	ADD 	RAX,16					
	JMP 	cicloProdVet_x_Scalare
prod_x_Scal_lavoromezzi:
	ADD RDX, 16

cicloProdVet_x_Scalare_mezzi:
	SUB RDX, 8
	CMP RDX,0
	JL fineProdVet_x_Scalare
	VMOVAPD YMM1, [RDI + RAX*8]			
	VMULPD	YMM1, YMM0				    
	
	VMOVAPD YMM2, [RDI + RAX*8+32]
	VMULPD	YMM2, YMM0
	
	VMOVAPD [RSI+RAX*8   ], YMM1
	VMOVAPD [RSI+RAX*8+32], YMM2
	ADD RAX,8
	JMP cicloProdVet_x_Scalare_mezzi
	
fineProdVet_x_Scalare:
	ADD 	RDX,7				
	
	
ciclofineProdVet_x_Scalare:
	CMP 	RDX,0					; dim == 0 ?
	JL 		e5						; se si fine
    VMOVSD 	XMM1, [RDI + RAX*8]		; XMM1 = [0, 0, 0, V1]
	VMULSD	XMM1, XMM0				; XMM1 = [0*S, 0*S, 0*S, V1*S]
	VMOVSD 	[RSI+RAX*8], XMM1			; LO RIMETTO A POSTO in memoria
	SUB 	RDX,1					; dim--
	ADD 	RAX,1					; i++
	JMP 		ciclofineProdVet_x_Scalare

e5:
	popaq_prod							; ripristina i registri generali
	mov		rsp, rbp					; ripristina lo Stack Pointer
	pop		rbp						; ripristina il Base Pointer
			
	ret
;--------------------------------------------------

global prodScalare

prodScalare:
	push		rbp
	mov			rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
	pushaq_prodSc					; salva i registri generali
	
	; RDI  = v1
	; RSI  = v2
	; RDX = dim
	
	XOR 	RAX,RAX
	VXORPD 	YMM1, YMM1   ;su YMM1 sommo i prodotti a multipli 
	VXORPD 	YMM2, YMM2	 ;così si azzera anche XMM2 su cui sommeremo singolarmente i restanti float		
	
cicloProdScalare:
	SUB 		RDX,16						; DIM-16
	CMP 		RDX,0						; DIM == 0?
	JL 		prodScal_lavoromezzi				; se si lavoriamo prendendo 8 float
	VMOVAPD YMM0, [RDI + RAX*8]				; YMM0 = [V1, V1, V1, V1]
	VMULPD	YMM0, [RSI + RAX*8]				; YMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	VADDPD	YMM1, YMM0
	
	VMOVAPD YMM2, [RDI + RAX*8+32]
	VMULPD	YMM2, [RSI + RAX*8+32]	
	VADDPD	YMM1, YMM2
	
	VMOVAPD YMM3, [RDI + RAX*8+64]
	VMULPD	YMM3, [RSI + RAX*8+64]	
	VADDPD	YMM1, YMM3
	
	VMOVAPD YMM4, [RDI + RAX*8+96]
	VMULPD	YMM4, [RSI + RAX*8+96]	
	VADDPD	YMM1, YMM4
	
	ADD 	RAX,16					
	JMP 		cicloProdScalare

prodScal_lavoromezzi:
	ADD RDX,16
	
cicloProdScalare_mezzi:	
	SUB RDX,8                               ;dim=dim-8
	CMP RDX,0
	JL fineProdScalare
	VMOVAPD YMM0, [RDI + RAX*8]				; YMM0 = [V1, V1, V1, V1]
	VMULPD	YMM0, [RSI + RAX*8]				; YMM0 = [V1*V2, V1*V2, V1*V2, V1*V2]
	VADDPD	YMM1, YMM0
	
	VMOVAPD YMM2, [RDI + RAX*8+32]
	VMULPD	YMM2, [RSI + RAX*8+32]	
	VADDPD	YMM1, YMM2
	ADD RAX, 8
	JMP cicloProdScalare_mezzi
	
fineProdScalare:
	ADD 	RDX,7			
	VXORPD 	YMM0, YMM0
	
ciclofineProdScalare:
	CMP 	RDX,0						; dim == 0 ?
	JL 		e6							; se si fine
    VMOVSD 	XMM0, [RDI + RAX]			; YMM0 = [0, 0, 0, V1]
	VMULSD	XMM0, [RSI + RAX]			; YMM0 = [0*V2, 0*V2, 0*V2, V1*V2]
	VADDSD	XMM2, XMM0					; XMM2 = [0, 0, 0, +=V1*V2]
	SUB 	RDX,1						; dim--
	ADD 	RAX,8						; i++
	JMP 		ciclofineProdScalare

e6:
	VHADDPD 	    YMM1, YMM1
	VEXTRACTF128 	XMM3, YMM1, 1b 
	
	VADDSD 		XMM1, XMM3
	VADDSD 		XMM2, XMM1
	
	VMOVSD 		[risSommaEu], XMM2

	popaq_prodScal								; ripristina i registri generali
	mov			rsp, 	rbp					; ripristina lo Stack Pointer
	pop			rbp						; ripristina il Base Pointer
	VMOVSD 		XMM0, [risSommaEu] 
	ret									; torna alla funzione C chiamante
;-----------------------------------------------------


global prodVet_x_ScalareUn

prodVet_x_ScalareUn:
	push		rbp
	mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
	pushaq_prodScUn							; salva i registri generali
	
	; RDI = v1
	; XMM0 = S
	; RSI = ris
	; RDX = dim
	
	XOR 		RAX,RAX
	
    VBROADCASTSD	YMM0, XMM0
		
cicloProdVet_x_ScalareUn:
	SUB 		RDX,16					; DIM-16
	CMP 		RDX,0					; DIM == 0?
	JL 		prod_x_Scal_lavoromezziUn		; se si lavora prendendo 8 double (anzicchè 16)
	VMOVUPD YMM1, [RDI + RAX*8]			; YMM1 = [V1, V1, V1, V1]
	VMULPD	YMM1, YMM0				    ; YMM1 = [V1*S, V1*S, V1*S, V1*S]
				
	VMOVUPD YMM2, [RDI + RAX*8+32]
	VMULPD	YMM2, YMM0
	
	VMOVUPD YMM3, [RDI + RAX*8+64]
	VMULPD	YMM3, YMM0
	
	VMOVUPD YMM4, [RDI + RAX*8+96]
	VMULPD	YMM4, YMM0
	
	VMOVAPD [RSI+RAX*8   ], YMM1        ; Rimetto a posto in memoria  
	VMOVAPD [RSI+RAX*8+32], YMM2
	VMOVAPD [RSI+RAX*8+64], YMM3
	VMOVAPD [RSI+RAX*8+96], YMM4
	ADD 	RAX,16					
	JMP 	cicloProdVet_x_ScalareUn
prod_x_Scal_lavoromezziUn:
	ADD RDX, 16

cicloProdVet_x_Scalare_mezziUn:
	SUB RDX, 8
	CMP RDX,0
	JL fineProdVet_x_Scalare
	VMOVUPD YMM1, [RDI + RAX*8]			
	VMULPD	YMM1, YMM0				    
	
	VMOVUPD YMM2, [RDI + RAX*8+32]
	VMULPD	YMM2, YMM0
	
	VMOVAPD [RSI+RAX*8   ], YMM1
	VMOVAPD [RSI+RAX*8+32], YMM2
	ADD RAX,8
	JMP cicloProdVet_x_Scalare_mezziUn
	
fineProdVet_x_ScalareUn:
	ADD 	RDX,7				
	
	
ciclofineProdVet_x_ScalareUn:
	CMP 	RDX,0					; dim == 0 ?
	JL 		e5						; se si fine
    VMOVSD 	XMM1, [RDI + RAX*8]		; XMM1 = [0, 0, 0, V1]
	VMULSD	XMM1, XMM0				; XMM1 = [0*S, 0*S, 0*S, V1*S]
	VMOVSD 	[RSI+RAX*8], XMM1			; LO RIMETTO A POSTO in memoria
	SUB 	RDX,1					; dim--
	ADD 	RAX,1					; i++
	JMP 		ciclofineProdVet_x_ScalareUn


;--------------------------------------------------

