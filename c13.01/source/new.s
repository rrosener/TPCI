	.file	"mole_solve.cpp"
	.globl	__powidf2
	.section	.text._ZSt3powdi,"axG",@progbits,_ZSt3powdi,comdat
	.weak	_ZSt3powdi
	.type	_ZSt3powdi, @function
_ZSt3powdi:
.LFB79:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movsd	%xmm0, -8(%rbp)
	movl	%edi, -12(%rbp)
	movl	-12(%rbp), %edi
	movsd	-8(%rbp), %xmm0
	call	__powidf2
	movsd	%xmm0, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, -24(%rbp)
	movsd	-24(%rbp), %xmm0
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE79:
	.size	_ZSt3powdi, .-_ZSt3powdi
	.section	.text._ZNSt14numeric_limitsIfE3minEv,"axG",@progbits,_ZNSt14numeric_limitsIfE3minEv,comdat
	.weak	_ZNSt14numeric_limitsIfE3minEv
	.type	_ZNSt14numeric_limitsIfE3minEv, @function
_ZNSt14numeric_limitsIfE3minEv:
.LFB238:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movl	.LC0(%rip), %eax
	movl	%eax, -4(%rbp)
	movss	-4(%rbp), %xmm0
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE238:
	.size	_ZNSt14numeric_limitsIfE3minEv, .-_ZNSt14numeric_limitsIfE3minEv
	.section	.text._ZNSt14numeric_limitsIfE3maxEv,"axG",@progbits,_ZNSt14numeric_limitsIfE3maxEv,comdat
	.weak	_ZNSt14numeric_limitsIfE3maxEv
	.type	_ZNSt14numeric_limitsIfE3maxEv, @function
_ZNSt14numeric_limitsIfE3maxEv:
.LFB239:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movl	.LC1(%rip), %eax
	movl	%eax, -4(%rbp)
	movss	-4(%rbp), %xmm0
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE239:
	.size	_ZNSt14numeric_limitsIfE3maxEv, .-_ZNSt14numeric_limitsIfE3maxEv
	.section	.text._ZNSt11char_traitsIcE6lengthEPKc,"axG",@progbits,_ZNSt11char_traitsIcE6lengthEPKc,comdat
	.weak	_ZNSt11char_traitsIcE6lengthEPKc
	.type	_ZNSt11char_traitsIcE6lengthEPKc, @function
_ZNSt11char_traitsIcE6lengthEPKc:
.LFB447:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	$-1, %rcx
	movq	%rax, %rdx
	movl	$0, %eax
	movq	%rdx, %rdi
	repnz scasb
	movq	%rcx, %rax
	notq	%rax
	subq	$1, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE447:
	.size	_ZNSt11char_traitsIcE6lengthEPKc, .-_ZNSt11char_traitsIcE6lengthEPKc
	.section	.text._ZnwmPv,"axG",@progbits,_ZnwmPv,comdat
	.weak	_ZnwmPv
	.type	_ZnwmPv, @function
_ZnwmPv:
.LFB480:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE480:
	.size	_ZnwmPv, .-_ZnwmPv
	.section	.text._ZStorSt13_Ios_OpenmodeS_,"axG",@progbits,_ZStorSt13_Ios_OpenmodeS_,comdat
	.weak	_ZStorSt13_Ios_OpenmodeS_
	.type	_ZStorSt13_Ios_OpenmodeS_, @function
_ZStorSt13_Ios_OpenmodeS_:
.LFB853:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movl	%edi, -4(%rbp)
	movl	%esi, -8(%rbp)
	movl	-8(%rbp), %eax
	movl	-4(%rbp), %edx
	orl	%edx, %eax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE853:
	.size	_ZStorSt13_Ios_OpenmodeS_, .-_ZStorSt13_Ios_OpenmodeS_
	.section	.text._ZSt21__valarray_get_memorym,"axG",@progbits,_ZSt21__valarray_get_memorym,comdat
	.weak	_ZSt21__valarray_get_memorym
	.type	_ZSt21__valarray_get_memorym, @function
_ZSt21__valarray_get_memorym:
.LFB1727:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_Znwm
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE1727:
	.size	_ZSt21__valarray_get_memorym, .-_ZSt21__valarray_get_memorym
	.section	.text._ZSt25__valarray_release_memoryPv,"axG",@progbits,_ZSt25__valarray_release_memoryPv,comdat
	.weak	_ZSt25__valarray_release_memoryPv
	.type	_ZSt25__valarray_release_memoryPv, @function
_ZSt25__valarray_release_memoryPv:
.LFB1729:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZdlPv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE1729:
	.size	_ZSt25__valarray_release_memoryPv, .-_ZSt25__valarray_release_memoryPv
	.local	_ZL8BIGFLOAT
	.comm	_ZL8BIGFLOAT,4,4
	.local	_ZL10SMALLFLOAT
	.comm	_ZL10SMALLFLOAT,4,4
	.local	_ZL6mode_w
	.comm	_ZL6mode_w,4,4
	.local	_ZL6mode_a
	.comm	_ZL6mode_a,4,4
	.local	_ZL7mode_rp
	.comm	_ZL7mode_rp,4,4
	.local	_ZL7mode_wp
	.comm	_ZL7mode_wp,4,4
	.local	_ZL7mode_ap
	.comm	_ZL7mode_ap,4,4
	.local	_ZL7mode_rb
	.comm	_ZL7mode_rb,4,4
	.local	_ZL7mode_wb
	.comm	_ZL7mode_wb,4,4
	.local	_ZL7mode_ab
	.comm	_ZL7mode_ab,4,4
	.local	_ZL8mode_rpb
	.comm	_ZL8mode_rpb,4,4
	.local	_ZL8mode_wpb
	.comm	_ZL8mode_wpb,4,4
	.local	_ZL8mode_apb
	.comm	_ZL8mode_apb,4,4
	.section	.text._ZNK7t_cpu_i13lgAssertAbortEv,"axG",@progbits,_ZNK7t_cpu_i13lgAssertAbortEv,comdat
	.align 2
	.weak	_ZNK7t_cpu_i13lgAssertAbortEv
	.type	_ZNK7t_cpu_i13lgAssertAbortEv, @function
_ZNK7t_cpu_i13lgAssertAbortEv:
.LFB3069:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movzbl	336(%rax), %eax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3069:
	.size	_ZNK7t_cpu_i13lgAssertAbortEv, .-_ZNK7t_cpu_i13lgAssertAbortEv
	.section	.text._ZN5t_cpu1iEv,"axG",@progbits,_ZN5t_cpu1iEv,comdat
	.align 2
	.weak	_ZN5t_cpu1iEv
	.type	_ZN5t_cpu1iEv, @function
_ZN5t_cpu1iEv:
.LFB3083:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	_ZN5t_cpu3m_iE(%rip), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3083:
	.size	_ZN5t_cpu1iEv, .-_ZN5t_cpu1iEv
	.local	_ZL3cpu
	.comm	_ZL3cpu,1,1
	.section	.rodata
	.align 8
.LC2:
	.string	"DISASTER Assertion failure at %s:%ld\n%s\n"
	.section	.text._ZNK10bad_assert5printEv,"axG",@progbits,_ZNK10bad_assert5printEv,comdat
	.align 2
	.weak	_ZNK10bad_assert5printEv
	.type	_ZNK10bad_assert5printEv, @function
_ZNK10bad_assert5printEv:
.LFB3105:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	24(%rax), %rsi
	movq	-8(%rbp), %rax
	movq	16(%rax), %rcx
	movq	-8(%rbp), %rax
	movq	8(%rax), %rdx
	movq	ioQQQ(%rip), %rax
	movq	%rsi, %r8
	movl	$.LC2, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	fprintf
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3105:
	.size	_ZNK10bad_assert5printEv, .-_ZNK10bad_assert5printEv
	.section	.text._ZN10bad_assertD2Ev,"axG",@progbits,_ZN10bad_assertD5Ev,comdat
	.align 2
	.weak	_ZN10bad_assertD2Ev
	.type	_ZN10bad_assertD2Ev, @function
_ZN10bad_assertD2Ev:
.LFB3107:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	$_ZTV10bad_assert+16, (%rax)
	movq	-8(%rbp), %rax
	movq	$0, 8(%rax)
	movl	$0, %eax
	testl	%eax, %eax
	je	.L21
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZdlPv
.L21:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3107:
	.size	_ZN10bad_assertD2Ev, .-_ZN10bad_assertD2Ev
	.weak	_ZN10bad_assertD1Ev
	.set	_ZN10bad_assertD1Ev,_ZN10bad_assertD2Ev
	.section	.text._ZN10bad_assertD0Ev,"axG",@progbits,_ZN10bad_assertD0Ev,comdat
	.align 2
	.weak	_ZN10bad_assertD0Ev
	.type	_ZN10bad_assertD0Ev, @function
_ZN10bad_assertD0Ev:
.LFB3109:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZdlPv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3109:
	.size	_ZN10bad_assertD0Ev, .-_ZN10bad_assertD0Ev
	.section	.text._ZN10bad_assertC2ERKS_,"axG",@progbits,_ZN10bad_assertC5ERKS_,comdat
	.align 2
	.weak	_ZN10bad_assertC2ERKS_
	.type	_ZN10bad_assertC2ERKS_, @function
_ZN10bad_assertC2ERKS_:
.LFB3162:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	$_ZTV10bad_assert+16, (%rax)
	movq	-16(%rbp), %rax
	movq	8(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, 8(%rax)
	movq	-16(%rbp), %rax
	movq	16(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, 16(%rax)
	movq	-16(%rbp), %rax
	movq	24(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, 24(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3162:
	.size	_ZN10bad_assertC2ERKS_, .-_ZN10bad_assertC2ERKS_
	.weak	_ZN10bad_assertC1ERKS_
	.set	_ZN10bad_assertC1ERKS_,_ZN10bad_assertC2ERKS_
	.section	.text._Z4SDIVf,"axG",@progbits,_Z4SDIVf,comdat
	.weak	_Z4SDIVf
	.type	_Z4SDIVf, @function
_Z4SDIVf:
.LFB3174:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movss	%xmm0, -4(%rbp)
	movss	-4(%rbp), %xmm1
	movss	.LC3(%rip), %xmm0
	andps	%xmm1, %xmm0
	movss	_ZL10SMALLFLOAT(%rip), %xmm1
	ucomiss	%xmm0, %xmm1
	jbe	.L33
	movl	_ZL10SMALLFLOAT(%rip), %eax
	jmp	.L30
.L33:
	movl	-4(%rbp), %eax
.L30:
	movl	%eax, -8(%rbp)
	movss	-8(%rbp), %xmm0
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3174:
	.size	_Z4SDIVf, .-_Z4SDIVf
	.section	.text._Z4SDIVd,"axG",@progbits,_Z4SDIVd,comdat
	.weak	_Z4SDIVd
	.type	_Z4SDIVd, @function
_Z4SDIVd:
.LFB3175:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movsd	%xmm0, -8(%rbp)
	movsd	-8(%rbp), %xmm1
	movsd	.LC4(%rip), %xmm0
	andpd	%xmm0, %xmm1
	movss	_ZL10SMALLFLOAT(%rip), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	ucomisd	%xmm1, %xmm0
	jbe	.L40
	movss	_ZL10SMALLFLOAT(%rip), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm2
	movsd	%xmm2, -16(%rbp)
	movq	-16(%rbp), %rax
	jmp	.L37
.L40:
	movq	-8(%rbp), %rax
.L37:
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3175:
	.size	_Z4SDIVd, .-_Z4SDIVd
	.section	.text._ZN8tree_vec8p_clear0Ev,"axG",@progbits,_ZN8tree_vec8p_clear0Ev,comdat
	.align 2
	.weak	_ZN8tree_vec8p_clear0Ev
	.type	_ZN8tree_vec8p_clear0Ev, @function
_ZN8tree_vec8p_clear0Ev:
.LFB3250:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$40, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -40(%rbp)
	movq	-40(%rbp), %rax
	movq	8(%rax), %rax
	testq	%rax, %rax
	je	.L41
	movq	$0, -24(%rbp)
	jmp	.L43
.L44:
	movq	-40(%rbp), %rax
	movq	8(%rax), %rax
	movq	-24(%rbp), %rdx
	salq	$4, %rdx
	addq	%rdx, %rax
	movq	%rax, %rdi
	call	_ZN8tree_vec5clearEv
	addq	$1, -24(%rbp)
.L43:
	movq	-40(%rbp), %rax
	movq	(%rax), %rax
	cmpq	-24(%rbp), %rax
	ja	.L44
	movq	-40(%rbp), %rax
	movq	8(%rax), %rax
	testq	%rax, %rax
	je	.L45
	movq	-40(%rbp), %rax
	movq	8(%rax), %rdx
	movq	-40(%rbp), %rax
	movq	8(%rax), %rax
	subq	$8, %rax
	movq	(%rax), %rax
	salq	$4, %rax
	leaq	(%rdx,%rax), %rbx
.L47:
	movq	-40(%rbp), %rax
	movq	8(%rax), %rax
	cmpq	%rax, %rbx
	je	.L46
	subq	$16, %rbx
	movq	%rbx, %rdi
	call	_ZN8tree_vecD1Ev
	jmp	.L47
.L46:
	movq	-40(%rbp), %rax
	movq	8(%rax), %rax
	subq	$8, %rax
	movq	%rax, %rdi
	call	_ZdaPv
.L45:
.L41:
	addq	$40, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3250:
	.size	_ZN8tree_vec8p_clear0Ev, .-_ZN8tree_vec8p_clear0Ev
	.section	.text._ZN8tree_vec8p_clear1Ev,"axG",@progbits,_ZN8tree_vec8p_clear1Ev,comdat
	.align 2
	.weak	_ZN8tree_vec8p_clear1Ev
	.type	_ZN8tree_vec8p_clear1Ev, @function
_ZN8tree_vec8p_clear1Ev:
.LFB3251:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	$0, (%rax)
	movq	-8(%rbp), %rax
	movq	$0, 8(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3251:
	.size	_ZN8tree_vec8p_clear1Ev, .-_ZN8tree_vec8p_clear1Ev
	.section	.text._ZN8tree_vecC2Ev,"axG",@progbits,_ZN8tree_vecC5Ev,comdat
	.align 2
	.weak	_ZN8tree_vecC2Ev
	.type	_ZN8tree_vecC2Ev, @function
_ZN8tree_vecC2Ev:
.LFB3253:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN8tree_vec8p_clear1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3253:
	.size	_ZN8tree_vecC2Ev, .-_ZN8tree_vecC2Ev
	.weak	_ZN8tree_vecC1Ev
	.set	_ZN8tree_vecC1Ev,_ZN8tree_vecC2Ev
	.section	.text._ZN8tree_vecD2Ev,"axG",@progbits,_ZN8tree_vecD5Ev,comdat
	.align 2
	.weak	_ZN8tree_vecD2Ev
	.type	_ZN8tree_vecD2Ev, @function
_ZN8tree_vecD2Ev:
.LFB3259:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN8tree_vec8p_clear0Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3259:
	.size	_ZN8tree_vecD2Ev, .-_ZN8tree_vecD2Ev
	.weak	_ZN8tree_vecD1Ev
	.set	_ZN8tree_vecD1Ev,_ZN8tree_vecD2Ev
	.section	.text._ZN8tree_vec5clearEv,"axG",@progbits,_ZN8tree_vec5clearEv,comdat
	.align 2
	.weak	_ZN8tree_vec5clearEv
	.type	_ZN8tree_vec5clearEv, @function
_ZN8tree_vec5clearEv:
.LFB3261:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN8tree_vec8p_clear0Ev
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN8tree_vec8p_clear1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3261:
	.size	_ZN8tree_vec5clearEv, .-_ZN8tree_vec5clearEv
	.section	.text._ZN8tree_vec6getvecEmPKm,"axG",@progbits,_ZN8tree_vec6getvecEmPKm,comdat
	.align 2
	.weak	_ZN8tree_vec6getvecEmPKm
	.type	_ZN8tree_vec6getvecEmPKm, @function
_ZN8tree_vec6getvecEmPKm:
.LFB3263:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	cmpq	$0, -16(%rbp)
	jne	.L54
	movq	-8(%rbp), %rax
	jmp	.L55
.L54:
	movq	-16(%rbp), %rax
	leaq	-1(%rax), %rcx
	movq	-24(%rbp), %rdx
	movq	-8(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN8tree_vec6getvecEmPKm
	movq	8(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$3, %rdx
	leaq	-8(%rdx), %rcx
	movq	-24(%rbp), %rdx
	addq	%rcx, %rdx
	movq	(%rdx), %rdx
	salq	$4, %rdx
	addq	%rdx, %rax
.L55:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3263:
	.size	_ZN8tree_vec6getvecEmPKm, .-_ZN8tree_vec6getvecEmPKm
	.section	.text._ZNK8tree_vec6getvecEmPKm,"axG",@progbits,_ZNK8tree_vec6getvecEmPKm,comdat
	.align 2
	.weak	_ZNK8tree_vec6getvecEmPKm
	.type	_ZNK8tree_vec6getvecEmPKm, @function
_ZNK8tree_vec6getvecEmPKm:
.LFB3264:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	cmpq	$0, -16(%rbp)
	jne	.L57
	movq	-8(%rbp), %rax
	jmp	.L58
.L57:
	movq	-16(%rbp), %rax
	leaq	-1(%rax), %rcx
	movq	-24(%rbp), %rdx
	movq	-8(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNK8tree_vec6getvecEmPKm
	movq	8(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$3, %rdx
	leaq	-8(%rdx), %rcx
	movq	-24(%rbp), %rdx
	addq	%rcx, %rdx
	movq	(%rdx), %rdx
	salq	$4, %rdx
	addq	%rdx, %rax
.L58:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3264:
	.size	_ZNK8tree_vec6getvecEmPKm, .-_ZNK8tree_vec6getvecEmPKm
	.local	_ZL7SQAS1SR
	.comm	_ZL7SQAS1SR,8,8
	.local	_ZL8SQAS_SKY
	.comm	_ZL8SQAS_SKY,8,8
	.local	_ZL14ELECTRIC_CONST
	.comm	_ZL14ELECTRIC_CONST,8,8
	.local	_ZL12HION_LTE_POP
	.comm	_ZL12HION_LTE_POP,8,8
	.local	_ZL4SAHA
	.comm	_ZL4SAHA,8,8
	.local	_ZL6HNU3C2
	.comm	_ZL6HNU3C2,8,8
	.local	_ZL12STEFAN_BOLTZ
	.comm	_ZL12STEFAN_BOLTZ,8,8
	.local	_ZL14FINE_STRUCTURE
	.comm	_ZL14FINE_STRUCTURE,8,8
	.local	_ZL15FINE_STRUCTURE2
	.comm	_ZL15FINE_STRUCTURE2,8,8
	.local	_ZL14BOHR_RADIUS_CM
	.comm	_ZL14BOHR_RADIUS_CM,8,8
	.local	_ZL14TWO_PHOT_CONST
	.comm	_ZL14TWO_PHOT_CONST,8,8
	.local	_ZL10COLL_CONST
	.comm	_ZL10COLL_CONST,8,8
	.local	_ZL11MILNE_CONST
	.comm	_ZL11MILNE_CONST,8,8
	.local	_ZL16TRANS_PROB_CONST
	.comm	_ZL16TRANS_PROB_CONST,8,8
	.section	.text._ZN6t_conv16setConvIonizFailEPKcdd,"axG",@progbits,_ZN6t_conv16setConvIonizFailEPKcdd,comdat
	.align 2
	.weak	_ZN6t_conv16setConvIonizFailEPKcdd
	.type	_ZN6t_conv16setConvIonizFailEPKcdd, @function
_ZN6t_conv16setConvIonizFailEPKcdd:
.LFB3872:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movsd	%xmm0, -24(%rbp)
	movsd	%xmm1, -32(%rbp)
	movq	-8(%rbp), %rax
	movb	$0, 6000(%rax)
	movq	-8(%rbp), %rax
	leaq	2000(%rax), %rcx
	movq	-16(%rbp), %rax
	movl	$1999, %edx
	movq	%rax, %rsi
	movq	%rcx, %rdi
	call	strncpy
	movq	-8(%rbp), %rax
	movb	$0, 3999(%rax)
	movq	-8(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rax, 6008(%rdx)
	movq	-8(%rbp), %rdx
	movq	-32(%rbp), %rax
	movq	%rax, 6016(%rdx)
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3872:
	.size	_ZN6t_conv16setConvIonizFailEPKcdd, .-_ZN6t_conv16setConvIonizFailEPKcdd
	.section	.text._ZN6t_conv16incrementCounterE12counter_type,"axG",@progbits,_ZN6t_conv16incrementCounterE12counter_type,comdat
	.align 2
	.weak	_ZN6t_conv16incrementCounterE12counter_type
	.type	_ZN6t_conv16incrementCounterE12counter_type, @function
_ZN6t_conv16incrementCounterE12counter_type:
.LFB3877:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movl	%esi, -12(%rbp)
	movl	-12(%rbp), %eax
	movq	-8(%rbp), %rdx
	movslq	%eax, %rcx
	addq	$818, %rcx
	movq	8(%rdx,%rcx,8), %rdx
	leaq	1(%rdx), %rcx
	movq	-8(%rbp), %rdx
	cltq
	addq	$818, %rax
	movq	%rcx, 8(%rdx,%rax,8)
	movl	-12(%rbp), %eax
	movq	-8(%rbp), %rdx
	movslq	%eax, %rcx
	addq	$830, %rcx
	movq	8(%rdx,%rcx,8), %rdx
	leaq	1(%rdx), %rcx
	movq	-8(%rbp), %rdx
	cltq
	addq	$830, %rax
	movq	%rcx, 8(%rdx,%rax,8)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3877:
	.size	_ZN6t_conv16incrementCounterE12counter_type, .-_ZN6t_conv16incrementCounterE12counter_type
	.section	.text._ZNK9chem_atom15lgMeanAbundanceEv,"axG",@progbits,_ZNK9chem_atom15lgMeanAbundanceEv,comdat
	.align 2
	.weak	_ZNK9chem_atom15lgMeanAbundanceEv
	.type	_ZNK9chem_atom15lgMeanAbundanceEv, @function
_ZNK9chem_atom15lgMeanAbundanceEv:
.LFB3894:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movl	8(%rax), %eax
	shrl	$31, %eax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3894:
	.size	_ZNK9chem_atom15lgMeanAbundanceEv, .-_ZNK9chem_atom15lgMeanAbundanceEv
	.section	.rodata
.LC5:
	.string	"D"
	.string	""
.LC6:
	.string	"^%d"
	.section	.text._ZNK9chem_atom5labelEv,"axG",@progbits,_ZNK9chem_atom5labelEv,comdat
	.align 2
	.weak	_ZNK9chem_atom5labelEv
	.type	_ZNK9chem_atom5labelEv, @function
_ZNK9chem_atom5labelEv:
.LFB3895:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA3895
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$56, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -56(%rbp)
	movq	%rsi, -64(%rbp)
	movq	%fs:40, %rax
	movq	%rax, -24(%rbp)
	xorl	%eax, %eax
	movq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNK9chem_atom15lgMeanAbundanceEv
	testb	%al, %al
	je	.L64
	movq	-64(%rbp), %rax
	movq	(%rax), %rax
	leaq	8(%rax), %rdx
	movq	-56(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB0:
	call	_ZNSsC1ERKSs
.LEHE0:
	jmp	.L63
.L64:
	movq	-64(%rbp), %rax
	movq	(%rax), %rax
	movl	(%rax), %eax
	cmpl	$1, %eax
	jne	.L66
	movq	-64(%rbp), %rax
	movl	8(%rax), %eax
	cmpl	$2, %eax
	jne	.L66
	leaq	-33(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSaIcEC1Ev
	leaq	-33(%rbp), %rdx
	movq	-56(%rbp), %rax
	movl	$.LC5, %esi
	movq	%rax, %rdi
.LEHB1:
	call	_ZNSsC1EPKcRKSaIcE
.LEHE1:
	leaq	-33(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSaIcED1Ev
	jmp	.L63
.L66:
	movq	-64(%rbp), %rax
	movl	8(%rax), %edx
	leaq	-32(%rbp), %rax
	movl	$.LC6, %esi
	movq	%rax, %rdi
	movl	$0, %eax
	call	sprintf
	movq	-64(%rbp), %rax
	movq	(%rax), %rax
	leaq	8(%rax), %rdx
	movq	-56(%rbp), %rax
	leaq	-32(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
.LEHB2:
	call	_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_
	jmp	.L63
.L69:
	movq	%rax, %rbx
	leaq	-33(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSaIcED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE2:
.L63:
	movq	-56(%rbp), %rax
	movq	-24(%rbp), %rcx
	xorq	%fs:40, %rcx
	je	.L68
	call	__stack_chk_fail
.L68:
	addq	$56, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3895:
	.globl	__gxx_personality_v0
	.section	.gcc_except_table._ZNK9chem_atom5labelEv,"aG",@progbits,_ZNK9chem_atom5labelEv,comdat
.LLSDA3895:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE3895-.LLSDACSB3895
.LLSDACSB3895:
	.uleb128 .LEHB0-.LFB3895
	.uleb128 .LEHE0-.LEHB0
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB1-.LFB3895
	.uleb128 .LEHE1-.LEHB1
	.uleb128 .L69-.LFB3895
	.uleb128 0
	.uleb128 .LEHB2-.LFB3895
	.uleb128 .LEHE2-.LEHB2
	.uleb128 0
	.uleb128 0
.LLSDACSE3895:
	.section	.text._ZNK9chem_atom5labelEv,"axG",@progbits,_ZNK9chem_atom5labelEv,comdat
	.size	_ZNK9chem_atom5labelEv, .-_ZNK9chem_atom5labelEv
	.section	.text._ZNK8molecule11isMonatomicEv,"axG",@progbits,_ZNK8molecule11isMonatomicEv,comdat
	.align 2
	.weak	_ZNK8molecule11isMonatomicEv
	.type	_ZNK8molecule11isMonatomicEv, @function
_ZNK8molecule11isMonatomicEv:
.LFB3905:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA3905
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$40, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -40(%rbp)
	movl	$0, %ebx
	movq	-40(%rbp), %rax
	addq	$24, %rax
	movq	%rax, %rdi
.LEHB3:
	call	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv
	cmpq	$1, %rax
	jne	.L71
	movq	-40(%rbp), %rax
	addq	$24, %rax
	movq	%rax, %rdi
	call	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	movq	%rax, -32(%rbp)
	movl	$1, %ebx
	leaq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
.LEHE3:
	movl	16(%rax), %eax
	cmpl	$1, %eax
	jne	.L71
	movl	$1, %eax
	jmp	.L72
.L71:
	movl	$0, %eax
.L72:
	testb	%bl, %bl
	testb	%al, %al
	je	.L74
	movl	$1, %eax
	jmp	.L79
.L74:
	movl	$0, %eax
	jmp	.L79
.L78:
	testb	%bl, %bl
	nop
	movq	%rax, %rdi
.LEHB4:
	call	_Unwind_Resume
.LEHE4:
.L79:
	addq	$40, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3905:
	.section	.gcc_except_table._ZNK8molecule11isMonatomicEv,"aG",@progbits,_ZNK8molecule11isMonatomicEv,comdat
.LLSDA3905:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE3905-.LLSDACSB3905
.LLSDACSB3905:
	.uleb128 .LEHB3-.LFB3905
	.uleb128 .LEHE3-.LEHB3
	.uleb128 .L78-.LFB3905
	.uleb128 0
	.uleb128 .LEHB4-.LFB3905
	.uleb128 .LEHE4-.LEHB4
	.uleb128 0
	.uleb128 0
.LLSDACSE3905:
	.section	.text._ZNK8molecule11isMonatomicEv,"axG",@progbits,_ZNK8molecule11isMonatomicEv,comdat
	.size	_ZNK8molecule11isMonatomicEv, .-_ZNK8molecule11isMonatomicEv
	.section	.text._ZN8GroupMapC2Em,"axG",@progbits,_ZN8GroupMapC5Em,comdat
	.align 2
	.weak	_ZN8GroupMapC2Em
	.type	_ZN8GroupMapC2Em, @function
_ZN8GroupMapC2Em:
.LFB3914:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA3914
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$40, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
.LEHB5:
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC1Ev
.LEHE5:
	movq	-40(%rbp), %rax
	addq	$144, %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEC1Ev
	movq	-40(%rbp), %rax
	leaq	144(%rax), %rdx
	movq	-48(%rbp), %rax
	xorpd	%xmm0, %xmm0
	movq	%rax, %rsi
	movq	%rdx, %rdi
.LEHB6:
	call	_ZNSt8valarrayIdE6resizeEmd
	movq	-40(%rbp), %rax
	movq	-48(%rbp), %rdx
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm
	movq	$0, -24(%rbp)
	jmp	.L81
.L82:
	movq	-48(%rbp), %rax
	leaq	1(%rax), %rdx
	movq	-40(%rbp), %rax
	movq	-24(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm
	addq	$1, -24(%rbp)
.L81:
	movq	-24(%rbp), %rax
	cmpq	-48(%rbp), %rax
	jb	.L82
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv
.LEHE6:
	jmp	.L85
.L84:
	movq	%rax, %rbx
	movq	-40(%rbp), %rax
	addq	$144, %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdED1Ev
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB7:
	call	_Unwind_Resume
.LEHE7:
.L85:
	addq	$40, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE3914:
	.section	.gcc_except_table._ZN8GroupMapC2Em,"aG",@progbits,_ZN8GroupMapC5Em,comdat
.LLSDA3914:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE3914-.LLSDACSB3914
.LLSDACSB3914:
	.uleb128 .LEHB5-.LFB3914
	.uleb128 .LEHE5-.LEHB5
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB6-.LFB3914
	.uleb128 .LEHE6-.LEHB6
	.uleb128 .L84-.LFB3914
	.uleb128 0
	.uleb128 .LEHB7-.LFB3914
	.uleb128 .LEHE7-.LEHB7
	.uleb128 0
	.uleb128 0
.LLSDACSE3914:
	.section	.text._ZN8GroupMapC2Em,"axG",@progbits,_ZN8GroupMapC5Em,comdat
	.size	_ZN8GroupMapC2Em, .-_ZN8GroupMapC2Em
	.weak	_ZN8GroupMapC1Em
	.set	_ZN8GroupMapC1Em,_ZN8GroupMapC2Em
	.section	.text._ZN8GroupMapD2Ev,"axG",@progbits,_ZN8GroupMapD5Ev,comdat
	.align 2
	.weak	_ZN8GroupMapD2Ev
	.type	_ZN8GroupMapD2Ev, @function
_ZN8GroupMapD2Ev:
.LFB4014:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4014
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	addq	$144, %rax
	movq	%rax, %rdi
.LEHB8:
	call	_ZNSt8valarrayIdED1Ev
.LEHE8:
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
.LEHB9:
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED1Ev
.LEHE9:
	jmp	.L90
.L89:
	movq	%rax, %rbx
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB10:
	call	_Unwind_Resume
.LEHE10:
.L90:
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4014:
	.section	.gcc_except_table._ZN8GroupMapD2Ev,"aG",@progbits,_ZN8GroupMapD5Ev,comdat
.LLSDA4014:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4014-.LLSDACSB4014
.LLSDACSB4014:
	.uleb128 .LEHB8-.LFB4014
	.uleb128 .LEHE8-.LEHB8
	.uleb128 .L89-.LFB4014
	.uleb128 0
	.uleb128 .LEHB9-.LFB4014
	.uleb128 .LEHE9-.LEHB9
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB10-.LFB4014
	.uleb128 .LEHE10-.LEHB10
	.uleb128 0
	.uleb128 0
.LLSDACSE4014:
	.section	.text._ZN8GroupMapD2Ev,"axG",@progbits,_ZN8GroupMapD5Ev,comdat
	.size	_ZN8GroupMapD2Ev, .-_ZN8GroupMapD2Ev
	.weak	_ZN8GroupMapD1Ev
	.set	_ZN8GroupMapD1Ev,_ZN8GroupMapD2Ev
	.section	.rodata
	.align 8
.LC8:
	.string	"Need to treat hmi.H2_frac_abund_set in revised network\n"
.LC9:
	.string	"%g\n"
.LC10:
	.string	"Failed: lgElemsConserved()"
.LC11:
	.string	"mole_solve.cpp"
.LC17:
	.string	"Failed: rlimit < BIGFLOAT"
.LC18:
	.string	"Failed: rlimit > 0.0"
.LC22:
	.string	"Chem conv"
.LC24:
	.string	"chem eft chg"
.LC25:
	.string	"Mole zn %3ld -- %7.2f\n"
	.align 8
.LC26:
	.string	"Internal error %15.8g nBad %4d loops %3d\n"
	.align 8
.LC27:
	.string	"Scaling delta on ions %15.8g; threshold %15.8g\n"
	.text
	.globl	_Z10mole_solvev
	.type	_Z10mole_solvev, @function
_Z10mole_solvev:
.LFB4012:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4012
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$392, %rsp
	.cfi_offset 3, -24
	movl	mole_global+64(%rip), %eax
	movslq	%eax, %rdx
	leaq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB11:
	call	_ZNSt8valarrayIdEC1Em
.LEHE11:
	movl	mole_global+64(%rip), %eax
	movslq	%eax, %rdx
	leaq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB12:
	call	_ZNSt8valarrayIdEC1Em
.LEHE12:
	movl	$atom_list, %edi
	call	_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv
	movq	%rax, %rdx
	leaq	-176(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB13:
	call	_ZN8GroupMapC1Em
.LEHE13:
	movsd	hmi+544(%rip), %xmm0
	xorpd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	jbe	.L143
.LEHB14:
	call	_ZL12mole_h_fixupv
	call	_Z5fixitv
	movq	stderr(%rip), %rax
	movq	%rax, %rcx
	movl	$55, %edx
	movl	$1, %esi
	movl	$.LC8, %edi
	call	fwrite
	movq	hmi+544(%rip), %rax
	movq	stderr(%rip), %rdx
	movq	%rax, -360(%rbp)
	movsd	-360(%rbp), %xmm0
	movl	$.LC9, %esi
	movq	%rdx, %rdi
	movl	$1, %eax
	call	fprintf
	movl	$-1, %edi
	call	exit
.L143:
	call	_Z16lgElemsConservedv
	xorl	$1, %eax
	movzbl	%al, %eax
	testq	%rax, %rax
	setne	%al
	testb	%al, %al
	je	.L94
	leaq	-208(%rbp), %rax
	movl	$.LC10, %ecx
	movl	$69, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE14:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L95
	leaq	-208(%rbp), %rax
	movq	%rax, %rdi
.LEHB15:
	call	_ZNK10bad_assert5printEv
	call	abort
.L95:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-208(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE15:
.L94:
	movabsq	$4487126258331716666, %rax
	movq	%rax, -272(%rbp)
	movl	$100, -328(%rbp)
	movl	$0, -324(%rbp)
	movl	-324(%rbp), %eax
	movl	%eax, -336(%rbp)
	movl	.LC13(%rip), %eax
	movl	%eax, -344(%rbp)
	movl	mole_global+64(%rip), %eax
	movslq	%eax, %rdx
	leaq	-224(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB16:
	call	_ZNSt8valarrayIdEC1Em
.LEHE16:
	movabsq	$-4616189618054758400, %rax
	movq	%rax, -312(%rbp)
	movl	$0, %eax
	movq	%rax, -304(%rbp)
	movl	$0, -332(%rbp)
	jmp	.L96
.L117:
	movl	-336(%rbp), %eax
	movl	%eax, -324(%rbp)
	leaq	-256(%rbp), %rax
	movq	%rax, %rdi
.LEHB17:
	call	_Z7get_ptrIdEPT_RSt8valarrayIS0_E
	movq	%rax, %rdx
	leaq	-176(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZN8GroupMap5setupEPd
	movl	$0, %esi
	movl	$conv, %edi
	call	_ZN6t_conv16incrementCounterE12counter_type
	movq	$0, -296(%rbp)
	jmp	.L97
.L111:
	movl	$1, %esi
	movl	$conv, %edi
	call	_ZN6t_conv16incrementCounterE12counter_type
	movl	mole_global+64(%rip), %eax
	movslq	%eax, %r9
	leaq	-340(%rbp), %r8
	leaq	-344(%rbp), %rcx
	leaq	-240(%rbp), %rdx
	leaq	-256(%rbp), %rsi
	leaq	-176(%rbp), %rax
	movq	$_ZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPb, 24(%rsp)
	leaq	-224(%rbp), %rdi
	movq	%rdi, 16(%rsp)
	leaq	-304(%rbp), %rdi
	movq	%rdi, 8(%rsp)
	leaq	-312(%rbp), %rdi
	movq	%rdi, (%rsp)
	movq	%rax, %rdi
	call	_Z11newton_stepR8GroupMapRKSt8valarrayIdERS2_PfS6_lPdS7_S5_PFvS0_S4_S7_S7_bPbE
	movb	%al, -345(%rbp)
	cmpb	$0, -345(%rbp)
	je	.L98
	movl	$0, -336(%rbp)
	movl	$0, %eax
	movq	%rax, -288(%rbp)
	movq	$-1, -264(%rbp)
	movq	$0, -280(%rbp)
	jmp	.L99
.L103:
	movq	-280(%rbp), %rdx
	leaq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movsd	(%rax), %xmm1
	xorpd	%xmm0, %xmm0
	ucomisd	%xmm1, %xmm0
	seta	%al
	testb	%al, %al
	je	.L100
	movq	-280(%rbp), %rdx
	leaq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movsd	(%rax), %xmm1
	movsd	.LC15(%rip), %xmm0
	xorpd	%xmm0, %xmm1
	movsd	%xmm1, -360(%rbp)
	movq	-280(%rbp), %rdx
	leaq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	(%rax), %rax
	movq	%rax, -368(%rbp)
	movsd	-368(%rbp), %xmm0
	call	_Z4SDIVd
	mulsd	-288(%rbp), %xmm0
	movsd	-360(%rbp), %xmm1
	ucomisd	%xmm0, %xmm1
	seta	%al
	testb	%al, %al
	je	.L101
	movq	-280(%rbp), %rdx
	leaq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movsd	(%rax), %xmm1
	movsd	.LC15(%rip), %xmm0
	xorpd	%xmm0, %xmm1
	movsd	%xmm1, -360(%rbp)
	movq	-280(%rbp), %rdx
	leaq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	(%rax), %rax
	movq	%rax, -368(%rbp)
	movsd	-368(%rbp), %xmm0
	call	_Z4SDIVd
	movsd	-360(%rbp), %xmm1
	divsd	%xmm0, %xmm1
	movapd	%xmm1, %xmm0
	movsd	%xmm0, -288(%rbp)
	movq	-280(%rbp), %rax
	movq	%rax, -264(%rbp)
.L101:
	movq	-280(%rbp), %rdx
	leaq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movsd	(%rax), %xmm1
	movsd	.LC16(%rip), %xmm0
	ucomisd	%xmm1, %xmm0
	seta	%al
	testb	%al, %al
	je	.L102
	addl	$1, -336(%rbp)
.L102:
	movq	-280(%rbp), %rdx
	leaq	-240(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	%rax, %rdx
	movl	$0, %eax
	movq	%rax, (%rdx)
.L100:
	addq	$1, -280(%rbp)
.L99:
	movl	mole_global+64(%rip), %eax
	cltq
	cmpq	-280(%rbp), %rax
	jg	.L103
	cmpl	$0, -336(%rbp)
	je	.L104
	movb	$0, -345(%rbp)
.L104:
.L98:
	cmpb	$0, -345(%rbp)
	je	.L105
	jmp	.L106
.L105:
	movss	_ZL8BIGFLOAT(%rip), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	-312(%rbp), %xmm1
	ucomisd	%xmm1, %xmm0
	seta	%al
	xorl	$1, %eax
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L107
	leaq	-208(%rbp), %rax
	movl	$.LC17, %ecx
	movl	$158, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE17:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L108
	leaq	-208(%rbp), %rax
	movq	%rax, %rdi
.LEHB18:
	call	_ZNK10bad_assert5printEv
	call	abort
.L108:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-208(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE18:
.L107:
	movsd	-312(%rbp), %xmm0
	xorpd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	seta	%al
	xorl	$1, %eax
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L109
	leaq	-208(%rbp), %rax
	movl	$.LC18, %ecx
	movl	$159, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
.LEHB19:
	call	_ZN10bad_assertC1EPKclS1_
.LEHE19:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L110
	leaq	-208(%rbp), %rax
	movq	%rax, %rdi
.LEHB20:
	call	_ZNK10bad_assert5printEv
	call	abort
.L110:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-208(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE20:
.L109:
	movsd	-312(%rbp), %xmm1
	movsd	.LC19(%rip), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, -312(%rbp)
	addq	$1, -296(%rbp)
.L97:
	cmpq	$29, -296(%rbp)
	jle	.L111
.L106:
	movss	-340(%rbp), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm1
	movss	-344(%rbp), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	.LC20(%rip), %xmm2
	mulsd	%xmm2, %xmm0
	ucomisd	%xmm1, %xmm0
	jbe	.L112
	movsd	-312(%rbp), %xmm1
	movsd	.LC21(%rip), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, -312(%rbp)
.L112:
	leaq	-240(%rbp), %rdx
	leaq	-176(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB21:
	call	_ZN8GroupMap15updateMoleculesERKSt8valarrayIdE
	movss	-344(%rbp), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	.LC12(%rip), %xmm1
	ucomisd	%xmm0, %xmm1
	jbe	.L114
	cmpl	$0, -336(%rbp)
	jne	.L114
	cmpl	$0, -324(%rbp)
	jne	.L114
	jmp	.L116
.L114:
	addl	$1, -332(%rbp)
.L96:
	cmpl	$99, -332(%rbp)
	jle	.L117
.L116:
	cmpl	$100, -332(%rbp)
	jne	.L118
	movss	-344(%rbp), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	ucomisd	.LC12(%rip), %xmm0
	ja	.L119
.L118:
	cmpl	$0, -336(%rbp)
	je	.L120
.L119:
	cvtsi2sd	-336(%rbp), %xmm1
	movss	-344(%rbp), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movl	$.LC22, %esi
	movl	$conv, %edi
	call	_ZN6t_conv16setConvIonizFailEPKcdd
.L120:
	leaq	-176(%rbp), %rax
	movq	%rax, %rdi
	call	_Z26mole_return_cached_speciesRK8GroupMap
	movss	%xmm0, -360(%rbp)
	movl	-360(%rbp), %eax
	movl	%eax, -320(%rbp)
	movss	conv+6364(%rip), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	.LC23(%rip), %xmm1
	mulsd	%xmm1, %xmm0
	unpcklpd	%xmm0, %xmm0
	cvtpd2ps	%xmm0, %xmm3
	movss	%xmm3, -316(%rbp)
	movss	-320(%rbp), %xmm0
	ucomiss	-316(%rbp), %xmm0
	jbe	.L121
	movss	-320(%rbp), %xmm0
	cvtps2pd	%xmm0, %xmm0
	xorpd	%xmm1, %xmm1
	movl	$.LC24, %esi
	movl	$conv, %edi
	call	_ZN6t_conv16setConvIonizFailEPKcdd
.L121:
	movzbl	trace(%rip), %eax
	testb	%al, %al
	je	.L123
	movq	fnzone(%rip), %rax
	movq	nzone(%rip), %rdx
	movq	ioQQQ(%rip), %rcx
	movq	%rax, -360(%rbp)
	movsd	-360(%rbp), %xmm0
	movl	$.LC25, %esi
	movq	%rcx, %rdi
	movl	$1, %eax
	call	fprintf
	movss	-344(%rbp), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movq	ioQQQ(%rip), %rax
	movl	-332(%rbp), %ecx
	movl	-336(%rbp), %edx
	movl	$.LC26, %esi
	movq	%rax, %rdi
	movl	$1, %eax
	call	fprintf
	movss	-316(%rbp), %xmm1
	cvtps2pd	%xmm1, %xmm1
	movss	-320(%rbp), %xmm0
	cvtps2pd	%xmm0, %xmm0
	movq	ioQQQ(%rip), %rax
	movl	$.LC27, %esi
	movq	%rax, %rdi
	movl	$2, %eax
	call	fprintf
.L123:
	call	_Z21check_co_ion_convergev
.LEHE21:
	movss	-344(%rbp), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm4
	movsd	%xmm4, -360(%rbp)
	leaq	-224(%rbp), %rax
	movq	%rax, %rdi
.LEHB22:
	call	_ZNSt8valarrayIdED1Ev
.LEHE22:
	leaq	-176(%rbp), %rax
	movq	%rax, %rdi
.LEHB23:
	call	_ZN8GroupMapD1Ev
.LEHE23:
	leaq	-240(%rbp), %rax
	movq	%rax, %rdi
.LEHB24:
	call	_ZNSt8valarrayIdED1Ev
.LEHE24:
	leaq	-256(%rbp), %rax
	movq	%rax, %rdi
.LEHB25:
	call	_ZNSt8valarrayIdED1Ev
.LEHE25:
	movq	-360(%rbp), %rax
	jmp	.L144
.L135:
	movq	%rax, %rbx
	leaq	-208(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	jmp	.L126
.L137:
	movq	%rax, %rbx
	leaq	-208(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	jmp	.L128
.L138:
	movq	%rax, %rbx
	leaq	-208(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	jmp	.L128
.L136:
	movq	%rax, %rbx
.L128:
	leaq	-224(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdED1Ev
	jmp	.L126
.L134:
	movq	%rax, %rbx
.L126:
	leaq	-176(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN8GroupMapD1Ev
	jmp	.L130
.L133:
	movq	%rax, %rbx
.L130:
	leaq	-240(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdED1Ev
	jmp	.L131
.L132:
	movq	%rax, %rbx
.L131:
	leaq	-256(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB26:
	call	_Unwind_Resume
.LEHE26:
.L144:
	movq	%rax, -360(%rbp)
	movsd	-360(%rbp), %xmm0
	addq	$392, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4012:
	.section	.gcc_except_table,"a",@progbits
.LLSDA4012:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4012-.LLSDACSB4012
.LLSDACSB4012:
	.uleb128 .LEHB11-.LFB4012
	.uleb128 .LEHE11-.LEHB11
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB12-.LFB4012
	.uleb128 .LEHE12-.LEHB12
	.uleb128 .L132-.LFB4012
	.uleb128 0
	.uleb128 .LEHB13-.LFB4012
	.uleb128 .LEHE13-.LEHB13
	.uleb128 .L133-.LFB4012
	.uleb128 0
	.uleb128 .LEHB14-.LFB4012
	.uleb128 .LEHE14-.LEHB14
	.uleb128 .L134-.LFB4012
	.uleb128 0
	.uleb128 .LEHB15-.LFB4012
	.uleb128 .LEHE15-.LEHB15
	.uleb128 .L135-.LFB4012
	.uleb128 0
	.uleb128 .LEHB16-.LFB4012
	.uleb128 .LEHE16-.LEHB16
	.uleb128 .L134-.LFB4012
	.uleb128 0
	.uleb128 .LEHB17-.LFB4012
	.uleb128 .LEHE17-.LEHB17
	.uleb128 .L136-.LFB4012
	.uleb128 0
	.uleb128 .LEHB18-.LFB4012
	.uleb128 .LEHE18-.LEHB18
	.uleb128 .L137-.LFB4012
	.uleb128 0
	.uleb128 .LEHB19-.LFB4012
	.uleb128 .LEHE19-.LEHB19
	.uleb128 .L136-.LFB4012
	.uleb128 0
	.uleb128 .LEHB20-.LFB4012
	.uleb128 .LEHE20-.LEHB20
	.uleb128 .L138-.LFB4012
	.uleb128 0
	.uleb128 .LEHB21-.LFB4012
	.uleb128 .LEHE21-.LEHB21
	.uleb128 .L136-.LFB4012
	.uleb128 0
	.uleb128 .LEHB22-.LFB4012
	.uleb128 .LEHE22-.LEHB22
	.uleb128 .L134-.LFB4012
	.uleb128 0
	.uleb128 .LEHB23-.LFB4012
	.uleb128 .LEHE23-.LEHB23
	.uleb128 .L133-.LFB4012
	.uleb128 0
	.uleb128 .LEHB24-.LFB4012
	.uleb128 .LEHE24-.LEHB24
	.uleb128 .L132-.LFB4012
	.uleb128 0
	.uleb128 .LEHB25-.LFB4012
	.uleb128 .LEHE25-.LEHB25
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB26-.LFB4012
	.uleb128 .LEHE26-.LEHB26
	.uleb128 0
	.uleb128 0
.LLSDACSE4012:
	.text
	.size	_Z10mole_solvev, .-_Z10mole_solvev
	.section	.rodata
.LC28:
	.string	"C"
.LC30:
	.string	"CO C0 con"
.LC31:
	.string	"C+"
.LC32:
	.string	"CO C1 con"
.LC33:
	.string	"O"
.LC34:
	.string	"CO O0 con"
.LC35:
	.string	"O+"
.LC36:
	.string	"CO O1 con"
	.text
	.globl	_Z21check_co_ion_convergev
	.type	_Z21check_co_ion_convergev, @function
_Z21check_co_ion_convergev:
.LFB4016:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movzbl	dense+17641(%rip), %eax
	testb	%al, %al
	je	.L146
	movsd	dense+2128(%rip), %xmm6
	movsd	%xmm6, -8(%rbp)
	movl	$.LC28, %edi
	call	_Z16findspecieslocalPKc
	movsd	40(%rax), %xmm0
	movsd	-8(%rbp), %xmm6
	subsd	%xmm0, %xmm6
	movapd	%xmm6, %xmm0
	movsd	.LC4(%rip), %xmm1
	andpd	%xmm1, %xmm0
	movsd	%xmm0, -8(%rbp)
	movl	dense+20(%rip), %eax
	movl	%eax, -12(%rbp)
	movss	-12(%rbp), %xmm0
	call	_Z4SDIVf
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	-8(%rbp), %xmm2
	divsd	%xmm0, %xmm2
	movapd	%xmm2, %xmm0
	ucomisd	.LC29(%rip), %xmm0
	jbe	.L146
	movl	$1, %eax
	jmp	.L148
.L146:
	movl	$0, %eax
.L148:
	testb	%al, %al
	je	.L149
	movl	$.LC28, %edi
	call	_Z16findspecieslocalPKc
	movq	40(%rax), %rdx
	movq	dense+2128(%rip), %rax
	movq	%rdx, -8(%rbp)
	movsd	-8(%rbp), %xmm1
	movq	%rax, -8(%rbp)
	movsd	-8(%rbp), %xmm0
	movl	$.LC30, %esi
	movl	$conv, %edi
	call	_ZN6t_conv16setConvIonizFailEPKcdd
	jmp	.L145
.L149:
	movzbl	dense+17641(%rip), %eax
	testb	%al, %al
	je	.L151
	movsd	dense+2136(%rip), %xmm7
	movsd	%xmm7, -8(%rbp)
	movl	$.LC31, %edi
	call	_Z16findspecieslocalPKc
	movsd	40(%rax), %xmm0
	movsd	-8(%rbp), %xmm7
	subsd	%xmm0, %xmm7
	movapd	%xmm7, %xmm0
	movsd	.LC4(%rip), %xmm1
	andpd	%xmm1, %xmm0
	movsd	%xmm0, -8(%rbp)
	movl	dense+20(%rip), %eax
	movl	%eax, -12(%rbp)
	movss	-12(%rbp), %xmm0
	call	_Z4SDIVf
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	-8(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movapd	%xmm3, %xmm0
	ucomisd	.LC29(%rip), %xmm0
	jbe	.L151
	movl	$1, %eax
	jmp	.L153
.L151:
	movl	$0, %eax
.L153:
	testb	%al, %al
	je	.L154
	movl	$.LC31, %edi
	call	_Z16findspecieslocalPKc
	movq	40(%rax), %rdx
	movq	dense+2136(%rip), %rax
	movq	%rdx, -8(%rbp)
	movsd	-8(%rbp), %xmm1
	movq	%rax, -8(%rbp)
	movsd	-8(%rbp), %xmm0
	movl	$.LC32, %esi
	movl	$conv, %edi
	call	_ZN6t_conv16setConvIonizFailEPKcdd
	jmp	.L145
.L154:
	movzbl	dense+17643(%rip), %eax
	testb	%al, %al
	je	.L155
	movsd	dense+2624(%rip), %xmm2
	movsd	%xmm2, -8(%rbp)
	movl	$.LC33, %edi
	call	_Z16findspecieslocalPKc
	movsd	40(%rax), %xmm0
	movsd	-8(%rbp), %xmm2
	subsd	%xmm0, %xmm2
	movapd	%xmm2, %xmm0
	movsd	.LC4(%rip), %xmm1
	andpd	%xmm1, %xmm0
	movsd	%xmm0, -8(%rbp)
	movl	dense+28(%rip), %eax
	movl	%eax, -12(%rbp)
	movss	-12(%rbp), %xmm0
	call	_Z4SDIVf
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	-8(%rbp), %xmm4
	divsd	%xmm0, %xmm4
	movapd	%xmm4, %xmm0
	ucomisd	.LC29(%rip), %xmm0
	jbe	.L155
	movl	$1, %eax
	jmp	.L157
.L155:
	movl	$0, %eax
.L157:
	testb	%al, %al
	je	.L158
	movl	$.LC33, %edi
	call	_Z16findspecieslocalPKc
	movq	40(%rax), %rdx
	movq	dense+2624(%rip), %rax
	movq	%rdx, -8(%rbp)
	movsd	-8(%rbp), %xmm1
	movq	%rax, -8(%rbp)
	movsd	-8(%rbp), %xmm0
	movl	$.LC34, %esi
	movl	$conv, %edi
	call	_ZN6t_conv16setConvIonizFailEPKcdd
	jmp	.L145
.L158:
	movzbl	dense+17643(%rip), %eax
	testb	%al, %al
	je	.L159
	movsd	dense+2632(%rip), %xmm3
	movsd	%xmm3, -8(%rbp)
	movl	$.LC35, %edi
	call	_Z16findspecieslocalPKc
	movsd	40(%rax), %xmm0
	movsd	-8(%rbp), %xmm3
	subsd	%xmm0, %xmm3
	movapd	%xmm3, %xmm0
	movsd	.LC4(%rip), %xmm1
	andpd	%xmm1, %xmm0
	movsd	%xmm0, -8(%rbp)
	movl	dense+28(%rip), %eax
	movl	%eax, -12(%rbp)
	movss	-12(%rbp), %xmm0
	call	_Z4SDIVf
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	-8(%rbp), %xmm5
	divsd	%xmm0, %xmm5
	movapd	%xmm5, %xmm0
	ucomisd	.LC29(%rip), %xmm0
	jbe	.L159
	movl	$1, %eax
	jmp	.L161
.L159:
	movl	$0, %eax
.L161:
	testb	%al, %al
	je	.L145
	movl	$.LC35, %edi
	call	_Z16findspecieslocalPKc
	movq	40(%rax), %rdx
	movq	dense+2632(%rip), %rax
	movq	%rdx, -8(%rbp)
	movsd	-8(%rbp), %xmm1
	movq	%rax, -8(%rbp)
	movsd	-8(%rbp), %xmm0
	movl	$.LC36, %esi
	movl	$conv, %edi
	call	_ZN6t_conv16setConvIonizFailEPKcdd
.L145:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4016:
	.size	_Z21check_co_ion_convergev, .-_Z21check_co_ion_convergev
	.section	.rodata
.LC37:
	.string	"H2"
.LC38:
	.string	"H2*"
	.text
	.type	_ZL12mole_h_fixupv, @function
_ZL12mole_h_fixupv:
.LFB4017:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movsd	hmi+544(%rip), %xmm0
	xorpd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	jbe	.L166
	movq	$0, -16(%rbp)
	jmp	.L169
.L170:
	movq	-16(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movq	%rax, %rdx
	movl	$0, %eax
	movq	%rax, 40(%rdx)
	addq	$1, -16(%rbp)
.L169:
	movl	mole_global+60(%rip), %eax
	cltq
	cmpq	-16(%rbp), %rax
	jg	.L170
	movss	_ZL10SMALLFLOAT(%rip), %xmm0
	addss	%xmm0, %xmm0
	movss	dense(%rip), %xmm1
	mulss	%xmm1, %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	%xmm0, dense+896(%rip)
	movq	dense+896(%rip), %rax
	movq	%rax, dense+888(%rip)
	movl	$.LC37, %edi
	call	_Z16findspecieslocalPKc
	movss	dense(%rip), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	hmi+544(%rip), %xmm1
	mulsd	%xmm1, %xmm0
	unpcklpd	%xmm0, %xmm0
	cvtpd2ps	%xmm0, %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	%xmm0, 40(%rax)
	movl	$.LC38, %edi
	call	_Z16findspecieslocalPKc
	movq	%rax, %rdx
	movl	$0, %eax
	movq	%rax, 40(%rdx)
	movl	$.LC37, %edi
	call	_Z16findspecieslocalPKc
	movsd	40(%rax), %xmm2
	movsd	%xmm2, -24(%rbp)
	movl	$.LC38, %edi
	call	_Z16findspecieslocalPKc
	movsd	40(%rax), %xmm0
	addsd	-24(%rbp), %xmm0
	movsd	%xmm0, hmi(%rip)
	movsd	hmi(%rip), %xmm1
	movsd	.LC39(%rip), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, h2+560(%rip)
	movsd	hmi(%rip), %xmm1
	movsd	.LC40(%rip), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, h2+568(%rip)
	movsd	hmi(%rip), %xmm0
	unpcklpd	%xmm0, %xmm0
	cvtpd2ps	%xmm0, %xmm0
	movss	%xmm0, hmi+8(%rip)
	movsd	h2+560(%rip), %xmm0
	unpcklpd	%xmm0, %xmm0
	cvtpd2ps	%xmm0, %xmm0
	movss	%xmm0, h2+576(%rip)
	movsd	h2+568(%rip), %xmm0
	unpcklpd	%xmm0, %xmm0
	cvtpd2ps	%xmm0, %xmm0
	movss	%xmm0, h2+580(%rip)
	movl	$0, %eax
	movq	%rax, hmi+32(%rip)
	movl	$0, %eax
	movq	%rax, hmi+96(%rip)
	movl	.LC41(%rip), %eax
	movl	%eax, hmi+56(%rip)
	movl	$0, %eax
	movq	%rax, hmi+48(%rip)
	movl	$0, %eax
	movq	%rax, hmi+392(%rip)
	movl	$0, %eax
	movq	%rax, hmi+432(%rip)
	movl	.LC41(%rip), %eax
	movl	%eax, hmi+468(%rip)
	movabsq	$4607182418800017408, %rax
	movq	%rax, hmi+64(%rip)
	movq	$0, -8(%rbp)
	jmp	.L171
.L172:
	movq	-8(%rbp), %rax
	movq	%rax, %rsi
	movl	$gv+117880, %edi
	call	_ZNSt6vectorIP8GrainBinSaIS1_EEixEm
	movq	(%rax), %rdx
	movl	$0, %eax
	movq	%rax, 10552(%rdx)
	addq	$1, -8(%rbp)
.L171:
	movl	$gv+117880, %edi
	call	_ZNKSt6vectorIP8GrainBinSaIS1_EE4sizeEv
	cmpq	-8(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L172
	nop
.L166:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4017:
	.size	_ZL12mole_h_fixupv, .-_ZL12mole_h_fixupv
	.local	_ZGVZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c
	.comm	_ZGVZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c,8,8
	.section	.rodata
.LC43:
	.string	"Failed: false"
	.text
	.type	_ZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPb, @function
_ZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPb:
.LFB4018:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4018
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r12
	pushq	%rbx
	subq	$112, %rsp
	.cfi_offset 12, -24
	.cfi_offset 3, -32
	movq	%rdi, -88(%rbp)
	movq	%rsi, -96(%rbp)
	movq	%rdx, -104(%rbp)
	movq	%rcx, -112(%rbp)
	movl	%r8d, %eax
	movq	%r9, -128(%rbp)
	movb	%al, -116(%rbp)
	movl	$_ZGVZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c, %eax
	movzbl	(%rax), %eax
	testb	%al, %al
	jne	.L175
	movl	$_ZGVZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c, %edi
	call	__cxa_guard_acquire
	testl	%eax, %eax
	setne	%al
	testb	%al, %al
	je	.L175
	movl	$0, %r12d
	movl	mole_global+56(%rip), %eax
	movslq	%eax, %rdx
	movl	mole_global+56(%rip), %eax
	cltq
	movq	%rax, %rsi
	movl	$_ZZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c, %edi
.LEHB27:
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC1Emm
.LEHE27:
	movl	$_ZGVZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c, %edi
	call	__cxa_guard_release
	movl	$__dso_handle, %edx
	movl	$_ZZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c, %esi
	movl	$_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED1Ev, %edi
	call	__cxa_atexit
.L175:
	movb	$0, -65(%rbp)
	movl	mole_global+56(%rip), %eax
	movslq	%eax, %rdx
	leaq	-64(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB28:
	call	_ZNSt8valarrayIdEC1Em
.LEHE28:
	leaq	-48(%rbp), %rax
	movl	$.LC43, %ecx
	movl	$348, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
.LEHB29:
	call	_ZN10bad_assertC1EPKclS1_
.LEHE29:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L176
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB30:
	call	_ZNK10bad_assert5printEv
	call	abort
.L176:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE30:
.L181:
	movq	%rax, %rbx
	testb	%r12b, %r12b
	jne	.L178
	movl	$_ZGVZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c, %edi
	call	__cxa_guard_abort
.L178:
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB31:
	call	_Unwind_Resume
.LEHE31:
.L182:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	jmp	.L180
.L183:
	movq	%rax, %rbx
.L180:
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB32:
	call	_Unwind_Resume
.LEHE32:
	.cfi_endproc
.LFE4018:
	.section	.gcc_except_table
.LLSDA4018:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4018-.LLSDACSB4018
.LLSDACSB4018:
	.uleb128 .LEHB27-.LFB4018
	.uleb128 .LEHE27-.LEHB27
	.uleb128 .L181-.LFB4018
	.uleb128 0
	.uleb128 .LEHB28-.LFB4018
	.uleb128 .LEHE28-.LEHB28
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB29-.LFB4018
	.uleb128 .LEHE29-.LEHB29
	.uleb128 .L183-.LFB4018
	.uleb128 0
	.uleb128 .LEHB30-.LFB4018
	.uleb128 .LEHE30-.LEHB30
	.uleb128 .L182-.LFB4018
	.uleb128 0
	.uleb128 .LEHB31-.LFB4018
	.uleb128 .LEHE31-.LEHB31
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB32-.LFB4018
	.uleb128 .LEHE32-.LEHB32
	.uleb128 0
	.uleb128 0
.LLSDACSE4018:
	.text
	.size	_ZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPb, .-_ZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPb
	.section	.text._ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED2Ev,"axG",@progbits,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED5Ev,comdat
	.align 2
	.weak	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED2Ev
	.type	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED2Ev, @function
_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED2Ev:
.LFB4021:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4021:
	.size	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED2Ev, .-_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED2Ev
	.weak	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED1Ev
	.set	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED1Ev,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED2Ev
	.text
	.type	_ZL13grouped_elemsPKdPd, @function
_ZL13grouped_elemsPKdPd:
.LFB4019:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4019
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$152, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -136(%rbp)
	movq	%rsi, -144(%rbp)
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
.LEHB33:
	call	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC1Ev
.LEHE33:
	movq	$0, -80(%rbp)
	jmp	.L187
.L188:
	movq	-80(%rbp), %rax
	leaq	0(,%rax,8), %rdx
	movq	-144(%rbp), %rax
	addq	%rax, %rdx
	movl	$0, %eax
	movq	%rax, (%rdx)
	movq	-80(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomE7get_ptrEv
	movq	%rax, -96(%rbp)
	leaq	-96(%rbp), %rdx
	leaq	-64(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB34:
	call	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_
	movq	-80(%rbp), %rdx
	movq	%rdx, (%rax)
	addq	$1, -80(%rbp)
.L187:
	movl	$atom_list, %edi
	call	_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv
	cmpq	-80(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L188
	movq	$0, -72(%rbp)
	jmp	.L189
.L194:
	movq	groupspecies(%rip), %rax
	movq	-72(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	_ZNKSs5emptyEv
	xorl	$1, %eax
	testb	%al, %al
	je	.L190
	jmp	.L191
.L190:
	movq	groupspecies(%rip), %rax
	movq	-72(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rax
	addq	$24, %rax
	movq	%rax, %rdi
	call	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	movq	%rax, -96(%rbp)
	leaq	-96(%rbp), %rdx
	leaq	-128(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1ERKSt17_Rb_tree_iteratorIS5_E
	jmp	.L192
.L193:
	movq	-72(%rbp), %rax
	leaq	0(,%rax,8), %rdx
	movq	-136(%rbp), %rax
	addq	%rdx, %rax
	movsd	(%rax), %xmm1
	movsd	%xmm1, -152(%rbp)
	leaq	-128(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	movl	16(%rax), %eax
	cvtsi2sd	%eax, %xmm0
	mulsd	-152(%rbp), %xmm0
	movsd	%xmm0, -152(%rbp)
	leaq	-128(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomE7get_ptrEv
	movq	%rax, -96(%rbp)
	leaq	-96(%rbp), %rdx
	leaq	-64(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_
	movq	(%rax), %rax
	movq	%rax, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-144(%rbp), %rdx
	addq	%rcx, %rdx
	leaq	0(,%rax,8), %rcx
	movq	-144(%rbp), %rax
	addq	%rcx, %rax
	movsd	(%rax), %xmm0
	addsd	-152(%rbp), %xmm0
	movsd	%xmm0, (%rdx)
	leaq	-128(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv
.L192:
	movq	groupspecies(%rip), %rax
	movq	-72(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rax
	addq	$24, %rax
	movq	%rax, %rdi
	call	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv
.LEHE34:
	movq	%rax, -112(%rbp)
	leaq	-112(%rbp), %rdx
	leaq	-96(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1ERKSt17_Rb_tree_iteratorIS5_E
	leaq	-96(%rbp), %rdx
	leaq	-128(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_
	testb	%al, %al
	jne	.L193
.L191:
	addq	$1, -72(%rbp)
.L189:
	movl	mole_global+64(%rip), %eax
	cltq
	cmpq	-72(%rbp), %rax
	jg	.L194
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
.LEHB35:
	call	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED1Ev
.LEHE35:
	jmp	.L197
.L196:
	movq	%rax, %rbx
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB36:
	call	_Unwind_Resume
.LEHE36:
.L197:
	addq	$152, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4019:
	.section	.gcc_except_table
.LLSDA4019:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4019-.LLSDACSB4019
.LLSDACSB4019:
	.uleb128 .LEHB33-.LFB4019
	.uleb128 .LEHE33-.LEHB33
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB34-.LFB4019
	.uleb128 .LEHE34-.LEHB34
	.uleb128 .L196-.LFB4019
	.uleb128 0
	.uleb128 .LEHB35-.LFB4019
	.uleb128 .LEHE35-.LEHB35
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB36-.LFB4019
	.uleb128 .LEHE36-.LEHB36
	.uleb128 0
	.uleb128 0
.LLSDACSE4019:
	.text
	.size	_ZL13grouped_elemsPKdPd, .-_ZL13grouped_elemsPKdPd
	.section	.rodata
.LC44:
	.string	"Failed: deut.lgElmtOn"
	.align 8
.LC46:
	.string	"PROBLEM molElems BAD  %li\t%s\t%.12e\t%.12e\t%.12e\n"
.LC47:
	.string	"Failed: lgTest"
	.text
	.align 2
	.globl	_ZN8GroupMap5setupEPd
	.type	_ZN8GroupMap5setupEPd, @function
_ZN8GroupMap5setupEPd:
.LFB4023:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4023
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r12
	pushq	%rbx
	subq	$288, %rsp
	.cfi_offset 12, -24
	.cfi_offset 3, -32
	movq	%rdi, -280(%rbp)
	movq	%rsi, -288(%rbp)
	movl	mole_global+56(%rip), %eax
	movslq	%eax, %rdx
	leaq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB37:
	call	_ZNSt8valarrayIdEC1Em
.LEHE37:
	movq	$0, -248(%rbp)
	jmp	.L199
.L200:
	movq	-248(%rbp), %rdx
	leaq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	%rax, %rbx
	movq	-248(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movq	40(%rax), %rax
	movq	%rax, (%rbx)
	addq	$1, -248(%rbp)
.L199:
	movl	mole_global+56(%rip), %eax
	cltq
	cmpq	-248(%rbp), %rax
	jg	.L200
	movq	$0, -240(%rbp)
	jmp	.L201
.L224:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cmpl	$-1, %eax
	setne	%al
	testb	%al, %al
	je	.L202
	movl	$0, %eax
	movq	%rax, -232(%rbp)
	movq	$0, -224(%rbp)
	jmp	.L203
.L205:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-224(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cmpl	$-1, %eax
	setne	%al
	testb	%al, %al
	je	.L204
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-224(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	movslq	%eax, %rdx
	leaq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movsd	(%rax), %xmm0
	movsd	-232(%rbp), %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -232(%rbp)
.L204:
	addq	$1, -224(%rbp)
.L203:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movq	%rax, %rdi
	call	_ZNKSt6vectorIiSaIiEE4sizeEv
	cmpq	-224(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L205
	movss	_ZL10SMALLFLOAT(%rip), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movsd	-232(%rbp), %xmm1
	ucomisd	%xmm0, %xmm1
	jbe	.L256
	movsd	.LC42(%rip), %xmm0
	divsd	-232(%rbp), %xmm0
	movsd	%xmm0, -168(%rbp)
	movq	$0, -216(%rbp)
	jmp	.L208
.L211:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-216(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cmpl	$-1, %eax
	setne	%al
	testb	%al, %al
	je	.L209
	movq	-280(%rbp), %rcx
	leaq	-144(%rbp), %rax
	movq	-240(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
.LEHB38:
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm
	movq	-216(%rbp), %rdx
	leaq	-144(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm
	movq	%rax, %rbx
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-216(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	movslq	%eax, %rdx
	leaq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movsd	(%rax), %xmm0
	mulsd	-168(%rbp), %xmm0
	movsd	%xmm0, (%rbx)
	jmp	.L210
.L209:
	movq	-280(%rbp), %rcx
	leaq	-112(%rbp), %rax
	movq	-240(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm
	movq	-216(%rbp), %rdx
	leaq	-112(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm
	movq	%rax, %rdx
	movl	$0, %eax
	movq	%rax, (%rdx)
.L210:
	addq	$1, -216(%rbp)
.L208:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movq	%rax, %rdi
	call	_ZNKSt6vectorIiSaIiEE4sizeEv
	cmpq	-216(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L211
	jmp	.L212
.L256:
	movb	$0, -262(%rbp)
	movq	$0, -208(%rbp)
	jmp	.L213
.L218:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-208(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cmpl	$-1, %eax
	je	.L214
	movzbl	-262(%rbp), %eax
	xorl	$1, %eax
	testb	%al, %al
	je	.L214
	movl	$1, %eax
	jmp	.L215
.L214:
	movl	$0, %eax
.L215:
	testb	%al, %al
	je	.L216
	movq	-280(%rbp), %rcx
	leaq	-80(%rbp), %rax
	movq	-240(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm
	movq	-208(%rbp), %rdx
	leaq	-80(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm
	movq	%rax, %rdx
	movabsq	$4607182418800017408, %rax
	movq	%rax, (%rdx)
	movb	$1, -262(%rbp)
	jmp	.L217
.L216:
	movq	-280(%rbp), %rcx
	leaq	-48(%rbp), %rax
	movq	-240(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm
	movq	-208(%rbp), %rdx
	leaq	-48(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm
	movq	%rax, %rdx
	movl	$0, %eax
	movq	%rax, (%rdx)
.L217:
	addq	$1, -208(%rbp)
.L213:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movq	%rax, %rdi
	call	_ZNKSt6vectorIiSaIiEE4sizeEv
	cmpq	-208(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L218
.L212:
	movb	$0, -262(%rbp)
	movq	$0, -200(%rbp)
	jmp	.L219
.L223:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cmpl	$-1, %eax
	setne	%al
	testb	%al, %al
	je	.L220
	movzbl	-262(%rbp), %eax
	xorl	$1, %eax
	testb	%al, %al
	je	.L221
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	movslq	%eax, %rdx
	leaq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	%rax, %rdx
	movq	-232(%rbp), %rax
	movq	%rax, (%rdx)
	jmp	.L222
.L221:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	movslq	%eax, %rdx
	leaq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	%rax, %rdx
	movl	$0, %eax
	movq	%rax, (%rdx)
.L222:
	movb	$1, -262(%rbp)
.L220:
	addq	$1, -200(%rbp)
.L219:
	movq	-240(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movq	%rax, %rdi
	call	_ZNKSt6vectorIiSaIiEE4sizeEv
	cmpq	-200(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L223
.L202:
	addq	$1, -240(%rbp)
.L201:
	movl	$atom_list, %edi
	call	_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv
	cmpq	-240(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L224
	movq	$0, -192(%rbp)
	jmp	.L225
.L226:
	movq	-192(%rbp), %rax
	leaq	0(,%rax,8), %rdx
	movq	-288(%rbp), %rax
	leaq	(%rdx,%rax), %rbx
	movq	groupspecies(%rip), %rax
	movq	-192(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rax
	movl	92(%rax), %eax
	movslq	%eax, %rdx
	leaq	-160(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	(%rax), %rax
	movq	%rax, (%rbx)
	addq	$1, -192(%rbp)
.L225:
	movl	mole_global+64(%rip), %eax
	cltq
	cmpq	-192(%rbp), %rax
	jg	.L226
	movq	-280(%rbp), %rax
	addq	$144, %rax
	movq	%rax, %rdi
	call	_Z7get_ptrIdEPT_RSt8valarrayIS0_E
	movq	%rax, %rbx
	movq	-288(%rbp), %rax
	movq	%rax, %rdi
	call	_Z7get_ptrIdEPT_S1_
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	_ZL13grouped_elemsPKdPd
	movq	$0, -184(%rbp)
	jmp	.L227
.L245:
	movl	$0, %eax
	movq	%rax, -176(%rbp)
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	(%rax), %rax
	movl	(%rax), %eax
	cmpl	$1, %eax
	jne	.L228
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movl	8(%rax), %eax
	cmpl	$2, %eax
	jne	.L228
	movl	$1, %eax
	jmp	.L229
.L228:
	movl	$0, %eax
.L229:
	testb	%al, %al
	je	.L230
	movzbl	deut(%rip), %eax
	xorl	$1, %eax
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L231
	leaq	-48(%rbp), %rax
	movl	$.LC44, %ecx
	movl	$643, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE38:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L232
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB39:
	call	_ZNK10bad_assert5printEv
	call	abort
.L232:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE39:
.L231:
	movss	deut+4(%rip), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm3
	movsd	%xmm3, -176(%rbp)
	jmp	.L233
.L230:
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	%rax, %rdi
	call	_ZNK9chem_atom15lgMeanAbundanceEv
	xorl	$1, %eax
	testb	%al, %al
	je	.L234
	jmp	.L235
.L234:
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	(%rax), %rax
	movl	(%rax), %eax
	subl	$1, %eax
	movl	%eax, -260(%rbp)
	movl	-260(%rbp), %eax
	cltq
	movss	dense(,%rax,4), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm4
	movsd	%xmm4, -176(%rbp)
.L233:
	movsd	.LC45(%rip), %xmm0
	ucomisd	-176(%rbp), %xmm0
	jbe	.L236
	movq	-280(%rbp), %rax
	leaq	144(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt8valarrayIdEixEm
	movsd	(%rax), %xmm1
	movsd	.LC45(%rip), %xmm0
	ucomisd	%xmm1, %xmm0
	ja	.L238
.L236:
	movq	-280(%rbp), %rax
	leaq	144(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt8valarrayIdEixEm
	movsd	(%rax), %xmm0
	subsd	-176(%rbp), %xmm0
	movsd	.LC4(%rip), %xmm1
	andpd	%xmm0, %xmm1
	movss	conv+6368(%rip), %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	mulsd	-176(%rbp), %xmm0
	ucomisd	%xmm1, %xmm0
	jb	.L257
.L238:
	movl	$1, %eax
	jmp	.L241
.L257:
	movl	$0, %eax
.L241:
	movb	%al, -261(%rbp)
	movzbl	-261(%rbp), %eax
	xorl	$1, %eax
	testb	%al, %al
	je	.L242
	movq	-280(%rbp), %rax
	leaq	144(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	(%rax), %r12
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	48(%rax), %rbx
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	%rax, %rdx
	leaq	-256(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB40:
	call	_ZNK9chem_atom5labelEv
.LEHE40:
	leaq	-256(%rbp), %rax
	movq	%rax, %rdi
.LEHB41:
	call	_ZNKSs5c_strEv
	movq	%rax, %rcx
	movq	ioQQQ(%rip), %rdi
	movq	-176(%rbp), %rax
	movq	-184(%rbp), %rdx
	movq	%r12, -296(%rbp)
	movsd	-296(%rbp), %xmm2
	movq	%rax, -296(%rbp)
	movsd	-296(%rbp), %xmm1
	movq	%rbx, -296(%rbp)
	movsd	-296(%rbp), %xmm0
	movl	$.LC46, %esi
	movl	$3, %eax
	call	fprintf
.LEHE41:
	leaq	-256(%rbp), %rax
	movq	%rax, %rdi
.LEHB42:
	call	_ZNSsD1Ev
.L242:
	movzbl	-261(%rbp), %eax
	xorl	$1, %eax
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L243
	leaq	-48(%rbp), %rax
	movl	$.LC47, %ecx
	movl	$661, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE42:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L244
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB43:
	call	_ZNK10bad_assert5printEv
	call	abort
.L244:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE43:
.L243:
	movq	-280(%rbp), %rax
	leaq	144(%rax), %rdx
	movq	-184(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt8valarrayIdEixEm
	movq	%rax, %rdx
	movq	-176(%rbp), %rax
	movq	%rax, (%rdx)
.L235:
	addq	$1, -184(%rbp)
.L227:
	movl	$atom_list, %edi
	call	_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv
	cmpq	-184(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L245
	leaq	-160(%rbp), %rax
	movq	%rax, %rdi
.LEHB44:
	call	_ZNSt8valarrayIdED1Ev
.LEHE44:
	jmp	.L258
.L251:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	jmp	.L247
.L252:
	movq	%rax, %rbx
	leaq	-256(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSsD1Ev
	jmp	.L247
.L253:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	jmp	.L247
.L250:
	movq	%rax, %rbx
.L247:
	leaq	-160(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB45:
	call	_Unwind_Resume
.LEHE45:
.L258:
	addq	$288, %rsp
	popq	%rbx
	popq	%r12
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4023:
	.section	.gcc_except_table
.LLSDA4023:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4023-.LLSDACSB4023
.LLSDACSB4023:
	.uleb128 .LEHB37-.LFB4023
	.uleb128 .LEHE37-.LEHB37
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB38-.LFB4023
	.uleb128 .LEHE38-.LEHB38
	.uleb128 .L250-.LFB4023
	.uleb128 0
	.uleb128 .LEHB39-.LFB4023
	.uleb128 .LEHE39-.LEHB39
	.uleb128 .L251-.LFB4023
	.uleb128 0
	.uleb128 .LEHB40-.LFB4023
	.uleb128 .LEHE40-.LEHB40
	.uleb128 .L250-.LFB4023
	.uleb128 0
	.uleb128 .LEHB41-.LFB4023
	.uleb128 .LEHE41-.LEHB41
	.uleb128 .L252-.LFB4023
	.uleb128 0
	.uleb128 .LEHB42-.LFB4023
	.uleb128 .LEHE42-.LEHB42
	.uleb128 .L250-.LFB4023
	.uleb128 0
	.uleb128 .LEHB43-.LFB4023
	.uleb128 .LEHE43-.LEHB43
	.uleb128 .L253-.LFB4023
	.uleb128 0
	.uleb128 .LEHB44-.LFB4023
	.uleb128 .LEHE44-.LEHB44
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB45-.LFB4023
	.uleb128 .LEHE45-.LEHB45
	.uleb128 0
	.uleb128 0
.LLSDACSE4023:
	.text
	.size	_ZN8GroupMap5setupEPd, .-_ZN8GroupMap5setupEPd
	.section	.rodata
	.align 8
.LC48:
	.string	"Failed: !mole_global.list[mol]->parentLabel.empty()"
	.align 8
.LC50:
	.string	"Failed: fabs(sum-grouptot) <= 1e-10 * fabs(grouptot)"
	.text
	.align 2
	.globl	_ZN8GroupMap15updateMoleculesERKSt8valarrayIdE
	.type	_ZN8GroupMap15updateMoleculesERKSt8valarrayIdE, @function
_ZN8GroupMap15updateMoleculesERKSt8valarrayIdE:
.LFB4024:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4024
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$152, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -136(%rbp)
	movq	%rsi, -144(%rbp)
	movq	$0, -112(%rbp)
	jmp	.L260
.L261:
	movq	-112(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movq	%rax, %rdx
	movl	$0, %eax
	movq	%rax, 40(%rdx)
	addq	$1, -112(%rbp)
.L260:
	movl	mole_global+60(%rip), %eax
	cltq
	cmpq	-112(%rbp), %rax
	jg	.L261
	movq	$0, -104(%rbp)
	jmp	.L262
.L263:
	movq	groupspecies(%rip), %rax
	movq	-104(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rax
	movl	92(%rax), %eax
	cltq
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movq	%rax, %rbx
	movq	-104(%rbp), %rdx
	movq	-144(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt8valarrayIdEixEm
	movq	(%rax), %rax
	movq	%rax, 40(%rbx)
	addq	$1, -104(%rbp)
.L262:
	movl	mole_global+64(%rip), %eax
	cltq
	cmpq	-104(%rbp), %rax
	jg	.L263
	movq	$0, -96(%rbp)
	jmp	.L264
.L271:
	movq	-96(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movl	8(%rax), %eax
	notl	%eax
	shrl	$31, %eax
	testb	%al, %al
	je	.L265
	movq	-96(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movq	%rax, %rdi
.LEHB46:
	call	_ZNKSs5emptyEv
	movzbl	%al, %eax
	testq	%rax, %rax
	setne	%al
	testb	%al, %al
	je	.L266
	leaq	-48(%rbp), %rax
	movl	$.LC48, %ecx
	movl	$685, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE46:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L267
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB47:
	call	_ZNK10bad_assert5printEv
	call	abort
.L267:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE47:
.L266:
	movq	-96(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movl	8(%rax), %eax
	cltq
	movq	%rax, -64(%rbp)
	movq	-96(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movq	%rax, %rbx
	movq	-64(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movq	40(%rax), %rax
	movq	%rax, 40(%rbx)
	movq	-96(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	addq	$24, %rax
	movq	%rax, %rdi
.LEHB48:
	call	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	movq	%rax, -128(%rbp)
	jmp	.L268
.L270:
	leaq	-128(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	%rax, %rdi
	call	_ZNK9chem_atom15lgMeanAbundanceEv
	xorl	$1, %eax
	testb	%al, %al
	je	.L269
	leaq	-128(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	movl	16(%rax), %ebx
	leaq	-128(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	48(%rax), %rax
	movl	%ebx, %edi
	movq	%rax, -152(%rbp)
	movsd	-152(%rbp), %xmm0
	call	_ZSt3powdi
	movsd	%xmm0, -152(%rbp)
	movq	-96(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movsd	40(%rax), %xmm0
	mulsd	-152(%rbp), %xmm0
	movsd	%xmm0, 40(%rax)
.L269:
	leaq	-128(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv
.L268:
	movq	-96(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	addq	$24, %rax
	movq	%rax, %rdi
	call	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv
	movq	%rax, -48(%rbp)
	leaq	-48(%rbp), %rdx
	leaq	-128(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_
	testb	%al, %al
	jne	.L270
.L265:
	addq	$1, -96(%rbp)
.L264:
	movl	mole_global+60(%rip), %eax
	cltq
	cmpq	-96(%rbp), %rax
	jg	.L271
	movq	$0, -88(%rbp)
	jmp	.L272
.L279:
	movq	-88(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cmpl	$-1, %eax
	setne	%al
	testb	%al, %al
	je	.L273
	movq	-88(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cltq
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movq	40(%rax), %rax
	movq	%rax, -56(%rbp)
	movl	$0, %eax
	movq	%rax, -80(%rbp)
	movq	$0, -72(%rbp)
	jmp	.L274
.L276:
	movq	-88(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cmpl	$-1, %eax
	setne	%al
	testb	%al, %al
	je	.L275
	movq	-88(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cltq
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movq	%rax, %rbx
	movq	-136(%rbp), %rcx
	leaq	-48(%rbp), %rax
	movq	-88(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm
	movq	-72(%rbp), %rdx
	leaq	-48(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm
	movsd	(%rax), %xmm0
	mulsd	-56(%rbp), %xmm0
	movsd	%xmm0, 40(%rbx)
	movq	-88(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	leaq	16(%rax), %rdx
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movq	%rdx, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cltq
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movsd	40(%rax), %xmm0
	movsd	-80(%rbp), %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -80(%rbp)
.L275:
	addq	$1, -72(%rbp)
.L274:
	movq	-88(%rbp), %rax
	movq	%rax, %rsi
	movl	$atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movq	%rax, %rdi
	call	_ZNKSt6vectorIiSaIiEE4sizeEv
	cmpq	-72(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L276
	movsd	-80(%rbp), %xmm0
	subsd	-56(%rbp), %xmm0
	movsd	.LC4(%rip), %xmm1
	andpd	%xmm0, %xmm1
	movsd	-56(%rbp), %xmm2
	movsd	.LC4(%rip), %xmm0
	andpd	%xmm2, %xmm0
	movsd	.LC49(%rip), %xmm2
	mulsd	%xmm2, %xmm0
	ucomisd	%xmm1, %xmm0
	setnb	%al
	xorl	$1, %eax
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L277
	leaq	-48(%rbp), %rax
	movl	$.LC50, %ecx
	movl	$713, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE48:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L278
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB49:
	call	_ZNK10bad_assert5printEv
	call	abort
.L278:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE49:
.L277:
.L273:
	addq	$1, -88(%rbp)
.L272:
	movl	$atom_list, %edi
	call	_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv
	cmpq	-88(%rbp), %rax
	seta	%al
	testb	%al, %al
	jne	.L279
	movl	$mole, %edi
.LEHB50:
	call	_ZN12t_mole_local22set_isotope_abundancesEv
	jmp	.L285
.L283:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.L284:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE50:
.L285:
	addq	$152, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4024:
	.section	.gcc_except_table
.LLSDA4024:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4024-.LLSDACSB4024
.LLSDACSB4024:
	.uleb128 .LEHB46-.LFB4024
	.uleb128 .LEHE46-.LEHB46
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB47-.LFB4024
	.uleb128 .LEHE47-.LEHB47
	.uleb128 .L283-.LFB4024
	.uleb128 0
	.uleb128 .LEHB48-.LFB4024
	.uleb128 .LEHE48-.LEHB48
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB49-.LFB4024
	.uleb128 .LEHE49-.LEHB49
	.uleb128 .L284-.LFB4024
	.uleb128 0
	.uleb128 .LEHB50-.LFB4024
	.uleb128 .LEHE50-.LEHB50
	.uleb128 0
	.uleb128 0
.LLSDACSE4024:
	.text
	.size	_ZN8GroupMap15updateMoleculesERKSt8valarrayIdE, .-_ZN8GroupMap15updateMoleculesERKSt8valarrayIdE
	.section	.rodata
	.align 8
.LC51:
	.string	"Failed: mole_global.list[i]->isMonatomic()"
	.align 8
.LC52:
	.string	"Failed: (int)mole_global.list[i]->nAtom.size() == 1"
	.text
	.type	_ZL25mole_eval_dynamic_balancelPdbR9multi_arrIdLi2EL10mem_layout0ELb0EE, @function
_ZL25mole_eval_dynamic_balancelPdbR9multi_arrIdLi2EL10mem_layout0ELb0EE:
.LFB4025:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4025
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$152, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -120(%rbp)
	movq	%rsi, -128(%rbp)
	movl	%edx, %eax
	movq	%rcx, -144(%rbp)
	movb	%al, -132(%rbp)
	movzbl	-132(%rbp), %edx
	movq	-144(%rbp), %rcx
	movq	-128(%rbp), %rsi
	movq	-120(%rbp), %rax
	movq	%rax, %rdi
.LEHB51:
	call	_Z17mole_eval_balancelPdbR9multi_arrIdLi2EL10mem_layout0ELb0EE
	movq	dynamics+160(%rip), %rax
	leaq	1(%rax), %rdx
	movq	iteration(%rip), %rax
	cmpq	%rax, %rdx
	jg	.L286
	movzbl	dynamics(%rip), %eax
	testb	%al, %al
	je	.L286
	movsd	dynamics+48(%rip), %xmm0
	xorpd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L317
	xorpd	%xmm1, %xmm1
	ucomisd	%xmm1, %xmm0
	jne	.L317
	jmp	.L286
.L317:
	movq	$0, -72(%rbp)
	jmp	.L289
.L309:
	cmpb	$0, -132(%rbp)
	je	.L290
	movq	-72(%rbp), %rbx
	movq	-72(%rbp), %rdx
	leaq	-48(%rbp), %rax
	movq	-144(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm
	leaq	-48(%rbp), %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm
	movsd	(%rax), %xmm0
	movsd	dynamics+48(%rip), %xmm1
	subsd	%xmm1, %xmm0
	movsd	%xmm0, (%rax)
.L290:
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movq	%rax, %rdi
	call	_ZNKSs5emptyEv
.LEHE51:
	xorl	$1, %eax
	testb	%al, %al
	je	.L291
	jmp	.L292
.L291:
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole+56, %edi
	call	_ZNSt8valarrayI8molezoneEixEm
	movsd	40(%rax), %xmm1
	movsd	dynamics+48(%rip), %xmm0
	mulsd	%xmm1, %xmm0
	movq	-72(%rbp), %rax
	leaq	0(,%rax,8), %rdx
	movq	-128(%rbp), %rax
	addq	%rdx, %rax
	movq	-72(%rbp), %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-128(%rbp), %rdx
	addq	%rcx, %rdx
	movsd	(%rdx), %xmm1
	subsd	%xmm0, %xmm1
	movapd	%xmm1, %xmm0
	movsd	%xmm0, (%rax)
	movl	$0, %ebx
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movq	%rax, %rdi
.LEHB52:
	call	_ZNK8molecule11isMonatomicEv
	xorl	$1, %eax
	testb	%al, %al
	jne	.L293
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movl	72(%rax), %eax
	testl	%eax, %eax
	js	.L293
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movzbl	77(%rax), %eax
	xorl	$1, %eax
	testb	%al, %al
	jne	.L293
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movq	%rax, %rdi
	call	_ZNK8molecule11isMonatomicEv
	testb	%al, %al
	je	.L294
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	addq	$24, %rax
	movq	%rax, %rdi
	call	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	movq	%rax, -112(%rbp)
	movl	$1, %ebx
	leaq	-112(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
.LEHE52:
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	%rax, %rdi
	call	_ZNK9chem_atom15lgMeanAbundanceEv
	xorl	$1, %eax
	testb	%al, %al
	je	.L294
.L293:
	movl	$1, %eax
	jmp	.L295
.L294:
	movl	$0, %eax
.L295:
	testb	%bl, %bl
	testb	%al, %al
	je	.L297
	movq	-72(%rbp), %rax
	leaq	0(,%rax,8), %rdx
	movq	-128(%rbp), %rax
	addq	%rdx, %rax
	movq	-72(%rbp), %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-128(%rbp), %rdx
	addq	%rcx, %rdx
	movsd	(%rdx), %xmm1
	movq	dynamics+72(%rip), %rdx
	movq	-72(%rbp), %rcx
	salq	$3, %rcx
	addq	%rcx, %rdx
	movsd	(%rdx), %xmm2
	movsd	dynamics+48(%rip), %xmm0
	mulsd	%xmm2, %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, (%rax)
	jmp	.L292
.L297:
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movl	72(%rax), %eax
	testl	%eax, %eax
	sete	%al
	testb	%al, %al
	je	.L292
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	movq	%rax, %rdi
.LEHB53:
	call	_ZNK8molecule11isMonatomicEv
	xorl	$1, %eax
	movzbl	%al, %eax
	testq	%rax, %rax
	setne	%al
	testb	%al, %al
	je	.L298
	leaq	-48(%rbp), %rax
	movl	$.LC51, %ecx
	movl	$753, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE53:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L299
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB54:
	call	_ZNK10bad_assert5printEv
	call	abort
.L299:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE54:
.L298:
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	addq	$24, %rax
	movq	%rax, %rdi
.LEHB55:
	call	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv
	cmpl	$1, %eax
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	setne	%al
	testb	%al, %al
	je	.L300
	leaq	-48(%rbp), %rax
	movl	$.LC52, %ecx
	movl	$754, %edx
	movl	$.LC11, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE55:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L301
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB56:
	call	_ZNK10bad_assert5printEv
	call	abort
.L301:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE56:
.L300:
	movq	-72(%rbp), %rax
	movq	%rax, %rsi
	movl	$mole_global+72, %edi
	call	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI8moleculeEptEv
	addq	$24, %rax
	movq	%rax, %rdi
.LEHB57:
	call	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	movq	%rax, -96(%rbp)
	leaq	-96(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	movq	%rax, %rdx
	leaq	-48(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZN9count_ptrI9chem_atomEC1ERKS1_
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movq	(%rax), %rax
	movl	(%rax), %eax
	subl	$1, %eax
	cltq
	movq	%rax, -56(%rbp)
	cmpq	$29, -56(%rbp)
	jle	.L302
	movl	$0, %ebx
	jmp	.L303
.L302:
	movl	$0, %eax
	movq	%rax, -80(%rbp)
	movq	-56(%rbp), %rax
	addq	$48, %rax
	movq	dense+8(,%rax,8), %rax
	movq	%rax, -64(%rbp)
	jmp	.L304
.L306:
	movq	dynamics+56(%rip), %rax
	movq	-56(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rax
	movq	-64(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movsd	(%rax), %xmm3
	movsd	%xmm3, -152(%rbp)
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movsd	48(%rax), %xmm0
	mulsd	-152(%rbp), %xmm0
	movsd	-80(%rbp), %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -80(%rbp)
	movq	-64(%rbp), %rbx
	movq	-56(%rbp), %rax
	movq	%rax, %rsi
	movl	$unresolved_atom_list, %edi
	call	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	addq	$16, %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	_ZNSt6vectorIiSaIiEEixEm
	movl	(%rax), %eax
	cmpl	$-1, %eax
	sete	%al
	testb	%al, %al
	je	.L305
	movq	-56(%rbp), %rdx
	movq	%rdx, %rax
	salq	$5, %rax
	subq	%rdx, %rax
	movq	-64(%rbp), %rdx
	addq	%rdx, %rax
	addq	$110, %rax
	movsd	dense+8(,%rax,8), %xmm1
	movsd	dynamics+48(%rip), %xmm0
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, -152(%rbp)
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNK9count_ptrI9chem_atomEptEv
	movsd	48(%rax), %xmm0
	mulsd	-152(%rbp), %xmm0
	movsd	-80(%rbp), %xmm1
	subsd	%xmm0, %xmm1
	movapd	%xmm1, %xmm0
	movsd	%xmm0, -80(%rbp)
.L305:
	addq	$1, -64(%rbp)
.L304:
	movq	-56(%rbp), %rax
	addq	$80, %rax
	movq	dense(,%rax,8), %rax
	cmpq	-64(%rbp), %rax
	jge	.L306
	movq	-72(%rbp), %rax
	leaq	0(,%rax,8), %rdx
	movq	-128(%rbp), %rax
	addq	%rdx, %rax
	movq	-72(%rbp), %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-128(%rbp), %rdx
	addq	%rcx, %rdx
	movsd	(%rdx), %xmm0
	addsd	-80(%rbp), %xmm0
	movsd	%xmm0, (%rax)
	movl	$1, %ebx
.L303:
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9count_ptrI9chem_atomED1Ev
	cmpl	$1, %ebx
	jne	.L292
	nop
.L292:
	addq	$1, -72(%rbp)
.L289:
	movl	mole_global+60(%rip), %eax
	cltq
	cmpq	-72(%rbp), %rax
	jg	.L309
	jmp	.L286
.L314:
	testb	%bl, %bl
	nop
	movq	%rax, %rdi
	call	_Unwind_Resume
.L315:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.L316:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE57:
.L286:
	addq	$152, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4025:
	.section	.gcc_except_table
.LLSDA4025:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4025-.LLSDACSB4025
.LLSDACSB4025:
	.uleb128 .LEHB51-.LFB4025
	.uleb128 .LEHE51-.LEHB51
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB52-.LFB4025
	.uleb128 .LEHE52-.LEHB52
	.uleb128 .L314-.LFB4025
	.uleb128 0
	.uleb128 .LEHB53-.LFB4025
	.uleb128 .LEHE53-.LEHB53
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB54-.LFB4025
	.uleb128 .LEHE54-.LEHB54
	.uleb128 .L315-.LFB4025
	.uleb128 0
	.uleb128 .LEHB55-.LFB4025
	.uleb128 .LEHE55-.LEHB55
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB56-.LFB4025
	.uleb128 .LEHE56-.LEHB56
	.uleb128 .L316-.LFB4025
	.uleb128 0
	.uleb128 .LEHB57-.LFB4025
	.uleb128 .LEHE57-.LEHB57
	.uleb128 0
	.uleb128 0
.LLSDACSE4025:
	.text
	.size	_ZL25mole_eval_dynamic_balancelPdbR9multi_arrIdLi2EL10mem_layout0ELb0EE, .-_ZL25mole_eval_dynamic_balancelPdbR9multi_arrIdLi2EL10mem_layout0ELb0EE
	.section	.text._Z4pow2IdET_S0_,"axG",@progbits,_Z4pow2IdET_S0_,comdat
	.weak	_Z4pow2IdET_S0_
	.type	_Z4pow2IdET_S0_, @function
_Z4pow2IdET_S0_:
.LFB4071:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movsd	%xmm0, -8(%rbp)
	movsd	-8(%rbp), %xmm0
	mulsd	-8(%rbp), %xmm0
	movsd	%xmm0, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4071:
	.size	_Z4pow2IdET_S0_, .-_Z4pow2IdET_S0_
	.section	.text._Z4pow3IdET_S0_,"axG",@progbits,_Z4pow3IdET_S0_,comdat
	.weak	_Z4pow3IdET_S0_
	.type	_Z4pow3IdET_S0_, @function
_Z4pow3IdET_S0_:
.LFB4072:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movsd	%xmm0, -8(%rbp)
	movsd	-8(%rbp), %xmm0
	mulsd	-8(%rbp), %xmm0
	mulsd	-8(%rbp), %xmm0
	movsd	%xmm0, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4072:
	.size	_Z4pow3IdET_S0_, .-_Z4pow3IdET_S0_
	.section	.text._ZNSt6vectorIiSaIiEED2Ev,"axG",@progbits,_ZNSt6vectorIiSaIiEED5Ev,comdat
	.align 2
	.weak	_ZNSt6vectorIiSaIiEED2Ev
	.type	_ZNSt6vectorIiSaIiEED2Ev, @function
_ZNSt6vectorIiSaIiEED2Ev:
.LFB4092:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4092
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt12_Vector_baseIiSaIiEE19_M_get_Tp_allocatorEv
	movq	%rax, %rdx
	movq	-24(%rbp), %rax
	movq	8(%rax), %rcx
	movq	-24(%rbp), %rax
	movq	(%rax), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
.LEHB58:
	call	_ZSt8_DestroyIPiiEvT_S1_RSaIT0_E
.LEHE58:
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
.LEHB59:
	call	_ZNSt12_Vector_baseIiSaIiEED2Ev
.LEHE59:
	jmp	.L327
.L326:
	movq	%rax, %rbx
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt12_Vector_baseIiSaIiEED2Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB60:
	call	_Unwind_Resume
.LEHE60:
.L327:
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4092:
	.section	.gcc_except_table
.LLSDA4092:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4092-.LLSDACSB4092
.LLSDACSB4092:
	.uleb128 .LEHB58-.LFB4092
	.uleb128 .LEHE58-.LEHB58
	.uleb128 .L326-.LFB4092
	.uleb128 0
	.uleb128 .LEHB59-.LFB4092
	.uleb128 .LEHE59-.LEHB59
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB60-.LFB4092
	.uleb128 .LEHE60-.LEHB60
	.uleb128 0
	.uleb128 0
.LLSDACSE4092:
	.section	.text._ZNSt6vectorIiSaIiEED2Ev,"axG",@progbits,_ZNSt6vectorIiSaIiEED5Ev,comdat
	.size	_ZNSt6vectorIiSaIiEED2Ev, .-_ZNSt6vectorIiSaIiEED2Ev
	.weak	_ZNSt6vectorIiSaIiEED1Ev
	.set	_ZNSt6vectorIiSaIiEED1Ev,_ZNSt6vectorIiSaIiEED2Ev
	.section	.text._ZNSt6vectorIiSaIiEEixEm,"axG",@progbits,_ZNSt6vectorIiSaIiEEixEm,comdat
	.align 2
	.weak	_ZNSt6vectorIiSaIiEEixEm
	.type	_ZNSt6vectorIiSaIiEEixEm, @function
_ZNSt6vectorIiSaIiEEixEm:
.LFB4127:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$2, %rdx
	addq	%rdx, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4127:
	.size	_ZNSt6vectorIiSaIiEEixEm, .-_ZNSt6vectorIiSaIiEEixEm
	.section	.text._ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_,"axG",@progbits,_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_,comdat
	.weak	_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_
	.type	_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_, @function
_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_:
.LFB4234:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4234
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$56, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	%rdx, -56(%rbp)
	movq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt11char_traitsIcE6lengthEPKc
	movq	%rax, -24(%rbp)
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
.LEHB61:
	call	_ZNSsC1Ev
.LEHE61:
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
.LEHB62:
	call	_ZNKSs4sizeEv
	movq	-24(%rbp), %rdx
	addq	%rax, %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSs7reserveEm
	movq	-24(%rbp), %rdx
	movq	-48(%rbp), %rcx
	movq	-40(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSs6appendEPKcm
	movq	-56(%rbp), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSs6appendERKSs
.LEHE62:
	jmp	.L334
.L333:
	movq	%rax, %rbx
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSsD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB63:
	call	_Unwind_Resume
.LEHE63:
.L334:
	movq	-40(%rbp), %rax
	addq	$56, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4234:
	.section	.gcc_except_table
.LLSDA4234:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4234-.LLSDACSB4234
.LLSDACSB4234:
	.uleb128 .LEHB61-.LFB4234
	.uleb128 .LEHE61-.LEHB61
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB62-.LFB4234
	.uleb128 .LEHE62-.LEHB62
	.uleb128 .L333-.LFB4234
	.uleb128 0
	.uleb128 .LEHB63-.LFB4234
	.uleb128 .LEHE63-.LEHB63
	.uleb128 0
	.uleb128 0
.LLSDACSE4234:
	.section	.text._ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_,"axG",@progbits,_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_,comdat
	.size	_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_, .-_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_EPKS3_RKS6_
	.section	.text._ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv,"axG",@progbits,_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv,comdat
	.align 2
	.weak	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	.type	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv, @function
_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv:
.LFB4236:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4236:
	.size	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv, .-_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	.section	.text._ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_,"axG",@progbits,_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_,comdat
	.align 2
	.weak	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_
	.type	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_, @function
_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_:
.LFB4238:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rdx
	movq	-16(%rbp), %rax
	movq	(%rax), %rax
	cmpq	%rax, %rdx
	setne	%al
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4238:
	.size	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_, .-_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_
	.section	.text._ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv,"axG",@progbits,_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv,comdat
	.align 2
	.weak	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv
	.type	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv, @function
_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv:
.LFB4239:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	_ZSt18_Rb_tree_incrementPKSt18_Rb_tree_node_base
	movq	-8(%rbp), %rdx
	movq	%rax, (%rdx)
	movq	-8(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4239:
	.size	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv, .-_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv
	.section	.text._ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv,"axG",@progbits,_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv,comdat
	.align 2
	.weak	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	.type	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv, @function
_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv:
.LFB4240:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	addq	$32, %rax
	movq	%rax, %rdi
	call	_ZSt11__addressofIKSt4pairIK9count_ptrI9chem_atomEiEEPT_RS7_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4240:
	.size	_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv, .-_ZNKSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	.section	.text._ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv,"axG",@progbits,_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv,comdat
	.align 2
	.weak	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv
	.type	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv, @function
_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv:
.LFB4241:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE4sizeEv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4241:
	.size	_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv, .-_ZNKSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE4sizeEv
	.section	.text._ZNK9count_ptrI9chem_atomE7get_ptrEv,"axG",@progbits,_ZNK9count_ptrI9chem_atomE7get_ptrEv,comdat
	.align 2
	.weak	_ZNK9count_ptrI9chem_atomE7get_ptrEv
	.type	_ZNK9count_ptrI9chem_atomE7get_ptrEv, @function
_ZNK9count_ptrI9chem_atomE7get_ptrEv:
.LFB4253:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4253:
	.size	_ZNK9count_ptrI9chem_atomE7get_ptrEv, .-_ZNK9count_ptrI9chem_atomE7get_ptrEv
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Ev,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC5Ev,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Ev
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Ev, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Ev:
.LFB4268:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EEC1Ev
	movq	-8(%rbp), %rax
	addq	$80, %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEC1Ev
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4268:
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Ev, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Ev
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC1Ev
	.set	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC1Ev,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Ev
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EED2Ev,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED5Ev,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED2Ev
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED2Ev, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED2Ev:
.LFB4271:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4271
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
.LEHB64:
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev
.LEHE64:
	movq	-24(%rbp), %rax
	addq	$80, %rax
	movq	%rax, %rdi
.LEHB65:
	call	_ZNSt8valarrayIdED1Ev
.LEHE65:
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
.LEHB66:
	call	_ZN10multi_geomILi2EL10mem_layout0EED1Ev
.LEHE66:
	jmp	.L354
.L352:
	movq	%rax, %rbx
	movq	-24(%rbp), %rax
	addq	$80, %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdED1Ev
	jmp	.L351
.L353:
	movq	%rax, %rbx
.L351:
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB67:
	call	_Unwind_Resume
.LEHE67:
.L354:
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4271:
	.section	.gcc_except_table
.LLSDA4271:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4271-.LLSDACSB4271
.LLSDACSB4271:
	.uleb128 .LEHB64-.LFB4271
	.uleb128 .LEHE64-.LEHB64
	.uleb128 .L352-.LFB4271
	.uleb128 0
	.uleb128 .LEHB65-.LFB4271
	.uleb128 .LEHE65-.LEHB65
	.uleb128 .L353-.LFB4271
	.uleb128 0
	.uleb128 .LEHB66-.LFB4271
	.uleb128 .LEHE66-.LEHB66
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB67-.LFB4271
	.uleb128 .LEHE67-.LEHB67
	.uleb128 0
	.uleb128 0
.LLSDACSE4271:
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EED2Ev,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED5Ev,comdat
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED2Ev, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED2Ev
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED1Ev
	.set	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED1Ev,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EED2Ev
	.section	.text._ZNSt8valarrayIdEC2Ev,"axG",@progbits,_ZNSt8valarrayIdEC5Ev,comdat
	.align 2
	.weak	_ZNSt8valarrayIdEC2Ev
	.type	_ZNSt8valarrayIdEC2Ev, @function
_ZNSt8valarrayIdEC2Ev:
.LFB4274:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	$0, (%rax)
	movq	-8(%rbp), %rax
	movq	$0, 8(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4274:
	.size	_ZNSt8valarrayIdEC2Ev, .-_ZNSt8valarrayIdEC2Ev
	.weak	_ZNSt8valarrayIdEC1Ev
	.set	_ZNSt8valarrayIdEC1Ev,_ZNSt8valarrayIdEC2Ev
	.section	.text._ZNSt8valarrayIdED2Ev,"axG",@progbits,_ZNSt8valarrayIdED5Ev,comdat
	.align 2
	.weak	_ZNSt8valarrayIdED2Ev
	.type	_ZNSt8valarrayIdED2Ev, @function
_ZNSt8valarrayIdED2Ev:
.LFB4277:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZSt27__valarray_destroy_elementsIdEvPT_S1_
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rax, %rdi
	call	_ZSt25__valarray_release_memoryPv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4277:
	.size	_ZNSt8valarrayIdED2Ev, .-_ZNSt8valarrayIdED2Ev
	.weak	_ZNSt8valarrayIdED1Ev
	.set	_ZNSt8valarrayIdED1Ev,_ZNSt8valarrayIdED2Ev
	.section	.text._ZNSt8valarrayIdE6resizeEmd,"axG",@progbits,_ZNSt8valarrayIdE6resizeEmd,comdat
	.align 2
	.weak	_ZNSt8valarrayIdE6resizeEmd
	.type	_ZNSt8valarrayIdE6resizeEmd, @function
_ZNSt8valarrayIdE6resizeEmd:
.LFB4279:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movsd	%xmm0, -24(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	salq	$3, %rax
	addq	%rax, %rdx
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZSt27__valarray_destroy_elementsIdEvPT_S1_
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	cmpq	-16(%rbp), %rax
	je	.L359
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rax, %rdi
	call	_ZSt25__valarray_release_memoryPv
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-16(%rbp), %rax
	movq	%rax, %rdi
	call	_ZSt22__valarray_get_storageIdErPT_m
	movq	-8(%rbp), %rdx
	movq	%rax, 8(%rdx)
.L359:
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rcx
	movq	-8(%rbp), %rax
	movq	8(%rax), %rdx
	movq	-24(%rbp), %rax
	movq	%rax, -32(%rbp)
	movsd	-32(%rbp), %xmm0
	movq	%rcx, %rsi
	movq	%rdx, %rdi
	call	_ZSt25__valarray_fill_constructIdEvPT_S1_S0_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4279:
	.size	_ZNSt8valarrayIdE6resizeEmd, .-_ZNSt8valarrayIdE6resizeEmd
	.section	.rodata
.LC53:
	.string	"Failed: vals().size() == 0"
.LC54:
	.string	"container_classes.h"
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm:
.LFB4280:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4280
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$56, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -56(%rbp)
	movq	%rsi, -64(%rbp)
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv
	movq	%rax, %rdi
	call	_ZNKSt8valarrayIdE4sizeEv
	testq	%rax, %rax
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	setne	%al
	testb	%al, %al
	je	.L361
	leaq	-48(%rbp), %rax
	movl	$.LC53, %ecx
	movl	$1082, %edx
	movl	$.LC54, %esi
	movq	%rax, %rdi
.LEHB68:
	call	_ZN10bad_assertC1EPKclS1_
.LEHE68:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L362
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB69:
	call	_ZNK10bad_assert5printEv
	call	abort
.L362:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE69:
.L361:
	movq	-64(%rbp), %rax
	movq	%rax, -48(%rbp)
	movq	-56(%rbp), %rax
	leaq	-48(%rbp), %rdx
	movl	$1, %esi
	movq	%rax, %rdi
.LEHB70:
	call	_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm
	jmp	.L365
.L364:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE70:
.L365:
	addq	$56, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4280:
	.section	.gcc_except_table
.LLSDA4280:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4280-.LLSDACSB4280
.LLSDACSB4280:
	.uleb128 .LEHB68-.LFB4280
	.uleb128 .LEHE68-.LEHB68
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB69-.LFB4280
	.uleb128 .LEHE69-.LEHB69
	.uleb128 .L364-.LFB4280
	.uleb128 0
	.uleb128 .LEHB70-.LFB4280
	.uleb128 .LEHE70-.LEHB70
	.uleb128 0
	.uleb128 0
.LLSDACSE4280:
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm,comdat
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEm
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm:
.LFB4281:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4281
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$72, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -56(%rbp)
	movq	%rsi, -64(%rbp)
	movq	%rdx, -72(%rbp)
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv
	movq	%rax, %rdi
	call	_ZNKSt8valarrayIdE4sizeEv
	testq	%rax, %rax
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	setne	%al
	testb	%al, %al
	je	.L367
	leaq	-48(%rbp), %rax
	movl	$.LC53, %ecx
	movl	$1088, %edx
	movl	$.LC54, %esi
	movq	%rax, %rdi
.LEHB71:
	call	_ZN10bad_assertC1EPKclS1_
.LEHE71:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L368
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB72:
	call	_ZNK10bad_assert5printEv
	call	abort
.L368:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE72:
.L367:
	movq	-64(%rbp), %rax
	movq	%rax, -48(%rbp)
	movq	-72(%rbp), %rax
	movq	%rax, -40(%rbp)
	movq	-56(%rbp), %rax
	leaq	-48(%rbp), %rdx
	movl	$2, %esi
	movq	%rax, %rdi
.LEHB73:
	call	_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm
	jmp	.L371
.L370:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE73:
.L371:
	addq	$72, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4281:
	.section	.gcc_except_table
.LLSDA4281:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4281-.LLSDACSB4281
.LLSDACSB4281:
	.uleb128 .LEHB71-.LFB4281
	.uleb128 .LEHE71-.LEHB71
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB72-.LFB4281
	.uleb128 .LEHE72-.LEHB72
	.uleb128 .L370-.LFB4281
	.uleb128 0
	.uleb128 .LEHB73-.LFB4281
	.uleb128 .LEHE73-.LEHB73
	.uleb128 0
	.uleb128 0
.LLSDACSE4281:
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm,comdat
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE7reserveEmm
	.section	.rodata
.LC55:
	.string	"Failed: p_psl[dim] == NULL"
.LC56:
	.string	"Failed: p_dsl.size() == 0"
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv:
.LFB4282:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4282
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$104, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -104(%rbp)
	movq	-104(%rbp), %rax
	movq	%rax, %rdi
.LEHB74:
	call	_ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv
	movl	$0, -84(%rbp)
	jmp	.L373
.L383:
	movl	-84(%rbp), %eax
	cltq
	movq	$0, -64(%rbp,%rax,8)
	movl	-84(%rbp), %eax
	cltq
	movq	-64(%rbp,%rax,8), %rdx
	movl	-84(%rbp), %eax
	cltq
	movq	%rdx, -80(%rbp,%rax,8)
	cmpl	$1, -84(%rbp)
	je	.L374
	movq	-104(%rbp), %rax
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	addq	$8, %rdx
	movq	8(%rax,%rdx,8), %rax
	testq	%rax, %rax
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L375
	leaq	-48(%rbp), %rax
	movl	$.LC55, %ecx
	movl	$1134, %edx
	movl	$.LC54, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE74:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L376
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB75:
	call	_ZNK10bad_assert5printEv
	call	abort
.L376:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE75:
.L375:
	movq	-104(%rbp), %rax
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	addq	$6, %rdx
	movq	8(%rax,%rdx,8), %rax
	testq	%rax, %rax
	je	.L377
	movq	-104(%rbp), %rax
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	addq	$6, %rdx
	movq	8(%rax,%rdx,8), %rax
	movabsq	$1143914305352105984, %rdx
	cmpq	%rdx, %rax
	ja	.L378
	salq	$3, %rax
	jmp	.L379
.L378:
	movq	$-1, %rax
.L379:
	movq	%rax, %rdi
.LEHB76:
	call	_Znam
	movq	-104(%rbp), %rdx
	movl	-84(%rbp), %ecx
	movslq	%ecx, %rcx
	addq	$8, %rcx
	movq	%rax, 8(%rdx,%rcx,8)
.L377:
	jmp	.L380
.L374:
	movq	-104(%rbp), %rax
	addq	$80, %rax
	movq	%rax, %rdi
	call	_ZNKSt8valarrayIdE4sizeEv
	testq	%rax, %rax
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	setne	%al
	testb	%al, %al
	je	.L381
	leaq	-48(%rbp), %rax
	movl	$.LC56, %ecx
	movl	$1140, %edx
	movl	$.LC54, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE76:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L382
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB77:
	call	_ZNK10bad_assert5printEv
	call	abort
.L382:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE77:
.L381:
	movq	-104(%rbp), %rax
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	addq	$6, %rdx
	movq	8(%rax,%rdx,8), %rax
	movq	-104(%rbp), %rdx
	addq	$80, %rdx
	xorpd	%xmm0, %xmm0
	movq	%rax, %rsi
	movq	%rdx, %rdi
.LEHB78:
	call	_ZNSt8valarrayIdE6resizeEmd
.L380:
	addl	$1, -84(%rbp)
.L373:
	cmpl	$1, -84(%rbp)
	jle	.L383
	movq	-104(%rbp), %rcx
	leaq	-64(%rbp), %rdx
	leaq	-80(%rbp), %rsi
	movq	-104(%rbp), %rax
	movl	$0, %r8d
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci
	movq	-104(%rbp), %rax
	movq	72(%rax), %rdx
	movq	-104(%rbp), %rax
	movq	%rdx, 96(%rax)
	movq	-104(%rbp), %rax
	movq	96(%rax), %rdx
	movq	-104(%rbp), %rax
	movq	%rdx, 104(%rax)
	movq	-104(%rbp), %rax
	movq	96(%rax), %rdx
	movq	-104(%rbp), %rax
	movq	%rdx, 112(%rax)
	movq	-104(%rbp), %rax
	movq	96(%rax), %rdx
	movq	-104(%rbp), %rax
	movq	%rdx, 120(%rax)
	movq	-104(%rbp), %rax
	movq	96(%rax), %rdx
	movq	-104(%rbp), %rax
	movq	%rdx, 128(%rax)
	movq	-104(%rbp), %rax
	movq	96(%rax), %rdx
	movq	-104(%rbp), %rax
	movq	%rdx, 136(%rax)
	jmp	.L388
.L386:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.L387:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE78:
.L388:
	addq	$104, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4282:
	.section	.gcc_except_table
.LLSDA4282:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4282-.LLSDACSB4282
.LLSDACSB4282:
	.uleb128 .LEHB74-.LFB4282
	.uleb128 .LEHE74-.LEHB74
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB75-.LFB4282
	.uleb128 .LEHE75-.LEHB75
	.uleb128 .L386-.LFB4282
	.uleb128 0
	.uleb128 .LEHB76-.LFB4282
	.uleb128 .LEHE76-.LEHB76
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB77-.LFB4282
	.uleb128 .LEHE77-.LEHB77
	.uleb128 .L387-.LFB4282
	.uleb128 0
	.uleb128 .LEHB78-.LFB4282
	.uleb128 .LEHE78-.LEHB78
	.uleb128 0
	.uleb128 0
.LLSDACSE4282:
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv,comdat
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv
	.section	.text._ZNKSt6vectorIP8GrainBinSaIS1_EE4sizeEv,"axG",@progbits,_ZNKSt6vectorIP8GrainBinSaIS1_EE4sizeEv,comdat
	.align 2
	.weak	_ZNKSt6vectorIP8GrainBinSaIS1_EE4sizeEv
	.type	_ZNKSt6vectorIP8GrainBinSaIS1_EE4sizeEv, @function
_ZNKSt6vectorIP8GrainBinSaIS1_EE4sizeEv:
.LFB4334:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rax, %rdx
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	subq	%rax, %rdx
	movq	%rdx, %rax
	sarq	$3, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4334:
	.size	_ZNKSt6vectorIP8GrainBinSaIS1_EE4sizeEv, .-_ZNKSt6vectorIP8GrainBinSaIS1_EE4sizeEv
	.section	.text._ZNSt8valarrayIdEC2Em,"axG",@progbits,_ZNSt8valarrayIdEC5Em,comdat
	.align 2
	.weak	_ZNSt8valarrayIdEC2Em
	.type	_ZNSt8valarrayIdEC2Em, @function
_ZNSt8valarrayIdEC2Em:
.LFB4433:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-16(%rbp), %rax
	movq	%rax, %rdi
	call	_ZSt22__valarray_get_storageIdErPT_m
	movq	-8(%rbp), %rdx
	movq	%rax, 8(%rdx)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$3, %rdx
	addq	%rax, %rdx
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZSt28__valarray_default_constructIdEvPT_S1_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4433:
	.size	_ZNSt8valarrayIdEC2Em, .-_ZNSt8valarrayIdEC2Em
	.weak	_ZNSt8valarrayIdEC1Em
	.set	_ZNSt8valarrayIdEC1Em,_ZNSt8valarrayIdEC2Em
	.section	.text._ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv,"axG",@progbits,_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv,comdat
	.align 2
	.weak	_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv
	.type	_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv, @function
_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv:
.LFB4435:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rax, %rdx
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	subq	%rax, %rdx
	movq	%rdx, %rax
	sarq	$4, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4435:
	.size	_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv, .-_ZNKSt6vectorI9count_ptrI9chem_atomESaIS2_EE4sizeEv
	.section	.text._ZNSt8valarrayI8molezoneEixEm,"axG",@progbits,_ZNSt8valarrayI8molezoneEixEm,comdat
	.align 2
	.weak	_ZNSt8valarrayI8molezoneEixEm
	.type	_ZNSt8valarrayI8molezoneEixEm, @function
_ZNSt8valarrayI8molezoneEixEm:
.LFB4436:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$6, %rdx
	addq	%rdx, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4436:
	.size	_ZNSt8valarrayI8molezoneEixEm, .-_ZNSt8valarrayI8molezoneEixEm
	.section	.text._Z7get_ptrIdEPT_RSt8valarrayIS0_E,"axG",@progbits,_Z7get_ptrIdEPT_RSt8valarrayIS0_E,comdat
	.weak	_Z7get_ptrIdEPT_RSt8valarrayIS0_E
	.type	_Z7get_ptrIdEPT_RSt8valarrayIS0_E, @function
_Z7get_ptrIdEPT_RSt8valarrayIS0_E:
.LFB4437:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4437:
	.size	_Z7get_ptrIdEPT_RSt8valarrayIS0_E, .-_Z7get_ptrIdEPT_RSt8valarrayIS0_E
	.section	.text._ZNSt8valarrayIdEixEm,"axG",@progbits,_ZNSt8valarrayIdEixEm,comdat
	.align 2
	.weak	_ZNSt8valarrayIdEixEm
	.type	_ZNSt8valarrayIdEixEm, @function
_ZNSt8valarrayIdEixEm:
.LFB4438:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4438:
	.size	_ZNSt8valarrayIdEixEm, .-_ZNSt8valarrayIdEixEm
	.section	.text._ZNSt6vectorIP8GrainBinSaIS1_EEixEm,"axG",@progbits,_ZNSt6vectorIP8GrainBinSaIS1_EEixEm,comdat
	.align 2
	.weak	_ZNSt6vectorIP8GrainBinSaIS1_EEixEm
	.type	_ZNSt6vectorIP8GrainBinSaIS1_EEixEm, @function
_ZNSt6vectorIP8GrainBinSaIS1_EEixEm:
.LFB4439:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4439:
	.size	_ZNSt6vectorIP8GrainBinSaIS1_EEixEm, .-_ZNSt6vectorIP8GrainBinSaIS1_EEixEm
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Emm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC5Emm,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Emm
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Emm, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Emm:
.LFB4441:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4441
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$56, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	%rdx, -56(%rbp)
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
.LEHB79:
	call	_ZN10multi_geomILi2EL10mem_layout0EEC1Ev
.LEHE79:
	movq	-40(%rbp), %rax
	addq	$80, %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEC1Ev
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev
	movq	-48(%rbp), %rax
	movq	%rax, -32(%rbp)
	movq	-56(%rbp), %rax
	movq	%rax, -24(%rbp)
	leaq	-32(%rbp), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB80:
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm
.LEHE80:
	jmp	.L405
.L404:
	movq	%rax, %rbx
	movq	-40(%rbp), %rax
	addq	$80, %rax
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdED1Ev
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB81:
	call	_Unwind_Resume
.LEHE81:
.L405:
	addq	$56, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4441:
	.section	.gcc_except_table
.LLSDA4441:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4441-.LLSDACSB4441
.LLSDACSB4441:
	.uleb128 .LEHB79-.LFB4441
	.uleb128 .LEHE79-.LEHB79
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB80-.LFB4441
	.uleb128 .LEHE80-.LEHB80
	.uleb128 .L404-.LFB4441
	.uleb128 0
	.uleb128 .LEHB81-.LFB4441
	.uleb128 .LEHE81-.LEHB81
	.uleb128 0
	.uleb128 0
.LLSDACSE4441:
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Emm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC5Emm,comdat
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Emm, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Emm
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC1Emm
	.set	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC1Emm,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEC2Emm
	.section	.text._ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm,"axG",@progbits,_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm,comdat
	.align 2
	.weak	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	.type	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm, @function
_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm:
.LFB4443:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$4, %rdx
	addq	%rdx, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4443:
	.size	_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm, .-_ZNSt6vectorI9count_ptrI8moleculeESaIS2_EEixEm
	.section	.text._ZNK9count_ptrI8moleculeEptEv,"axG",@progbits,_ZNK9count_ptrI8moleculeEptEv,comdat
	.align 2
	.weak	_ZNK9count_ptrI8moleculeEptEv
	.type	_ZNK9count_ptrI8moleculeEptEv, @function
_ZNK9count_ptrI8moleculeEptEv:
.LFB4444:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4444:
	.size	_ZNK9count_ptrI8moleculeEptEv, .-_ZNK9count_ptrI8moleculeEptEv
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm:
.LFB4446:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$64, %rsp
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	%rdx, -56(%rbp)
	leaq	-32(%rbp), %rax
	movq	-48(%rbp), %rdx
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5n_ptrEv
	movq	-40(%rbp), %rax
	movq	-56(%rbp), %rdx
	leaq	-32(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNK9n_pointerIdLi2EL10mem_layout0ELb0EEixEm
	movq	-40(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4446:
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EEixEm
	.section	.text._ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm,"axG",@progbits,_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm,comdat
	.align 2
	.weak	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm
	.type	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm, @function
_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm:
.LFB4447:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4447:
	.size	_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm, .-_ZNK9n_pointerIdLi1EL10mem_layout0ELb0EEixEm
	.section	.text._ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm,"axG",@progbits,_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm,comdat
	.align 2
	.weak	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	.type	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm, @function
_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm:
.LFB4448:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$4, %rdx
	addq	%rdx, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4448:
	.size	_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm, .-_ZNSt6vectorI9count_ptrI9chem_atomESaIS2_EEixEm
	.section	.text._ZNK9count_ptrI9chem_atomEptEv,"axG",@progbits,_ZNK9count_ptrI9chem_atomEptEv,comdat
	.align 2
	.weak	_ZNK9count_ptrI9chem_atomEptEv
	.type	_ZNK9count_ptrI9chem_atomEptEv, @function
_ZNK9count_ptrI9chem_atomEptEv:
.LFB4449:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4449:
	.size	_ZNK9count_ptrI9chem_atomEptEv, .-_ZNK9count_ptrI9chem_atomEptEv
	.section	.text._ZNKSt6vectorIiSaIiEE4sizeEv,"axG",@progbits,_ZNKSt6vectorIiSaIiEE4sizeEv,comdat
	.align 2
	.weak	_ZNKSt6vectorIiSaIiEE4sizeEv
	.type	_ZNKSt6vectorIiSaIiEE4sizeEv, @function
_ZNKSt6vectorIiSaIiEE4sizeEv:
.LFB4450:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rax, %rdx
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	subq	%rax, %rdx
	movq	%rdx, %rax
	sarq	$2, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4450:
	.size	_ZNKSt6vectorIiSaIiEE4sizeEv, .-_ZNKSt6vectorIiSaIiEE4sizeEv
	.section	.text._ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv,"axG",@progbits,_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv,comdat
	.align 2
	.weak	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv
	.type	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv, @function
_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv:
.LFB4452:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE3endEv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4452:
	.size	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv, .-_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE3endEv
	.section	.text._ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_,"axG",@progbits,_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_,comdat
	.align 2
	.weak	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_
	.type	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_, @function
_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_:
.LFB4453:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rdx
	movq	-16(%rbp), %rax
	movq	(%rax), %rax
	cmpq	%rax, %rdx
	setne	%al
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4453:
	.size	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_, .-_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEneERKS6_
	.section	.text._ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC2Ev,"axG",@progbits,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC5Ev,comdat
	.align 2
	.weak	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC2Ev
	.type	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC2Ev, @function
_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC2Ev:
.LFB4460:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4460:
	.size	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC2Ev, .-_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC2Ev
	.weak	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC1Ev
	.set	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC1Ev,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEC2Ev
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED2Ev,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED5Ev,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED2Ev
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED2Ev, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED2Ev:
.LFB4464:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4464:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED2Ev, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED2Ev
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED1Ev
	.set	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED1Ev,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED2Ev
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED2Ev,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED5Ev,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED2Ev
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED2Ev, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED2Ev:
.LFB4466:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4466
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv
	movq	%rax, %rdx
	movq	-24(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB82:
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E
.LEHE82:
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED1Ev
	jmp	.L431
.L430:
	movq	%rax, %rbx
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EED1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB83:
	call	_Unwind_Resume
.LEHE83:
.L431:
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4466:
	.section	.gcc_except_table
.LLSDA4466:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4466-.LLSDACSB4466
.LLSDACSB4466:
	.uleb128 .LEHB82-.LFB4466
	.uleb128 .LEHE82-.LEHB82
	.uleb128 .L430-.LFB4466
	.uleb128 0
	.uleb128 .LEHB83-.LFB4466
	.uleb128 .LEHE83-.LEHB83
	.uleb128 0
	.uleb128 0
.LLSDACSE4466:
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED2Ev,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED5Ev,comdat
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED2Ev, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED2Ev
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED1Ev
	.set	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED1Ev,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EED2Ev
	.section	.text._ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_,"axG",@progbits,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_,comdat
	.align 2
	.weak	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_
	.type	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_, @function
_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_:
.LFB4468:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4468
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$88, %rsp
	.cfi_offset 13, -24
	.cfi_offset 12, -32
	.cfi_offset 3, -40
	movq	%rdi, -104(%rbp)
	movq	%rsi, -112(%rbp)
	movq	-112(%rbp), %rdx
	movq	-104(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB84:
	call	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE11lower_boundERS5_
.LEHE84:
	movq	%rax, -80(%rbp)
	movl	$0, %ebx
	movl	$0, %r12d
	movq	-104(%rbp), %rax
	movq	%rax, %rdi
.LEHB85:
	call	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE3endEv
	movq	%rax, -64(%rbp)
	movl	$1, %ebx
	leaq	-64(%rbp), %rdx
	leaq	-80(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_
	testb	%al, %al
	jne	.L433
	leaq	-80(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv
	movq	%rax, %r13
	movq	-104(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE8key_compEv
.LEHE85:
	movl	$1, %r12d
	movq	-112(%rbp), %rcx
	leaq	-81(%rbp), %rax
	movq	%r13, %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	testb	%al, %al
	je	.L434
.L433:
	movl	$1, %eax
	jmp	.L435
.L434:
	movl	$0, %eax
.L435:
	testb	%r12b, %r12b
	testb	%bl, %bl
	testb	%al, %al
	je	.L438
	movq	$0, -56(%rbp)
	leaq	-56(%rbp), %rdx
	movq	-112(%rbp), %rcx
	leaq	-48(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIKP9chem_atomlEC1ERS2_RKl
	leaq	-48(%rbp), %rdx
	movq	-80(%rbp), %rcx
	movq	-104(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
.LEHB86:
	call	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE6insertESt17_Rb_tree_iteratorIS6_ERKS6_
	movq	%rax, -80(%rbp)
.L438:
	leaq	-80(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv
	addq	$8, %rax
	jmp	.L444
.L443:
	testb	%r12b, %r12b
	testb	%bl, %bl
	nop
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE86:
.L444:
	addq	$88, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4468:
	.section	.gcc_except_table
.LLSDA4468:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4468-.LLSDACSB4468
.LLSDACSB4468:
	.uleb128 .LEHB84-.LFB4468
	.uleb128 .LEHE84-.LEHB84
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB85-.LFB4468
	.uleb128 .LEHE85-.LEHB85
	.uleb128 .L443-.LFB4468
	.uleb128 0
	.uleb128 .LEHB86-.LFB4468
	.uleb128 .LEHE86-.LEHB86
	.uleb128 0
	.uleb128 0
.LLSDACSE4468:
	.section	.text._ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_,"axG",@progbits,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_,comdat
	.size	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_, .-_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEEixERS5_
	.section	.text._ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv,"axG",@progbits,_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv,comdat
	.align 2
	.weak	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	.type	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv, @function
_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv:
.LFB4470:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4470:
	.size	_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv, .-_ZNSt3mapIK9count_ptrI9chem_atomEi26element_pointer_value_lessSaISt4pairIS3_iEEE5beginEv
	.section	.text._ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2ERKSt17_Rb_tree_iteratorIS5_E,"axG",@progbits,_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC5ERKSt17_Rb_tree_iteratorIS5_E,comdat
	.align 2
	.weak	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2ERKSt17_Rb_tree_iteratorIS5_E
	.type	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2ERKSt17_Rb_tree_iteratorIS5_E, @function
_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2ERKSt17_Rb_tree_iteratorIS5_E:
.LFB4472:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, (%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4472:
	.size	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2ERKSt17_Rb_tree_iteratorIS5_E, .-_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2ERKSt17_Rb_tree_iteratorIS5_E
	.weak	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1ERKSt17_Rb_tree_iteratorIS5_E
	.set	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1ERKSt17_Rb_tree_iteratorIS5_E,_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2ERKSt17_Rb_tree_iteratorIS5_E
	.section	.text._Z7get_ptrIdEPT_S1_,"axG",@progbits,_Z7get_ptrIdEPT_S1_,comdat
	.weak	_Z7get_ptrIdEPT_S1_
	.type	_Z7get_ptrIdEPT_S1_, @function
_Z7get_ptrIdEPT_S1_:
.LFB4474:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4474:
	.size	_Z7get_ptrIdEPT_S1_, .-_Z7get_ptrIdEPT_S1_
	.section	.text._ZNKSt8valarrayIdEixEm,"axG",@progbits,_ZNKSt8valarrayIdEixEm,comdat
	.align 2
	.weak	_ZNKSt8valarrayIdEixEm
	.type	_ZNKSt8valarrayIdEixEm, @function
_ZNKSt8valarrayIdEixEm:
.LFB4475:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	-16(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4475:
	.size	_ZNKSt8valarrayIdEixEm, .-_ZNKSt8valarrayIdEixEm
	.section	.text._ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv,"axG",@progbits,_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv,comdat
	.align 2
	.weak	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv
	.type	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv, @function
_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv:
.LFB4476:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	_ZSt18_Rb_tree_incrementPSt18_Rb_tree_node_base
	movq	-8(%rbp), %rdx
	movq	%rax, (%rdx)
	movq	-8(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4476:
	.size	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv, .-_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEppEv
	.section	.text._ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv,"axG",@progbits,_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv,comdat
	.align 2
	.weak	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	.type	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv, @function
_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv:
.LFB4477:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	addq	$32, %rax
	movq	%rax, %rdi
	call	_ZSt11__addressofISt4pairIK9count_ptrI9chem_atomEiEEPT_RS6_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4477:
	.size	_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv, .-_ZNKSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEptEv
	.section	.text._ZN9count_ptrI9chem_atomEC2ERKS1_,"axG",@progbits,_ZN9count_ptrI9chem_atomEC5ERKS1_,comdat
	.align 2
	.weak	_ZN9count_ptrI9chem_atomEC2ERKS1_
	.type	_ZN9count_ptrI9chem_atomEC2ERKS1_, @function
_ZN9count_ptrI9chem_atomEC2ERKS1_:
.LFB4479:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, (%rax)
	movq	-16(%rbp), %rax
	movq	8(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, 8(%rax)
	movq	-8(%rbp), %rax
	movq	8(%rax), %rax
	movq	(%rax), %rdx
	addq	$1, %rdx
	movq	%rdx, (%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4479:
	.size	_ZN9count_ptrI9chem_atomEC2ERKS1_, .-_ZN9count_ptrI9chem_atomEC2ERKS1_
	.weak	_ZN9count_ptrI9chem_atomEC1ERKS1_
	.set	_ZN9count_ptrI9chem_atomEC1ERKS1_,_ZN9count_ptrI9chem_atomEC2ERKS1_
	.section	.text._ZN9count_ptrI9chem_atomED2Ev,"axG",@progbits,_ZN9count_ptrI9chem_atomED5Ev,comdat
	.align 2
	.weak	_ZN9count_ptrI9chem_atomED2Ev
	.type	_ZN9count_ptrI9chem_atomED2Ev, @function
_ZN9count_ptrI9chem_atomED2Ev:
.LFB4482:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4482
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
.LEHB87:
	call	_ZN9count_ptrI9chem_atomE6cancelEv
.LEHE87:
	jmp	.L462
.L461:
	cmpq	$-1, %rdx
	je	.L460
	movq	%rax, %rdi
.LEHB88:
	call	_Unwind_Resume
.L460:
	movq	%rax, %rdi
	call	__cxa_call_unexpected
.LEHE88:
.L462:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4482:
	.section	.gcc_except_table
	.align 4
.LLSDA4482:
	.byte	0xff
	.byte	0x3
	.uleb128 .LLSDATT4482-.LLSDATTD4482
.LLSDATTD4482:
	.byte	0x1
	.uleb128 .LLSDACSE4482-.LLSDACSB4482
.LLSDACSB4482:
	.uleb128 .LEHB87-.LFB4482
	.uleb128 .LEHE87-.LEHB87
	.uleb128 .L461-.LFB4482
	.uleb128 0x1
	.uleb128 .LEHB88-.LFB4482
	.uleb128 .LEHE88-.LEHB88
	.uleb128 0
	.uleb128 0
.LLSDACSE4482:
	.byte	0x7f
	.byte	0
	.align 4
.LLSDATT4482:
	.byte	0
	.section	.text._ZN9count_ptrI9chem_atomED2Ev,"axG",@progbits,_ZN9count_ptrI9chem_atomED5Ev,comdat
	.size	_ZN9count_ptrI9chem_atomED2Ev, .-_ZN9count_ptrI9chem_atomED2Ev
	.weak	_ZN9count_ptrI9chem_atomED1Ev
	.set	_ZN9count_ptrI9chem_atomED1Ev,_ZN9count_ptrI9chem_atomED2Ev
	.section	.text._ZNSt12_Vector_baseIiSaIiEE12_Vector_implD2Ev,"axG",@progbits,_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD5Ev,comdat
	.align 2
	.weak	_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD2Ev
	.type	_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD2Ev, @function
_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD2Ev:
.LFB4545:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSaIiED2Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4545:
	.size	_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD2Ev, .-_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD2Ev
	.weak	_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD1Ev
	.set	_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD1Ev,_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD2Ev
	.section	.text._ZNSt12_Vector_baseIiSaIiEED2Ev,"axG",@progbits,_ZNSt12_Vector_baseIiSaIiEED5Ev,comdat
	.align 2
	.weak	_ZNSt12_Vector_baseIiSaIiEED2Ev
	.type	_ZNSt12_Vector_baseIiSaIiEED2Ev, @function
_ZNSt12_Vector_baseIiSaIiEED2Ev:
.LFB4550:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4550
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	16(%rax), %rax
	movq	%rax, %rdx
	movq	-24(%rbp), %rax
	movq	(%rax), %rax
	subq	%rax, %rdx
	movq	%rdx, %rax
	sarq	$2, %rax
	movq	%rax, %rdx
	movq	-24(%rbp), %rax
	movq	(%rax), %rcx
	movq	-24(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
.LEHB89:
	call	_ZNSt12_Vector_baseIiSaIiEE13_M_deallocateEPim
.LEHE89:
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD1Ev
	jmp	.L469
.L468:
	movq	%rax, %rbx
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt12_Vector_baseIiSaIiEE12_Vector_implD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB90:
	call	_Unwind_Resume
.LEHE90:
.L469:
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4550:
	.section	.gcc_except_table
.LLSDA4550:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4550-.LLSDACSB4550
.LLSDACSB4550:
	.uleb128 .LEHB89-.LFB4550
	.uleb128 .LEHE89-.LEHB89
	.uleb128 .L468-.LFB4550
	.uleb128 0
	.uleb128 .LEHB90-.LFB4550
	.uleb128 .LEHE90-.LEHB90
	.uleb128 0
	.uleb128 0
.LLSDACSE4550:
	.section	.text._ZNSt12_Vector_baseIiSaIiEED2Ev,"axG",@progbits,_ZNSt12_Vector_baseIiSaIiEED5Ev,comdat
	.size	_ZNSt12_Vector_baseIiSaIiEED2Ev, .-_ZNSt12_Vector_baseIiSaIiEED2Ev
	.weak	_ZNSt12_Vector_baseIiSaIiEED1Ev
	.set	_ZNSt12_Vector_baseIiSaIiEED1Ev,_ZNSt12_Vector_baseIiSaIiEED2Ev
	.section	.text._ZNSt12_Vector_baseIiSaIiEE19_M_get_Tp_allocatorEv,"axG",@progbits,_ZNSt12_Vector_baseIiSaIiEE19_M_get_Tp_allocatorEv,comdat
	.align 2
	.weak	_ZNSt12_Vector_baseIiSaIiEE19_M_get_Tp_allocatorEv
	.type	_ZNSt12_Vector_baseIiSaIiEE19_M_get_Tp_allocatorEv, @function
_ZNSt12_Vector_baseIiSaIiEE19_M_get_Tp_allocatorEv:
.LFB4552:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4552:
	.size	_ZNSt12_Vector_baseIiSaIiEE19_M_get_Tp_allocatorEv, .-_ZNSt12_Vector_baseIiSaIiEE19_M_get_Tp_allocatorEv
	.section	.text._ZSt8_DestroyIPiiEvT_S1_RSaIT0_E,"axG",@progbits,_ZSt8_DestroyIPiiEvT_S1_RSaIT0_E,comdat
	.weak	_ZSt8_DestroyIPiiEvT_S1_RSaIT0_E
	.type	_ZSt8_DestroyIPiiEvT_S1_RSaIT0_E, @function
_ZSt8_DestroyIPiiEvT_S1_RSaIT0_E:
.LFB4553:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-16(%rbp), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZSt8_DestroyIPiEvT_S1_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4553:
	.size	_ZSt8_DestroyIPiiEvT_S1_RSaIT0_E, .-_ZSt8_DestroyIPiiEvT_S1_RSaIT0_E
	.section	.text._ZNSt12_Vector_baseIiSaIiEE13_M_deallocateEPim,"axG",@progbits,_ZNSt12_Vector_baseIiSaIiEE13_M_deallocateEPim,comdat
	.align 2
	.weak	_ZNSt12_Vector_baseIiSaIiEE13_M_deallocateEPim
	.type	_ZNSt12_Vector_baseIiSaIiEE13_M_deallocateEPim, @function
_ZNSt12_Vector_baseIiSaIiEE13_M_deallocateEPim:
.LFB4657:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	cmpq	$0, -16(%rbp)
	je	.L473
	movq	-8(%rbp), %rax
	movq	-24(%rbp), %rdx
	movq	-16(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorIiE10deallocateEPim
.L473:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4657:
	.size	_ZNSt12_Vector_baseIiSaIiEE13_M_deallocateEPim, .-_ZNSt12_Vector_baseIiSaIiEE13_M_deallocateEPim
	.section	.text._ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv,"axG",@progbits,_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv,comdat
	.align 2
	.weak	_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv
	.type	_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv, @function
_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv:
.LFB4681:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	24(%rax), %rdx
	leaq	-16(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1EPKSt13_Rb_tree_nodeIS5_E
	movq	-16(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4681:
	.size	_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv, .-_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv
	.section	.text._ZSt11__addressofIKSt4pairIK9count_ptrI9chem_atomEiEEPT_RS7_,"axG",@progbits,_ZSt11__addressofIKSt4pairIK9count_ptrI9chem_atomEiEEPT_RS7_,comdat
	.weak	_ZSt11__addressofIKSt4pairIK9count_ptrI9chem_atomEiEEPT_RS7_
	.type	_ZSt11__addressofIKSt4pairIK9count_ptrI9chem_atomEiEEPT_RS7_, @function
_ZSt11__addressofIKSt4pairIK9count_ptrI9chem_atomEiEEPT_RS7_:
.LFB4683:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4683:
	.size	_ZSt11__addressofIKSt4pairIK9count_ptrI9chem_atomEiEEPT_RS7_, .-_ZSt11__addressofIKSt4pairIK9count_ptrI9chem_atomEiEEPT_RS7_
	.section	.text._ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE4sizeEv,"axG",@progbits,_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE4sizeEv,comdat
	.align 2
	.weak	_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE4sizeEv
	.type	_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE4sizeEv, @function
_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE4sizeEv:
.LFB4684:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	40(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4684:
	.size	_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE4sizeEv, .-_ZNKSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE4sizeEv
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EEC2Ev,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EEC5Ev,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EEC2Ev
	.type	_ZN10multi_geomILi2EL10mem_layout0EEC2Ev, @function
_ZN10multi_geomILi2EL10mem_layout0EEC2Ev:
.LFB4705:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN8tree_vecC1Ev
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4705:
	.size	_ZN10multi_geomILi2EL10mem_layout0EEC2Ev, .-_ZN10multi_geomILi2EL10mem_layout0EEC2Ev
	.weak	_ZN10multi_geomILi2EL10mem_layout0EEC1Ev
	.set	_ZN10multi_geomILi2EL10mem_layout0EEC1Ev,_ZN10multi_geomILi2EL10mem_layout0EEC2Ev
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EED2Ev,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EED5Ev,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EED2Ev
	.type	_ZN10multi_geomILi2EL10mem_layout0EED2Ev, @function
_ZN10multi_geomILi2EL10mem_layout0EED2Ev:
.LFB4708:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4708
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
.LEHB91:
	call	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev
.LEHE91:
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
.LEHB92:
	call	_ZN8tree_vecD1Ev
.LEHE92:
	jmp	.L486
.L485:
	movq	%rax, %rbx
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN8tree_vecD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB93:
	call	_Unwind_Resume
.LEHE93:
.L486:
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4708:
	.section	.gcc_except_table
.LLSDA4708:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4708-.LLSDACSB4708
.LLSDACSB4708:
	.uleb128 .LEHB91-.LFB4708
	.uleb128 .LEHE91-.LEHB91
	.uleb128 .L485-.LFB4708
	.uleb128 0
	.uleb128 .LEHB92-.LFB4708
	.uleb128 .LEHE92-.LEHB92
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB93-.LFB4708
	.uleb128 .LEHE93-.LEHB93
	.uleb128 0
	.uleb128 0
.LLSDACSE4708:
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EED2Ev,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EED5Ev,comdat
	.size	_ZN10multi_geomILi2EL10mem_layout0EED2Ev, .-_ZN10multi_geomILi2EL10mem_layout0EED2Ev
	.weak	_ZN10multi_geomILi2EL10mem_layout0EED1Ev
	.set	_ZN10multi_geomILi2EL10mem_layout0EED1Ev,_ZN10multi_geomILi2EL10mem_layout0EED2Ev
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev:
.LFB4710:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -24(%rbp)
	movl	$0, -4(%rbp)
	jmp	.L488
.L489:
	movq	-24(%rbp), %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	addq	$8, %rdx
	movq	$0, 8(%rax,%rdx,8)
	addl	$1, -4(%rbp)
.L488:
	cmpl	$0, -4(%rbp)
	jle	.L489
	movq	-24(%rbp), %rax
	movq	$0, 96(%rax)
	movq	-24(%rbp), %rax
	movq	$0, 104(%rax)
	movq	-24(%rbp), %rax
	movq	$0, 112(%rax)
	movq	-24(%rbp), %rax
	movq	$0, 120(%rax)
	movq	-24(%rbp), %rax
	movq	$0, 128(%rax)
	movq	-24(%rbp), %rax
	movq	$0, 136(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4710:
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev:
.LFB4711:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE5clearEv
	movl	$0, -4(%rbp)
	jmp	.L491
.L493:
	movq	-24(%rbp), %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	addq	$8, %rdx
	movq	8(%rax,%rdx,8), %rax
	testq	%rax, %rax
	je	.L492
	movq	-24(%rbp), %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	addq	$8, %rdx
	movq	8(%rax,%rdx,8), %rax
	movq	%rax, %rdi
	call	_ZdaPv
.L492:
	addl	$1, -4(%rbp)
.L491:
	cmpl	$0, -4(%rbp)
	jle	.L493
	movq	-24(%rbp), %rax
	addq	$80, %rax
	xorpd	%xmm0, %xmm0
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdE6resizeEmd
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4711:
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev
	.section	.text._ZSt27__valarray_destroy_elementsIdEvPT_S1_,"axG",@progbits,_ZSt27__valarray_destroy_elementsIdEvPT_S1_,comdat
	.weak	_ZSt27__valarray_destroy_elementsIdEvPT_S1_
	.type	_ZSt27__valarray_destroy_elementsIdEvPT_S1_, @function
_ZSt27__valarray_destroy_elementsIdEvPT_S1_:
.LFB4712:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4712:
	.size	_ZSt27__valarray_destroy_elementsIdEvPT_S1_, .-_ZSt27__valarray_destroy_elementsIdEvPT_S1_
	.section	.text._ZSt22__valarray_get_storageIdErPT_m,"axG",@progbits,_ZSt22__valarray_get_storageIdErPT_m,comdat
	.weak	_ZSt22__valarray_get_storageIdErPT_m
	.type	_ZSt22__valarray_get_storageIdErPT_m, @function
_ZSt22__valarray_get_storageIdErPT_m:
.LFB4713:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	salq	$3, %rax
	movq	%rax, %rdi
	call	_ZSt21__valarray_get_memorym
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4713:
	.size	_ZSt22__valarray_get_storageIdErPT_m, .-_ZSt22__valarray_get_storageIdErPT_m
	.section	.text._ZSt25__valarray_fill_constructIdEvPT_S1_S0_,"axG",@progbits,_ZSt25__valarray_fill_constructIdEvPT_S1_S0_,comdat
	.weak	_ZSt25__valarray_fill_constructIdEvPT_S1_S0_
	.type	_ZSt25__valarray_fill_constructIdEvPT_S1_S0_, @function
_ZSt25__valarray_fill_constructIdEvPT_S1_S0_:
.LFB4714:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movsd	%xmm0, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	-16(%rbp), %rcx
	movq	-8(%rbp), %rdx
	movq	%rax, -32(%rbp)
	movsd	-32(%rbp), %xmm0
	movq	%rcx, %rsi
	movq	%rdx, %rdi
	call	_ZNSt16_Array_init_ctorIdLb1EE8_S_do_itEPdS1_d
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4714:
	.size	_ZSt25__valarray_fill_constructIdEvPT_S1_S0_, .-_ZSt25__valarray_fill_constructIdEvPT_S1_S0_
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv:
.LFB4715:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$80, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4715:
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE4valsEv
	.section	.text._ZNKSt8valarrayIdE4sizeEv,"axG",@progbits,_ZNKSt8valarrayIdE4sizeEv,comdat
	.align 2
	.weak	_ZNKSt8valarrayIdE4sizeEv
	.type	_ZNKSt8valarrayIdE4sizeEv, @function
_ZNKSt8valarrayIdE4sizeEv:
.LFB4716:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4716:
	.size	_ZNKSt8valarrayIdE4sizeEv, .-_ZNKSt8valarrayIdE4sizeEv
	.section	.rodata
	.align 8
.LC57:
	.string	"Failed: n <= d && index[n-1] > 0 && lgInbounds( n-1, index )"
.LC58:
	.string	"Failed: w.d == NULL"
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm
	.type	_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm, @function
_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm:
.LFB4717:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4717
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$88, %rsp
	.cfi_offset 13, -24
	.cfi_offset 12, -32
	.cfi_offset 3, -40
	movq	%rdi, -88(%rbp)
	movq	%rsi, -96(%rbp)
	movq	%rdx, -104(%rbp)
	cmpq	$2, -96(%rbp)
	ja	.L503
	movq	-96(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-104(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rax
	testq	%rax, %rax
	je	.L503
	movq	-96(%rbp), %rax
	leaq	-1(%rax), %rcx
	movq	-104(%rbp), %rdx
	movq	-88(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
.LEHB94:
	call	_ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm
	xorl	$1, %eax
	testb	%al, %al
	je	.L504
.L503:
	movl	$1, %eax
	jmp	.L505
.L504:
	movl	$0, %eax
.L505:
	movzbl	%al, %eax
	testq	%rax, %rax
	setne	%al
	testb	%al, %al
	je	.L506
	leaq	-64(%rbp), %rax
	movl	$.LC57, %ecx
	movl	$366, %edx
	movl	$.LC54, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE94:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L507
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
.LEHB95:
	call	_ZNK10bad_assert5printEv
	call	abort
.L507:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-64(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE95:
.L506:
	movq	-96(%rbp), %rax
	leaq	-1(%rax), %rcx
	movq	-88(%rbp), %rax
	movq	-104(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
.LEHB96:
	call	_ZN8tree_vec6getvecEmPKm
	movq	%rax, -72(%rbp)
	cmpq	$1, -96(%rbp)
	ja	.L508
	movq	-72(%rbp), %rax
	movq	8(%rax), %rax
	testq	%rax, %rax
	setne	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L509
	leaq	-64(%rbp), %rax
	movl	$.LC58, %ecx
	movl	$371, %edx
	movl	$.LC54, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE96:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L510
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
.LEHB97:
	call	_ZNK10bad_assert5printEv
	call	abort
.L510:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-64(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE97:
.L509:
	movq	-96(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-104(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rdx
	movabsq	$571957152676052992, %rax
	cmpq	%rax, %rdx
	ja	.L511
	movq	-96(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-104(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rax
	salq	$4, %rax
	addq	$8, %rax
	jmp	.L512
.L511:
	movq	$-1, %rax
.L512:
	movq	%rax, %rdi
.LEHB98:
	call	_Znam
	movq	%rax, %rbx
	movq	-96(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-104(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rax
	movq	%rax, (%rbx)
	leaq	8(%rbx), %rax
	movq	-96(%rbp), %rdx
	salq	$3, %rdx
	leaq	-8(%rdx), %rcx
	movq	-104(%rbp), %rdx
	addq	%rcx, %rdx
	movq	(%rdx), %rdx
	subq	$1, %rdx
	movq	%rdx, %r12
	movq	%rax, %r13
	jmp	.L513
.L514:
	movq	%r13, %rdi
	call	_ZN8tree_vecC1Ev
	addq	$16, %r13
	subq	$1, %r12
.L513:
	cmpq	$-1, %r12
	jne	.L514
	leaq	8(%rbx), %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, 8(%rax)
.L508:
	movq	-96(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-104(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, (%rax)
	movq	-96(%rbp), %rax
	leaq	-1(%rax), %rbx
	movq	-96(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-104(%rbp), %rax
	addq	%rax, %rdx
	movq	-96(%rbp), %rax
	subq	$1, %rax
	addq	$2, %rax
	leaq	0(,%rax,8), %rcx
	movq	-88(%rbp), %rax
	addq	%rcx, %rax
	addq	$8, %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZSt3maxImERKT_S2_S2_
	movq	(%rax), %rdx
	movq	-88(%rbp), %rax
	leaq	2(%rbx), %rcx
	movq	%rdx, 8(%rax,%rcx,8)
	movq	-96(%rbp), %rax
	leaq	-1(%rax), %rsi
	movq	-96(%rbp), %rax
	leaq	-1(%rax), %rdx
	movq	-88(%rbp), %rax
	addq	$6, %rdx
	movq	8(%rax,%rdx,8), %rdx
	movq	-96(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rcx
	movq	-104(%rbp), %rax
	addq	%rcx, %rax
	movq	(%rax), %rax
	leaq	(%rdx,%rax), %rcx
	movq	-88(%rbp), %rax
	leaq	6(%rsi), %rdx
	movq	%rcx, 8(%rax,%rdx,8)
	jmp	.L519
.L517:
	movq	%rax, %rbx
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.L518:
	movq	%rax, %rbx
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE98:
.L519:
	addq	$88, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4717:
	.section	.gcc_except_table
.LLSDA4717:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4717-.LLSDACSB4717
.LLSDACSB4717:
	.uleb128 .LEHB94-.LFB4717
	.uleb128 .LEHE94-.LEHB94
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB95-.LFB4717
	.uleb128 .LEHE95-.LEHB95
	.uleb128 .L517-.LFB4717
	.uleb128 0
	.uleb128 .LEHB96-.LFB4717
	.uleb128 .LEHE96-.LEHB96
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB97-.LFB4717
	.uleb128 .LEHE97-.LEHB97
	.uleb128 .L518-.LFB4717
	.uleb128 0
	.uleb128 .LEHB98-.LFB4717
	.uleb128 .LEHE98-.LEHB98
	.uleb128 0
	.uleb128 0
.LLSDACSE4717:
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm,comdat
	.size	_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm, .-_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm
	.section	.rodata
	.align 8
.LC59:
	.string	"Failed: n1[dim] == nsl[dim] && n2[dim] == nsl[dim+1]"
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv
	.type	_ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv, @function
_ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv:
.LFB4718:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4718
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$104, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -104(%rbp)
	movl	$0, -88(%rbp)
	jmp	.L521
.L522:
	movl	-88(%rbp), %eax
	cltq
	movq	$0, -64(%rbp,%rax,8)
	movl	-88(%rbp), %eax
	cltq
	movq	-64(%rbp,%rax,8), %rdx
	movl	-88(%rbp), %eax
	cltq
	movq	%rdx, -80(%rbp,%rax,8)
	addl	$1, -88(%rbp)
.L521:
	cmpl	$1, -88(%rbp)
	jle	.L522
	movq	-104(%rbp), %rcx
	leaq	-64(%rbp), %rdx
	leaq	-80(%rbp), %rsi
	movq	-104(%rbp), %rax
	movl	$0, %r8d
	movq	%rax, %rdi
.LEHB99:
	call	_ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci
	movl	$0, -84(%rbp)
	jmp	.L523
.L529:
	movl	-84(%rbp), %eax
	cltq
	movq	-80(%rbp,%rax,8), %rdx
	movq	-104(%rbp), %rax
	movl	-84(%rbp), %ecx
	movslq	%ecx, %rcx
	addq	$6, %rcx
	movq	8(%rax,%rcx,8), %rax
	cmpq	%rax, %rdx
	jne	.L524
	movl	-84(%rbp), %eax
	cltq
	movq	-64(%rbp,%rax,8), %rdx
	movl	-84(%rbp), %eax
	leal	1(%rax), %ecx
	movq	-104(%rbp), %rax
	movslq	%ecx, %rcx
	addq	$6, %rcx
	movq	8(%rax,%rcx,8), %rax
	cmpq	%rax, %rdx
	je	.L525
.L524:
	movl	$1, %eax
	jmp	.L526
.L525:
	movl	$0, %eax
.L526:
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L527
	leaq	-48(%rbp), %rax
	movl	$.LC59, %ecx
	movl	$416, %edx
	movl	$.LC54, %esi
	movq	%rax, %rdi
	call	_ZN10bad_assertC1EPKclS1_
.LEHE99:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L528
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB100:
	call	_ZNK10bad_assert5printEv
	call	abort
.L528:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE100:
.L527:
	addl	$1, -84(%rbp)
.L523:
	cmpl	$0, -84(%rbp)
	jle	.L529
	movq	-104(%rbp), %rax
	movq	64(%rax), %rdx
	movq	-104(%rbp), %rax
	movq	%rdx, 16(%rax)
	jmp	.L532
.L531:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
.LEHB101:
	call	_Unwind_Resume
.LEHE101:
.L532:
	addq	$104, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4718:
	.section	.gcc_except_table
.LLSDA4718:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4718-.LLSDACSB4718
.LLSDACSB4718:
	.uleb128 .LEHB99-.LFB4718
	.uleb128 .LEHE99-.LEHB99
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB100-.LFB4718
	.uleb128 .LEHE100-.LEHB100
	.uleb128 .L531-.LFB4718
	.uleb128 0
	.uleb128 .LEHB101-.LFB4718
	.uleb128 .LEHE101-.LEHB101
	.uleb128 0
	.uleb128 0
.LLSDACSE4718:
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv,comdat
	.size	_ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv, .-_ZN10multi_geomILi2EL10mem_layout0EE8finalizeEv
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci:
.LFB4719:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$72, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	%rdx, -56(%rbp)
	movq	%rcx, -64(%rbp)
	movl	%r8d, -68(%rbp)
	cmpl	$0, -68(%rbp)
	jns	.L534
	call	_Z13TotalInsanityv
.L534:
	movq	$0, -24(%rbp)
	jmp	.L535
.L538:
	cmpl	$0, -68(%rbp)
	jns	.L536
	movq	-40(%rbp), %rax
	movl	-68(%rbp), %edx
	movslq	%edx, %rdx
	addq	$8, %rdx
	movq	8(%rax,%rdx,8), %rsi
	movl	-68(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-48(%rbp), %rax
	addq	%rax, %rdx
	movq	(%rdx), %rax
	leaq	1(%rax), %rcx
	movq	%rcx, (%rdx)
	salq	$3, %rax
	leaq	(%rsi,%rax), %rdx
	movl	-68(%rbp), %eax
	leal	1(%rax), %ecx
	movq	-40(%rbp), %rax
	movslq	%ecx, %rcx
	addq	$8, %rcx
	movq	8(%rax,%rcx,8), %rax
	movl	-68(%rbp), %ecx
	movslq	%ecx, %rcx
	leaq	0(,%rcx,8), %rsi
	movq	-56(%rbp), %rcx
	addq	%rsi, %rcx
	movq	(%rcx), %rcx
	salq	$3, %rcx
	addq	%rcx, %rax
	movq	%rax, (%rdx)
	movl	-68(%rbp), %eax
	leal	1(%rax), %edi
	movq	-64(%rbp), %rax
	movq	8(%rax), %rax
	movq	-24(%rbp), %rdx
	salq	$4, %rdx
	leaq	(%rax,%rdx), %rcx
	movq	-56(%rbp), %rdx
	movq	-48(%rbp), %rsi
	movq	-40(%rbp), %rax
	movl	%edi, %r8d
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci
	jmp	.L537
.L536:
	movq	-40(%rbp), %rax
	movl	-68(%rbp), %edx
	movslq	%edx, %rdx
	addq	$8, %rdx
	movq	8(%rax,%rdx,8), %rsi
	movl	-68(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-48(%rbp), %rax
	addq	%rax, %rdx
	movq	(%rdx), %rax
	leaq	1(%rax), %rcx
	movq	%rcx, (%rdx)
	salq	$3, %rax
	leaq	(%rsi,%rax), %rbx
	movq	-40(%rbp), %rax
	addq	$80, %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZNSt8valarrayIdEixEm
	movl	-68(%rbp), %edx
	movslq	%edx, %rdx
	leaq	0(,%rdx,8), %rcx
	movq	-56(%rbp), %rdx
	addq	%rcx, %rdx
	movq	(%rdx), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	%rax, (%rbx)
.L537:
	movl	-68(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-56(%rbp), %rax
	addq	%rax, %rdx
	movl	-68(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rcx
	movq	-56(%rbp), %rax
	addq	%rcx, %rax
	movq	(%rax), %rcx
	movq	-64(%rbp), %rax
	movq	8(%rax), %rax
	movq	-24(%rbp), %rsi
	salq	$4, %rsi
	addq	%rsi, %rax
	movq	(%rax), %rax
	addq	%rcx, %rax
	movq	%rax, (%rdx)
	addq	$1, -24(%rbp)
.L535:
	movq	-64(%rbp), %rax
	movq	(%rax), %rax
	cmpq	-24(%rbp), %rax
	ja	.L538
	addq	$72, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4719:
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE12p_setupArrayEPmS2_PK8tree_veci
	.section	.text._ZSt28__valarray_default_constructIdEvPT_S1_,"axG",@progbits,_ZSt28__valarray_default_constructIdEvPT_S1_,comdat
	.weak	_ZSt28__valarray_default_constructIdEvPT_S1_
	.type	_ZSt28__valarray_default_constructIdEvPT_S1_, @function
_ZSt28__valarray_default_constructIdEvPT_S1_:
.LFB4873:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt19_Array_default_ctorIdLb1EE8_S_do_itEPdS1_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4873:
	.size	_ZSt28__valarray_default_constructIdEvPT_S1_, .-_ZSt28__valarray_default_constructIdEvPT_S1_
	.section	.rodata
.LC60:
	.string	"Failed: index[n] > 0"
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm:
.LFB4874:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA4874
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$72, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -72(%rbp)
	movq	%rsi, -80(%rbp)
	movl	$0, -52(%rbp)
	jmp	.L541
.L544:
	movl	-52(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-80(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rax
	testq	%rax, %rax
	sete	%al
	movzbl	%al, %eax
	testq	%rax, %rax
	je	.L542
	leaq	-48(%rbp), %rax
	movl	$.LC60, %ecx
	movl	$1205, %edx
	movl	$.LC54, %esi
	movq	%rax, %rdi
.LEHB102:
	call	_ZN10bad_assertC1EPKclS1_
.LEHE102:
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpu1iEv
	movq	%rax, %rdi
	call	_ZNK7t_cpu_i13lgAssertAbortEv
	testb	%al, %al
	je	.L543
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
.LEHB103:
	call	_ZNK10bad_assert5printEv
	call	abort
.L543:
	movl	$32, %edi
	call	__cxa_allocate_exception
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZN10bad_assertC1ERKS_
	movl	$_ZN10bad_assertD1Ev, %edx
	movl	$_ZTI10bad_assert, %esi
	movq	%rbx, %rdi
	call	__cxa_throw
.LEHE103:
.L542:
	addl	$1, -52(%rbp)
.L541:
	cmpl	$1, -52(%rbp)
	jle	.L544
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
.LEHB104:
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5clearEv
	movq	-72(%rbp), %rax
	movq	-80(%rbp), %rdx
	movl	$0, %esi
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEv
	jmp	.L547
.L546:
	movq	%rax, %rbx
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10bad_assertD1Ev
	movq	%rbx, %rax
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE104:
.L547:
	addq	$72, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4874:
	.section	.gcc_except_table
.LLSDA4874:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE4874-.LLSDACSB4874
.LLSDACSB4874:
	.uleb128 .LEHB102-.LFB4874
	.uleb128 .LEHE102-.LEHB102
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB103-.LFB4874
	.uleb128 .LEHE103-.LEHB103
	.uleb128 .L546-.LFB4874
	.uleb128 0
	.uleb128 .LEHB104-.LFB4874
	.uleb128 .LEHE104-.LEHB104
	.uleb128 0
	.uleb128 0
.LLSDACSE4874:
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm,comdat
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5allocEPm
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5n_ptrEv,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5n_ptrEv,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5n_ptrEv
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5n_ptrEv, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5n_ptrEv:
.LFB4875:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rdx
	movq	-16(%rbp), %rax
	leaq	48(%rax), %rdi
	movq	-16(%rbp), %rax
	movq	96(%rax), %rsi
	movq	-8(%rbp), %rax
	movq	%rdx, %rcx
	movq	%rdi, %rdx
	movq	%rax, %rdi
	call	_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC1EPdPKmPK8tree_vec
	movq	-8(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4875:
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5n_ptrEv, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5n_ptrEv
	.section	.text._ZNK9n_pointerIdLi2EL10mem_layout0ELb0EEixEm,"axG",@progbits,_ZNK9n_pointerIdLi2EL10mem_layout0ELb0EEixEm,comdat
	.align 2
	.weak	_ZNK9n_pointerIdLi2EL10mem_layout0ELb0EEixEm
	.type	_ZNK9n_pointerIdLi2EL10mem_layout0ELb0EEixEm, @function
_ZNK9n_pointerIdLi2EL10mem_layout0ELb0EEixEm:
.LFB4876:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-16(%rbp), %rax
	movq	(%rax), %rax
	movq	-24(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, %rax
	movq	(%rax), %rsi
	movq	-8(%rbp), %rax
	movl	$0, %ecx
	movl	$0, %edx
	movq	%rax, %rdi
	call	_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC1EPdPKmPK8tree_vec
	movq	-8(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4876:
	.size	_ZNK9n_pointerIdLi2EL10mem_layout0ELb0EEixEm, .-_ZNK9n_pointerIdLi2EL10mem_layout0ELb0EEixEm
	.section	.text._ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE3endEv,"axG",@progbits,_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE3endEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE3endEv
	.type	_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE3endEv, @function
_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE3endEv:
.LFB4878:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	leaq	8(%rax), %rdx
	leaq	-16(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1EPSt13_Rb_tree_nodeIS5_E
	movq	-16(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4878:
	.size	_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE3endEv, .-_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE3endEv
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC2Ev,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC5Ev,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC2Ev
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC2Ev, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC2Ev:
.LFB4888:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4888:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC2Ev, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC2Ev
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC1Ev
	.set	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC1Ev,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EEC2Ev
	.section	.text._ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev,"axG",@progbits,_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED5Ev,comdat
	.align 2
	.weak	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev
	.type	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev, @function
_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev:
.LFB4891:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4891:
	.size	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev, .-_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev
	.weak	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED1Ev
	.set	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED1Ev,_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E:
.LFB4893:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	jmp	.L558
.L559:
	movq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base
	movq	%rax, %rdx
	movq	-24(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E
	movq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base
	movq	%rax, -8(%rbp)
	movq	-32(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE15_M_destroy_nodeEPSt13_Rb_tree_nodeIS4_E
	movq	-8(%rbp), %rax
	movq	%rax, -32(%rbp)
.L558:
	cmpq	$0, -32(%rbp)
	jne	.L559
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4893:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_eraseEPSt13_Rb_tree_nodeIS4_E
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv:
.LFB4894:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	16(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4894:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv
	.section	.text._ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE11lower_boundERS5_,"axG",@progbits,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE11lower_boundERS5_,comdat
	.align 2
	.weak	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE11lower_boundERS5_
	.type	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE11lower_boundERS5_, @function
_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE11lower_boundERS5_:
.LFB4895:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rdx
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11lower_boundERS3_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4895:
	.size	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE11lower_boundERS5_, .-_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE11lower_boundERS5_
	.section	.text._ZNKSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE8key_compEv,"axG",@progbits,_ZNKSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE8key_compEv,comdat
	.align 2
	.weak	_ZNKSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE8key_compEv
	.type	_ZNKSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE8key_compEv, @function
_ZNKSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE8key_compEv:
.LFB4896:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8key_compEv
	movl	%ebx, %eax
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4896:
	.size	_ZNKSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE8key_compEv, .-_ZNKSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE8key_compEv
	.section	.text._ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv,"axG",@progbits,_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv,comdat
	.align 2
	.weak	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv
	.type	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv, @function
_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv:
.LFB4897:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	addq	$32, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4897:
	.size	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv, .-_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEdeEv
	.section	.text._ZNKSt4lessIP9chem_atomEclERKS1_S4_,"axG",@progbits,_ZNKSt4lessIP9chem_atomEclERKS1_S4_,comdat
	.align 2
	.weak	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	.type	_ZNKSt4lessIP9chem_atomEclERKS1_S4_, @function
_ZNKSt4lessIP9chem_atomEclERKS1_S4_:
.LFB4898:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-16(%rbp), %rax
	movq	(%rax), %rdx
	movq	-24(%rbp), %rax
	movq	(%rax), %rax
	cmpq	%rax, %rdx
	setb	%al
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4898:
	.size	_ZNKSt4lessIP9chem_atomEclERKS1_S4_, .-_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	.section	.text._ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE3endEv,"axG",@progbits,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE3endEv,comdat
	.align 2
	.weak	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE3endEv
	.type	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE3endEv, @function
_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE3endEv:
.LFB4899:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE3endEv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4899:
	.size	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE3endEv, .-_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE3endEv
	.section	.text._ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_,"axG",@progbits,_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_,comdat
	.align 2
	.weak	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_
	.type	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_, @function
_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_:
.LFB4900:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rdx
	movq	-16(%rbp), %rax
	movq	(%rax), %rax
	cmpq	%rax, %rdx
	sete	%al
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4900:
	.size	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_, .-_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_
	.section	.text._ZNSt4pairIKP9chem_atomlEC2ERS2_RKl,"axG",@progbits,_ZNSt4pairIKP9chem_atomlEC5ERS2_RKl,comdat
	.align 2
	.weak	_ZNSt4pairIKP9chem_atomlEC2ERS2_RKl
	.type	_ZNSt4pairIKP9chem_atomlEC2ERS2_RKl, @function
_ZNSt4pairIKP9chem_atomlEC2ERS2_RKl:
.LFB4902:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-16(%rbp), %rax
	movq	(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, (%rax)
	movq	-24(%rbp), %rax
	movq	(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, 8(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4902:
	.size	_ZNSt4pairIKP9chem_atomlEC2ERS2_RKl, .-_ZNSt4pairIKP9chem_atomlEC2ERS2_RKl
	.weak	_ZNSt4pairIKP9chem_atomlEC1ERS2_RKl
	.set	_ZNSt4pairIKP9chem_atomlEC1ERS2_RKl,_ZNSt4pairIKP9chem_atomlEC2ERS2_RKl
	.section	.text._ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE6insertESt17_Rb_tree_iteratorIS6_ERKS6_,"axG",@progbits,_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE6insertESt17_Rb_tree_iteratorIS6_ERKS6_,comdat
	.align 2
	.weak	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE6insertESt17_Rb_tree_iteratorIS6_ERKS6_
	.type	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE6insertESt17_Rb_tree_iteratorIS6_ERKS6_, @function
_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE6insertESt17_Rb_tree_iteratorIS6_ERKS6_:
.LFB4904:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	leaq	-32(%rbp), %rdx
	leaq	-16(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC1ERKSt17_Rb_tree_iteratorIS4_E
	movq	-24(%rbp), %rax
	movq	-40(%rbp), %rdx
	movq	-16(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE17_M_insert_unique_ESt23_Rb_tree_const_iteratorIS4_ERKS4_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4904:
	.size	_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE6insertESt17_Rb_tree_iteratorIS6_ERKS6_, .-_ZNSt3mapIP9chem_atomlSt4lessIS1_ESaISt4pairIKS1_lEEE6insertESt17_Rb_tree_iteratorIS6_ERKS6_
	.section	.text._ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv,"axG",@progbits,_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv
	.type	_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv, @function
_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv:
.LFB4905:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	24(%rax), %rdx
	leaq	-16(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1EPSt13_Rb_tree_nodeIS5_E
	movq	-16(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4905:
	.size	_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv, .-_ZNSt8_Rb_treeIK9count_ptrI9chem_atomESt4pairIS3_iESt10_Select1stIS5_E26element_pointer_value_lessSaIS5_EE5beginEv
	.section	.text._ZSt11__addressofISt4pairIK9count_ptrI9chem_atomEiEEPT_RS6_,"axG",@progbits,_ZSt11__addressofISt4pairIK9count_ptrI9chem_atomEiEEPT_RS6_,comdat
	.weak	_ZSt11__addressofISt4pairIK9count_ptrI9chem_atomEiEEPT_RS6_
	.type	_ZSt11__addressofISt4pairIK9count_ptrI9chem_atomEiEEPT_RS6_, @function
_ZSt11__addressofISt4pairIK9count_ptrI9chem_atomEiEEPT_RS6_:
.LFB4906:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4906:
	.size	_ZSt11__addressofISt4pairIK9count_ptrI9chem_atomEiEEPT_RS6_, .-_ZSt11__addressofISt4pairIK9count_ptrI9chem_atomEiEEPT_RS6_
	.section	.text._ZN9chem_atomD2Ev,"axG",@progbits,_ZN9chem_atomD5Ev,comdat
	.align 2
	.weak	_ZN9chem_atomD2Ev
	.type	_ZN9chem_atomD2Ev, @function
_ZN9chem_atomD2Ev:
.LFB4909:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$16, %rax
	movq	%rax, %rdi
	call	_ZNSt6vectorIiSaIiEED1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4909:
	.size	_ZN9chem_atomD2Ev, .-_ZN9chem_atomD2Ev
	.weak	_ZN9chem_atomD1Ev
	.set	_ZN9chem_atomD1Ev,_ZN9chem_atomD2Ev
	.section	.text._ZN9count_ptrI9chem_atomE6cancelEv,"axG",@progbits,_ZN9count_ptrI9chem_atomE6cancelEv,comdat
	.align 2
	.weak	_ZN9count_ptrI9chem_atomE6cancelEv
	.type	_ZN9count_ptrI9chem_atomE6cancelEv, @function
_ZN9count_ptrI9chem_atomE6cancelEv:
.LFB4907:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	8(%rax), %rax
	movq	(%rax), %rdx
	subq	$1, %rdx
	movq	%rdx, (%rax)
	movq	(%rax), %rax
	testq	%rax, %rax
	sete	%al
	testb	%al, %al
	je	.L583
	movq	-24(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rax, %rdi
	call	_ZdlPv
	movq	-24(%rbp), %rax
	movq	(%rax), %rbx
	testq	%rbx, %rbx
	je	.L583
	movq	%rbx, %rdi
	call	_ZN9chem_atomD1Ev
	movq	%rbx, %rdi
	call	_ZdlPv
.L583:
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4907:
	.size	_ZN9count_ptrI9chem_atomE6cancelEv, .-_ZN9count_ptrI9chem_atomE6cancelEv
	.section	.text._ZNSaIiED2Ev,"axG",@progbits,_ZNSaIiED5Ev,comdat
	.align 2
	.weak	_ZNSaIiED2Ev
	.type	_ZNSaIiED2Ev, @function
_ZNSaIiED2Ev:
.LFB4958:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorIiED2Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4958:
	.size	_ZNSaIiED2Ev, .-_ZNSaIiED2Ev
	.weak	_ZNSaIiED1Ev
	.set	_ZNSaIiED1Ev,_ZNSaIiED2Ev
	.section	.text._ZSt8_DestroyIPiEvT_S1_,"axG",@progbits,_ZSt8_DestroyIPiEvT_S1_,comdat
	.weak	_ZSt8_DestroyIPiEvT_S1_
	.type	_ZSt8_DestroyIPiEvT_S1_, @function
_ZSt8_DestroyIPiEvT_S1_:
.LFB4960:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt12_Destroy_auxILb1EE9__destroyIPiEEvT_S3_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4960:
	.size	_ZSt8_DestroyIPiEvT_S1_, .-_ZSt8_DestroyIPiEvT_S1_
	.section	.text._ZN9__gnu_cxx13new_allocatorIiE10deallocateEPim,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorIiE10deallocateEPim,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorIiE10deallocateEPim
	.type	_ZN9__gnu_cxx13new_allocatorIiE10deallocateEPim, @function
_ZN9__gnu_cxx13new_allocatorIiE10deallocateEPim:
.LFB5050:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, %rdi
	call	_ZdlPv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5050:
	.size	_ZN9__gnu_cxx13new_allocatorIiE10deallocateEPim, .-_ZN9__gnu_cxx13new_allocatorIiE10deallocateEPim
	.section	.text._ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPKSt13_Rb_tree_nodeIS5_E,"axG",@progbits,_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC5EPKSt13_Rb_tree_nodeIS5_E,comdat
	.align 2
	.weak	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPKSt13_Rb_tree_nodeIS5_E
	.type	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPKSt13_Rb_tree_nodeIS5_E, @function
_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPKSt13_Rb_tree_nodeIS5_E:
.LFB5072:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rdx
	movq	%rdx, (%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5072:
	.size	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPKSt13_Rb_tree_nodeIS5_E, .-_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPKSt13_Rb_tree_nodeIS5_E
	.weak	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1EPKSt13_Rb_tree_nodeIS5_E
	.set	_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1EPKSt13_Rb_tree_nodeIS5_E,_ZNSt23_Rb_tree_const_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPKSt13_Rb_tree_nodeIS5_E
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev
	.type	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev, @function
_ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev:
.LFB5082:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	$0, 16(%rax)
	movl	$0, -4(%rbp)
	jmp	.L591
.L592:
	movq	-24(%rbp), %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	addq	$2, %rdx
	movq	$0, 8(%rax,%rdx,8)
	movq	-24(%rbp), %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	addq	$4, %rdx
	movq	$0, 8(%rax,%rdx,8)
	movq	-24(%rbp), %rax
	movl	-4(%rbp), %edx
	movslq	%edx, %rdx
	addq	$6, %rdx
	movq	$0, 8(%rax,%rdx,8)
	addl	$1, -4(%rbp)
.L591:
	cmpl	$1, -4(%rbp)
	jle	.L592
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5082:
	.size	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev, .-_ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev
	.type	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev, @function
_ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev:
.LFB5083:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN8tree_vec5clearEv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5083:
	.size	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev, .-_ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE5clearEv,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE5clearEv,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EE5clearEv
	.type	_ZN10multi_geomILi2EL10mem_layout0EE5clearEv, @function
_ZN10multi_geomILi2EL10mem_layout0EE5clearEv:
.LFB5084:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear0Ev
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE8p_clear1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5084:
	.size	_ZN10multi_geomILi2EL10mem_layout0EE5clearEv, .-_ZN10multi_geomILi2EL10mem_layout0EE5clearEv
	.section	.text._ZNSt16_Array_init_ctorIdLb1EE8_S_do_itEPdS1_d,"axG",@progbits,_ZNSt16_Array_init_ctorIdLb1EE8_S_do_itEPdS1_d,comdat
	.weak	_ZNSt16_Array_init_ctorIdLb1EE8_S_do_itEPdS1_d
	.type	_ZNSt16_Array_init_ctorIdLb1EE8_S_do_itEPdS1_d, @function
_ZNSt16_Array_init_ctorIdLb1EE8_S_do_itEPdS1_d:
.LFB5085:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movsd	%xmm0, -24(%rbp)
	jmp	.L596
.L597:
	movq	-8(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	%rdx, -8(%rbp)
	movq	-24(%rbp), %rdx
	movq	%rdx, (%rax)
.L596:
	movq	-8(%rbp), %rax
	cmpq	-16(%rbp), %rax
	jne	.L597
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5085:
	.size	_ZNSt16_Array_init_ctorIdLb1EE8_S_do_itEPdS1_d, .-_ZNSt16_Array_init_ctorIdLb1EE8_S_do_itEPdS1_d
	.section	.text._ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm,"axG",@progbits,_ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm,comdat
	.align 2
	.weak	_ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm
	.type	_ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm, @function
_ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm:
.LFB5086:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$40, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	cmpq	$0, -32(%rbp)
	je	.L599
	movq	-32(%rbp), %rax
	leaq	-1(%rax), %rcx
	movq	-40(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm
	testb	%al, %al
	je	.L600
	movq	-32(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-40(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rbx
	movq	-32(%rbp), %rax
	leaq	-1(%rax), %rcx
	movq	-24(%rbp), %rax
	movq	-40(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNK8tree_vec6getvecEmPKm
	movq	(%rax), %rax
	cmpq	%rax, %rbx
	jnb	.L600
	movl	$1, %eax
	jmp	.L602
.L600:
	movl	$0, %eax
	jmp	.L602
.L599:
	movl	$1, %eax
.L602:
	addq	$40, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5086:
	.size	_ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm, .-_ZNK10multi_geomILi2EL10mem_layout0EE10lgInboundsEmPKm
	.section	.text._ZSt3maxImERKT_S2_S2_,"axG",@progbits,_ZSt3maxImERKT_S2_S2_,comdat
	.weak	_ZSt3maxImERKT_S2_S2_
	.type	_ZSt3maxImERKT_S2_S2_, @function
_ZSt3maxImERKT_S2_S2_:
.LFB5087:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rdx
	movq	-16(%rbp), %rax
	movq	(%rax), %rax
	cmpq	%rax, %rdx
	jnb	.L604
	movq	-16(%rbp), %rax
	jmp	.L605
.L604:
	movq	-8(%rbp), %rax
.L605:
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5087:
	.size	_ZSt3maxImERKT_S2_S2_, .-_ZSt3maxImERKT_S2_S2_
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci
	.type	_ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci, @function
_ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci:
.LFB5088:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$64, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	movq	%rcx, -48(%rbp)
	movl	%r8d, -52(%rbp)
	movq	$0, -8(%rbp)
	jmp	.L607
.L609:
	movl	-52(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-32(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rdx
	addq	$1, %rdx
	movq	%rdx, (%rax)
	cmpl	$0, -52(%rbp)
	jns	.L608
	movl	-52(%rbp), %eax
	leal	1(%rax), %edi
	movq	-48(%rbp), %rax
	movq	8(%rax), %rax
	movq	-8(%rbp), %rdx
	salq	$4, %rdx
	leaq	(%rax,%rdx), %rcx
	movq	-40(%rbp), %rdx
	movq	-32(%rbp), %rsi
	movq	-24(%rbp), %rax
	movl	%edi, %r8d
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci
.L608:
	movl	-52(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-40(%rbp), %rax
	addq	%rax, %rdx
	movl	-52(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rcx
	movq	-40(%rbp), %rax
	addq	%rcx, %rax
	movq	(%rax), %rcx
	movq	-48(%rbp), %rax
	movq	8(%rax), %rax
	movq	-8(%rbp), %rsi
	salq	$4, %rsi
	addq	%rsi, %rax
	movq	(%rax), %rax
	addq	%rcx, %rax
	movq	%rax, (%rdx)
	addq	$1, -8(%rbp)
.L607:
	movq	-48(%rbp), %rax
	movq	(%rax), %rax
	cmpq	-8(%rbp), %rax
	ja	.L609
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5088:
	.size	_ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci, .-_ZN10multi_geomILi2EL10mem_layout0EE12p_setupArrayEPmS2_PK8tree_veci
	.section	.text._ZNSt19_Array_default_ctorIdLb1EE8_S_do_itEPdS1_,"axG",@progbits,_ZNSt19_Array_default_ctorIdLb1EE8_S_do_itEPdS1_,comdat
	.weak	_ZNSt19_Array_default_ctorIdLb1EE8_S_do_itEPdS1_
	.type	_ZNSt19_Array_default_ctorIdLb1EE8_S_do_itEPdS1_, @function
_ZNSt19_Array_default_ctorIdLb1EE8_S_do_itEPdS1_:
.LFB5178:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rdx
	movq	-8(%rbp), %rax
	subq	%rax, %rdx
	movq	%rdx, %rax
	sarq	$3, %rax
	leaq	0(,%rax,8), %rdx
	movq	-8(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	memset
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5178:
	.size	_ZNSt19_Array_default_ctorIdLb1EE8_S_do_itEPdS1_, .-_ZNSt19_Array_default_ctorIdLb1EE8_S_do_itEPdS1_
	.section	.text._ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5clearEv,"axG",@progbits,_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5clearEv,comdat
	.align 2
	.weak	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5clearEv
	.type	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5clearEv, @function
_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5clearEv:
.LFB5179:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear0Ev
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE8p_clear1Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5179:
	.size	_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5clearEv, .-_ZN9multi_arrIdLi2EL10mem_layout0ELb0EE5clearEv
	.section	.text._ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm,"axG",@progbits,_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm,comdat
	.align 2
	.weak	_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm
	.type	_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm, @function
_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm:
.LFB5180:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	cmpq	$0, -32(%rbp)
	jne	.L613
	movq	-32(%rbp), %rax
	leaq	1(%rax), %rcx
	movq	-40(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm
	movq	-32(%rbp), %rax
	addq	$1, %rax
	cmpq	$1, %rax
	ja	.L612
	movq	-32(%rbp), %rax
	leaq	1(%rax), %rcx
	movq	-40(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm
	jmp	.L612
.L613:
	movq	-32(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-40(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rax
	movq	%rax, -8(%rbp)
	movq	$0, -16(%rbp)
	jmp	.L615
.L617:
	movq	-32(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-40(%rbp), %rax
	addq	%rax, %rdx
	movq	-16(%rbp), %rax
	movq	%rax, (%rdx)
	movq	-32(%rbp), %rax
	leaq	1(%rax), %rcx
	movq	-40(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE7reserveEmPKm
	movq	-32(%rbp), %rax
	addq	$1, %rax
	cmpq	$1, %rax
	ja	.L616
	movq	-32(%rbp), %rax
	leaq	1(%rax), %rcx
	movq	-40(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm
.L616:
	addq	$1, -16(%rbp)
.L615:
	movq	-16(%rbp), %rax
	cmpq	-8(%rbp), %rax
	jb	.L617
	movq	-32(%rbp), %rax
	salq	$3, %rax
	leaq	-8(%rax), %rdx
	movq	-40(%rbp), %rax
	addq	%rax, %rdx
	movq	-8(%rbp), %rax
	movq	%rax, (%rdx)
.L612:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5180:
	.size	_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm, .-_ZN10multi_geomILi2EL10mem_layout0EE17reserve_recursiveEmPm
	.section	.text._ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec,"axG",@progbits,_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC5EPdPKmPK8tree_vec,comdat
	.align 2
	.weak	_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec
	.type	_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec, @function
_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec:
.LFB5182:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	%rcx, -32(%rbp)
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-8(%rbp), %rax
	movq	-24(%rbp), %rdx
	movq	%rdx, 8(%rax)
	movq	-8(%rbp), %rax
	movq	-32(%rbp), %rdx
	movq	%rdx, 16(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5182:
	.size	_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec, .-_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec
	.weak	_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC1EPdPKmPK8tree_vec
	.set	_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC1EPdPKmPK8tree_vec,_ZN9n_pointerIdLi2EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec
	.section	.text._ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec,"axG",@progbits,_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC5EPdPKmPK8tree_vec,comdat
	.align 2
	.weak	_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec
	.type	_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec, @function
_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec:
.LFB5185:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	%rcx, -32(%rbp)
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rdx
	movq	%rdx, (%rax)
	movq	-8(%rbp), %rax
	movq	-24(%rbp), %rdx
	movq	%rdx, 8(%rax)
	movq	-8(%rbp), %rax
	movq	-32(%rbp), %rdx
	movq	%rdx, 16(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5185:
	.size	_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec, .-_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec
	.weak	_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC1EPdPKmPK8tree_vec
	.set	_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC1EPdPKmPK8tree_vec,_ZN9n_pointerIdLi1EL10mem_layout0ELb0EEC2EPdPKmPK8tree_vec
	.section	.text._ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPSt13_Rb_tree_nodeIS5_E,"axG",@progbits,_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC5EPSt13_Rb_tree_nodeIS5_E,comdat
	.align 2
	.weak	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPSt13_Rb_tree_nodeIS5_E
	.type	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPSt13_Rb_tree_nodeIS5_E, @function
_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPSt13_Rb_tree_nodeIS5_E:
.LFB5192:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rdx
	movq	%rdx, (%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5192:
	.size	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPSt13_Rb_tree_nodeIS5_E, .-_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPSt13_Rb_tree_nodeIS5_E
	.weak	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1EPSt13_Rb_tree_nodeIS5_E
	.set	_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC1EPSt13_Rb_tree_nodeIS5_E,_ZNSt17_Rb_tree_iteratorISt4pairIK9count_ptrI9chem_atomEiEEC2EPSt13_Rb_tree_nodeIS5_E
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC2Ev,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC5Ev,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC2Ev
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC2Ev, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC2Ev:
.LFB5198:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev
	movq	-8(%rbp), %rax
	movl	$0, 8(%rax)
	movq	-8(%rbp), %rax
	movq	$0, 16(%rax)
	movq	-8(%rbp), %rax
	movq	$0, 24(%rax)
	movq	-8(%rbp), %rax
	movq	$0, 32(%rax)
	movq	-8(%rbp), %rax
	movq	$0, 40(%rax)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EE13_M_initializeEv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5198:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC2Ev, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC2Ev
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC1Ev
	.set	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC1Ev,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EEC2Ev
	.section	.text._ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED5Ev,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev
	.type	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev, @function
_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev:
.LFB5201:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5201:
	.size	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev, .-_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev
	.weak	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED1Ev
	.set	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED1Ev,_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEED2Ev
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base,comdat
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base:
.LFB5203:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	24(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5203:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base,comdat
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base:
.LFB5204:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	16(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5204:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE15_M_destroy_nodeEPSt13_Rb_tree_nodeIS4_E,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE15_M_destroy_nodeEPSt13_Rb_tree_nodeIS4_E,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE15_M_destroy_nodeEPSt13_Rb_tree_nodeIS4_E
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE15_M_destroy_nodeEPSt13_Rb_tree_nodeIS4_E, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE15_M_destroy_nodeEPSt13_Rb_tree_nodeIS4_E:
.LFB5205:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$40, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	-48(%rbp), %rax
	addq	$32, %rax
	movq	%rax, %rdi
	call	_ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_
	movq	%rax, %rbx
	leaq	-17(%rbp), %rax
	movq	-40(%rbp), %rdx
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv
	leaq	-17(%rbp), %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE7destroyEPS5_
	leaq	-17(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSaISt4pairIKP9chem_atomlEED1Ev
	movq	-48(%rbp), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_put_nodeEPSt13_Rb_tree_nodeIS4_E
	addq	$40, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5205:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE15_M_destroy_nodeEPSt13_Rb_tree_nodeIS4_E, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE15_M_destroy_nodeEPSt13_Rb_tree_nodeIS4_E
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11lower_boundERS3_,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11lower_boundERS3_,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11lower_boundERS3_
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11lower_boundERS3_, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11lower_boundERS3_:
.LFB5206:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$24, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv
	movq	%rax, %rbx
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv
	movq	%rax, %rsi
	movq	-32(%rbp), %rdx
	movq	-24(%rbp), %rax
	movq	%rdx, %rcx
	movq	%rbx, %rdx
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_lower_boundEPSt13_Rb_tree_nodeIS4_ESD_RS3_
	addq	$24, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5206:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11lower_boundERS3_, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11lower_boundERS3_
	.section	.text._ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8key_compEv,"axG",@progbits,_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8key_compEv,comdat
	.align 2
	.weak	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8key_compEv
	.type	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8key_compEv, @function
_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8key_compEv:
.LFB5207:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5207:
	.size	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8key_compEv, .-_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8key_compEv
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE3endEv,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE3endEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE3endEv
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE3endEv, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE3endEv:
.LFB5208:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	leaq	8(%rax), %rdx
	leaq	-16(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E
	movq	-16(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5208:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE3endEv, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE3endEv
	.section	.text._ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC2ERKSt17_Rb_tree_iteratorIS4_E,"axG",@progbits,_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC5ERKSt17_Rb_tree_iteratorIS4_E,comdat
	.align 2
	.weak	_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC2ERKSt17_Rb_tree_iteratorIS4_E
	.type	_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC2ERKSt17_Rb_tree_iteratorIS4_E, @function
_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC2ERKSt17_Rb_tree_iteratorIS4_E:
.LFB5210:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, (%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5210:
	.size	_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC2ERKSt17_Rb_tree_iteratorIS4_E, .-_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC2ERKSt17_Rb_tree_iteratorIS4_E
	.weak	_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC1ERKSt17_Rb_tree_iteratorIS4_E
	.set	_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC1ERKSt17_Rb_tree_iteratorIS4_E,_ZNSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEEC2ERKSt17_Rb_tree_iteratorIS4_E
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE17_M_insert_unique_ESt23_Rb_tree_const_iteratorIS4_ERKS4_,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE17_M_insert_unique_ESt23_Rb_tree_const_iteratorIS4_ERKS4_,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE17_M_insert_unique_ESt23_Rb_tree_const_iteratorIS4_ERKS4_
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE17_M_insert_unique_ESt23_Rb_tree_const_iteratorIS4_ERKS4_, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE17_M_insert_unique_ESt23_Rb_tree_const_iteratorIS4_ERKS4_:
.LFB5212:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$64, %rsp
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	%rdx, -56(%rbp)
	movq	-56(%rbp), %rdx
	leaq	-32(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_
	movq	%rax, %rdx
	movq	-48(%rbp), %rcx
	movq	-40(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE29_M_get_insert_hint_unique_posESt23_Rb_tree_const_iteratorIS4_ERS3_
	movq	%rax, -16(%rbp)
	movq	%rdx, -8(%rbp)
	movq	-8(%rbp), %rax
	testq	%rax, %rax
	je	.L637
	movq	-8(%rbp), %rdx
	movq	-16(%rbp), %rsi
	movq	-56(%rbp), %rcx
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_
	jmp	.L639
.L637:
	movq	-16(%rbp), %rdx
	leaq	-32(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E
	movq	-32(%rbp), %rax
.L639:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5212:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE17_M_insert_unique_ESt23_Rb_tree_const_iteratorIS4_ERKS4_, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE17_M_insert_unique_ESt23_Rb_tree_const_iteratorIS4_ERKS4_
	.section	.text._ZN9__gnu_cxx13new_allocatorIiED2Ev,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorIiED5Ev,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorIiED2Ev
	.type	_ZN9__gnu_cxx13new_allocatorIiED2Ev, @function
_ZN9__gnu_cxx13new_allocatorIiED2Ev:
.LFB5247:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5247:
	.size	_ZN9__gnu_cxx13new_allocatorIiED2Ev, .-_ZN9__gnu_cxx13new_allocatorIiED2Ev
	.weak	_ZN9__gnu_cxx13new_allocatorIiED1Ev
	.set	_ZN9__gnu_cxx13new_allocatorIiED1Ev,_ZN9__gnu_cxx13new_allocatorIiED2Ev
	.section	.text._ZNSt12_Destroy_auxILb1EE9__destroyIPiEEvT_S3_,"axG",@progbits,_ZNSt12_Destroy_auxILb1EE9__destroyIPiEEvT_S3_,comdat
	.weak	_ZNSt12_Destroy_auxILb1EE9__destroyIPiEEvT_S3_
	.type	_ZNSt12_Destroy_auxILb1EE9__destroyIPiEEvT_S3_, @function
_ZNSt12_Destroy_auxILb1EE9__destroyIPiEEvT_S3_:
.LFB5249:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5249:
	.size	_ZNSt12_Destroy_auxILb1EE9__destroyIPiEEvT_S3_, .-_ZNSt12_Destroy_auxILb1EE9__destroyIPiEEvT_S3_
	.section	.text._ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev,"axG",@progbits,_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC5Ev,comdat
	.align 2
	.weak	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev
	.type	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev, @function
_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev:
.LFB5487:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5487:
	.size	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev, .-_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev
	.weak	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC1Ev
	.set	_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC1Ev,_ZNSaISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EE13_M_initializeEv,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EE13_M_initializeEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EE13_M_initializeEv
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EE13_M_initializeEv, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EE13_M_initializeEv:
.LFB5489:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movl	$0, 8(%rax)
	movq	-8(%rbp), %rax
	movq	$0, 16(%rax)
	movq	-8(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, 24(%rax)
	movq	-8(%rbp), %rax
	leaq	8(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, 32(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5489:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EE13_M_initializeEv, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13_Rb_tree_implIS8_Lb0EE13_M_initializeEv
	.section	.text._ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv,"axG",@progbits,_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv,comdat
	.align 2
	.weak	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv
	.type	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv, @function
_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv:
.LFB5490:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE21_M_get_Node_allocatorEv
	movq	%rax, %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSaISt4pairIKP9chem_atomlEEC1ISt13_Rb_tree_nodeIS3_EEERKSaIT_E
	movq	-8(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5490:
	.size	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv, .-_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv
	.section	.text._ZNSaISt4pairIKP9chem_atomlEED2Ev,"axG",@progbits,_ZNSaISt4pairIKP9chem_atomlEED5Ev,comdat
	.align 2
	.weak	_ZNSaISt4pairIKP9chem_atomlEED2Ev
	.type	_ZNSaISt4pairIKP9chem_atomlEED2Ev, @function
_ZNSaISt4pairIKP9chem_atomlEED2Ev:
.LFB5492:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED2Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5492:
	.size	_ZNSaISt4pairIKP9chem_atomlEED2Ev, .-_ZNSaISt4pairIKP9chem_atomlEED2Ev
	.weak	_ZNSaISt4pairIKP9chem_atomlEED1Ev
	.set	_ZNSaISt4pairIKP9chem_atomlEED1Ev,_ZNSaISt4pairIKP9chem_atomlEED2Ev
	.section	.text._ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_,"axG",@progbits,_ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_,comdat
	.weak	_ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_
	.type	_ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_, @function
_ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_:
.LFB5494:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5494:
	.size	_ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_, .-_ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_
	.section	.text._ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE7destroyEPS5_,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE7destroyEPS5_,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE7destroyEPS5_
	.type	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE7destroyEPS5_, @function
_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE7destroyEPS5_:
.LFB5495:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5495:
	.size	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE7destroyEPS5_, .-_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE7destroyEPS5_
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_put_nodeEPSt13_Rb_tree_nodeIS4_E,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_put_nodeEPSt13_Rb_tree_nodeIS4_E,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_put_nodeEPSt13_Rb_tree_nodeIS4_E
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_put_nodeEPSt13_Rb_tree_nodeIS4_E, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_put_nodeEPSt13_Rb_tree_nodeIS4_E:
.LFB5496:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rcx
	movl	$1, %edx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE10deallocateEPS7_m
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5496:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_put_nodeEPSt13_Rb_tree_nodeIS4_E, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_put_nodeEPSt13_Rb_tree_nodeIS4_E
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv:
.LFB5497:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$8, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5497:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_lower_boundEPSt13_Rb_tree_nodeIS4_ESD_RS3_,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_lower_boundEPSt13_Rb_tree_nodeIS4_ESD_RS3_,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_lower_boundEPSt13_Rb_tree_nodeIS4_ESD_RS3_
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_lower_boundEPSt13_Rb_tree_nodeIS4_ESD_RS3_, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_lower_boundEPSt13_Rb_tree_nodeIS4_ESD_RS3_:
.LFB5498:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$48, %rsp
	movq	%rdi, -24(%rbp)
	movq	%rsi, -32(%rbp)
	movq	%rdx, -40(%rbp)
	movq	%rcx, -48(%rbp)
	jmp	.L656
.L658:
	movq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E
	movq	%rax, %rcx
	movq	-24(%rbp), %rax
	movq	-48(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	xorl	$1, %eax
	testb	%al, %al
	je	.L657
	movq	-32(%rbp), %rax
	movq	%rax, -40(%rbp)
	movq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base
	movq	%rax, -32(%rbp)
	jmp	.L656
.L657:
	movq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base
	movq	%rax, -32(%rbp)
.L656:
	cmpq	$0, -32(%rbp)
	jne	.L658
	movq	-40(%rbp), %rdx
	leaq	-16(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E
	movq	-16(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5498:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_lower_boundEPSt13_Rb_tree_nodeIS4_ESD_RS3_, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_lower_boundEPSt13_Rb_tree_nodeIS4_ESD_RS3_
	.section	.text._ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC2EPSt13_Rb_tree_nodeIS4_E,"axG",@progbits,_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC5EPSt13_Rb_tree_nodeIS4_E,comdat
	.align 2
	.weak	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC2EPSt13_Rb_tree_nodeIS4_E
	.type	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC2EPSt13_Rb_tree_nodeIS4_E, @function
_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC2EPSt13_Rb_tree_nodeIS4_E:
.LFB5500:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	-16(%rbp), %rdx
	movq	%rdx, (%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5500:
	.size	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC2EPSt13_Rb_tree_nodeIS4_E, .-_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC2EPSt13_Rb_tree_nodeIS4_E
	.weak	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E
	.set	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E,_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC2EPSt13_Rb_tree_nodeIS4_E
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE29_M_get_insert_hint_unique_posESt23_Rb_tree_const_iteratorIS4_ERS3_,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE29_M_get_insert_hint_unique_posESt23_Rb_tree_const_iteratorIS4_ERS3_,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE29_M_get_insert_hint_unique_posESt23_Rb_tree_const_iteratorIS4_ERS3_
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE29_M_get_insert_hint_unique_posESt23_Rb_tree_const_iteratorIS4_ERS3_, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE29_M_get_insert_hint_unique_posESt23_Rb_tree_const_iteratorIS4_ERS3_:
.LFB5502:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$88, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -72(%rbp)
	movq	%rsi, -80(%rbp)
	movq	%rdx, -88(%rbp)
	leaq	-80(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEE13_M_const_castEv
	movq	%rax, -64(%rbp)
	movq	-64(%rbp), %rbx
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv
	cmpq	%rax, %rbx
	sete	%al
	testb	%al, %al
	je	.L662
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE4sizeEv
	testq	%rax, %rax
	je	.L663
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
	movq	%rax, %rcx
	movq	-72(%rbp), %rax
	movq	-88(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	testb	%al, %al
	je	.L663
	movl	$1, %eax
	jmp	.L664
.L663:
	movl	$0, %eax
.L664:
	testb	%al, %al
	je	.L665
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv
	movq	%rax, %rdx
	movq	$0, -40(%rbp)
	leaq	-40(%rbp), %rcx
	leaq	-32(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rdx
	jmp	.L677
.L665:
	movq	-88(%rbp), %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_
	jmp	.L677
.L662:
	movq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
	movq	%rax, %rdx
	movq	-72(%rbp), %rax
	movq	-88(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	testb	%al, %al
	je	.L667
	movq	-64(%rbp), %rax
	movq	%rax, -48(%rbp)
	movq	-64(%rbp), %rbx
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv
	movq	(%rax), %rax
	cmpq	%rax, %rbx
	sete	%al
	testb	%al, %al
	je	.L668
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv
	movq	%rax, %rbx
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv
	movq	%rax, %rcx
	leaq	-32(%rbp), %rax
	movq	%rbx, %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rdx
	jmp	.L677
.L668:
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
	movq	%rax, %rcx
	movq	-72(%rbp), %rax
	movq	-88(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	testb	%al, %al
	je	.L670
	movq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base
	testq	%rax, %rax
	sete	%al
	testb	%al, %al
	je	.L671
	movq	$0, -40(%rbp)
	leaq	-48(%rbp), %rdx
	leaq	-40(%rbp), %rcx
	leaq	-32(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rdx
	jmp	.L677
.L671:
	leaq	-64(%rbp), %rdx
	leaq	-64(%rbp), %rcx
	leaq	-32(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rdx
	jmp	.L677
.L670:
	movq	-88(%rbp), %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_
	jmp	.L677
.L667:
	movq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
	movq	%rax, %rcx
	movq	-72(%rbp), %rax
	movq	-88(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	testb	%al, %al
	je	.L672
	movq	-64(%rbp), %rax
	movq	%rax, -48(%rbp)
	movq	-64(%rbp), %rbx
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv
	movq	(%rax), %rax
	cmpq	%rax, %rbx
	sete	%al
	testb	%al, %al
	je	.L673
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv
	movq	%rax, %rdx
	movq	$0, -40(%rbp)
	leaq	-40(%rbp), %rcx
	leaq	-32(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rdx
	jmp	.L677
.L673:
	leaq	-48(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEppEv
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
	movq	%rax, %rdx
	movq	-72(%rbp), %rax
	movq	-88(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	testb	%al, %al
	je	.L675
	movq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base
	testq	%rax, %rax
	sete	%al
	testb	%al, %al
	je	.L676
	movq	$0, -40(%rbp)
	leaq	-64(%rbp), %rdx
	leaq	-40(%rbp), %rcx
	leaq	-32(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rdx
	jmp	.L677
.L676:
	leaq	-48(%rbp), %rdx
	leaq	-48(%rbp), %rcx
	leaq	-32(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rdx
	jmp	.L677
.L675:
	movq	-88(%rbp), %rdx
	movq	-72(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_
	jmp	.L677
.L672:
	movq	$0, -40(%rbp)
	leaq	-40(%rbp), %rdx
	leaq	-64(%rbp), %rcx
	leaq	-32(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-32(%rbp), %rax
	movq	-24(%rbp), %rdx
.L677:
	addq	$88, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5502:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE29_M_get_insert_hint_unique_posESt23_Rb_tree_const_iteratorIS4_ERS3_, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE29_M_get_insert_hint_unique_posESt23_Rb_tree_const_iteratorIS4_ERS3_
	.section	.text._ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_,"axG",@progbits,_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_,comdat
	.align 2
	.weak	_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_
	.type	_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_, @function
_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_:
.LFB5503:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-16(%rbp), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5503:
	.size	_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_, .-_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_:
.LFB5504:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	.cfi_lsda 0x3,.LLSDA5504
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r12
	pushq	%rbx
	subq	$64, %rsp
	.cfi_offset 12, -24
	.cfi_offset 3, -32
	movq	%rdi, -56(%rbp)
	movq	%rsi, -64(%rbp)
	movq	%rdx, -72(%rbp)
	movq	%rcx, -80(%rbp)
	movl	$0, %ebx
	cmpq	$0, -64(%rbp)
	jne	.L681
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv
	cmpq	-72(%rbp), %rax
	je	.L681
	movq	-72(%rbp), %rax
	movq	%rax, %rdi
.LEHB105:
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
.LEHE105:
	movq	%rax, %r12
	movl	$1, %ebx
	movq	-80(%rbp), %rdx
	leaq	-34(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_
	movq	%rax, %rcx
	movq	-56(%rbp), %rax
	movq	%r12, %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	testb	%al, %al
	je	.L682
.L681:
	movl	$1, %eax
	jmp	.L683
.L682:
	movl	$0, %eax
.L683:
	movb	%al, -33(%rbp)
	testb	%bl, %bl
	nop
	movq	-80(%rbp), %rdx
	movq	-56(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
.LEHB106:
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_create_nodeERKS4_
	movq	%rax, -24(%rbp)
	movq	-56(%rbp), %rax
	leaq	8(%rax), %rcx
	movzbl	-33(%rbp), %eax
	movq	-72(%rbp), %rdx
	movq	-24(%rbp), %rsi
	movl	%eax, %edi
	call	_ZSt29_Rb_tree_insert_and_rebalancebPSt18_Rb_tree_node_baseS0_RS_
	movq	-56(%rbp), %rax
	movq	40(%rax), %rax
	leaq	1(%rax), %rdx
	movq	-56(%rbp), %rax
	movq	%rdx, 40(%rax)
	movq	-24(%rbp), %rdx
	leaq	-32(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E
	movq	-32(%rbp), %rax
	jmp	.L689
.L688:
	testb	%bl, %bl
	nop
	movq	%rax, %rdi
	call	_Unwind_Resume
.LEHE106:
.L689:
	addq	$64, %rsp
	popq	%rbx
	popq	%r12
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5504:
	.section	.gcc_except_table
.LLSDA5504:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE5504-.LLSDACSB5504
.LLSDACSB5504:
	.uleb128 .LEHB105-.LFB5504
	.uleb128 .LEHE105-.LEHB105
	.uleb128 .L688-.LFB5504
	.uleb128 0
	.uleb128 .LEHB106-.LFB5504
	.uleb128 .LEHE106-.LEHB106
	.uleb128 0
	.uleb128 0
.LLSDACSE5504:
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_,comdat
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE10_M_insert_EPSt18_Rb_tree_node_baseSC_RKS4_
	.section	.text._ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC2ERKS1_S4_,"axG",@progbits,_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC5ERKS1_S4_,comdat
	.align 2
	.weak	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC2ERKS1_S4_
	.type	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC2ERKS1_S4_, @function
_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC2ERKS1_S4_:
.LFB5666:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-16(%rbp), %rax
	movq	(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, (%rax)
	movq	-24(%rbp), %rax
	movq	(%rax), %rdx
	movq	-8(%rbp), %rax
	movq	%rdx, 8(%rax)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5666:
	.size	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC2ERKS1_S4_, .-_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC2ERKS1_S4_
	.weak	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	.set	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_,_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC2ERKS1_S4_
	.section	.text._ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC5Ev,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev
	.type	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev, @function
_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev:
.LFB5672:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5672:
	.size	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev, .-_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev
	.weak	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC1Ev
	.set	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC1Ev,_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEEC2Ev
	.section	.text._ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE21_M_get_Node_allocatorEv,"axG",@progbits,_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE21_M_get_Node_allocatorEv,comdat
	.align 2
	.weak	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE21_M_get_Node_allocatorEv
	.type	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE21_M_get_Node_allocatorEv, @function
_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE21_M_get_Node_allocatorEv:
.LFB5674:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5674:
	.size	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE21_M_get_Node_allocatorEv, .-_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE21_M_get_Node_allocatorEv
	.section	.text._ZNSaISt4pairIKP9chem_atomlEEC2ISt13_Rb_tree_nodeIS3_EEERKSaIT_E,"axG",@progbits,_ZNSaISt4pairIKP9chem_atomlEEC5ISt13_Rb_tree_nodeIS3_EEERKSaIT_E,comdat
	.align 2
	.weak	_ZNSaISt4pairIKP9chem_atomlEEC2ISt13_Rb_tree_nodeIS3_EEERKSaIT_E
	.type	_ZNSaISt4pairIKP9chem_atomlEEC2ISt13_Rb_tree_nodeIS3_EEERKSaIT_E, @function
_ZNSaISt4pairIKP9chem_atomlEEC2ISt13_Rb_tree_nodeIS3_EEERKSaIT_E:
.LFB5676:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC2Ev
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5676:
	.size	_ZNSaISt4pairIKP9chem_atomlEEC2ISt13_Rb_tree_nodeIS3_EEERKSaIT_E, .-_ZNSaISt4pairIKP9chem_atomlEEC2ISt13_Rb_tree_nodeIS3_EEERKSaIT_E
	.weak	_ZNSaISt4pairIKP9chem_atomlEEC1ISt13_Rb_tree_nodeIS3_EEERKSaIT_E
	.set	_ZNSaISt4pairIKP9chem_atomlEEC1ISt13_Rb_tree_nodeIS3_EEERKSaIT_E,_ZNSaISt4pairIKP9chem_atomlEEC2ISt13_Rb_tree_nodeIS3_EEERKSaIT_E
	.section	.text._ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED2Ev,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED5Ev,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED2Ev
	.type	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED2Ev, @function
_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED2Ev:
.LFB5682:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5682:
	.size	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED2Ev, .-_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED2Ev
	.weak	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED1Ev
	.set	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED1Ev,_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEED2Ev
	.section	.text._ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE10deallocateEPS7_m,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE10deallocateEPS7_m,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE10deallocateEPS7_m
	.type	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE10deallocateEPS7_m, @function
_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE10deallocateEPS7_m:
.LFB5684:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, %rdi
	call	_ZdlPv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5684:
	.size	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE10deallocateEPS7_m, .-_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE10deallocateEPS7_m
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E,comdat
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E:
.LFB5685:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt13_Rb_tree_nodeIS4_E
	movq	%rax, %rdx
	leaq	-1(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5685:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E
	.section	.text._ZNKSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEE13_M_const_castEv,"axG",@progbits,_ZNKSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEE13_M_const_castEv,comdat
	.align 2
	.weak	_ZNKSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEE13_M_const_castEv
	.type	_ZNKSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEE13_M_const_castEv, @function
_ZNKSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEE13_M_const_castEv:
.LFB5686:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	(%rax), %rdx
	leaq	-16(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E
	movq	-16(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5686:
	.size	_ZNKSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEE13_M_const_castEv, .-_ZNKSt23_Rb_tree_const_iteratorISt4pairIKP9chem_atomlEE13_M_const_castEv
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv:
.LFB5687:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$32, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5687:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE12_M_rightmostEv
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base,comdat
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base:
.LFB5688:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt18_Rb_tree_node_base
	movq	%rax, %rdx
	leaq	-1(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt10_Select1stISt4pairIKP9chem_atomlEEclERKS4_
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5688:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
	.section	.text._ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE4sizeEv,"axG",@progbits,_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE4sizeEv,comdat
	.align 2
	.weak	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE4sizeEv
	.type	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE4sizeEv, @function
_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE4sizeEv:
.LFB5689:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	40(%rax), %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5689:
	.size	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE4sizeEv, .-_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE4sizeEv
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_:
.LFB5690:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$96, %rsp
	movq	%rdi, -88(%rbp)
	movq	%rsi, -96(%rbp)
	movq	-88(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_M_beginEv
	movq	%rax, -32(%rbp)
	movq	-88(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_M_endEv
	movq	%rax, -24(%rbp)
	movb	$1, -65(%rbp)
	jmp	.L709
.L712:
	movq	-32(%rbp), %rax
	movq	%rax, -24(%rbp)
	movq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt13_Rb_tree_nodeIS4_E
	movq	%rax, %rdx
	movq	-88(%rbp), %rax
	movq	-96(%rbp), %rcx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	movb	%al, -65(%rbp)
	cmpb	$0, -65(%rbp)
	je	.L710
	movq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE7_S_leftEPSt18_Rb_tree_node_base
	jmp	.L711
.L710:
	movq	-32(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_rightEPSt18_Rb_tree_node_base
.L711:
	movq	%rax, -32(%rbp)
.L709:
	cmpq	$0, -32(%rbp)
	jne	.L712
	movq	-24(%rbp), %rdx
	leaq	-64(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E
	cmpb	$0, -65(%rbp)
	je	.L713
	movq	-88(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE5beginEv
	movq	%rax, -16(%rbp)
	leaq	-16(%rbp), %rdx
	leaq	-64(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEeqERKS5_
	testb	%al, %al
	je	.L714
	movq	-24(%rbp), %rax
	movq	%rax, -40(%rbp)
	movq	-32(%rbp), %rax
	movq	%rax, -48(%rbp)
	leaq	-40(%rbp), %rdx
	leaq	-48(%rbp), %rcx
	leaq	-16(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-16(%rbp), %rax
	movq	-8(%rbp), %rdx
	jmp	.L717
.L714:
	leaq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv
.L713:
	movq	-64(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE6_S_keyEPKSt18_Rb_tree_node_base
	movq	%rax, %rcx
	movq	-88(%rbp), %rax
	movq	-96(%rbp), %rdx
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt4lessIP9chem_atomEclERKS1_S4_
	testb	%al, %al
	je	.L716
	movq	-24(%rbp), %rax
	movq	%rax, -40(%rbp)
	movq	-32(%rbp), %rax
	movq	%rax, -48(%rbp)
	leaq	-40(%rbp), %rdx
	leaq	-48(%rbp), %rcx
	leaq	-16(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-16(%rbp), %rax
	movq	-8(%rbp), %rdx
	jmp	.L717
.L716:
	movq	$0, -40(%rbp)
	leaq	-40(%rbp), %rdx
	leaq	-64(%rbp), %rcx
	leaq	-16(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	_ZNSt4pairIPSt18_Rb_tree_node_baseS1_EC1ERKS1_S4_
	movq	-16(%rbp), %rax
	movq	-8(%rbp), %rdx
.L717:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5690:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE24_M_get_insert_unique_posERS3_
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv:
.LFB5691:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$24, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5691:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_leftmostEv
	.section	.text._ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv,"axG",@progbits,_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv,comdat
	.align 2
	.weak	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv
	.type	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv, @function
_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv:
.LFB5692:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	_ZSt18_Rb_tree_decrementPSt18_Rb_tree_node_base
	movq	-8(%rbp), %rdx
	movq	%rax, (%rdx)
	movq	-8(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5692:
	.size	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv, .-_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEmmEv
	.section	.text._ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEppEv,"axG",@progbits,_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEppEv,comdat
	.align 2
	.weak	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEppEv
	.type	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEppEv, @function
_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEppEv:
.LFB5693:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, %rdi
	call	_ZSt18_Rb_tree_incrementPSt18_Rb_tree_node_base
	movq	-8(%rbp), %rdx
	movq	%rax, (%rdx)
	movq	-8(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5693:
	.size	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEppEv, .-_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEppEv
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_create_nodeERKS4_,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_create_nodeERKS4_,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_create_nodeERKS4_
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_create_nodeERKS4_, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_create_nodeERKS4_:
.LFB5694:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$40, %rsp
	.cfi_offset 3, -24
	movq	%rdi, -40(%rbp)
	movq	%rsi, -48(%rbp)
	movq	-40(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_get_nodeEv
	movq	%rax, -24(%rbp)
	movq	-24(%rbp), %rax
	addq	$32, %rax
	movq	%rax, %rdi
	call	_ZSt11__addressofISt4pairIKP9chem_atomlEEPT_RS5_
	movq	%rax, %rbx
	leaq	-25(%rbp), %rax
	movq	-40(%rbp), %rdx
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNKSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE13get_allocatorEv
	movq	-48(%rbp), %rdx
	leaq	-25(%rbp), %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE9constructEPS5_RKS5_
	leaq	-25(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNSaISt4pairIKP9chem_atomlEED1Ev
	movq	-24(%rbp), %rax
	addq	$40, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5694:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_create_nodeERKS4_, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE14_M_create_nodeERKS4_
	.section	.text._ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC2Ev,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC5Ev,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC2Ev
	.type	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC2Ev, @function
_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC2Ev:
.LFB5763:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5763:
	.size	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC2Ev, .-_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC2Ev
	.weak	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC1Ev
	.set	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC1Ev,_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEEC2Ev
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt13_Rb_tree_nodeIS4_E,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt13_Rb_tree_nodeIS4_E,comdat
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt13_Rb_tree_nodeIS4_E
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt13_Rb_tree_nodeIS4_E, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt13_Rb_tree_nodeIS4_E:
.LFB5768:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$32, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5768:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt13_Rb_tree_nodeIS4_E, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt13_Rb_tree_nodeIS4_E
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt18_Rb_tree_node_base,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt18_Rb_tree_node_base,comdat
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt18_Rb_tree_node_base
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt18_Rb_tree_node_base, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt18_Rb_tree_node_base:
.LFB5769:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	addq	$32, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5769:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt18_Rb_tree_node_base, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE8_S_valueEPKSt18_Rb_tree_node_base
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE5beginEv,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE5beginEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE5beginEv
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE5beginEv, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE5beginEv:
.LFB5770:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movq	24(%rax), %rdx
	leaq	-16(%rbp), %rax
	movq	%rdx, %rsi
	movq	%rax, %rdi
	call	_ZNSt17_Rb_tree_iteratorISt4pairIKP9chem_atomlEEC1EPSt13_Rb_tree_nodeIS4_E
	movq	-16(%rbp), %rax
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5770:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE5beginEv, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE5beginEv
	.section	.text._ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_get_nodeEv,"axG",@progbits,_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_get_nodeEv,comdat
	.align 2
	.weak	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_get_nodeEv
	.type	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_get_nodeEv, @function
_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_get_nodeEv:
.LFB5771:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$16, %rsp
	movq	%rdi, -8(%rbp)
	movq	-8(%rbp), %rax
	movl	$0, %edx
	movl	$1, %esi
	movq	%rax, %rdi
	call	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8allocateEmPKv
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5771:
	.size	_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_get_nodeEv, .-_ZNSt8_Rb_treeIP9chem_atomSt4pairIKS1_lESt10_Select1stIS4_ESt4lessIS1_ESaIS4_EE11_M_get_nodeEv
	.section	.text._ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE9constructEPS5_RKS5_,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE9constructEPS5_RKS5_,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE9constructEPS5_RKS5_
	.type	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE9constructEPS5_RKS5_, @function
_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE9constructEPS5_RKS5_:
.LFB5772:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, %rsi
	movl	$16, %edi
	call	_ZnwmPv
	movq	%rax, %rcx
	testq	%rcx, %rcx
	je	.L735
	movq	-24(%rbp), %rax
	movq	8(%rax), %rdx
	movq	(%rax), %rax
	movq	%rax, (%rcx)
	movq	%rdx, 8(%rcx)
.L735:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5772:
	.size	_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE9constructEPS5_RKS5_, .-_ZN9__gnu_cxx13new_allocatorISt4pairIKP9chem_atomlEE9constructEPS5_RKS5_
	.section	.text._ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8allocateEmPKv,"axG",@progbits,_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8allocateEmPKv,comdat
	.align 2
	.weak	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8allocateEmPKv
	.type	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8allocateEmPKv, @function
_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8allocateEmPKv:
.LFB5809:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -8(%rbp)
	movq	%rsi, -16(%rbp)
	movq	%rdx, -24(%rbp)
	movq	-8(%rbp), %rax
	movq	%rax, %rdi
	call	_ZNK9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8max_sizeEv
	cmpq	-16(%rbp), %rax
	setb	%al
	testb	%al, %al
	je	.L739
	call	_ZSt17__throw_bad_allocv
.L739:
	movq	-16(%rbp), %rdx
	movq	%rdx, %rax
	addq	%rax, %rax
	addq	%rdx, %rax
	salq	$4, %rax
	movq	%rax, %rdi
	call	_Znwm
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5809:
	.size	_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8allocateEmPKv, .-_ZN9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8allocateEmPKv
	.section	.text._ZNK9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8max_sizeEv,"axG",@progbits,_ZNK9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8max_sizeEv,comdat
	.align 2
	.weak	_ZNK9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8max_sizeEv
	.type	_ZNK9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8max_sizeEv, @function
_ZNK9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8max_sizeEv:
.LFB5827:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movq	%rdi, -8(%rbp)
	movabsq	$384307168202282325, %rax
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5827:
	.size	_ZNK9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8max_sizeEv, .-_ZNK9__gnu_cxx13new_allocatorISt13_Rb_tree_nodeISt4pairIKP9chem_atomlEEE8max_sizeEv
	.weak	_ZTV10bad_assert
	.section	.rodata._ZTV10bad_assert,"aG",@progbits,_ZTV10bad_assert,comdat
	.align 32
	.type	_ZTV10bad_assert, @object
	.size	_ZTV10bad_assert, 32
_ZTV10bad_assert:
	.quad	0
	.quad	_ZTI10bad_assert
	.quad	_ZN10bad_assertD1Ev
	.quad	_ZN10bad_assertD0Ev
	.weak	_ZTS10bad_assert
	.section	.rodata._ZTS10bad_assert,"aG",@progbits,_ZTS10bad_assert,comdat
	.type	_ZTS10bad_assert, @object
	.size	_ZTS10bad_assert, 13
_ZTS10bad_assert:
	.string	"10bad_assert"
	.weak	_ZTI10bad_assert
	.section	.rodata._ZTI10bad_assert,"aG",@progbits,_ZTI10bad_assert,comdat
	.align 16
	.type	_ZTI10bad_assert, @object
	.size	_ZTI10bad_assert, 16
_ZTI10bad_assert:
	.quad	_ZTVN10__cxxabiv117__class_type_infoE+16
	.quad	_ZTS10bad_assert
	.text
	.type	_Z41__static_initialization_and_destruction_0ii, @function
_Z41__static_initialization_and_destruction_0ii:
.LFB5835:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movl	%edi, -4(%rbp)
	movl	%esi, -8(%rbp)
	cmpl	$1, -4(%rbp)
	jne	.L743
	cmpl	$65535, -8(%rbp)
	jne	.L743
	call	_ZNSt14numeric_limitsIfE3maxEv
	movss	.LC61(%rip), %xmm1
	divss	%xmm1, %xmm0
	movss	%xmm0, _ZL8BIGFLOAT(%rip)
	call	_ZNSt14numeric_limitsIfE3minEv
	movss	.LC61(%rip), %xmm1
	mulss	%xmm1, %xmm0
	movss	%xmm0, _ZL10SMALLFLOAT(%rip)
	movl	$32, %esi
	movl	$16, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL6mode_w(%rip)
	movl	$1, %esi
	movl	$16, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL6mode_a(%rip)
	movl	$16, %esi
	movl	$8, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL7mode_rp(%rip)
	movl	$16, %esi
	movl	$8, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	$32, %esi
	movl	%eax, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL7mode_wp(%rip)
	movl	$16, %esi
	movl	$8, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	$1, %esi
	movl	%eax, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL7mode_ap(%rip)
	movl	$4, %esi
	movl	$8, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL7mode_rb(%rip)
	movl	_ZL6mode_w(%rip), %eax
	movl	$4, %esi
	movl	%eax, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL7mode_wb(%rip)
	movl	_ZL6mode_a(%rip), %eax
	movl	$4, %esi
	movl	%eax, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL7mode_ab(%rip)
	movl	_ZL7mode_rp(%rip), %eax
	movl	$4, %esi
	movl	%eax, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL8mode_rpb(%rip)
	movl	_ZL7mode_wp(%rip), %eax
	movl	$4, %esi
	movl	%eax, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL8mode_wpb(%rip)
	movl	_ZL7mode_ap(%rip), %eax
	movl	$4, %esi
	movl	%eax, %edi
	call	_ZStorSt13_Ios_OpenmodeS_
	movl	%eax, _ZL8mode_apb(%rip)
	movl	$_ZL3cpu, %edi
	call	_ZN5t_cpuC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZL3cpu, %esi
	movl	$_ZN5t_cpuD1Ev, %edi
	call	__cxa_atexit
	movabsq	$4686327217615121787, %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow2IdET_S0_
	movsd	%xmm0, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, _ZL7SQAS1SR(%rip)
	movsd	_ZL7SQAS1SR(%rip), %xmm1
	movsd	.LC63(%rip), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL8SQAS_SKY(%rip)
	movabsq	$4763660085914238976, %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow2IdET_S0_
	movsd	.LC63(%rip), %xmm1
	mulsd	%xmm1, %xmm0
	movsd	.LC65(%rip), %xmm1
	divsd	%xmm0, %xmm1
	movapd	%xmm1, %xmm0
	movsd	%xmm0, _ZL14ELECTRIC_CONST(%rip)
	movabsq	$4215483348127298365, %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow2IdET_S0_
	movsd	.LC67(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL12HION_LTE_POP(%rip)
	movq	_ZL12HION_LTE_POP(%rip), %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow3IdET_S0_
	movsd	%xmm0, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	sqrt
	movsd	%xmm0, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, _ZL4SAHA(%rip)
	movabsq	$4682277504122888094, %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow3IdET_S0_
	movsd	.LC69(%rip), %xmm1
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL6HNU3C2(%rip)
	movabsq	$4369588615818088655, %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow2IdET_S0_
	movsd	.LC71(%rip), %xmm1
	mulsd	%xmm1, %xmm0
	call	_Z4pow2IdET_S0_
	movsd	%xmm0, -16(%rbp)
	movabsq	$4203234286544377125, %rax
	movq	%rax, -24(%rbp)
	movsd	-24(%rbp), %xmm0
	call	_Z4pow3IdET_S0_
	movsd	.LC73(%rip), %xmm1
	movapd	%xmm0, %xmm2
	mulsd	%xmm1, %xmm2
	movsd	%xmm2, -24(%rbp)
	movabsq	$4763660085914238976, %rax
	movq	%rax, -32(%rbp)
	movsd	-32(%rbp), %xmm0
	call	_Z4pow2IdET_S0_
	mulsd	-24(%rbp), %xmm0
	movsd	-16(%rbp), %xmm3
	divsd	%xmm0, %xmm3
	movapd	%xmm3, %xmm0
	movsd	%xmm0, _ZL12STEFAN_BOLTZ(%rip)
	movabsq	$4467712605079238603, %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow2IdET_S0_
	movsd	.LC64(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movsd	.LC72(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL14FINE_STRUCTURE(%rip)
	movq	_ZL14FINE_STRUCTURE(%rip), %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow2IdET_S0_
	movsd	%xmm0, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, _ZL15FINE_STRUCTURE2(%rip)
	movsd	_ZL14FINE_STRUCTURE(%rip), %xmm0
	movsd	.LC75(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL14BOHR_RADIUS_CM(%rip)
	movq	_ZL15FINE_STRUCTURE2(%rip), %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow3IdET_S0_
	movsd	.LC76(%rip), %xmm1
	mulsd	%xmm1, %xmm0
	movsd	.LC77(%rip), %xmm1
	mulsd	%xmm1, %xmm0
	movsd	.LC78(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL14TWO_PHOT_CONST(%rip)
	movsd	_ZL4SAHA(%rip), %xmm1
	movsd	.LC70(%rip), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	.LC66(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL10COLL_CONST(%rip)
	movq	_ZL15FINE_STRUCTURE2(%rip), %rax
	movq	%rax, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	call	_Z4pow3IdET_S0_
	movsd	%xmm0, -16(%rbp)
	movabsq	$4684664986766835159, %rax
	movq	%rax, -24(%rbp)
	movsd	-24(%rbp), %xmm0
	call	_Z4pow3IdET_S0_
	mulsd	-16(%rbp), %xmm0
	movsd	.LC71(%rip), %xmm1
	divsd	%xmm1, %xmm0
	call	sqrt
	movsd	.LC64(%rip), %xmm1
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL11MILNE_CONST(%rip)
	movsd	_ZL14FINE_STRUCTURE(%rip), %xmm1
	movsd	.LC80(%rip), %xmm0
	mulsd	%xmm1, %xmm0
	movsd	.LC81(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movsd	%xmm0, _ZL16TRANS_PROB_CONST(%rip)
.L743:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5835:
	.size	_Z41__static_initialization_and_destruction_0ii, .-_Z41__static_initialization_and_destruction_0ii
	.section	.rodata
	.align 8
	.type	_ZL9BIGDOUBLE, @object
	.size	_ZL9BIGDOUBLE, 8
_ZL9BIGDOUBLE:
	.long	1202590842
	.long	2139388641
	.align 8
	.type	_ZL11SMALLDOUBLE, @object
	.size	_ZL11SMALLDOUBLE, 8
_ZL11SMALLDOUBLE:
	.long	0
	.long	7929856
	.align 4
	.type	_ZL6STDLEN, @object
	.size	_ZL6STDLEN, 4
_ZL6STDLEN:
	.long	32
	.align 4
	.type	_ZL6mode_r, @object
	.size	_ZL6mode_r, 4
_ZL6mode_r:
	.long	8
	.align 4
	.type	_ZL20FILENAME_PATH_LENGTH, @object
	.size	_ZL20FILENAME_PATH_LENGTH, 4
_ZL20FILENAME_PATH_LENGTH:
	.long	200
	.align 4
	.type	_ZL22FILENAME_PATH_LENGTH_2, @object
	.size	_ZL22FILENAME_PATH_LENGTH_2, 4
_ZL22FILENAME_PATH_LENGTH_2:
	.long	400
	.align 4
	.type	_ZL17INPUT_LINE_LENGTH, @object
	.size	_ZL17INPUT_LINE_LENGTH, 4
_ZL17INPUT_LINE_LENGTH:
	.long	2000
	.align 4
	.type	_ZL6LIMELM, @object
	.size	_ZL6LIMELM, 4
_ZL6LIMELM:
	.long	30
	.align 4
	.type	_ZL4NISO, @object
	.size	_ZL4NISO, 4
_ZL4NISO:
	.long	2
	.align 4
	.type	_ZL16NHYDRO_MAX_LEVEL, @object
	.size	_ZL16NHYDRO_MAX_LEVEL, 4
_ZL16NHYDRO_MAX_LEVEL:
	.long	401
	.align 8
	.type	_ZL11MAX_DENSITY, @object
	.size	_ZL11MAX_DENSITY, 8
_ZL11MAX_DENSITY:
	.long	2044304820
	.long	1156216899
	.align 8
	.type	_ZL12DEPTH_OFFSET, @object
	.size	_ZL12DEPTH_OFFSET, 8
_ZL12DEPTH_OFFSET:
	.long	4276863648
	.long	968116299
	.align 4
	.type	_ZL8ipRecEsc, @object
	.size	_ZL8ipRecEsc, 4
_ZL8ipRecEsc:
	.long	2
	.align 4
	.type	_ZL11ipRecNetEsc, @object
	.size	_ZL11ipRecNetEsc, 4
_ZL11ipRecNetEsc:
	.long	1
	.align 4
	.type	_ZL8ipRecRad, @object
	.size	_ZL8ipRecRad, 4
_ZL8ipRecRad:
	.zero	4
	.align 4
	.type	_ZL5ipPRD, @object
	.size	_ZL5ipPRD, 4
_ZL5ipPRD:
	.long	1
	.align 4
	.type	_ZL5ipCRD, @object
	.size	_ZL5ipCRD, 4
_ZL5ipCRD:
	.long	-1
	.align 4
	.type	_ZL6ipCRDW, @object
	.size	_ZL6ipCRDW, 4
_ZL6ipCRDW:
	.long	2
	.align 4
	.type	_ZL6ipLY_A, @object
	.size	_ZL6ipLY_A, 4
_ZL6ipLY_A:
	.long	-2
	.align 4
	.type	_ZL9ipDEST_K2, @object
	.size	_ZL9ipDEST_K2, 4
_ZL9ipDEST_K2:
	.long	1
	.align 4
	.type	_ZL12ipDEST_INCOM, @object
	.size	_ZL12ipDEST_INCOM, 4
_ZL12ipDEST_INCOM:
	.long	2
	.align 4
	.type	_ZL12ipDEST_SIMPL, @object
	.size	_ZL12ipDEST_SIMPL, 4
_ZL12ipDEST_SIMPL:
	.long	3
	.align 4
	.type	_ZL10ipHYDROGEN, @object
	.size	_ZL10ipHYDROGEN, 4
_ZL10ipHYDROGEN:
	.zero	4
	.align 4
	.type	_ZL8ipHELIUM, @object
	.size	_ZL8ipHELIUM, 4
_ZL8ipHELIUM:
	.long	1
	.align 4
	.type	_ZL9ipLITHIUM, @object
	.size	_ZL9ipLITHIUM, 4
_ZL9ipLITHIUM:
	.long	2
	.align 4
	.type	_ZL11ipBERYLLIUM, @object
	.size	_ZL11ipBERYLLIUM, 4
_ZL11ipBERYLLIUM:
	.long	3
	.align 4
	.type	_ZL7ipBORON, @object
	.size	_ZL7ipBORON, 4
_ZL7ipBORON:
	.long	4
	.align 4
	.type	_ZL8ipCARBON, @object
	.size	_ZL8ipCARBON, 4
_ZL8ipCARBON:
	.long	5
	.align 4
	.type	_ZL10ipNITROGEN, @object
	.size	_ZL10ipNITROGEN, 4
_ZL10ipNITROGEN:
	.long	6
	.align 4
	.type	_ZL8ipOXYGEN, @object
	.size	_ZL8ipOXYGEN, 4
_ZL8ipOXYGEN:
	.long	7
	.align 4
	.type	_ZL10ipFLUORINE, @object
	.size	_ZL10ipFLUORINE, 4
_ZL10ipFLUORINE:
	.long	8
	.align 4
	.type	_ZL6ipNEON, @object
	.size	_ZL6ipNEON, 4
_ZL6ipNEON:
	.long	9
	.align 4
	.type	_ZL8ipSODIUM, @object
	.size	_ZL8ipSODIUM, 4
_ZL8ipSODIUM:
	.long	10
	.align 4
	.type	_ZL11ipMAGNESIUM, @object
	.size	_ZL11ipMAGNESIUM, 4
_ZL11ipMAGNESIUM:
	.long	11
	.align 4
	.type	_ZL11ipALUMINIUM, @object
	.size	_ZL11ipALUMINIUM, 4
_ZL11ipALUMINIUM:
	.long	12
	.align 4
	.type	_ZL9ipSILICON, @object
	.size	_ZL9ipSILICON, 4
_ZL9ipSILICON:
	.long	13
	.align 4
	.type	_ZL12ipPHOSPHORUS, @object
	.size	_ZL12ipPHOSPHORUS, 4
_ZL12ipPHOSPHORUS:
	.long	14
	.align 4
	.type	_ZL9ipSULPHUR, @object
	.size	_ZL9ipSULPHUR, 4
_ZL9ipSULPHUR:
	.long	15
	.align 4
	.type	_ZL10ipCHLORINE, @object
	.size	_ZL10ipCHLORINE, 4
_ZL10ipCHLORINE:
	.long	16
	.align 4
	.type	_ZL7ipARGON, @object
	.size	_ZL7ipARGON, 4
_ZL7ipARGON:
	.long	17
	.align 4
	.type	_ZL11ipPOTASSIUM, @object
	.size	_ZL11ipPOTASSIUM, 4
_ZL11ipPOTASSIUM:
	.long	18
	.align 4
	.type	_ZL9ipCALCIUM, @object
	.size	_ZL9ipCALCIUM, 4
_ZL9ipCALCIUM:
	.long	19
	.align 4
	.type	_ZL10ipSCANDIUM, @object
	.size	_ZL10ipSCANDIUM, 4
_ZL10ipSCANDIUM:
	.long	20
	.align 4
	.type	_ZL10ipTITANIUM, @object
	.size	_ZL10ipTITANIUM, 4
_ZL10ipTITANIUM:
	.long	21
	.align 4
	.type	_ZL10ipVANADIUM, @object
	.size	_ZL10ipVANADIUM, 4
_ZL10ipVANADIUM:
	.long	22
	.align 4
	.type	_ZL10ipCHROMIUM, @object
	.size	_ZL10ipCHROMIUM, 4
_ZL10ipCHROMIUM:
	.long	23
	.align 4
	.type	_ZL11ipMANGANESE, @object
	.size	_ZL11ipMANGANESE, 4
_ZL11ipMANGANESE:
	.long	24
	.align 4
	.type	_ZL6ipIRON, @object
	.size	_ZL6ipIRON, 4
_ZL6ipIRON:
	.long	25
	.align 4
	.type	_ZL8ipCOBALT, @object
	.size	_ZL8ipCOBALT, 4
_ZL8ipCOBALT:
	.long	26
	.align 4
	.type	_ZL8ipNICKEL, @object
	.size	_ZL8ipNICKEL, 4
_ZL8ipNICKEL:
	.long	27
	.align 4
	.type	_ZL8ipCOPPER, @object
	.size	_ZL8ipCOPPER, 4
_ZL8ipCOPPER:
	.long	28
	.align 4
	.type	_ZL6ipZINC, @object
	.size	_ZL6ipZINC, 4
_ZL6ipZINC:
	.long	29
	.align 4
	.type	_ZL9ipKRYPTON, @object
	.size	_ZL9ipKRYPTON, 4
_ZL9ipKRYPTON:
	.long	35
	.align 4
	.type	_ZL7MA_VERS, @object
	.size	_ZL7MA_VERS, 12
_ZL7MA_VERS:
	.long	120070905
	.long	220070803
	.long	320071126
	.align 8
	.type	_ZL2EE, @object
	.size	_ZL2EE, 8
_ZL2EE:
	.long	2333366121
	.long	1074118410
	.align 8
	.type	_ZL5EULER, @object
	.size	_ZL5EULER, 8
_ZL5EULER:
	.long	4235179545
	.long	1071806604
	.align 8
	.type	_ZL2PI, @object
	.size	_ZL2PI, 8
_ZL2PI:
	.long	1413754136
	.long	1074340347
	.align 8
	.type	_ZL3PI2, @object
	.size	_ZL3PI2, 8
_ZL3PI2:
	.long	1413754136
	.long	1075388923
	.align 8
	.type	_ZL3PI4, @object
	.size	_ZL3PI4, 8
_ZL3PI4:
	.long	1413754136
	.long	1076437499
	.align 8
	.type	_ZL3PI8, @object
	.size	_ZL3PI8, 8
_ZL3PI8:
	.long	1413754136
	.long	1077486075
	.align 8
	.type	_ZL5SQRT2, @object
	.size	_ZL5SQRT2, 8
_ZL5SQRT2:
	.long	1719614413
	.long	1073127582
	.align 8
	.type	_ZL6SQRTPI, @object
	.size	_ZL6SQRTPI, 8
_ZL6SQRTPI:
	.long	2444554091
	.long	1073503224
	.align 8
	.type	_ZL9SQRTPIBY2, @object
	.size	_ZL9SQRTPIBY2, 8
_ZL9SQRTPIBY2:
	.long	536225542
	.long	1072958867
	.align 8
	.type	_ZL6LN_TWO, @object
	.size	_ZL6LN_TWO, 8
_ZL6LN_TWO:
	.long	4277811695
	.long	1072049730
	.align 8
	.type	_ZL6LN_TEN, @object
	.size	_ZL6LN_TEN, 8
_ZL6LN_TEN:
	.long	3149223190
	.long	1073900465
	.align 8
	.type	_ZL7LOG10_E, @object
	.size	_ZL7LOG10_E, 8
_ZL7LOG10_E:
	.long	354870542
	.long	1071369083
	.align 8
	.type	_ZL12OPTDEP2EXTIN, @object
	.size	_ZL12OPTDEP2EXTIN, 8
_ZL12OPTDEP2EXTIN:
	.long	3979890473
	.long	1072783148
	.align 8
	.type	_ZL6RADIAN, @object
	.size	_ZL6RADIAN, 8
_ZL6RADIAN:
	.long	442745336
	.long	1078765020
	.align 8
	.type	_ZL10SOLAR_MASS, @object
	.size	_ZL10SOLAR_MASS, 8
_ZL10SOLAR_MASS:
	.long	3128952630
	.long	1188594248
	.align 8
	.type	_ZL16SOLAR_LUMINOSITY, @object
	.size	_ZL16SOLAR_LUMINOSITY, 8
_ZL16SOLAR_LUMINOSITY:
	.long	729510093
	.long	1189588662
	.align 8
	.type	_ZL2AU, @object
	.size	_ZL2AU, 8
_ZL2AU:
	.long	2527354880
	.long	1118516785
	.align 8
	.type	_ZL16ATOMIC_MASS_UNIT, @object
	.size	_ZL16ATOMIC_MASS_UNIT, 8
_ZL16ATOMIC_MASS_UNIT:
	.long	2490883434
	.long	989859659
	.align 8
	.type	_ZL13ELECTRON_MASS, @object
	.size	_ZL13ELECTRON_MASS, 8
_ZL13ELECTRON_MASS:
	.long	2176129345
	.long	978455297
	.align 8
	.type	_ZL11PROTON_MASS, @object
	.size	_ZL11PROTON_MASS, 8
_ZL11PROTON_MASS:
	.long	4294041425
	.long	989867317
	.align 8
	.type	_ZL9BOLTZMANN, @object
	.size	_ZL9BOLTZMANN, 8
_ZL9BOLTZMANN:
	.long	3966603471
	.long	1017374129
	.align 8
	.type	_ZL10SPEEDLIGHT, @object
	.size	_ZL10SPEEDLIGHT, 8
_ZL10SPEEDLIGHT:
	.long	4087349248
	.long	1109126043
	.align 8
	.type	_ZL7HPLANCK, @object
	.size	_ZL7HPLANCK, 8
_ZL7HPLANCK:
	.long	1102295869
	.long	981493701
	.align 8
	.type	_ZL8AVOGADRO, @object
	.size	_ZL8AVOGADRO, 8
_ZL8AVOGADRO:
	.long	202468352
	.long	1155522950
	.align 8
	.type	_ZL10GRAV_CONST, @object
	.size	_ZL10GRAV_CONST, 8
_ZL10GRAV_CONST:
	.long	1499401656
	.long	1047652922
	.align 8
	.type	_ZL11ELEM_CHARGE, @object
	.size	_ZL11ELEM_CHARGE, 8
_ZL11ELEM_CHARGE:
	.long	402043217
	.long	1007133914
	.align 8
	.type	_ZL7RYD_INF, @object
	.size	_ZL7RYD_INF, 8
_ZL7RYD_INF:
	.long	218898334
	.long	1090177685
	.align 8
	.type	_ZL7HIONPOT, @object
	.size	_ZL7HIONPOT, 8
_ZL7HIONPOT:
	.long	802766897
	.long	1072692129
	.align 8
	.type	_ZL6AS1RAD, @object
	.size	_ZL6AS1RAD, 8
_ZL6AS1RAD:
	.long	1932635515
	.long	1091120582
	.align 8
	.type	_ZL6PARSEC, @object
	.size	_ZL6PARSEC, 8
_ZL6PARSEC:
	.long	3031701239
	.long	1137011011
	.align 8
	.type	_ZL10MEGAPARSEC, @object
	.size	_ZL10MEGAPARSEC, 8
_ZL10MEGAPARSEC:
	.long	1438715703
	.long	1157917527
	.align 8
	.type	_ZL5H_BAR, @object
	.size	_ZL5H_BAR, 8
_ZL5H_BAR:
	.long	1563972901
	.long	978641744
	.align 8
	.type	_ZL15ELEM_CHARGE_ESU, @object
	.size	_ZL15ELEM_CHARGE_ESU, 8
_ZL15ELEM_CHARGE_ESU:
	.long	2152232907
	.long	1040220401
	.align 8
	.type	_ZL6ERG1CM, @object
	.size	_ZL6ERG1CM, 8
_ZL6ERG1CM:
	.long	3878108119
	.long	1017946288
	.align 8
	.type	_ZL4T1CM, @object
	.size	_ZL4T1CM, 8
_ZL4T1CM:
	.long	4250265231
	.long	1073153338
	.align 8
	.type	_ZL8KJMOL1CM, @object
	.size	_ZL8KJMOL1CM, 8
_ZL8KJMOL1CM:
	.long	2755650162
	.long	1065910240
	.align 8
	.type	_ZL7WAVNRYD, @object
	.size	_ZL7WAVNRYD, 8
_ZL7WAVNRYD:
	.long	3855215639
	.long	1055071315
	.align 8
	.type	_ZL6RYDLAM, @object
	.size	_ZL6RYDLAM, 8
_ZL6RYDLAM:
	.long	3949343864
	.long	1082948130
	.align 8
	.type	_ZL6EN1RYD, @object
	.size	_ZL6EN1RYD, 8
_ZL6EN1RYD:
	.long	1885789440
	.long	1035466699
	.align 8
	.type	_ZL6TE1RYD, @object
	.size	_ZL6TE1RYD, 8
_ZL6TE1RYD:
	.long	737312215
	.long	1090733564
	.align 8
	.type	_ZL6EVDEGK, @object
	.size	_ZL6EVDEGK, 8
_ZL6EVDEGK:
	.long	2021896848
	.long	1086761538
	.align 8
	.type	_ZL5EVRYD, @object
	.size	_ZL5EVRYD, 8
_ZL5EVRYD:
	.long	1416395336
	.long	1076573725
	.align 8
	.type	_ZL5EN1EV, @object
	.size	_ZL5EN1EV, 8
_ZL5EN1EV:
	.long	3504545695
	.long	1031548815
	.align 8
	.type	_ZL6FR1RYD, @object
	.size	_ZL6FR1RYD, 8
_ZL6FR1RYD:
	.long	181687213
	.long	1126654000
	.align 8
	.type	_ZL9FR1RYDHYD, @object
	.size	_ZL9FR1RYDHYD, 8
_ZL9FR1RYDHYD:
	.long	3258471094
	.long	1126653182
	.align 8
	.type	_ZL6HBAReV, @object
	.size	_ZL6HBAReV, 8
_ZL6HBAReV:
	.long	4178815490
	.long	1019721454
	.align 8
	.type	_ZL9RYDLAMHYD, @object
	.size	_ZL9RYDLAMHYD, 8
_ZL9RYDLAMHYD:
	.long	399750008
	.long	1082949127
	.align 8
	.type	_ZL8FREQ_1EV, @object
	.size	_ZL8FREQ_1EV, 8
_ZL8FREQ_1EV:
	.long	1946364709
	.long	1122729286
	.align 8
	.type	_ZL10SEXP_LIMIT, @object
	.size	_ZL10SEXP_LIMIT, 8
_ZL10SEXP_LIMIT:
	.long	0
	.long	1079312384
	.align 8
	.type	_ZL11DSEXP_LIMIT, @object
	.size	_ZL11DSEXP_LIMIT, 8
_ZL11DSEXP_LIMIT:
	.long	0
	.long	1082474496
	.align 4
	.type	_ZL10LIMTABDLAW, @object
	.size	_ZL10LIMTABDLAW, 4
_ZL10LIMTABDLAW:
	.long	500
	.align 4
	.type	_ZL5ipH1s, @object
	.size	_ZL5ipH1s, 4
_ZL5ipH1s:
	.zero	4
	.align 4
	.type	_ZL5ipH2s, @object
	.size	_ZL5ipH2s, 4
_ZL5ipH2s:
	.long	1
	.align 4
	.type	_ZL5ipH2p, @object
	.size	_ZL5ipH2p, 4
_ZL5ipH2p:
	.long	2
	.align 4
	.type	_ZL5ipH3s, @object
	.size	_ZL5ipH3s, 4
_ZL5ipH3s:
	.long	3
	.align 4
	.type	_ZL5ipH3p, @object
	.size	_ZL5ipH3p, 4
_ZL5ipH3p:
	.long	4
	.align 4
	.type	_ZL5ipH3d, @object
	.size	_ZL5ipH3d, 4
_ZL5ipH3d:
	.long	5
	.align 4
	.type	_ZL5ipH4s, @object
	.size	_ZL5ipH4s, 4
_ZL5ipH4s:
	.long	6
	.align 4
	.type	_ZL5ipH4p, @object
	.size	_ZL5ipH4p, 4
_ZL5ipH4p:
	.long	7
	.align 4
	.type	_ZL5ipH4d, @object
	.size	_ZL5ipH4d, 4
_ZL5ipH4d:
	.long	8
	.align 4
	.type	_ZL5ipH4f, @object
	.size	_ZL5ipH4f, 4
_ZL5ipH4f:
	.long	9
	.align 4
	.type	_ZL8ipHe1s1S, @object
	.size	_ZL8ipHe1s1S, 4
_ZL8ipHe1s1S:
	.zero	4
	.align 4
	.type	_ZL8ipHe2s3S, @object
	.size	_ZL8ipHe2s3S, 4
_ZL8ipHe2s3S:
	.long	1
	.align 4
	.type	_ZL8ipHe2s1S, @object
	.size	_ZL8ipHe2s1S, 4
_ZL8ipHe2s1S:
	.long	2
	.align 4
	.type	_ZL9ipHe2p3P0, @object
	.size	_ZL9ipHe2p3P0, 4
_ZL9ipHe2p3P0:
	.long	3
	.align 4
	.type	_ZL9ipHe2p3P1, @object
	.size	_ZL9ipHe2p3P1, 4
_ZL9ipHe2p3P1:
	.long	4
	.align 4
	.type	_ZL9ipHe2p3P2, @object
	.size	_ZL9ipHe2p3P2, 4
_ZL9ipHe2p3P2:
	.long	5
	.align 4
	.type	_ZL8ipHe2p1P, @object
	.size	_ZL8ipHe2p1P, 4
_ZL8ipHe2p1P:
	.long	6
	.align 4
	.type	_ZL8ipHe3s3S, @object
	.size	_ZL8ipHe3s3S, 4
_ZL8ipHe3s3S:
	.long	7
	.align 4
	.type	_ZL8ipHe3s1S, @object
	.size	_ZL8ipHe3s1S, 4
_ZL8ipHe3s1S:
	.long	8
	.align 4
	.type	_ZL8ipHe3p3P, @object
	.size	_ZL8ipHe3p3P, 4
_ZL8ipHe3p3P:
	.long	9
	.align 4
	.type	_ZL8ipHe3d3D, @object
	.size	_ZL8ipHe3d3D, 4
_ZL8ipHe3d3D:
	.long	10
	.align 4
	.type	_ZL8ipHe3d1D, @object
	.size	_ZL8ipHe3d1D, 4
_ZL8ipHe3d1D:
	.long	11
	.align 4
	.type	_ZL8ipHe3p1P, @object
	.size	_ZL8ipHe3p1P, 4
_ZL8ipHe3p1P:
	.long	12
	.align 4
	.type	_ZL8ipH_LIKE, @object
	.size	_ZL8ipH_LIKE, 4
_ZL8ipH_LIKE:
	.zero	4
	.align 4
	.type	_ZL9ipHE_LIKE, @object
	.size	_ZL9ipHE_LIKE, 4
_ZL9ipHE_LIKE:
	.long	1
	.align 4
	.type	_ZL9ipLI_LIKE, @object
	.size	_ZL9ipLI_LIKE, 4
_ZL9ipLI_LIKE:
	.long	2
	.align 4
	.type	_ZL9ipBE_LIKE, @object
	.size	_ZL9ipBE_LIKE, 4
_ZL9ipBE_LIKE:
	.long	3
	.align 4
	.type	_ZL8ipB_LIKE, @object
	.size	_ZL8ipB_LIKE, 4
_ZL8ipB_LIKE:
	.long	4
	.align 4
	.type	_ZL8ipC_LIKE, @object
	.size	_ZL8ipC_LIKE, 4
_ZL8ipC_LIKE:
	.long	5
	.align 4
	.type	_ZL8ipN_LIKE, @object
	.size	_ZL8ipN_LIKE, 4
_ZL8ipN_LIKE:
	.long	6
	.align 4
	.type	_ZL8ipO_LIKE, @object
	.size	_ZL8ipO_LIKE, 4
_ZL8ipO_LIKE:
	.long	7
	.align 4
	.type	_ZL8ipF_LIKE, @object
	.size	_ZL8ipF_LIKE, 4
_ZL8ipF_LIKE:
	.long	8
	.align 4
	.type	_ZL9ipNE_LIKE, @object
	.size	_ZL9ipNE_LIKE, 4
_ZL9ipNE_LIKE:
	.long	9
	.align 4
	.type	_ZL9ipNA_LIKE, @object
	.size	_ZL9ipNA_LIKE, 4
_ZL9ipNA_LIKE:
	.long	10
	.align 4
	.type	_ZL9ipMG_LIKE, @object
	.size	_ZL9ipMG_LIKE, 4
_ZL9ipMG_LIKE:
	.long	11
	.align 4
	.type	_ZL9ipAL_LIKE, @object
	.size	_ZL9ipAL_LIKE, 4
_ZL9ipAL_LIKE:
	.long	12
	.align 4
	.type	_ZL9ipSI_LIKE, @object
	.size	_ZL9ipSI_LIKE, 4
_ZL9ipSI_LIKE:
	.long	13
	.align 4
	.type	_ZL8ipP_LIKE, @object
	.size	_ZL8ipP_LIKE, 4
_ZL8ipP_LIKE:
	.long	14
	.align 4
	.type	_ZL8ipS_LIKE, @object
	.size	_ZL8ipS_LIKE, 4
_ZL8ipS_LIKE:
	.long	15
	.align 4
	.type	_ZL9ipCL_LIKE, @object
	.size	_ZL9ipCL_LIKE, 4
_ZL9ipCL_LIKE:
	.long	16
	.align 4
	.type	_ZL9ipAR_LIKE, @object
	.size	_ZL9ipAR_LIKE, 4
_ZL9ipAR_LIKE:
	.long	17
	.align 4
	.type	_ZL10LIMTEMPTAB, @object
	.size	_ZL10LIMTEMPTAB, 4
_ZL10LIMTEMPTAB:
	.long	500
	.align 4
	.type	_ZL5NDEMS, @object
	.size	_ZL5NDEMS, 4
_ZL5NDEMS:
	.long	200
	.align 8
	.type	_ZL10GRAIN_TMIN, @object
	.size	_ZL10GRAIN_TMIN, 8
_ZL10GRAIN_TMIN:
	.long	3539053052
	.long	1062232653
	.align 8
	.type	_ZL10GRAIN_TMID, @object
	.size	_ZL10GRAIN_TMID, 8
_ZL10GRAIN_TMID:
	.long	0
	.long	1085507584
	.align 8
	.type	_ZL10GRAIN_TMAX, @object
	.size	_ZL10GRAIN_TMAX, 8
_ZL10GRAIN_TMAX:
	.long	0
	.long	1104273827
	.align 4
	.type	_ZL4NCHS, @object
	.size	_ZL4NCHS, 4
_ZL4NCHS:
	.long	30
	.align 4
	.type	_ZL4NCHU, @object
	.size	_ZL4NCHU, 4
_ZL4NCHU:
	.long	10
	.align 4
	.type	_ZL13NCHRG_DEFAULT, @object
	.size	_ZL13NCHRG_DEFAULT, 4
_ZL13NCHRG_DEFAULT:
	.long	2
	.align 4
	.type	_ZL6NQGRID, @object
	.size	_ZL6NQGRID, 4
_ZL6NQGRID:
	.long	10000
	.align 8
	.type	_ZL11CONSERV_TOL, @object
	.size	_ZL11CONSERV_TOL, 8
_ZL11CONSERV_TOL:
	.long	3539053052
	.long	1062232653
	.align 4
	.type	_ZL12N_X_COLLIDER, @object
	.size	_ZL12N_X_COLLIDER, 4
_ZL12N_X_COLLIDER:
	.long	5
	.align 4
	.type	_ZL14chN_X_COLLIDER, @object
	.size	_ZL14chN_X_COLLIDER, 4
_ZL14chN_X_COLLIDER:
	.long	10
	.align 4
	.type	_ZL10nTE_HMINUS, @object
	.size	_ZL10nTE_HMINUS, 4
_ZL10nTE_HMINUS:
	.long	7
	.align 4
	.type	_ZL6N_ELEC, @object
	.size	_ZL6N_ELEC, 4
_ZL6N_ELEC:
	.long	7
	.align 16
	.type	_ZL15H2_logte_hminus, @object
	.size	_ZL15H2_logte_hminus, 28
_ZL15H2_logte_hminus:
	.long	1065353216
	.long	1069355589
	.long	1073741824
	.long	1075743010
	.long	1077936128
	.long	1079937314
	.long	1082130432
	.text
	.type	_GLOBAL__sub_I__Z10mole_solvev, @function
_GLOBAL__sub_I__Z10mole_solvev:
.LFB5836:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	movl	$65535, %esi
	movl	$1, %edi
	call	_Z41__static_initialization_and_destruction_0ii
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE5836:
	.size	_GLOBAL__sub_I__Z10mole_solvev, .-_GLOBAL__sub_I__Z10mole_solvev
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I__Z10mole_solvev
	.local	_ZZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c
	.comm	_ZZL6funjacR8GroupMapRKSt8valarrayIdEPdS5_bPbE1c,144,32
	.section	.rodata
	.align 8
	.type	_ZZL18__gthread_active_pvE20__gthread_active_ptr, @object
	.size	_ZZL18__gthread_active_pvE20__gthread_active_ptr, 8
_ZZL18__gthread_active_pvE20__gthread_active_ptr:
	.quad	_ZL28__gthrw___pthread_key_createPjPFvPvE
	.weakref	_ZL28__gthrw___pthread_key_createPjPFvPvE,__pthread_key_create
	.align 4
.LC0:
	.long	8388608
	.align 4
.LC1:
	.long	2139095039
	.align 16
.LC3:
	.long	2147483647
	.long	0
	.long	0
	.long	0
	.align 16
.LC4:
	.long	4294967295
	.long	2147483647
	.long	0
	.long	0
	.align 8
.LC12:
	.long	3794832442
	.long	1044740494
	.align 4
.LC13:
	.long	3212836864
	.align 16
.LC15:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.align 8
.LC16:
	.long	2577960615
	.long	-1158457406
	.align 8
.LC19:
	.long	0
	.long	1076101120
	.align 8
.LC20:
	.long	1202590843
	.long	1065646817
	.align 8
.LC21:
	.long	2576980378
	.long	1069128089
	.align 8
.LC23:
	.long	2576980378
	.long	1070176665
	.align 8
.LC29:
	.long	3539053052
	.long	1062232653
	.align 8
.LC39:
	.long	0
	.long	1072168960
	.align 8
.LC40:
	.long	0
	.long	1070596096
	.align 4
.LC41:
	.long	0
	.align 8
.LC42:
	.long	0
	.long	1072693248
	.align 8
.LC45:
	.long	2577960615
	.long	989026242
	.align 8
.LC49:
	.long	3654794683
	.long	1037794527
	.align 4
.LC61:
	.long	1120403456
	.align 8
.LC63:
	.long	1413754136
	.long	1076437499
	.align 8
.LC64:
	.long	4087349248
	.long	1109126043
	.align 8
.LC65:
	.long	3892314112
	.long	1110919286
	.align 8
.LC66:
	.long	1102295869
	.long	981493701
	.align 8
.LC67:
	.long	2308816496
	.long	925998950
	.align 8
.LC69:
	.long	3878108119
	.long	1018994864
	.align 8
.LC70:
	.long	3966603471
	.long	1017374129
	.align 8
.LC71:
	.long	1413754136
	.long	1074340347
	.align 8
.LC72:
	.long	1563972901
	.long	978641744
	.align 8
.LC73:
	.long	0
	.long	1078853632
	.align 8
.LC75:
	.long	3346327307
	.long	1093995191
	.align 8
.LC76:
	.long	0
	.long	1075970048
	.align 8
.LC77:
	.long	181687213
	.long	1126654000
	.align 8
.LC78:
	.long	0
	.long	1084227584
	.align 8
.LC80:
	.long	4023219264
	.long	985253115
	.align 8
.LC81:
	.long	2176129345
	.long	978455297
	.hidden	__dso_handle
	.ident	"GCC: (Ubuntu 4.8.5-4ubuntu8) 4.8.5"
	.section	.note.GNU-stack,"",@progbits
