/*
 *  linux/arch/arm/mm/minicache.c
 *
 *  Copyright (C) 2001 Russell King
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation.
 *
 * This handles the mini data cache, as found on SA11x0 and XScale
 * processors.  When we copy a user page page, we map it in such a way
 * that accesses to this page will not touch the main data cache, but
 * will be cached in the mini data cache.  This prevents us thrashing
 * the main data cache on page faults.
 */
#include <linux/init.h>
#include <linux/mm.h>
#include <asm/page.h>
#include <asm/pgtable.h>

#define minicache_address (0xffff2000)
#define minicache_pgprot __pgprot(L_PTE_PRESENT | L_PTE_YOUNG | \
				  L_PTE_CACHEABLE)

static pte_t *minicache_pte;

/*
 * Note that this is intended to be called only from the copy_user_page
 * asm code; anything else will require special locking to prevent the
 * mini-cache space being re-used.  (Note: probably preempt unsafe).
 *
 * We rely on the fact that the minicache is 2K, and we'll be pushing
 * 4K of data through it, so we don't actually have to specifically
 * flush the minicache when we change the mapping.
 *
 * Note also: assert(PAGE_OFFSET <= virt < high_memory).
 * Unsafe: preempt, kmap.
 */
unsigned long map_page_minicache(unsigned long virt)
{
	set_pte(minicache_pte, mk_pte_phys(__pa(virt), minicache_pgprot));
	cpu_tlb_invalidate_page(minicache_address, 0);

	return minicache_address;
}

static int __init minicache_init(void)
{
	pgd_t *pgd;
	pmd_t *pmd;

	pgd = pgd_offset_k(minicache_address);
	pmd = pmd_alloc(&init_mm, pgd, minicache_address);
	if (!pmd)
		BUG();
	minicache_pte = pte_alloc(&init_mm, pmd, minicache_address);
	if (!minicache_pte)
		BUG();

	return 0;
}

__initcall(minicache_init);
