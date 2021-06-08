//
//  l0_mem.h
//  osc
//
//  Created by Futao Zhang on 12/02/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __osc__l0_mem__
#define __osc__l0_mem__

#include "l0_com.h"

extern size_t getMemorySize( );
extern uint64_t getMemSize_Plink( );
extern uint64_t getAllocMB_Plink(uint64_t llxx);
extern bool allocReserved(char** buf, uint64_t size2alloc, string msg); // malloc, should be with dealloc
extern uint64_t allocReserved(char** buf, uint64_t size2need); //
extern void deallocReserved(char** buf, uint64_t size2alloc);

extern bool sudoAllocReserved(uint64_t size2alloc,string msg); // before resize

#endif /* defined(__osc__l0_mem__) */
