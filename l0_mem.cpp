
#include "l0_mem.h"

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <Windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#else
#error "Unable to define getMemorySize( ) for an unknown OS."
#endif



/**
 * Returns the size of physical memory (RAM) in bytes.
 */
size_t getMemorySize( )
{
#if defined(_WIN32) && (defined(__CYGWIN__) || defined(__CYGWIN32__))
    /* Cygwin under Windows. ------------------------------------ */
    /* New 64-bit MEMORYSTATUSEX isn't available.  Use old 32.bit */
    MEMORYSTATUS status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatus( &status );
    return (size_t)status.dwTotalPhys;
    
#elif defined(_WIN32)
    /* Windows. ------------------------------------------------- */
    /* Use new 64-bit MEMORYSTATUSEX, not old 32-bit MEMORYSTATUS */
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx( &status );
    return (size_t)status.ullTotalPhys;
    
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* UNIX variants. ------------------------------------------- */
    /* Prefer sysctl() over sysconf() except sysctl() HW_REALMEM and HW_PHYSMEM */
    
#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))
    int mib[2];
    mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
    mib[1] = HW_MEMSIZE;            /* OSX. --------------------- */
#elif defined(HW_PHYSMEM64)
    mib[1] = HW_PHYSMEM64;          /* NetBSD, OpenBSD. --------- */
#endif
    int64_t size = 0;               /* 64-bit */
    size_t len = sizeof( size );
    if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
        return (size_t)size;
    return 0L;			/* Failed? */
    
#elif defined(_SC_AIX_REALMEM)
    /* AIX. ----------------------------------------------------- */
    return (size_t)sysconf( _SC_AIX_REALMEM ) * (size_t)1024L;
    
#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
    /* FreeBSD, Linux, OpenBSD, and Solaris. -------------------- */
    return (size_t)sysconf( _SC_PHYS_PAGES ) *
    (size_t)sysconf( _SC_PAGESIZE );
    
#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
    /* Legacy. -------------------------------------------------- */
    return (size_t)sysconf( _SC_PHYS_PAGES ) *
    (size_t)sysconf( _SC_PAGE_SIZE );
    
#elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))
    /* DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX. -------- */
    int mib[2];
    mib[0] = CTL_HW;
#if defined(HW_REALMEM)
    mib[1] = HW_REALMEM;		/* FreeBSD. ----------------- */
#elif defined(HW_PYSMEM)
    mib[1] = HW_PHYSMEM;		/* Others. ------------------ */
#endif
    unsigned int size = 0;		/* 32-bit */
    size_t len = sizeof( size );
    if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
        return (size_t)size;
    return 0L;			/* Failed? */
#endif /* sysctl and sysconf variants */
    
#else
    return 0L;			/* Unknown OS. */
#endif
}


uint64_t getMemSize_Plink()
{
    
#ifdef __APPLE__
    int32_t mib[2];
    size_t sztmp;
#endif
    int64_t llxx=0;
    
    // see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    llxx = 0;
    sztmp = sizeof(int64_t);
    sysctl(mib, 2, &llxx, &sztmp, NULL, 0);
    llxx /= 1048576;
#else
#ifdef _WIN32
    memstatus.dwLength = sizeof(memstatus);
    GlobalMemoryStatusEx(&memstatus);
    llxx = memstatus.ullTotalPhys / 1048576;
#else
    llxx = ((uint64_t)sysconf(_SC_PHYS_PAGES)) * ((size_t)sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
      return llxx;
}
uint64_t getAllocMB_Plink(uint64_t llxx)
{
    intptr_t default_alloc_mb=0;
    intptr_t malloc_size_mb = 0;

    if (!llxx) {
        default_alloc_mb = WKSPACE_DEFAULT_MB;
    } else if (llxx < (WKSPACE_MIN_MB * 2)) {
        default_alloc_mb = WKSPACE_MIN_MB;
    } else {
        default_alloc_mb = llxx*3/4;
    }
    if (!malloc_size_mb) {
        malloc_size_mb = default_alloc_mb;
    } else if (malloc_size_mb < WKSPACE_MIN_MB) {
        malloc_size_mb = WKSPACE_MIN_MB;
    }
#ifndef __LP64__
    if (malloc_size_mb > 2047) {
        malloc_size_mb = 2047;
    }
#endif
    return malloc_size_mb;
}
uint64_t allocReserved(char** buf, uint64_t size2need)
{
    uint64_t base=0x7ffe0000;
    while(mem_left<=base/1048576) base>>=1;
    if(base<0x1000000)
    {
        sprintf(logbuf, "Error: Not enough memory for writing buffer.\n");
        logprintb();
        TERMINATE();
    }
    if(base>size2need) base=size2need;
    mem_left-=base/1048576;
    *buf=(char*) malloc (sizeof(char)*(base));
    memset(*buf,0,sizeof(char)*(base));
    sprintf(logbuf, "Reserving %llu MB memory for I/O buffer.\n", base/1048576);
    logprintb();
    return base;
}

bool allocReserved(char** buf, uint64_t size2alloc,string msg)
{
    mem_left-=size2alloc/1048576;
    if(mem_left<0)
    {
        sprintf(logbuf, "Error: can not allocate more memory for %s than the reserved.\n",msg.c_str());
        logprintb();
        return 0;
    }
   *buf=(char*) malloc (sizeof(char)*(size2alloc));
    memset(*buf,0,sizeof(char)*(size2alloc));
    if(size2alloc>>20) sprintf(logbuf, "Reserving %llu MB memory for %s.\n", size2alloc/1048576,msg.c_str());
    else sprintf(logbuf, "Reserving %llu KB memory for %s.\n", size2alloc/1024,msg.c_str());
    logprintb();
   return 1;
}
bool sudoAllocReserved(uint64_t size2alloc,string msg)
{
    mem_left-=size2alloc/1048576;
    if(mem_left<0)
    {
        sprintf(logbuf, "Error: can not allocate more memory for %s than the reserved.\n",msg.c_str());
        logprintb();
        return 0;
    }
    return 1;
}

void deallocReserved(char** buf, uint64_t size2alloc)
{
    if(*buf!=NULL)
    {
        free(*buf);
        mem_left+=size2alloc/1048576;
    }
}
