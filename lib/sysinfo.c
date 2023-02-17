#define SYSINFO_SRC

#include "sysinfo.h"

#if defined __amd64__ || __x86_64__
#if defined __linux__
#include <unistd.h>
void
get_sysinfo(SYSINFO_ptr data_in)
{
    data_in->cpu_arch = CPU_ARCH_X86;
    data_in->os = OS_LINUX;
    long phypz = sysconf(_SC_PHYS_PAGES);
    long psize = sysconf(_SC_PAGE_SIZE);
    data_in->mem_size_byte = phypz * psize;
    unsigned int cpu, cpu_available;
    cpu = sysconf(_SC_NPROCESSORS_CONF);
    cpu_available = sysconf(_SC_NPROCESSORS_ONLN);
    data_in->cpu_total = cpu;
    data_in->cpu_available = cpu_available;

    return;
}

#endif
#endif


#ifdef SYSINFO_TEST
#include <stdio.h>
int
main(void)
{
    SYSINFO data_in;
    get_sysinfo(&data_in);
    printf("%lu %u\n", data_in.mem_size_byte, data_in.cpu_total);
    return 0;
}
#endif

