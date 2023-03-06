#define SYSINFO_SRC

#include "sysinfo.h"

#if defined __amd64__ || __x86_64__
#if defined __linux__



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


/*borrowed from https://gist.github.com/ChisholmKyle/0cbedcd3e64132243a39*/
/* recursive mkdir */
int
mkdir_p(const char *dir, const mode_t mode)
{
    char tmp[PATH_MAX_STRING_SIZE];
    char *p = NULL;
    struct stat sb;
    size_t len;
    
    /* copy path */
    len = strnlen (dir, PATH_MAX_STRING_SIZE);
    if (len == 0 || len == PATH_MAX_STRING_SIZE) {
        return -1;
    }
    memcpy (tmp, dir, len);
    tmp[len] = '\0';

    /* remove trailing slash */
    if(tmp[len - 1] == '/') {
        tmp[len - 1] = '\0';
    }

    /* check if path exists and is a directory */
    if (stat (tmp, &sb) == 0) {
        if (S_ISDIR (sb.st_mode)) {
            return 0;
        }
    }
    
    /* recursive mkdir */
    for(p = tmp + 1; *p; p++) {
        if(*p == '/') {
            *p = 0;
            /* test path */
            if (stat(tmp, &sb) != 0) {
                /* path does not exist - create directory */
                if (mkdir(tmp, mode) < 0) {
                    return -1;
                }
            } else if (!S_ISDIR(sb.st_mode)) {
                /* not a directory */
                return -1;
            }
            *p = '/';
        }
    }
    /* test path */
    if (stat(tmp, &sb) != 0) {
        /* path does not exist - create directory */
        if (mkdir(tmp, mode) < 0) {
            return -1;
        }
    } else if (!S_ISDIR(sb.st_mode)) {
        /* not a directory */
        return -1;
    }
    return 0;
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

