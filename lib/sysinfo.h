#ifndef SYSINFO_HEAD
#define SYSINFO_HEAD

#if !defined SYSINFO_SRC
#define SYSINFO_EXTERN extern
#else
#define SYSINFO_EXTERN
#endif

#ifdef __linux__

#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define CPU_ARCH_X86 1
#define CPU_ARCH_ARM 2
#define CPU_ARCH_RISCV 3
#define OS_LINUX 1
#define OS_MACOS 2
#define OS_WINDOWS 3
#define PATH_MAX_STRING_SIZE 1024

typedef struct sysinfo {
    char cpu_arch;
    char os;
    uint64_t mem_size_byte;
    uint32_t cpu_total;
    uint32_t cpu_available; 
} SYSINFO, *SYSINFO_ptr;

#ifdef __cplusplus
extern "C" {
#endif

SYSINFO_EXTERN void get_sysinfo(SYSINFO_ptr sysinfo_data);

SYSINFO_EXTERN int mkdir_p(const char *dir, const mode_t mode);

#ifdef __cplusplus
}
#endif

#endif

#endif

