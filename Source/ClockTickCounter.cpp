#include "ClockTickCounter.h"

unsigned long long int clockTickCounter(void)
{
#if defined (__x86_64__)                      // 64 bit mode
     unsigned hi, lo;
     __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
     return ( (unsigned long long int)lo)|( ((unsigned long long int)hi)<<32 );
#else                                    // 32 bit mode
     unsigned long long int val;
    __asm__ __volatile__("rdtsc" : "=A" (val) : );
     return(val);
#endif
}

/*
#include <sys/time.h>

double getSecond()
{
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
*/

#include <time.h>

#define NANOS_PER_SECF 1000000000.0
#define USECS_PER_SEC 1000000


double monotonic_seconds() {
  struct timespec time;
  // Note: Make sure to link with -lrt to define clock_gettime.
  clock_gettime(CLOCK_MONOTONIC, &time);
  return ((double) time.tv_sec) + ((double) time.tv_nsec / (NANOS_PER_SECF));
}
