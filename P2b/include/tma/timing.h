#ifndef TMA_TIMING_H_
#define TMA_TIMING_H_

#include <tma/types.h>

#include <algorithm>
#include <ctime>
#include <unistd.h>
#include <cstdio>

namespace tma
{

#define BILLION  1000000000L;

/******************************************************************************
 * Timing functions
 ******************************************************************************/

//XXX: not thread-safe.
static inline int tic(real &t)
{
  static struct timespec ts0; static struct timespec* tp0 = &ts0;
  static struct timespec ts1; static struct timespec* tp1 = &ts1;
  int ret = clock_gettime(CLOCK_MONOTONIC, tp1);
  t = static_cast<real>(tp1->tv_sec - tp0->tv_sec)
      + static_cast<real>(tp1->tv_nsec - tp0->tv_nsec) / BILLION;
  std::swap(tp0, tp1);
  return ret;
}

//-----------------------------------------------------------------------------

static inline uint64_t rdtsc()
{
  uint32_t lo, hi;
  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return (uint64_t) hi << 32 | lo;
}

typedef uint64_t cpu_ticks;

//-----------------------------------------------------------------------------

} /* namespace */

#endif /* TMA_TIMING_H_ */
