#ifndef TMA_CARRAY_H_
#define TMA_CARRAY_H_

#include <tma/types.h>

namespace tma
{

template<class T>
static inline T * ncarray(uidx n) { return (n ? new T[n] : NULL); }

template<class T>
static inline T * zcarray(uidx n) { return (n ? new T[n]() : NULL); }

} /* namespace */

#endif /* TMA_CARRAY_H_ */
