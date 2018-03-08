#if NDIM == 2
#define KJI_DO(lo,hi) j = lo, hi; do i = lo, hi
#define CLOSE_DO end do
#define IJK i, j
#define DTIMES(TXT) TXT, TXT
#elif NDIM == 3
#define KJI_DO(lo,hi) j = lo, hi; do i = lo, hi; do k = lo, hi
#define CLOSE_DO end do; end do
#define IJK i, j, k
#define DTIMES(TXT) TXT, TXT, TXT
#endif
