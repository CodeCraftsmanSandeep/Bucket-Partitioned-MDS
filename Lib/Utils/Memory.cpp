#include <sys/resource.h>

#include "Utils.h"


double get_peak_mb() 
{
    struct rusage u;
    getrusage(RUSAGE_SELF, &u);
    return u.ru_maxrss / 1024.0;  // Convert KB to MB
}
