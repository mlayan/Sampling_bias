from __future__ import print_function

cimport cython

from libc.stdlib cimport malloc, free

# Turn off Cython checks (for speed up)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

def ESS_cython1(double[:] trace, int sampleInterval, int MAX_LAG):
    # Wrapped function to compute the ESS of trace

    cdef int traceLength = trace.shape[0]
    cdef int maxLag = min(traceLength-1, MAX_LAG)
    cdef double varStat = 0.0
    cdef double *gammaStat = <double *> malloc(maxLag * sizeof(double))
    if not gammaStat:
        raise MemoryError()

    cdef int lag
    cdef int j
    cdef int i
    cdef double del1
    cdef double del2

    cdef double mean = 0.0
    for i in range(traceLength):
        mean += trace[i]
    mean /= traceLength

    print(sampleInterval)
    print(traceLength)
    print(maxLag)
    print(varStat)
    print(mean)

    for lag in range(maxLag):
        for j in range(traceLength - lag):
            del1 = trace[j] - mean
            del2 = trace[j + lag] - mean
            gammaStat[lag] += (del1 * del2)

        gammaStat[lag] /= (traceLength - lag)

        if lag == 0:
            varStat = gammaStat[0]
        elif lag % 2 == 0:
            if gammaStat[lag-1] + gammaStat[lag] > 0:
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag])
            else:
                maxLag = lag


    cdef double ACT = 0.0 
    cdef double ESS = 0.0 

    if gammaStat[0] == 0:
        ACT = 0.0
    else :
        ACT = sampleInterval * varStat / gammaStat[0]

    if ACT == 0.0:
        ESS = 1.0
    else:
        ESS = (sampleInterval * traceLength) / ACT

    free(gammaStat)

    return(ESS)



def ESS_cython2(double[:] trace, int sampleInterval, int MAX_LAG):
    # Wrapped function to compute the ESS of trace

    # sum of trace, excluding burn-in
    cdef double sum0 = 0.0    

    # Trace length
    cdef int traceLength = trace.shape[0]

    # keep track of sums of trace(i)*trace(i_+ lag) for all lags, excluding burn-in
    cdef double *squareLaggedSums = <double *> malloc(MAX_LAG * sizeof(double))
    if not squareLaggedSums:
        raise MemoryError()
    cdef double *autoCorrelation = <double *> malloc(MAX_LAG * sizeof(double))
    if not autoCorrelation:
        raise MemoryError()

    cdef double sum1
    cdef double sum2
    cdef double mean
    cdef int i
    cdef int lagIndex

    for i in range(MAX_LAG):
        squareLaggedSums[i] = 0.0
        autoCorrelation[i] = 0.0

    for i in range(traceLength):
        sum0 += trace[i]

        # calculate mean
        mean = sum0 / (i + 1)

        # calculate auto correlation for selected lag times
        sum1 = sum0
        sum2 = sum0
        for lagIndex in range(min(i + 1, MAX_LAG)):
            squareLaggedSums[lagIndex] += trace[i - lagIndex] * trace[i]
            #The following line is the same approximation as in Tracer
            autoCorrelation[lagIndex] = squareLaggedSums[lagIndex] - (sum1 + sum2) * mean + mean * mean * (i + 1 - lagIndex)
            autoCorrelation[lagIndex] /= (i + 1 - lagIndex)
            sum1 -= trace[i - lagIndex]
            sum2 -= trace[lagIndex]

    cdef int maxLag = min(traceLength, MAX_LAG)
    cdef double integralOfACFunctionTimes2 = 0.0
    for lagIndex in range(maxLag):
        if lagIndex == 0:
            integralOfACFunctionTimes2 = autoCorrelation[0]
        elif lagIndex % 2 == 0:
            if autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex] > 0:
                integralOfACFunctionTimes2 += 2.0 * (autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex])
            else:
                break

    # ACT
    cdef double ACT = sampleInterval * integralOfACFunctionTimes2 / autoCorrelation[0]

    # ESS
    cdef double ESS = traceLength / (ACT / sampleInterval)

    # Release memory
    free(squareLaggedSums)
    free(autoCorrelation)

    return(ESS)
