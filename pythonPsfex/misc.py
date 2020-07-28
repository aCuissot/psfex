import numpy as np

def fast_median(array, n):
    return np.percentile(array, 50)

def fast_quantile(array, n, quantile):
        return np.percentile(array, int(quantile*100))

def dqmedian(array, n):
    return np.percentile(array, 50)
    
def fqmedian(array, n):
    return np.percentile(array, 50)
