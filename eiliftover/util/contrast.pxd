import numpy as np
cimport numpy as np

DTYPE = np.int
ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_f

cpdef np.ndarray[DTYPE_f, ndim=6] array_compare(np.ndarray[DTYPE_t] ref, np.ndarray[DTYPE_t] pred, double identity)