import cython
cimport numpy as np
import numpy as np

@cython.profile(True)
@cython.cdivision(True)
cdef double calc_f1(double recall, double precision):
    """
    Static method to calculate the F1 statistic given precision
    and recall (order is unimportant). Definition:
    F1 = (2 * precision * recall) / (precision + recall)
    """
    cdef double result, summa, product

    if precision < 0 or recall < 0:
        raise ValueError("Negative values are an invalid input! ({0}, {1})".format(
            recall, precision))

    elif precision == 0 or recall == 0:
        return 0
    else:
        product = 2 * precision * recall
        summa = precision + recall
        result = product / summa
        return result


@cython.cdivision(True)
cpdef np.ndarray array_compare(np.ndarray[DTYPE_t] ref, np.ndarray[DTYPE_t] pred):

    # Calculate the junction intersection
    cdef:
        np.ndarray[DTYPE_t] intersect, intersect_junc
        np.ndarray result
        np.ndarray[DTYPE_t] ref_junc, pred_junc
        double junction_recall, junction_precision, junction_f1
        double exon_recall, exon_precision, exon_f1


    intersect = np.intersect1d(ref, pred)
    exon_recall = <double>intersect.shape[0] / <double>ref.shape[0]
    exon_precision = <double>intersect.shape[0] / <double>pred.shape[0]
    exon_f1 = calc_f1(exon_recall, exon_precision)

    if ref.shape[0] > 2 and pred.shape[0] > 2:
        ref_junc = ref[1:-1]
        pred_junc = pred[1:-1]
        intersect_junc = np.intersect1d(ref_junc, pred_junc)
        junction_recall = <double>intersect_junc.shape[0] / <double>ref_junc.shape[0]
        junction_precision = <double>intersect_junc.shape[0] / <double>pred_junc.shape[0]
        junction_f1 = calc_f1(junction_recall, junction_precision)
        pass
    else:
        junction_recall, junction_precision, junction_f1 = 0, 0, 0

    # TODO: Find a way to calculate faux nucleotide stats as well


    result = np.array([exon_recall, exon_precision, exon_f1,
                        junction_recall, junction_precision, junction_f1])

    return result
