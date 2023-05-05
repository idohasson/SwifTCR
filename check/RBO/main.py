from typing import Union, List, Optional
import numpy as np
x=Union[List, np.ndarray]

def rbo(
        S: Union[List, np.ndarray],
        T: Union[List, np.ndarray],
        k: Optional[float] = None,
        p: float = 1.0,
        ext: bool = False,
        verbose=False):
    """
    This the weighted non-conjoint measures, namely, rank-biased overlap.
    Unlike Kendall tau which is correlation based, this is intersection
    based.
    The implementation if from Eq. (4) or Eq. (7) (for p != 1) from the
    RBO paper: http://www.williamwebber.com/research/papers/wmz10_tois.pdf
    If p=1, it returns to the un-bounded set-intersection overlap,
    according to Fagin et al.
    https://researcher.watson.ibm.com/researcher/files/us-fagin/topk.pdf
    The fig. 5 in that RBO paper can be used as test case.
    Note there the choice of p is of great importance, since it
    essentically control the "top-weightness". Simply put, to an extreme,
    a small p value will only consider first few items, whereas a larger p
    value will consider more itmes. See Eq. (21) for quantitative measure.
    Args:
        S, T (list or numpy array): lists with alphanumeric elements. They
            could be of different lengths. Both of the them should be
            ranked, i.e., each element"s position reflects its respective
            ranking in the list. Also we will require that there is no
            duplicate element in each list.
        k (int), default None: The depth of evaluation.
        p (float), default 1.0: Weight of each agreement at depth d:
            p**(d-1). When set to 1.0, there is no weight, the rbo returns
            to average overlap.
        ext (Boolean) default False: If True, we will extrapolate the rbo,
            as in Eq. (23)
        verbose (bool). If True, print out intermediate results.
            Default to False.
    Returns:
        The rbo at depth k (or extrapolated beyond)
    """

    assert type(S) in [list, np.ndarray]
    assert type(T) in [list, np.ndarray]

    assert len(S) == len(set(S))
    assert len(T) == len(set(T))

    assert 0.0 < p < 1.0, "p must be between (0, 1)"

    N_S, N_T = len(S), len(T)

    if not N_S and not N_T:
        return 1  # both lists are empty

    if not N_S or not N_T:
        return 0  # one list empty, one non-empty

    if k is None:
        k = float("inf")
    k = min(N_S, N_T, k)

    # initialize the agreement and average overlap arrays
    A, AO = [0] * k, [0] * k

    if p == 1.0:
        weights = [1.0 for _ in range(k)]
    else:

        weights = [1.0 * (1 - p) * p ** d for d in range(k)]

    for i in range(k):
        i_S = i_T = 0
        while i_S < N_S and i_T < N_T:
            if S[i_S] == T[i_T]:
                A[i] += 1
                i_S += 1
                i_T += 1
            elif S[i_S] in T:
                i_T += 1
            else:
                i_S += 1
        AO[i] = (A[i] + AO[i - 1] * i) / (i + 1) if i > 0 else A[i]
    if ext:
        return 1 - (1 - AO[-1]) * (1 - p) / p
    return AO[-1]


def rbo_ext(S: Union[List, np.ndarray], T: Union[List, np.ndarray], p=0.98, verbose=False):
    """
    This is the ultimate implementation of the rbo, namely, the
    extrapolated version. The corresponding formula is Eq. (32) in the rbo
    paper.
    Args:
        S, T (list or numpy array): lists with alphanumeric elements. They
            could be of different lengths. Both of the them should be
            ranked, i.e., each element's position reflects its respective
            ranking in the list. Also we will require that there is no
            duplicate element in each list.
        p (float): The value p.
        verbose (bool). If True, print out intermediate results.
            Default to False.
    """

    assert type(S) in [list, np.ndarray]
    assert type(T) in [list, np.ndarray]

    assert len(S) == len(set(S))
    assert len(T) == len(set(T))

    N_S, N_T = len(S), len(T)

    if not N_S and not N_T:
        return 1  # both lists are empty

    if not N_S or not N_T:
        return 0  # one list empty, one non-empty

    # since we are dealing with un-even lists, we need to figure out the
    # long (L) and short (S) list first. The name S might be confusing
    # but in this function, S refers to short list, L refers to long list
    if len(S) > len(T):
        L, S = S, T
    else:
        S, L = S, T

    s, l = len(S), len(L)  # noqa

    # initialize the overlap and rbo arrays
    # the agreement can be simply calculated from the overlap
    X, A, rbo = [0] * l, [0] * l, [0] * l

    # first item
    S_running, L_running = {S[0]}, {L[0]}  # for O(1) look up
    X[0] = 1 if S[0] == L[0] else 0
    A[0] = X[0]
    rbo[0] = 1.0 * (1 - p) * A[0]

    # start the calculation
    disjoint = 0
    ext_term = A[0] * p

    for d in tqdm(range(1, l), disable=~verbose):
        if d < s:  # still overlapping in length

            S_running.add(S[d])
            L_running.add(L[d])

            # again I will revoke the DP-like step
            overlap_incr = 0  # overlap increment at step d

            # if the new items are the same
            if S[d] == L[d]:
                overlap_incr += 1
            else:
                # if the new item from S is in L already
                if S[d] in L_running:
                    overlap_incr += 1
                # if the new item from L is in M
                if L[d] in M:
                    A[d] = weights[d]
                    M.remove(L[d])

                # if the new item from M is in L
                if M and M[0] in L:
                    A[d] = weights[d]
                    M.pop(0)

                # update average overlap
                AO[d] = sum(A[:d + 1]) / (d + 1)

    # return RBO value
    if ext:
        if k < min(self.N_S, self.N_T):
            return self._bound_range(sum(A[:-1]) + AO[k-1] * (p**k))
        else:
            return self._bound_range(sum(A))
    else:
        return self._bound_range(sum(A[:k]))



def top_weightness(
        S: Union[List, np.ndarray],
        T: Union[List, np.ndarray],
        p: Optional[float] = None,
        d: Optional[int] = None,
        verbose=False):
    """
    This function will evaluate the degree of the top-weightness of the
    rbo. It is the implementation of Eq. (21) of the rbo paper.
    As a sanity check (per the rbo paper),
    top_weightness(p=0.9, d=10) should be 86%
    top_weightness(p=0.98, d=50) should be 86% too
    Args:
        S (list or numpy array): lists with alphanumeric elements. They
            could be of different lengths. Both of the them should be
            ranked, i.e., each element"s position reflects its respective
            ranking in the list. Also we will require that there is no
            duplicate element in each list.
        T (list or numpy array): lists with alphanumeric elements. They
            could be of different lengths. Both of the them should be
            ranked, i.e., each element"s position reflects its respective
            ranking in the list. Also we will require that there is no
            duplicate element in each list.
        p (float), default None: A value between zero and one.
        d (int), default None: Evaluation depth of the list.
        verbose (bool). If True, print out intermediate results.
            Default to False.
    Returns:
        A float between [0, 1], that indicates the top-weightness.
    """

    N_S, N_T = len(S), len(T)

    assert type(S) in [list, np.ndarray]
    assert type(T) in [list, np.ndarray]

    assert len(S) == len(set(S))
    assert len(T) == len(set(T))

    if p is None:
        p = 0.5

    assert 0.0 < p < 1.0, "p must be between (0, 1)"

    if d is None:
        d = min(N_S, N_T)
    else:
        d = min(N_S, N_T, int(d))

    if d == 0:
        top_w = 1
    elif d == 1:
        top_w = 1 - 1 + 1.0 * (1 - p) / p * (np.log(1.0 / (1 - p)))
    else:
        sum_1 = 0
        for i in range(1, d):
            sum_1 += 1.0 * p ** (i) / i
        top_w = 1 - p ** (i) + 1.0 * (1 - p) / p * (i + 1) * \
                (np.log(1.0 / (1 - p)) - sum_1)  # here i == d-1

    if verbose:
        print("The first {} ranks have {:6.3%} of the weight of "
              "the evaluation.".format(d, top_w))

    return min(1, self._bound_range(sum(AO[:k])))
