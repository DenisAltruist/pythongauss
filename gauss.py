def gauss(A, b):
    B = A.copy()
    c = b.copy()
    B = np.column_stack((B, c))
    m = B.shape[0]
    n = B.shape[1]
    k = 0
    freeVariables = np.zeros(n - 1, dtype = np.int)
    notFreeVariables = np.zeros(n - 1, dtype = np.int)
    cntOfFreeVariables = 0
    cntOfNotFreeVariables = 0
    for j in range(n - 1):
        idxOfNonZero = -1
        for i in range(k, m):
            if (B[i, j] != 0):
                idxOfNonZero = i
        if idxOfNonZero != -1:
            B[k, :], B[idxOfNonZero, :] = B[idxOfNonZero, :].copy(), B[k, :].copy()
            value = B[k, j]
            for t in range(j, n):
                B[k, t] /= value
            for t in range(m):
                if t == k:
                    continue
                B[t, :] -= B[t, j] * B[k, :]
            notFreeVariables[cntOfNotFreeVariables] = j
            cntOfNotFreeVariables += 1
            k += 1
        else:
            freeVariables[cntOfFreeVariables] = j
            cntOfFreeVariables += 1
    X = np.zeros(n - 1)
    for i in range(k, m):
        if B[i, n - 1] != 0:
            return "Has not solutions", X
    if (cntOfFreeVariables == 0):
        for j in range(n - 1):
            X[j] = B[j,n - 1]
        return "Has one solution", X
    ShiftVector = np.zeros(n - 1)
    for i in range(cntOfNotFreeVariables):
        idx = int(notFreeVariables[i])
        ShiftVector[idx] = B[i,n-1]
    
    numOfStrInBasis = np.zeros(n, dtype = np.int )
    BasisMatrix = np.zeros((cntOfFreeVariables, n - 1))
    for i in range (cntOfFreeVariables):
        idx = freeVariables[i]
        BasisMatrix[i][idx] = 1
        numOfStrInBasis[idx] = i
    for i in range (cntOfNotFreeVariables):
        j = notFreeVariables[i]
        for t in range (j + 1, n - 1):
            if B[i][t] != 0:
                BasisMatrix[numOfStrInBasis[t]][j] = -B[i][t]    
    return "Has infinitely many solutions", ShiftVector, BasisMatrix