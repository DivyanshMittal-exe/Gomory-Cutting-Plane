import numpy as np

tol = 1e-7


def dual_simplex(T,basis):
    m,n = T.shape
    n = n-1
    m = m-1
    while True:
        if np.all(T[:-1,-1] >= -tol):
            x = np.zeros(n)
            x[basis] = T[:-1, -1]
            return T,basis,x
        
        for it in range(m):
            if T[it,-1] < 0:
                l = it
                
                
                if np.all(T[l,:-1] >= -tol):
                    return None, None, None

                min_index = -1
                min_ratio = 1e9
                for rit in range(n):
                    if(T[l,rit] < 0):
                        if T[-1,rit]/abs(T[l,rit]) < min_ratio:
                            min_index = rit
                            min_ratio = T[-1,rit]/abs(T[l,rit])
                            
                j = min_index
        
                if min_index == -1:
                    # unbounded solution
                    return None, None,None
                
                
                # step 3: perform pivot
                basis[l] = j
                
                T[l,:] /= T[l,j]
                for k in range(m+1):
                    if k != l:
                        T[k,:] -= T[k,j] * T[l,:]
                                   
                break        
        # j = np.argmin(T[-1,:-1])
    

def tableau_phase1(T,n,m,basis):
    while True:
        j = np.argmin(T[-1,:-1])
        if T[-1,j] >= -tol:
            # optimal solution found
            x = np.zeros(n)
            # TODO
            # x[basis < n] = T[basis < n, -1]
            return T,basis, x
        # step 2: select leaving variable
        # T[:-1,-1]/T[:-1,j]
        
        min_index = -1
        min_ratio = 1e9
        for it in range(m):
            if T[it,j] > 0:
                if T[it,-1]/T[it,j] < min_ratio:
                    min_index = it
                    min_ratio = T[it,-1]/T[it,j] 
        
        # i = np.argmin(T[:-1,-1]/T[:-1,j]) #Doubt about correctness of this
        
        i = min_index
        
        if min_index == -1:
            # unbounded solution
            return None, None,None
        
        
        # step 3: perform pivot
        basis[i] = j
        
        T[i,:] /= T[i,j]
        for k in range(m+1):
            if k != i:
                T[k,:] -= T[k,j] * T[i,:]


def tableau_phase2(T,n,m,basis):
    while True:
        j = np.argmin(T[-1,:-1])
        if T[-1,j] >= -tol:
            # optimal solution found
            x = np.zeros(n)
            x[basis] = T[:-1, -1]
            return T,basis, x
        # step 2: select leaving variable
        # T[:-1,-1]/T[:-1,j]
        
        min_index = -1
        min_ratio = 1e9
        for it in range(m):
            if T[it,j] > 0:
                if T[it,-1]/T[it,j] < min_ratio:
                    min_index = it
                    min_ratio = T[it,-1]/T[it,j] 
        
        # i = np.argmin(T[:-1,-1]/T[:-1,j]) #Doubt about correctness of this
        
        i = min_index
        
        if min_index == -1:
            # unbounded solution
            return None, None,None
        
        
        # step 3: perform pivot
        basis[i] = j
        
        T[i,:] /= T[i,j]
        for k in range(m+1):
            if k != i:
                T[k,:] -= T[k,j] * T[i,:]


def gomory(c, A, b):
    m, n_old = A.shape
    n = n_old + m
    # Converting to standard form Ax = b
    A = np.hstack((A, np.eye(A.shape[0])))
    c = np.hstack((c, np.zeros(A.shape[0])))
    # Minimising c
    c = -c

    # making b positive
    for i in range(m):
        if(b[i] < 0):
            b[i] = -b[i]
            A[i,:] = -A[i,:]

    # Phase1
    # Ax + Iy = b
    c_phase1 = np.concatenate((np.zeros(n), np.ones(m)))
    A_phase1 = np.hstack((A, np.eye(m)))
    basis_indices = np.arange(n, n+m)

    # Constructing initial tableau
    T = np.zeros((m+1, n+m+1))
    # T[1:, 1:] = A
    # T[1:,0] = b
    # T[0, 1:] = c_phase1 - np.ones(m).T @ A
    # T[0,0] = -np.ones(m) @ b
    T[:-1, :-1] = A_phase1
    T[0:-1:,-1] = b
    T[-1, :-1] = c_phase1 - np.ones(m).T @ A_phase1
    T[-1,-1] = -np.ones(m) @ b

    # Solving the initial tableau
    T, basis_indices, x = tableau_phase1(T, n, m, basis_indices)

    if(abs(T[-1,-1]) > tol):
        
        print("F Cost not 0 ")
        return None

    T = np.hstack((T[:,:n],T[:,n+m:]))
    c_b = c[basis_indices]
    T[-1, :-1] = c - c_b @ T[:-1, :-1]
    T[-1,-1] = -c_b @ T[0:-1:,-1]


    T, basis_indices, x = tableau_phase2(T, n, m, basis_indices)



if __name__ == "__main__":

    n, m = map(int, input().split())
    b = np.array(list(map(int, input().split())))
    c = np.array(list(map(int, input().split())))

    A = np.empty((m, n))
    for i in range(m):
        row = np.array(list(map(int, input().split())))
        A[i, :] = row

    print(gomory(c,A,b))
