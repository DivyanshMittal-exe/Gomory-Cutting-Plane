import numpy as np

tol = 1e-7


def dual_simplex(T,basis):
    print("Entered Dual Simplex")
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
    # while True:
    while True:
        j = np.argmin(T[-1,:-1])
        if T[-1,j] >= -tol:
            # optimal solution found
            x = np.zeros(n)
            
            to_remove_rows = []
            
            # Removing slack variables
            for base_id,basis_it in enumerate(basis):
                if basis_it >= n:
                    for remover in range(n):
                        if abs(T[base_id,remover]) > tol:
                            basis[base_id] = remover
                            T[base_id,:] /= T[base_id,remover]
                            for k in range(m+1):
                                if k != base_id:
                                    T[k,:] -= T[k,remover] * T[base_id,:]
                            break
                    else:
                        to_remove_rows.append(base_id)
                             
            T = np.delete(T,to_remove_rows,0)
            basis = np.delete(basis,to_remove_rows)

            if len(to_remove_rows):
                print(f"Removing redundant rows {to_remove_rows}")
                pass
                

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
    iters = -1
    while True:
        iters+=1
        
        print("Us")
        print(np.round(T,2))
        print(basis)
        
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
        # print("Basis", basis)
        # print("Ghapla original", "number of iters",iters,'\n', T)
        if min_index == -1:
            print("Basis", basis)
            print("Ghapla original", "number of iters",iters,'\n', T)
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
    
    x_only_basis = T[:-1,-1]
    
    x_correct = np.zeros(n)
    x_correct[basis_indices] = x_only_basis
    
    if np.allclose(np.dot(A, x_correct), b) and np.all(x >= 0):  
        print("IBFS is Correct")

    
    T_p = np.round(T,2)
    np.set_printoptions(suppress=True)
    np.set_printoptions(linewidth=np.inf)


    print("Initial Tableu\n", T_p)
    
    
    print(basis_indices)
    # return None
    m = len(basis_indices)    

    if(abs(T[-1,-1]) > tol):
        
        print("F Cost not 0 ")
        pass
    
    T = np.hstack((T[:,:n],T[:,-1:]))
    c_b = c[basis_indices]
    T[-1, :-1] = c - c_b @ T[:-1, :-1]
    T[-1,-1] = -c_b @ T[0:-1:,-1]

    print(n,m)
    print("Pre Phase 2 Tableau\n", T)

    T, basis_indices, x = tableau_phase2(T, n, m, basis_indices)
    if(T is None):
        print("Ghapla")
    x_only_basis = T[:-1,-1]
    
    x_correct = np.zeros(n)
    x_correct[basis_indices] = x_only_basis

    return x_correct[:n_old]

    return None
    print(T)

    while True:
        # Check if solution is all integers
        print(x)
        x = np.round(x, 5)
        T = np.round(T, 5)
        flag = 1
        if  x is None:
            print("fat gaya")
            return None
        for it in range(len(x)):
            if abs(int(x[it]) - x[it]) > tol:
                flag = 0
                break
        if flag:
            print (x[:n_old].astype(int))
            print(T[-1,-1])
            return x[:n_old].astype(int)
        
        # Add another constraint
        for it in range(m):
            if abs(int(T[it,-1]) - T[it,-1]) > tol:
                j = it

        print("Taking j ", j)

        basis_indices = np.hstack((basis_indices, len(x)))
        b_gomory = T[:-1,-1]
        b_gomory = np.append(b_gomory, -(T[j,-1] - np.floor(T[j,-1])))
        A_gomory = T[:-1, :-1]
        a = np.zeros(n)
        tmp = T[j, :-1]
        for it in range(n):
            if it not in basis_indices and abs(int(tmp[it]) - tmp[it]) > tol:
                a[it] = -(tmp[it] - np.floor(tmp[it]))
        A_gomory = np.vstack([A_gomory, a])

        a = np.zeros(m+1)
        a[-1] = 1
        A_gomory = np.hstack((A_gomory, np.array([a]).T, np.array([b_gomory]).T))

        temp = np.hstack((T[-1, :-1], 0 ,T[-1, -1]))
        T = np.vstack((A_gomory, temp))
        print(T)
        print("Basis",basis_indices)

        m += 1
        n += 1

        T, basis_indices, x = dual_simplex(T, basis_indices)
        print("After Dual")
        print(basis_indices)
        print(T)
        print(x)


if __name__ == "__main__":

    n, m = map(int, input().split())
    b = np.array(list(map(int, input().split())))
    c = np.array(list(map(int, input().split())))

    A = np.empty((m, n))
    for i in range(m):
        row = np.array(list(map(int, input().split())))
        A[i, :] = row

    print(gomory(c,A,b))
