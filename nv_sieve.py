import numpy as np
from numpy.linalg import norm
from utils import remove_zeros

def nguyen_vidick_sieve(S, gamma):
    """
        Runs the NV sieve.

        Parameters:
            S:      The initial set of lattice points
            gamma:  The norm reduction factor

        Returns:
            v:      The shortest vector found by the sieve
    """
    
    S_0 = S
    while len(S) > 0:
        S_0 = S
        S = lattice_sieve(S, gamma)
        S = remove_zeros(S)
        if len(S):
            print("\r Min norm in S: " + str(norm(min(S, key=lambda v: norm(v)))), end="\r")
    
    # Return the vector with smallest vector norm
    return min(S_0, key=lambda v: norm(v))

def lattice_sieve(S, gamma):
    """
        Helper method for the main sieving loop. Builds the next set of the sieve
        by checking to see if a vector is small enough or if there is a 'center' in the 
        list of centers that can be used to reduce the vector.

        Parameters:
            S:      The current sieve set
            gamma:  Norm reduction factor

        Returns:
            S_p:    Set for the next step of the sieve
    """
    
    # Start with setting R as the max norm in our set
    #R = norm(max(S, key=lambda v: norm(v)))
    R = sum(norm(v) for v in S)/len(S)
    # Make an empty list of centers and a list for the next step of the sieve
    C = []
    S_p = []
    gR = gamma*R

    # Run on each vector in S
    for v in S:
        # If the vector is small enough add it to S_p
        if norm(v) <= gR:
            S_p.append(v)
        # Check if there is a close center
        else:
            res, c = exists_close_center(C, v, gR)

            # If a close center is found, reduce the vector, else add it as a center
            S_p.append(v - c) if res else C.append(v)
    
    return S_p

def exists_close_center(C, v, gR):
    """
        Given a vector v, check the list of centers to see if we can 
        reduce v by one of the centers in C

        Parameters:
            C:          The list of centers
            v:          The vector to reduce
            gR:         Norm reduction factor (gamma * R)

        Returns:
            True, c:    The center if it is found
            False, 0:   If no center is found

    """
    
    for c in C:
        if norm(v - c) <= gR:
            return (True, c)
    return (False, 0)
