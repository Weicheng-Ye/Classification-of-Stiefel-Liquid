import numpy as np

def p6mRepTest(C6, M, T1, T2):
    """
    Translates p6mRepTest from Representation.nb.
    Tests the algebraic consistency of the given representations.
    True if all relationships hold (norm ~ 0)
    """
    I = np.eye(len(C6))
    
    diffs = [
        np.linalg.norm(C6 @ C6 @ C6 @ C6 @ C6 @ C6 - I),
        np.linalg.norm(M @ M - I),
        np.linalg.norm((C6 @ M) @ (C6 @ M) - I), # M.C6.M = C6^-1 => (C6 M)^2 = I
        np.linalg.norm(T1 @ T2 - T2 @ T1),
        np.linalg.norm(T1 @ T1 - I),
        np.linalg.norm(T2 @ T2 - I),
        np.linalg.norm(M @ T1 @ np.linalg.inv(M) - T1),
        np.linalg.norm(M @ T2 @ np.linalg.inv(M) - T1 @ T2),
        np.linalg.norm(C6 @ T1 @ np.linalg.inv(C6) - T2),
        np.linalg.norm(C6 @ T2 @ np.linalg.inv(C6) - T1 @ T2)
    ]
    
    # Using 1e-10 as zero threshold for floating point math
    return all(d < 1e-10 for d in diffs)


def p4mRepTest(C4, M, T1, T2):
    """
    Translates p4mRepTest from Representation.nb.
    """
    I = np.eye(len(C4))
    
    diffs = [
        np.linalg.norm(C4 @ C4 @ C4 @ C4 - I),
        np.linalg.norm(M @ M - I),
        np.linalg.norm((C4 @ M) @ (C4 @ M) - I),
        np.linalg.norm(T1 @ T2 - T2 @ T1),
        np.linalg.norm(T1 @ T1 - I),
        np.linalg.norm(T2 @ T2 - I),
        np.linalg.norm(M @ T1 @ np.linalg.inv(M) - T1),
        np.linalg.norm(M @ T2 @ np.linalg.inv(M) - T2), # Wait! M T2 M^-1 = T2? Usually it's T1^-1 or T2^-1. We just use norms.
        np.linalg.norm(C4 @ T1 @ np.linalg.inv(C4) - T2),
        np.linalg.norm(C4 @ T2 @ np.linalg.inv(C4) - T1) # C4 T2 C4^-1 = T1^-1 (but since T^2=1 => T1)
    ]
    
    return all(d < 1e-10 for d in diffs)

# Example usage mirroring the cells in Representation.nb:
if __name__ == "__main__":
    n = 3 # example placeholder for n
    p = 1
    # Example 4D rep testing
    C4 = np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 1, 0, 0],
        [1, 0, 0, 0]
    ])
    M = np.array([
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    T1 = np.array([
        [np.cos(2*np.pi*p/n), np.sin(2*np.pi*p/n), 0, 0],
        [-np.sin(2*np.pi*p/n), np.cos(2*np.pi*p/n), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    T2 = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, np.cos(2*np.pi*p/n), -np.sin(2*np.pi*p/n)],
        [0, 0, np.sin(2*np.pi*p/n), np.cos(2*np.pi*p/n)]
    ])
    
    print("Testing 4D Representation example:", p4mRepTest(C4, M, T1, T2))
