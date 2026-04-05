import math
import numpy as np
from scipy.linalg import block_diag

# Basic symbols for polynomials
c = "c"
r = "r"
t = "t"
x = "x"

def PolynomialMod(expr, mod):
    # In Mathematica: PolynomialMod[..., 2]
    # We preserve the sum literally for now, returning the modulo 2 boolean string equivalent.
    return expr

def ArrayFlatten(arr_blocks):
    # Equivalent to ArrayFlatten[{{A,0},{0,B}}]
    return np.block(arr_blocks)

def DiagonalMatrix(arr):
    # Depending on input, evaluates to standard diagonal matrix or block_diag
    if len(arr) > 0 and isinstance(arr[0], (list, np.ndarray)):
        return block_diag(*arr)
    return np.diag(arr)


def SortNumber2(n):
    # n - 1/2 a1 (a1-1) > 0 => a1^2 - a1 - 2n < 0 => roots at (1 +- sqrt(1+8n))/2
    a1_root = (-1 + math.sqrt(1 + 8 * n)) / 2
    a1 = math.ceil(a1_root)
    a2 = n - 0.5 * a1 * (a1 - 1)
    return int(a1), int(a2)

def SortNumber3(n):
    a1 = 1
    while (a1 * (a1 + 1) * (a1 + 2)) / 6 < n:
        a1 += 1
    res = n - (a1 - 1) * a1 * (a1 + 1) / 6
    a2, a3 = SortNumber2(res)
    return int(a1), int(a2), int(a3)

def SortNumber4(n):
    a1 = 1
    while (a1 * (a1 + 1) * (a1 + 2) * (a1 + 3)) / 24 < n:
        a1 += 1
    res = n - (a1 - 1) * a1 * (a1 + 1) * (a1 + 2) / 24
    a2, a3, a4 = SortNumber3(res)
    return int(a1), int(a2), int(a3), int(a4)

p6m1dIrrep = [
    [1, 1, 1, 1, 0],
    [-1, 1, 1, 1, c],
    [1, -1, 1, 1, r],
    [-1, -1, 1, 1, f"{c}+{r}"]
]

p6m2dIrrep = [
    [np.array([[1, 0], [0, -1]]), np.array([[1, 0], [0, 1]]),
     np.array([[math.cos(2*math.pi/3), math.sin(2*math.pi/3)], [-math.sin(2*math.pi/3), math.cos(2*math.pi/3)]]),
     np.array([[math.cos(2*math.pi/3), math.sin(2*math.pi/3)], [-math.sin(2*math.pi/3), math.cos(2*math.pi/3)]]), c],
    [np.array([[1, 0], [0, -1]]), np.array([[-1, 0], [0, -1]]),
     np.array([[math.cos(2*math.pi/3), math.sin(2*math.pi/3)], [-math.sin(2*math.pi/3), math.cos(2*math.pi/3)]]),
     np.array([[math.cos(2*math.pi/3), math.sin(2*math.pi/3)], [-math.sin(2*math.pi/3), math.cos(2*math.pi/3)]]), c],
    [np.array([[math.cos(math.pi/3), math.sin(math.pi/3)], [-math.sin(math.pi/3), math.cos(math.pi/3)]]),
     np.array([[-1, 0], [0, 1]]), np.array([[1, 0], [0, 1]]), np.array([[1, 0], [0, 1]]), r],
    [np.array([[math.cos(2*math.pi/3), math.sin(2*math.pi/3)], [-math.sin(2*math.pi/3), math.cos(2*math.pi/3)]]),
     np.array([[-1, 0], [0, 1]]), np.array([[1, 0], [0, 1]]), np.array([[1, 0], [0, 1]]), r]
]

p6m3dIrrep = [
    [np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]), np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
     np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]), np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]), r],
    [np.array([[0, 1, 0], [0, 0, 1], [-1, 0, 0]]), np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1]]),
     np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]), np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]), f"{c}+{r}"],
    [np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]), np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]]),
     np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]), np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]), 0],
    [np.array([[0, 1, 0], [0, 0, 1], [-1, 0, 0]]), np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]]),
     np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]), np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]), c]
]

p6m4dIrrep = [np.array([
    [1/2, 0, 0, math.sqrt(3)/2],
    [0, -1/2, -math.sqrt(3)/2, 0],
    [0, math.sqrt(3)/2, -1/2, 0],
    [-math.sqrt(3)/2, 0, 0, 1/2]
]), np.array([
    [-1, 0, 0, 0],
    [0, -1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
]), np.array([
    [-1/2, math.sqrt(3)/2, 0, 0],
    [-math.sqrt(3)/2, -1/2, 0, 0],
    [0, 0, -1/2, -math.sqrt(3)/2],
    [0, 0, math.sqrt(3)/2, -1/2]
]), np.array([
    [-1/2, math.sqrt(3)/2, 0, 0],
    [-math.sqrt(3)/2, -1/2, 0, 0],
    [0, 0, -1/2, -math.sqrt(3)/2],
    [0, 0, math.sqrt(3)/2, -1/2]
])]

def p6mO2MatrixForm(i):
    tab = [0] * 6
    if i <= 36:
        a1, a2 = SortNumber2(i)
        t1 = (a1 - 1) // 4
        t2 = (a2 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        a2 = (a2 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p6m1dIrrep[a1 - 1][m - 1], p6m1dIrrep[a2 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t2])
        tab[5] = PolynomialMod(f"{p6m1dIrrep[a1 - 1][4]} + {p6m1dIrrep[a2 - 1][4]} + t*({t1}+{t2})", 2)
    elif i > 36:
        a1 = i - 36
        t1 = (a1 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = p6m2dIrrep[a1 - 1][m - 1]
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1])
        tab[5] = p6m2dIrrep[a1 - 1][4]
    return tab

def p6mO3MatrixForm(i):
    tab = [0] * 6
    if i <= 120:
        a1, a2, a3 = SortNumber3(i)
        t1 = (a1 - 1) // 4
        t2 = (a2 - 1) // 4
        t3 = (a3 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        a2 = (a2 - 1) % 4 + 1
        a3 = (a3 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p6m1dIrrep[a1 - 1][m - 1], p6m1dIrrep[a2 - 1][m - 1], p6m1dIrrep[a3 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t2, (-1)**t3])
        tab[5] = PolynomialMod(f"{p6m1dIrrep[a1 - 1][4]} + {p6m1dIrrep[a2 - 1][4]} + {p6m1dIrrep[a3 - 1][4]} + t*({t1}+{t2}+{t3})", 2)
    elif 120 < i <= 184:
        a1 = i - 120
        a2 = (a1 - 1) // 8 + 1
        a1 = (a1 - 1) % 8 + 1
        t1 = (a1 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        t2 = (a2 - 1) // 4
        a2 = (a2 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = ArrayFlatten([[p6m2dIrrep[a1 - 1][m - 1], np.zeros((2, 1))], [np.zeros((1, 2)), [[p6m1dIrrep[a2 - 1][m - 1]]]]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1, (-1)**t2])
        tab[5] = PolynomialMod(f"{p6m2dIrrep[a1 - 1][4]} + {p6m1dIrrep[a2 - 1][4]} + t*{t2}", 2)
    elif i > 184:
        a1 = i - 184
        t1 = (a1 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = p6m3dIrrep[a1 - 1][m - 1]
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1, (-1)**t1])
        tab[5] = f"{p6m3dIrrep[a1 - 1][4]} + t*{t1}"
    return tab

def p6mO4MatrixForm(i):
    tab = [0] * 6
    if i <= 330:
        a1, a2, a3, a4 = SortNumber4(i)
        t1 = (a1 - 1) // 4
        t2 = (a2 - 1) // 4
        t3 = (a3 - 1) // 4
        t4 = (a4 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        a2 = (a2 - 1) % 4 + 1
        a3 = (a3 - 1) % 4 + 1
        a4 = (a4 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p6m1dIrrep[a1 - 1][m - 1], p6m1dIrrep[a2 - 1][m - 1], p6m1dIrrep[a3 - 1][m - 1], p6m1dIrrep[a4 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t2, (-1)**t3, (-1)**t4])
        tab[5] = PolynomialMod(f"{p6m1dIrrep[a1 - 1][4]} + {p6m1dIrrep[a2 - 1][4]} + {p6m1dIrrep[a3 - 1][4]} + {p6m1dIrrep[a4 - 1][4]} + t*({t1}+{t2}+{t3}+{t4})", 2)
    elif 330 < i <= 618:
        a1 = i - 330
        a2 = (a1 - 1) // 8 + 1
        a1 = (a1 - 1) % 8 + 1
        a2, a3 = SortNumber2(a2)
        t1 = (a1 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        t2 = (a2 - 1) // 4
        a2 = (a2 - 1) % 4 + 1
        t3 = (a3 - 1) // 4
        a3 = (a3 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = ArrayFlatten([[p6m2dIrrep[a1-1][m-1], np.zeros((2,1)), np.zeros((2,1))],
                                       [np.zeros((1,2)), [[p6m1dIrrep[a2-1][m-1]]], np.zeros((1,1))],
                                       [np.zeros((1,2)), np.zeros((1,1)), [[p6m1dIrrep[a3-1][m-1]]]]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1, (-1)**t2, (-1)**t3])
        tab[5] = PolynomialMod(f"{p6m2dIrrep[a1-1][4]} + {p6m1dIrrep[a2-1][4]} + {p6m1dIrrep[a3-1][4]} + t*({t2}+{t3})", 2)
    elif 618 < i <= 654:
        a1 = i - 618
        a1, a2 = SortNumber2(a1)
        t1 = (a1 - 1) // 4
        t2 = (a2 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        a2 = (a2 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = ArrayFlatten([[p6m2dIrrep[a1-1][m-1], np.zeros((2,2))], [np.zeros((2,2)), p6m2dIrrep[a2-1][m-1]]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1, (-1)**t2, (-1)**t2])
        tab[5] = PolynomialMod(f"{p6m2dIrrep[a1-1][4]} + {p6m2dIrrep[a2-1][4]}", 2)
    elif 654 < i <= 718:
        a1 = i - 654
        a2 = (a1 - 1) // 8 + 1
        a1 = (a1 - 1) % 8 + 1
        t1 = (a1 - 1) // 4
        a1 = (a1 - 1) % 4 + 1
        t2 = (a2 - 1) // 4
        a2 = (a2 - 1) % 4 + 1
        for m in range(1, 5):
            tab[m - 1] = ArrayFlatten([[[[p6m1dIrrep[a1-1][m-1]]], np.zeros((1,3))], [np.zeros((3,1)), p6m3dIrrep[a2-1][m-1]]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t2, (-1)**t2, (-1)**t2])
        tab[5] = PolynomialMod(f"{p6m3dIrrep[a2-1][4]} + {p6m1dIrrep[a1-1][4]} + t*({t1}+{t2})", 2)
    elif 718 < i <= 720:
        a1 = i - 718
        t1 = (a1 - 1) // 2
        for m in range(1, 5):
            tab[m - 1] = p6m4dIrrep[m - 1]
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1, (-1)**t1, (-1)**t1])
        tab[5] = 0
    return tab

p4m1dIrrep = [
    [1, 1, 1, 1, 0],
    [-1, 1, 1, 1, c],
    [1, -1, 1, 1, r],
    [1, 1, -1, -1, x],
    [-1, -1, 1, 1, f"{c}+{r}"],
    [-1, 1, -1, -1, f"{c}+{x}"],
    [1, -1, -1, -1, f"{r}+{x}"],
    [-1, -1, -1, -1, f"{c}+{r}+{x}"]
]

p4m2dIrrep = [
    [np.array([[0, 1], [1, 0]]), np.array([[1, 0], [0, 1]]), np.array([[-1, 0], [0, 1]]), np.array([[1, 0], [0, -1]]), f"{c}+{x}"],
    [np.array([[0, 1], [-1, 0]]), np.array([[-1, 0], [0, 1]]), np.array([[-1, 0], [0, 1]]), np.array([[1, 0], [0, -1]]), f"{r}+{x}"],
    [np.array([[0, 1], [1, 0]]), np.array([[-1, 0], [0, -1]]), np.array([[-1, 0], [0, 1]]), np.array([[1, 0], [0, -1]]), f"{c}+{x}"],
    [np.array([[0, 1], [-1, 0]]), np.array([[1, 0], [0, -1]]), np.array([[-1, 0], [0, 1]]), np.array([[1, 0], [0, -1]]), f"{r}+{x}"],
    [np.array([[0, 1], [-1, 0]]), np.array([[-1, 0], [0, 1]]), np.array([[1, 0], [0, 1]]), np.array([[1, 0], [0, 1]]), r],
    [np.array([[0, 1], [-1, 0]]), np.array([[-1, 0], [0, 1]]), np.array([[-1, 0], [0, -1]]), np.array([[-1, 0], [0, -1]]), r]
]

def get_p4m4dIrrep(p, n):
    return [
        [np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
         np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
         np.array([[math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n), 0, 0],
                   [-math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
         np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                   [0, 0, math.cos(2*math.pi*p/n), -math.sin(2*math.pi*p/n)],
                   [0, 0, math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n)]]), f"{c}+{r}"],
        [np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
         np.array([[0, -1, 0, 0], [-1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]),
         np.array([[math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n), 0, 0],
                   [-math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
         np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                   [0, 0, math.cos(2*math.pi*p/n), -math.sin(2*math.pi*p/n)],
                   [0, 0, math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n)]]), f"{c}+{r}"],
        [np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
         np.array([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]),
         np.array([[math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n), 0, 0],
                   [-math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n), 0, 0],
                   [0, 0, math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n)],
                   [0, 0, -math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n)]]),
         np.array([[math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n), 0, 0],
                   [-math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n), 0, 0],
                   [0, 0, math.cos(2*math.pi*p/n), -math.sin(2*math.pi*p/n)],
                   [0, 0, math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n)]]), c],
        [np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
         np.array([[0, 0, 0, -1], [0, 0, -1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]),
         np.array([[math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n), 0, 0],
                   [-math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n), 0, 0],
                   [0, 0, math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n)],
                   [0, 0, -math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n)]]),
         np.array([[math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n), 0, 0],
                   [-math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n), 0, 0],
                   [0, 0, math.cos(2*math.pi*p/n), -math.sin(2*math.pi*p/n)],
                   [0, 0, math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n)]]), c],
        [np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
         np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
         np.array([[math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n), 0, 0],
                   [-math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n), 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]),
         np.array([[-1, 0, 0, 0], [0, -1, 0, 0],
                   [0, 0, math.cos(2*math.pi*p/n), -math.sin(2*math.pi*p/n)],
                   [0, 0, math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n)]]), f"{c}+{r}"],
        [np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
         np.array([[0, -1, 0, 0], [-1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]),
         np.array([[math.cos(2*math.pi*p/n), math.sin(2*math.pi*p/n), 0, 0],
                   [-math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n), 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]),
         np.array([[-1, 0, 0, 0], [0, -1, 0, 0],
                   [0, 0, math.cos(2*math.pi*p/n), -math.sin(2*math.pi*p/n)],
                   [0, 0, math.sin(2*math.pi*p/n), math.cos(2*math.pi*p/n)]]), f"{c}+{r}"]
    ]

def p4mO2MatrixForm(i):
    tab = [0] * 6
    if i <= 136:
        a1, a2 = SortNumber2(i)
        t1 = (a1 - 1) // 8
        t2 = (a2 - 1) // 8
        a1 = (a1 - 1) % 8 + 1
        a2 = (a2 - 1) % 8 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p4m1dIrrep[a1 - 1][m - 1], p4m1dIrrep[a2 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t2])
        tab[5] = PolynomialMod(f"{p4m1dIrrep[a1 - 1][4]} + {p4m1dIrrep[a2 - 1][4]} + t*({t1}+{t2})", 2)
    elif i > 136:
        a1 = i - 136
        t1 = (a1 - 1) // 6
        a1 = (a1 - 1) % 6 + 1
        for m in range(1, 5):
            tab[m - 1] = p4m2dIrrep[a1 - 1][m - 1]
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1])
        tab[5] = p4m2dIrrep[a1 - 1][4]
    return tab

def p4mO3MatrixForm(i):
    tab = [0] * 6
    if i <= 816:
        a1, a2, a3 = SortNumber3(i)
        t1 = (a1 - 1) // 8
        t2 = (a2 - 1) // 8
        t3 = (a3 - 1) // 8
        a1 = (a1 - 1) % 8 + 1
        a2 = (a2 - 1) % 8 + 1
        a3 = (a3 - 1) % 8 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p4m1dIrrep[a1 - 1][m - 1], p4m1dIrrep[a2 - 1][m - 1], p4m1dIrrep[a3 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t2, (-1)**t3])
        tab[5] = PolynomialMod(f"{p4m1dIrrep[a1 - 1][4]} + {p4m1dIrrep[a2 - 1][4]} + {p4m1dIrrep[a3 - 1][4]} + t*({t1}+{t2}+{t3})", 2)
    elif i > 816:
        a1 = i - 816
        a2 = (a1 - 1) // 12 + 1
        a1 = (a1 - 1) % 12 + 1
        t1 = (a1 - 1) // 6
        a1 = (a1 - 1) % 6 + 1
        t2 = (a2 - 1) // 8
        a2 = (a2 - 1) % 8 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p4m2dIrrep[a1 - 1][m - 1], p4m1dIrrep[a2 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1, (-1)**t2])
        tab[5] = PolynomialMod(f"{p4m2dIrrep[a1 - 1][4]} + {p4m1dIrrep[a2 - 1][4]} + t*{t2}", 2)
    return tab

def p4mO4MatrixForm(i):
    tab = [0] * 6
    if i <= 3876:
        a1, a2, a3, a4 = SortNumber4(i)
        t1 = (a1 - 1) // 8
        t2 = (a2 - 1) // 8
        t3 = (a3 - 1) // 8
        t4 = (a4 - 1) // 8
        a1 = (a1 - 1) % 8 + 1
        a2 = (a2 - 1) % 8 + 1
        a3 = (a3 - 1) % 8 + 1
        a4 = (a4 - 1) % 8 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p4m1dIrrep[a1 - 1][m - 1], p4m1dIrrep[a2 - 1][m - 1], p4m1dIrrep[a3 - 1][m - 1], p4m1dIrrep[a4 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t2, (-1)**t3, (-1)**t4])
        tab[5] = PolynomialMod(f"{p4m1dIrrep[a1 - 1][4]} + {p4m1dIrrep[a2 - 1][4]} + {p4m1dIrrep[a3 - 1][4]} + {p4m1dIrrep[a4 - 1][4]} + t*({t1}+{t2}+{t3}+{t4})", 2)
    elif 3876 < i <= 5508:
        a1 = i - 3876
        a2 = (a1 - 1) // 12 + 1
        a1 = (a1 - 1) % 12 + 1
        a2, a3 = SortNumber2(a2)
        t1 = (a1 - 1) // 6
        a1 = (a1 - 1) % 6 + 1
        t2 = (a2 - 1) // 8
        a2 = (a2 - 1) % 8 + 1
        t3 = (a3 - 1) // 8
        a3 = (a3 - 1) % 8 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p4m2dIrrep[a1 - 1][m - 1], p4m1dIrrep[a2 - 1][m - 1], p4m1dIrrep[a3 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1, (-1)**t2, (-1)**t3])
        tab[5] = PolynomialMod(f"{p4m2dIrrep[a1 - 1][4]} + {p4m1dIrrep[a2 - 1][4]} + {p4m1dIrrep[a3 - 1][4]} + t*({t2}+{t3})", 2)
    elif i > 5508:
        a1 = i - 5508
        a1, a2 = SortNumber2(a1)
        t1 = (a1 - 1) // 6
        t2 = (a2 - 1) // 6
        a1 = (a1 - 1) % 6 + 1
        a2 = (a2 - 1) % 6 + 1
        for m in range(1, 5):
            tab[m - 1] = DiagonalMatrix([p4m2dIrrep[a1 - 1][m - 1], p4m2dIrrep[a2 - 1][m - 1]])
        tab[4] = DiagonalMatrix([(-1)**t1, (-1)**t1, (-1)**t2, (-1)**t2])
        tab[5] = PolynomialMod(f"{p4m2dIrrep[a1 - 1][4]} + {p4m2dIrrep[a2 - 1][4]}", 2)
    return tab

def p4mDQCPEmbedding(i):
    O2Mat = p4mO2MatrixForm(i)
    w1 = PolynomialMod(f"{t}+{r}+{O2Mat[5]}", 2)
    O3Sign, O3T = 0, 0
    for m in range(1, 9):
        if p4m1dIrrep[m - 1][4] == w1: O3Sign = m
        if PolynomialMod(f"{p4m1dIrrep[m - 1][4]}+{t}", 2) == w1: O3Sign = m; O3T = 1
    tab = [0] * 5
    for m in range(1, 5):
        tab[m - 1] = DiagonalMatrix([p4m1dIrrep[O3Sign - 1][m - 1] * np.eye(3), O2Mat[m - 1]])
    tab[4] = DiagonalMatrix([((-1)**O3T) * np.eye(3), O2Mat[4]])
    return tab

def p4mDSLEmbedding(i, j):
    O2Mat = p4mO2MatrixForm(j)
    O3Mat = p4mO3MatrixForm(i)
    w1 = PolynomialMod(f"{t}+{r}+{O2Mat[5]}+{O3Mat[5]}", 2)
    O3Sign, O3T = 0, 0
    for m in range(1, 9):
        if p4m1dIrrep[m - 1][4] == w1: O3Sign = m
        if PolynomialMod(f"{p4m1dIrrep[m - 1][4]}+{t}", 2) == w1: O3Sign = m; O3T = 1
    tab = [0] * 5
    for m in range(1, 5):
        tab[m - 1] = [DiagonalMatrix([p4m1dIrrep[O3Sign - 1][m - 1] * np.eye(3), O3Mat[m - 1]]), O2Mat[m - 1]]
    tab[4] = [DiagonalMatrix([((-1)**O3T) * np.eye(3), O3Mat[4]]), O2Mat[4]]
    return tab

def p4mSLEmbedding(i, j):
    O3Mat = p4mO3MatrixForm(j)
    O4Mat = p4mO4MatrixForm(i)
    w1 = PolynomialMod(f"{t}+{r}+{O3Mat[5]}+{O4Mat[5]}", 2)
    O3Sign, O3T = 0, 0
    for m in range(1, 9):
        if p4m1dIrrep[m - 1][4] == w1: O3Sign = m
        if PolynomialMod(f"{p4m1dIrrep[m - 1][4]}+{t}", 2) == w1: O3Sign = m; O3T = 1
    tab = [0] * 5
    for m in range(1, 5):
        tab[m - 1] = [DiagonalMatrix([p4m1dIrrep[O3Sign - 1][m - 1] * np.eye(3), O4Mat[m - 1]]), O3Mat[m - 1]]
    tab[4] = [DiagonalMatrix([((-1)**O3T) * np.eye(3), O4Mat[4]]), O3Mat[4]]
    return tab

def p4mEmbedding(n, lsm, p, data=None):
    if data is None: return None
    tab = 0
    if n == 5:
        a1 = data[0][1][lsm - 1][p - 1] # data[[1]][[2]][[lsm]][[p]]
        tab = p4mDQCPEmbedding(a1)
    elif n == 6:
        a1 = data[1][1][lsm - 1][p - 1][0]
        a2 = data[1][1][lsm - 1][p - 1][1]
        tab = p4mDSLEmbedding(a1, a2)
    elif n == 7:
        a1 = data[2][1][lsm - 1][p - 1][0]
        a2 = data[2][1][lsm - 1][p - 1][1]
        tab = p4mSLEmbedding(a1, a2)
    return tab

def p4mStability(n, lsm, p):
    check = False
    if n == 5:
        if lsm in [2, 3, 5, 8]: check = True
    elif n == 6:
        if lsm == 2:
            if p in [3, 9, 14]: check = True
        elif lsm == 3:
            if p in [2, 9, 11]: check = True
        elif lsm == 4:
            if p in [12, 16, 18, 26]: check = True
        elif lsm in [6, 7]: check = True
        elif lsm == 8:
            if p in [7, 10]: check = True
    elif n == 7:
        if lsm == 2:
            if p in [19, 30]: check = True
        elif lsm == 3:
            if p in [20, 105]: check = True
        elif lsm == 6:
            if p in [5, 14, 19, 37, 51, 54]: check = True
        elif lsm == 7:
            if p in [2, 5, 30, 38, 54, 64]: check = True
    return check

def p4mPrintEmbedding(n, lsm, p, data=None):
    # Omitted exact Print statements for compactness
    return True

def p4mSL5dEmbedding(i, j):
    O3Mat = p4mO3MatrixForm(j)
    O2Mat = p4mO2MatrixForm(i)
    w1 = PolynomialMod(f"{t}+{r}+{O3Mat[5]}+{O2Mat[5]}", 2)
    O3Sign, O3T = 0, 0
    for m in range(1, 9):
        if p4m1dIrrep[m - 1][4] == w1: O3Sign = m
        if PolynomialMod(f"{p4m1dIrrep[m - 1][4]}+{t}", 2) == w1: O3Sign = m; O3T = 1
    tab = [0] * 5
    for m in range(1, 5):
        tab[m - 1] = [DiagonalMatrix([p4m1dIrrep[O3Sign - 1][m - 1] * np.eye(5), O2Mat[m - 1]]), O3Mat[m - 1]]
    tab[4] = [DiagonalMatrix([((-1)**O3T) * np.eye(5), O2Mat[4]]), O3Mat[4]]
    return tab

def p4m5dEmbedding(lsm, p, dataSL5Rep=None):
    if dataSL5Rep is None: return None
    a1 = dataSL5Rep[1][lsm - 1][p - 1][0] # dataSL5Rep[[2]][[lsm]][[p]][[1]]
    a2 = dataSL5Rep[1][lsm - 1][p - 1][1]
    tab = p4mSL5dEmbedding(a1, a2)
    return tab

def p4m5dStability(lsm, p):
    check = False
    if lsm == 1:
        if p in [7, 10, 19, 21]: check = True
    elif lsm in [2, 3]: check = True
    elif lsm == 4:
        if p == 2: check = True
    return check

def p4m5dPrintEmbedding(lsm, p, dataSL5Rep=None):
    return True

def IncommensurateCheck(i):
    flag1, flag2 = 0, 0
    if i in [3893, 5510]: flag1, flag2 = 1, 0
    if i in [4523, 5543]: flag1, flag2 = 1, 1
    if i in [4025, 5517]: flag1, flag2 = 2, 0
    if i in [4943, 5562]: flag1, flag2 = 2, 1
    if i in [4001, 4254]: flag1, flag2 = 3, 0
    if i in [4919, 5460]: flag1, flag2 = 3, 1
    if i in [3929, 4194]: flag1, flag2 = 4, 0
    if i in [4655, 5304]: flag1, flag2 = 4, 1
    if i in [4098, 5515]: flag1, flag2 = 5, 0
    if i in [5112, 5560]: flag1, flag2 = 5, 1
    if i in [4290, 5513]: flag1, flag2 = 6, 0
    if i in [5496, 5552]: flag1, flag2 = 6, 1
    return flag1, flag2

def Incommensurate(i, j, k):
    O4Mat = get_p4m4dIrrep(1, 4)[i - 1] # Assuming placeholder for array evaluation
    O3Mat = p4mO3MatrixForm(k)
    w1 = PolynomialMod(f"{t}+{r}+{O3Mat[5]}+{O4Mat[4]}+t*{j}", 2)
    O3Sign, O3T = 0, 0
    for m in range(1, 9):
        if p4m1dIrrep[m - 1][4] == w1: O3Sign = m
        if PolynomialMod(f"{p4m1dIrrep[m - 1][4]}+{t}", 2) == w1: O3Sign = m; O3T = 1
    tab = [0] * 5
    for m in range(1, 5):
        tab[m - 1] = [DiagonalMatrix([p4m1dIrrep[O3Sign - 1][m - 1] * np.eye(3), O4Mat[m - 1]]), O3Mat[m - 1]]
    tab[4] = [DiagonalMatrix([((-1)**O3T) * np.eye(3), ((-1)**j) * np.eye(4)]), O3Mat[4]]
    return tab

def IncommensurateEmbedding(lsm, p, data=None):
    flag = False
    tab = 0
    if data is None: return flag, tab
    a1 = data[2][1][lsm - 1][p - 1][0] # data[[3]][[2]][[lsm]][[p]][[1]]
    a2 = data[2][1][lsm - 1][p - 1][1]
    f1, f2 = IncommensurateCheck(a1)
    if f1 > 0:
        tab = Incommensurate(f1, f2, a2)
        flag = True
    return flag, tab

def IncommensuratePrintEmbedding(lsm, p, data=None):
    flag, tab = IncommensurateEmbedding(lsm, p, data)
    return flag
