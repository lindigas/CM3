import math
import numpy as np

def omega(h1, h2, l1, l2):
    lambdamax = 2. / (h1 * h1 + h2 * h2) * (h2 * h2 * math.pow(math.sin(math.pi * h1 / (2. * l1)), 2.0)+ h1 * h1 * math.pow(math.sin(math.pi * h2 / (2. * l2)), 2.0)) 
    return 2. / (1 + math.sqrt(lambdamax * (2 - lambdamax)))

def TopRelaxation(N, M, a, b, c, d, X, y, V, eps, Nmax):
    stepError = 0
    stepNumber = 1
    stop = False

    h = (b - a) / N
    k = (d - c) / M
    w = omega(h, k, b-a, d-c)
    while(stop == False):
        stepError = 0; 
        H = 1.0 / (h * h)
        K = 1.0 / (k * k)
        D = 1.0 / (2.0 * (H + K))
        for  j in range(1,M):
            for i in range(1, N):
                old = V[i][j] 
                V[i][j] = (1 - w) * old + D * w * (H * V[i - 1][j] + K * V[i][j - 1] + H * V[i + 1][j] + K * V[i][j + 1] + f(X[i], y[j]))
                Error = math.fabs(V[i, j] - old)
                if (Error > stepError): stepError = Error
        stepNumber += 1 
        if ((stepNumber > Nmax) or (stepError <= eps)): stop = True
        S = stepNumber - 1 
        E = stepError
    return S, E, V 

def discrepancy(V, N, M, X, y, h, k):
    H = 1.0/(h*h)
    K = 1.0/(k*k)
    max_nev = 0
    for j in range(1, M):
        for i in range(1, N):
            nev = math.fabs(-f(X[i],y[j]) - (H*(V[i-1][j]-2*V[i][j]+V[i+1][j]) + K*(V[i][j-1]-2*V[i][j]+V[i][j+1])))
            if nev > max_nev:
                max_nev = nev
    return max_nev


def zeidel(n, m, a, b, c, d, x, y, v, eps, nmax):
    s = 0
    eps_max = 0
    max_nev = 0
    h2 = - np.float64((n * n)/((b - a) * (b - a)))
    k2 = - np.float64((m * m) / ((d - c) * (d - c)))
    a2 = -2.0 * (h2 + k2)
    flag = 0
    while flag != 1:
        eps_max = 0
        for j in range(1, m):
            for i in range(1, n):
                v_old = v[i, j]
                v_new = -(h2 * (v[i + 1, j] + v[i - 1, j]) + k2 * (v[i, j + 1] + v[i, j - 1]))
                v_new = v_new + f(x[i], y[j])
                v_new = v_new / a2
                eps_cur = math.fabs(v_old - v_new)
                if eps_cur > eps_max:
                    eps_max = eps_cur
                v[i, j] = v_new
        s = s + 1
        if (eps_max < eps) or (s >= nmax):
            flag = 1
    ss = s
    ee = eps_max
    return ss, ee, v


def grid(a, b, c, d, n, m):
    h = (b - a) / n
    k = (d - c) / m

    x = [0] * (n + 1)
    y = [0] * (m + 1)

    x[0] = a
    x[n] = b
    y[0] = c
    y[m] = d

    for i in range(1, n):
        x[i] = a + i * h
    for j in range(1, m):
        y[j] = c + j * k

    q = [[0.] * (m + 1)] * (n + 1)
    v = np.array(q, dtype=np.float64)

    for i in range(0, n + 1):
        v[i][0] = mu3(x[i], c)
        v[i][m] = mu4(x[i], d)
    for j in range(0, m + 1):
        v[0][j] = mu1(y[j], a)
        v[n][j] = mu2(y[j], b)

    return x, y, v


def grid_2(a, b, c, d, n, m):
    h = (b - a) / n
    k = (d - c) / m

    x = [0] * (n + 1)
    y = [0] * (m + 1)

    x[0] = a
    x[n] = b
    y[0] = c
    y[m] = d

    for i in range(1, n):
        x[i] = a + i * h
    for j in range(1, m):
        y[j] = c + j * k

    q = [[0.] * (m + 1)] * (n + 1)
    v = np.array(q, dtype=np.float64)

    for i in range(0, math.ceil(n / 2)):
        v[i][0] = u(x[i], c)
        v[i][m] = u(x[i], d)
    for i in range(math.ceil(n / 2), n + 1):
        v[i][0] = u(x[i], c)
        v[i][np.int(m / 2)] = u(x[i], (c + d) / 2.0)

    for j in range(0, math.ceil(m / 2)):
        v[0][j] = u(a, y[j])
        v[n][j] = u(b, y[j])
    for j in range(math.ceil(m / 2) + 1, m + 1):
        v[0][j] = u(a, y[j])
        v[np.int(n / 2)][j] = u((a + b) / 2.0, y[j])
    return x, y, v


def f(x, y):
    return 4


def u(x, y):
    return 1-(x-1)*(x-1) - (y-0.5)*(y-0.5)


def mu1(y, a):
    return 1-(a-1)*(a-1) - (y-0.5)*(y-0.5)


def mu2(y, b):
    return 1-(b-1)*(b-1) - (y-0.5)*(y-0.5)


def mu3(x, c):
    return 1-(x-1)*(x-1) - (c-0.5)*(c-0.5)
	

def mu4(x, d):
    return 1-(x-1)*(x-1) - (d-0.5)*(d-0.5)


def accuracy(n, m, x, y, v):
    max = 0.0
    MaxX = 0.0
    MaxY = 0.0
    for i in range(0, n + 1):
        for j in range(0, m + 1):
            tmp = math.fabs(u(x[i], y[j]) - v[i][j])
            if tmp > max:
                max = tmp
                MaxX = x[i]
                MaxY = y[j]
    Max = max
    return Max, MaxX, MaxY


def accuracy_2(n, m, x, y, v):
    max = 0.0
    max1 = 0.0
    max2 = 0.0
    MaxX1 = 0.0
    MaxX2 = 0.0
    MaxY1 = 0.0
    MaxY2 = 0.0

    for i in range(0, math.ceil(n / 2)):
        for j in range(0, m + 1):
            tmp = math.fabs(u(x[i], y[j]) - v[i][j])
            if tmp > max1:
                max1 = tmp
                MaxX1 = x[i]
                MaxY1 = y[j]
    for i in range(math.ceil(m / 2), n + 1):
        for j in range(0, math.ceil(m / 2)):
            if tmp > max2:
                max2 = tmp
                MaxX2 = x[i]
                MaxY2 = y[j]
    if max1 > max2:
        Max = max1
        MaxX = MaxX1
        MaxY = MaxY1
    else:
        Max = max2
        MaxX = MaxX2
        MaxY = MaxY2
    return Max, MaxX, MaxY
