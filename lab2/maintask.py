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


def gridd(a, b, c, d, n, m):
    h = (b - a) / n
    k = (d - c) / m

    x = [0] * (n + 1)
    y = [0] * (m + 1)

    x2 = [0] * (2 * n + 1)
    y2 = [0] * (2 * m + 1)

    x[0] = a
    x[n] = b
    y[0] = c
    y[m] = d

    x2[0] = a
    x2[2 * n] = b
    y2[0] = c
    y2[2 * m] = d

    for i in range(1, n):
        x[i] = a + i * h
    for j in range(1, m):
        y[j] = c + j * k

    for i in range(1, 2 * n):
        x2[i] = a + i * h / 2
    for j in range(1, 2 * m):
        y2[j] = c + j * k / 2

    q = [[0.] * (m + 1)] * (n + 1)
    v = np.array(q, dtype=np.float64)

    q = [[0.] * (2 * m + 1)] * (2 * n + 1)
    v2 = np.array(q, dtype=np.float64)

    for i in range(0, n + 1):
        v[i][0] = mu3(x[i])
        v[i][m] = mu4(x[i])
    for j in range(0, m + 1):
        v[0][j] = mu1(y[j])
        v[n][j] = mu2(y[j])

    for i in range(0, 2 * n + 1):
        v2[i][0] = mu3(x2[i])
        v2[i][2 * m] = mu4(x2[i])
    for j in range(0, 2 * m + 1):
        v2[0][j] = mu1(y2[j])
        v2[2 * n][j] = mu2(y2[j])

    return x, y, v, x2, y2, v2


def f(x, y):
    return math.fabs(x - y)


def mu1(y):
    return -y * (y - 1)


def mu2(y):
    return -y * (y - 1)


def mu3(x):
    return math.fabs(math.sin(math.pi*x))


def mu4(x):
    return math.fabs(math.sin(math.pi*x)) * math.exp(x)


def accuracy(n, m, x, y, v, v2):
    max = 0.0
    MaxX = 0.0
    MaxY = 0.0
    for i in range(0, n + 1):
        for j in range(0, m + 1):
            tmp = math.fabs(v[i][j] - v2[2 * i][2 * j])
            if tmp > max:
                max = tmp
                MaxX = x[i]
                MaxY = y[j]
    Max = max
    return Max, MaxX, MaxY
