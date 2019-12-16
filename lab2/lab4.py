def TopRelaxation(N, M, X, y, V):
    stepError = 0
    stepNumber = 1
    stop = false

    h = (b - a) / N
    k = (d - c) / M
    w = omega(h, k, b-a, d-c)
    while(stop == false):
        stepError = 0; 
        H = 1.0 / (h * h)
        K = 1.0 / (k * k)
        D = 1.0 / (2.0 * (H + K))
        for  j in range(1,M-1):
            for i in range(1, N-1):
                old = V[i][j] 
                V[i][j] = (1 - w) * old + D * w * (H * V[i - 1][j] + K * V[i][j - 1] + H * V[i + 1][j] + K * V[i][j + 1] + f(X[i], Y[j]))
                Error = Math.Abs(V[i, j] - old)
                if (Error > stepError): stepError = Error
        stepNumber+=1; 
        if ((stepNumber > Nmax) or (stepError <= eps)): stop = true
        S = stepNumber - 1 
        E = stepError
     


TopRelaxation(int N, int M, double[] X, double[] Y, double[,] V) 
{ 
double stepError = 0; //Точность на шаге stepNumber 
int stepNumber = 1; //Номер текущей итерации 

bool stop = false; 

h = (b - a) / N; 
k = (d - c) / M; 

w = omega(h, k, b-a, d-c); 

while (!stop) 
{ 
stepError = 0; 
double H = 1.0 / (h * h); 
double K = 1.0 / (k * k); 
double D = 1.0 / (2.0 * (H + K)); 
for (int j = 1; j <= M - 1; j++) 
for (int i = 1; i <= N - 1; i++) 
{ 
double old = V[i, j]; 
V[i, j] = (1 - w) * old + D * w * (H * V[i - 1, j] + K * V[i, j - 1] + H * V[i + 1, j] + K * V[i, j + 1] + f(X[i], Y[j])); 
double Error = Math.Abs(V[i, j] - old); 
if (Error > stepError) stepError = Error; 
} 
stepNumber++; 
if ((stepNumber > Nmax) || (stepError <= eps)) stop = true; 
} 
S = stepNumber - 1; 
E = stepError; 
}


def omega(h1, h2, l1, l2) 
    lambdamin = 2 / (h1 * h1 + h2 * h2) * (h2 * h2 * Math.Pow(Math.Sin(Math.PI * h1 / ((double)2 * l1)), 2.0) + h1 * h1 * Math.Pow(Math.Sin(Math.PI * h2 / ((double)2 * l2)), 2.0)) 
    return 2 / (1 + Math.Sqrt(lambdamin * (2 - lambdamin)))
