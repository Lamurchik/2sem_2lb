const int n = 10;
decimal h = 1m / n;
const decimal e = 0.0000001m;
#region итерационый метод

decimal Fx(decimal x)
{
	decimal fx = (-14 * x * x * x * x + 15 * x * x * x - 4 * x * x - 5 * x + 4) * h * h;
	return fx;
}
decimal U(decimal x)
{
    decimal u = x * (1 - x) * (1 - x);
    return u;
}
decimal G(decimal x)
{
    return 1 + x;
}
decimal P(decimal x)
{
    return 1 + x * x * x;
}
decimal ai(int i)
{
    return P(i * h);
}

decimal A(int i)//ai
{
    return -ai(i);
}
decimal B(int i)
{
    return -ai(i + 1) - ai(i) - G(i * h) * h * h;
}
decimal C(int i)
{
    return -ai(i + 1);
}

decimal Alfaip1(int i)
{

    if (i < 0) Console.WriteLine("накасячмил в альфа");
    if (i == 0)
        return 0;
    decimal alfa = C(i) / (B(i) - A(i) * Alfaip1(i - 1));
    return alfa;
    /*
    int I = i + 1;
    if(I==1) return 0;
    decimal alfaIp1 = C(i) / (B(i) - A(i) * Alfaip1(i));
    return alfaIp1;
    */
}


void alfa1(decimal[] A, decimal[] B, decimal[] C, ref decimal[] alfa)
{
    alfa[1] = 0;
    for (int j = 1; j < n; j++)
    {
        alfa[j + 1] = C[j] / (B[j] - A[j] * alfa[j]);
    }
}

void beta1(decimal[] A, decimal[] B, decimal[] C, decimal[] alfa, ref decimal[] beta)
{

    beta[1] = 0;
    for (int j = 1; j < n; j++)
    {
        beta[j + 1] = (A[j] * beta[j] - Fx(j * h)) / (B[j] - A[j] * alfa[j]);
    }
}

decimal[] Yi1(decimal[] A, decimal[] B, decimal[] C)
{
    decimal[] alfa = new decimal[n + 1];
    decimal[] beta = new decimal[n + 1];
    alfa1(A, B, C, ref alfa);
    beta1(A, B, C, alfa, ref beta);

    decimal[] yi = new decimal[n + 1];
    yi[0] = 0;
    yi[n] = 0;

    for (int i = n - 1; i >= 1; i--)
    {
        yi[i] = alfa[i + 1] * yi[i + 1] + beta[i + 1];          //yi[i + 1] * (C[i]) / (B[i] - A[i] * alfa[i]) + (A[i] * beta[i]-Fx(i*h)) / (B[i] - A[i] * alfa[i]);
    }
    return yi;
}

decimal Betaip1(int i)
{
    if (i < 0)
    {
        Console.WriteLine("накасячмил в beta");
    }
    if (i == 0)
        return 0;
    return (A(i) * Betaip1(i - 1) - Fx(i * h)) / (B(i) - A(i) * Alfaip1(i - 1));
}

decimal Yi(int i)
{
    if (i == 0 || i == n)
        return 0;

    decimal y = Yi(i + 1) * C(i) / (-A(i) * Alfaip1(i - 1) + B(i)) + (A(i) * Betaip1(i - 1) - Fx(i * h)) / (B(i) - A(i) * Alfaip1(i - 1));
    return y;
}
#endregion


static decimal EuclideanNorm(decimal[] vector)
{
    decimal sum = 0;
    foreach (var element in vector)
    {
        sum += element * element;
    }
    return (decimal)Math.Sqrt((double)sum);
}

 decimal[] Solve( decimal tolerance = 1e-10m)
{
  
    decimal[] x = new decimal[n];
    decimal[] xPrev = new decimal[n];


    bool t = true;
    var count = 0;
    while(t)
    {
        for (int i = 1; i < n-1; i++)
        {
            x[i] = (ai(i) * x[i - 1] / (-B(i))) + ai(i + 1) * xPrev[i + 1] / (-B(i)) + Fx(i)/ (-B(i));
        }
        decimal normDiff = EuclideanNorm(xPrev.Zip(x, (prev, cur) => cur - prev).ToArray());
        if (normDiff < tolerance)
            t=false;
        Array.Copy(x, xPrev, n);
        count++;
    }
    return x;
}
static void ShowArr(decimal[] a)
{
    foreach (decimal c in a)
    {
        Console.WriteLine(c);
    }
}

ShowArr(Solve());
