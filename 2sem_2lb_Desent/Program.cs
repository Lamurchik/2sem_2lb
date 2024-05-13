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

static decimal DotProduct(decimal[] a, decimal[] b)
{
    decimal result = 0;
    for (int i = 0; i < a.Length; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

// Метод для умножения матрицы на вектор
static decimal[] Multiply(decimal[,] matrix, decimal[] vector)
{
    int n = matrix.GetLength(0);
    int m = matrix.GetLength(1);
    if (vector.Length != m)
    {
        throw new ArgumentException("Matrix and vector dimensions don't match.");
    }

    decimal[] result = new decimal[n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            result[i] += matrix[i, j] * vector[j];
        }
    }
    return result;
}

 decimal[] ConjugateGradient(decimal[,] A, decimal[] b, int maxIterations, decimal tolerance)
{
    int n = b.Length;
    decimal[] x = new decimal[n]; // Начальное приближение
    decimal[] r = new decimal[n]; // Остаток
    decimal[] p = new decimal[n]; // Направление
    decimal[] Ap = new decimal[n]; // Результат умножения матрицы A на вектор p

    // Инициализация начального приближения и остатка
    Array.Copy(b, r, n);
    Array.Copy(r, p, n);

    decimal rDotR = DotProduct(r, r); // Скалярное произведение r на r
    decimal initialRDotR = rDotR;
    // Итерационный процесс
    for (int iteration = 0; iteration < maxIterations && rDotR > tolerance * tolerance * initialRDotR; iteration++)
    {
        // Вычисление Ap
        Ap = Multiply(A, p);
        // Вычисление параметра alpha
        decimal alpha = rDotR / DotProduct(p, Ap);
        // Обновление x и r
        for (int i = 0; i < n; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        // Вычисление нового значения rDotR
        decimal newRDotR = DotProduct(r, r);
        // Вычисление параметра beta
        decimal beta = newRDotR / rDotR;
        // Обновление направления p
        for (int i = 0; i < n; i++)
        {
            p[i] = r[i] + beta * p[i];
        }
        rDotR = newRDotR; // Обновление значения rDotR для следующей итерации
    }
    return x; // Возвращаем найденное решение
}
static void ShowArr(decimal[] a)
{
	foreach (decimal c in a)
	{
		Console.WriteLine(c);
	}
}


// A - матрица коэффициентов системы уравнений
Decimal[,] Ap = new decimal[n, n];
for (int i = 0; i < n; i++)
{
	for (int j = 0; j < n; j++)
	{ Ap[i, j] = 0; }
}
Ap[0, 0] = 1;
Ap[n - 1, n - 1] = 1;


for (int i = 1; i < n - 1; i++)
{
	//кладем занчения в диоганаль
	Ap[i, i] = B(i);
	//в поддиоганаль
	Ap[i, i - 1] = A(i);
	//в наддиоганаль
	Ap[i, i + 1] = C(i);
}
// b - вектор свободных членов системы уравнений
decimal[] b = new decimal[n];
b[0] = 0;
b[n - 1] = 0;
for (int i = 1; i < n - 1; i++)
{
	b[i] = Fx(i * h);
}

var des = ConjugateGradient(Ap, b, 10000, e);

ShowArr(des);