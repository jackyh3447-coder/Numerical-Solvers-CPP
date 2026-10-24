#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

vector<vector<double>> getMatrixQ();
vector<vector<double>> getMatrixR(double epsilon);
vector<vector<double>> constant_multiply(vector<vector<double>>& matrix, double constant);
vector<vector<double>> transpose(vector<vector<double>> matrix);
vector<vector<double>> subtractMatrices(
        vector<vector<double>> A,
        vector<vector<double>> B);
double colNorm(vector<vector<double>> A, int col);
double dotProduct(vector<vector<double>> A, int col1, int col2);
void printMatrix(const vector<vector<double>> &matrix); 
double FrobeniusNorm(vector<vector<double>> matrix);
vector<vector<double>> MGS(vector<vector<double>> &A);
vector<vector<double>> ErrorMatrix(vector<vector<double>> &Q);
vector<vector<double>> matrixMultiply(vector<vector<double>> matrix1, vector<vector<double>> matrix2);
vector<vector<double>> CGS(vector<vector<double>>& A);
vector<vector<double>> CGS2(vector<vector<double>> &A);
vector<vector<double>> Q;
vector<vector<double>> R;
vector<vector<double>> E;
vector<vector<double>> A;

double epsilon = 1e-5;

int main()
{
     Q = getMatrixQ();
     R = getMatrixR(epsilon);
     A = matrixMultiply(Q, R);

     // A = CGS(A);
     // A = MGS(A);
     A = CGS2(A);
     A = ErrorMatrix(A);

     printMatrix(A);
     cout << endl << "Frobenius: " << FrobeniusNorm(A) << " and Epsilon: " << epsilon << endl;
     return 0;
}

vector<vector<double>> getMatrixQ()
{
     vector<vector<double>> A = {
         {-9.0/2.0, sqrt(2), sqrt(3), sqrt(5)},
         {sqrt(2), -7.0/2.0, sqrt(6), sqrt(10)},
         {sqrt(3), sqrt(6), -5.0/2.0, sqrt(15)},
         {sqrt(5), sqrt(10), sqrt(15), -0.5}};

     return constant_multiply(A, 2.0/11.0);
}

vector<vector<double>> getMatrixR(double epsilon)
{
     vector<vector<double>> B = {
         {1.0, 1.0, 1.0, 1.0},
         {0.0, epsilon, 1.0, 1.0},
         {0.0, 0.0, epsilon, 1.0},
         {0.0, 0.0, 0.0, epsilon}};

     return B;
}

vector<vector<double>> CGS(vector<vector<double>>& A)
{
     int n = A.size();
     double columnNorm = colNorm(A, 0);
     vector<double> dotProducts(n);

     for (int col = 0; col < n; col++)
          {
               for (int i = 0; i < col; i++)
                    dotProducts[i] = dotProduct(A, i, col);

               for (int row = 0; row < n; row++)
                    {
                         double total_subtract = 0.0;
                         for (int j = 0; j < col; j++) 
                              total_subtract += dotProducts[j] * A[row][j];
                         A[row][col] -= total_subtract;
                    }
               columnNorm = colNorm(A, col);
               for (int row = 0; row < n; row++)
                    A[row][col] /= columnNorm;
          }
     return A;
}

vector<vector<double>> CGS2(vector<vector<double>> &A)
{
     int n = A.size();
     double columnNorm = colNorm(A, 0);
     vector<double> dotProducts(n);

     for (int col = 0; col < n; col++)
     {
          for (int pass = 0; pass < 2; pass++)
          {
               for (int i = 0; i < col; i++)
                    dotProducts[i] = dotProduct(A, i, col);

               for (int row = 0; row < n; row++)
               {
                    double total_subtract = 0.0;
                    for (int j = 0; j < col; j++)
                         total_subtract += dotProducts[j] * A[row][j];
                    A[row][col] -= total_subtract;
               }
          }
          columnNorm = colNorm(A, col);
          for (int row = 0; row < n; row++)
               A[row][col] /= columnNorm;
     }
     return A;
}

vector<vector<double>> MGS(vector<vector<double>> &A)
{
     int n = A.size();
     vector<double> dotProducts(n);

     for (int step = 0; step < n; step++)
          {
               double columnNorm = colNorm(A, step);
               for (int row = 0; row < n; row++)
               {
                    A[row][step] /= columnNorm;
               }

               for (int col = step + 1; col < n; col++)
                    {
                         dotProducts[col] = dotProduct(A, col, step);
                    }
               for (int row = 0; row < n; row++)
                    {
                         for(int col = step + 1; col < n; col++)   
                              {
                                   A[row][col] -= dotProducts[col] * A[row][step];
                              }
                    }
          }
     return A;
}

vector<vector<double>> constant_multiply(vector<vector<double>> &matrix, double constant)
{
     int n = matrix.size();
     for (int i = 0; i < n; i++)
          {
               for (int j = 0; j < n; j++)
                    {
                         matrix[i][j] = constant * matrix[i][j];
                    }
          }
     return matrix;
}

void printMatrix(const vector<vector<double>> &matrix)
{
     int n = matrix.size();

     cout << scientific << setprecision(5);
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               cout << setw(12) << matrix[i][j];
          }
          cout << endl;
     }
}

vector<vector<double>> matrixMultiply(vector<vector<double>> matrix1, vector<vector<double>> matrix2)
{
     int n = matrix1.size();
     vector<vector<double>> result(n, vector<double>(n, 0.0));
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               for (int k = 0; k < n; k++)
               {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
               }
          }
     }
     return result;
}

double colNorm(const vector<vector<double>> A, int col)
{
     int n = A.size();
     double norm = 0.0;
     for (int row = 0; row < n; row++)
          {
               norm += A[row][col] * A[row][col];
          }
     return sqrt(norm);
}

double dotProduct(const vector<vector<double>> A, int col1, int col2)
{
     int n = A.size();
     double dot_product = 0.0;
     for (int i = 0; i < n; i++)
          {
               dot_product += A[i][col1] * A[i][col2];
          }
     return dot_product;
}

vector<vector<double>> transpose(vector<vector<double>> matrix)
{
     int n = matrix.size();
     vector<vector<double>> result(n, vector<double>(n, 0.0));
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               result[i][j] = matrix[j][i];
          }
     }
     return result;
}

vector<vector<double>> subtractMatrices(
    vector<vector<double>> A,
    vector<vector<double>> B)
{
     int n = A.size();
     vector<vector<double>> result(n, vector<double>(n, 0.0));

     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               result[i][j] = A[i][j] - B[i][j];
          }
     }
     return result;
}

vector<vector<double>> ErrorMatrix(vector<vector<double>>& Q)
{
     int n = Q.size();
     vector<vector<double>> I(n, vector<double>(n, 0.0));

     for (int i = 0; i < n; i++)
          I[i][i] = 1;

     Q = subtractMatrices(matrixMultiply(transpose(Q), Q), I);

     for (int i = 0; i < n; i++)
          {
          for (int j = 0; j < n; j++)
               {
                    Q[i][j] = fabs(Q[i][j]);
               }
          }

     return Q;
}

double FrobeniusNorm(vector<vector<double>> matrix)
{
     int n = matrix.size();
     double result = 0;
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               result += fabs(matrix[i][j]) * fabs(matrix[i][j]);
          }
     }
     return pow(result, 0.5);
}