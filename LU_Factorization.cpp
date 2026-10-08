#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<vector<double>> getMatrixA();
vector<vector<double>> getMatrixB();
void printMatrix(vector<vector<double>> matrix);
vector<vector<double>> LU_no_pivoting(vector<vector<double>> matrix);
vector<vector<double>> LU_partial_pivoting(vector<vector<double>> matrix);
vector<vector<double>> LU_complete_pivoting(vector<vector<double>> matrix);
vector<vector<double>> subtractMatrices(
    vector<vector<double>> matrix,
    vector<vector<double>> matrixB);
double FrobeniusNorm(vector<vector<double>> matrix);
void swapRows(vector<vector<double>> &matrix, int a, int b);
void swapColumns(vector<vector<double>> &matrix, int a, int b);
vector<vector<double>> matrixMultiply(vector<vector<double>> matrix1, vector<vector<double>> matrix2);
vector<vector<double>> transpose(vector<vector<double>> matrix);
int main()

{

vector<vector<double>> A = getMatrixA();
vector<vector<double>> originalA = getMatrixA();

vector<vector<double>> B = getMatrixB();
vector<vector<double>> originalB = getMatrixB();

cout << "r1A: " << FrobeniusNorm(subtractMatrices(originalA, LU_no_pivoting(A)))/FrobeniusNorm(originalA) << endl;
cout << "r2A: " << FrobeniusNorm(subtractMatrices(originalA, LU_partial_pivoting(A))) / FrobeniusNorm(originalA) << endl;
cout << "r3A: " << FrobeniusNorm(subtractMatrices(originalA, LU_complete_pivoting(A))) / FrobeniusNorm(originalA) << endl << endl << endl;

cout << "r1B: " << FrobeniusNorm(subtractMatrices(originalB, LU_no_pivoting(B))) / FrobeniusNorm(originalB) << endl;
cout << "r2B: " << FrobeniusNorm(subtractMatrices(originalB, LU_partial_pivoting(B))) / FrobeniusNorm(originalB) << endl;
cout << "r3B: " << FrobeniusNorm(subtractMatrices(originalB, LU_complete_pivoting(B))) / FrobeniusNorm(originalB) << endl;

return 0;
};



vector<vector<double>> getMatrixA()
{
     vector<vector<double>> A = {
         {1.000e-8, 2.291e0, 1.328e5, 6.550e10},
         {2.291e0, 5.126e-5, 6.255e0, 1.606e1},
         {1.328e5, 6.255e0, 7.836e2, 3.068e3},
         {6.550e10, 1.606e1, 3.068e3, 1.738e5}};

     return A;
}

vector<vector<double>> getMatrixB()
{
     vector<vector<double>> B = {
         {1.000e0, 2.291e0, 1.328e0, 6.550e0},
         {2.291e0, 5.126e0, 6.255e0, 1.606e0},
         {1.328e0, 6.255e0, 7.836e0, 3.068e0},
         {6.550e0, 1.606e0, 3.068e0, 1.738e0}};

     return B;
}

void printMatrix(vector<vector<double>> matrix)
{
     int n = matrix.size();

     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               cout << matrix[i][j] << " ";
          }

          cout << endl;
     }
}

vector<vector<double>> LU_no_pivoting(vector<vector<double>> matrix)
{
int n = matrix.size();
vector<vector<double>> LowerTriangle(n, vector<double>(n, 0.0));
double mult_ik;

for (int i = 0; i < n; i++)
     {
     LowerTriangle[i][i] = 1.0;
     }

for (int k = 0; k < n - 1; k++)
{
     for (int j = k + 1; j < n; j++)
     {
          mult_ik = matrix[j][k] / matrix[k][k];
          matrix[j][k] = 0;
          LowerTriangle[j][k] = mult_ik;
          for (int i = k + 1; i < n; i++)
          {
               matrix[j][i] = matrix[j][i] - mult_ik * matrix[k][i];
          }
     }
}
return matrixMultiply(LowerTriangle, matrix);
}

vector<vector<double>> LU_partial_pivoting(vector<vector<double>> matrix)
{
     int n = matrix.size();
     vector<vector<double>> LowerTriangle(n, vector<double>(n, 0.0));
     vector<vector<double>> Identity(n, vector<double>(n, 0.0));
     double mult_ik;
     int pivot_row;
     double pivot;
     for (int i = 0; i < n; i++)
     {
          LowerTriangle[i][i] = 1.0;
          Identity[i][i] = 1.0;
     }

     for (int k = 0; k < n - 1; k++)
     {
          pivot_row = k;
          pivot = fabs(matrix[k][k]);
          for (int j = k + 1; j < n; j++)
          {
               if (fabs(matrix[j][k]) > pivot)
                    {
                         pivot_row = j;
                         pivot = fabs(matrix[j][k]);
                    }
          }

          swapRows(matrix, pivot_row, k);
          swapRows(Identity, pivot_row, k);

          for (int i = 0; i < k; i ++)
          {
               double temp = LowerTriangle[pivot_row][i];
               LowerTriangle[pivot_row][i] = LowerTriangle[k][i];
               LowerTriangle[k][i] = temp;
          }

          for (int j = k + 1; j < n; j++)
          {
               mult_ik = matrix[j][k] / matrix[k][k];
               matrix[j][k] = 0;
               LowerTriangle[j][k] = mult_ik;
               for (int i = k + 1; i < n; i++)
               {
                    matrix[j][i] = matrix[j][i] - mult_ik * matrix[k][i];
               }
          }
     }
     return matrixMultiply(transpose(Identity), matrixMultiply(LowerTriangle, matrix));
}

vector<vector<double>> LU_complete_pivoting(vector<vector<double>> matrix)
{
     int n = matrix.size();
     vector<vector<double>> LowerTriangle(n, vector<double>(n, 0.0));
     vector<vector<double>> IdentityRow(n, vector<double>(n, 0.0));
     vector<vector<double>> IdentityColumn(n, vector<double>(n, 0.0));

     for (int i = 0; i < n; i++)
     {
          LowerTriangle[i][i] = 1.0;
          IdentityRow[i][i] = 1.0;
          IdentityColumn[i][i] = 1.0;
     }

     for (int k = 0; k < n - 1; k++)
     {
          // Search the whole remaining submatrix for the largest entry
          int pivot_row = k, pivot_column = k;
          double max_val = std::abs(matrix[k][k]);
          for (int i = k; i < n; i++)
          {
               for (int j = k; j < n; j++)
               {
                    if (std::abs(matrix[i][j]) > max_val)
                    {
                         max_val = std::abs(matrix[i][j]);
                         pivot_row = i;
                         pivot_column = j;
                    }
               }
          }

          if (max_val == 0.0) // remaining submatrix is all zeros -> singular
               break;

          swapRows(matrix, pivot_row, k);
          swapRows(IdentityRow, pivot_row, k);
          for (int i = 0; i < k; i++)
               std::swap(LowerTriangle[pivot_row][i], LowerTriangle[k][i]);

          swapColumns(matrix, pivot_column, k);
          swapColumns(IdentityColumn, pivot_column, k);

          for (int j = k + 1; j < n; j++)
          {
               double mult = matrix[j][k] / matrix[k][k];
               matrix[j][k] = 0;
               LowerTriangle[j][k] = mult;
               for (int i = k + 1; i < n; i++)
                    matrix[j][i] -= mult * matrix[k][i];
          }
     }

     return matrixMultiply(
         matrixMultiply(transpose(IdentityRow), matrixMultiply(LowerTriangle, matrix)),
         transpose(IdentityColumn));
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

void swapRows(vector<vector<double>> &matrix, int row1, int row2)
{
     vector<double> temp = matrix[row1];
     matrix[row1] = matrix[row2];
     matrix[row2] = temp;
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

void swapColumns(vector<vector<double>> &matrix, int col1, int col2)
{
     int n = matrix.size();

     for (int i = 0; i < n; i++)
     {
          double temp = matrix[i][col1];
          matrix[i][col1] = matrix[i][col2];
          matrix[i][col2] = temp;
     }
}
