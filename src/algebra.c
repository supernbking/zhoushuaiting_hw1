#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            m.data[i][j] = 0.0; // 将矩阵的元素初始化为0
        }
    }
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols) // 如果无法加
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }

    Matrix result;
    result.rows = a.rows;
    result.cols = a.cols;

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }

    return result;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols) // 如果无法减
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }

    Matrix result;
    result.rows = a.rows;
    result.cols = a.cols;

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }

    return result;
}
Matrix mul_matrix(Matrix a, Matrix b)
{
    if (a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }

    Matrix result;
    result.rows = a.rows;
    result.cols = b.cols;

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[i][j] = 0.0;
            for (int k = 0; k < a.cols; k++)
            {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix result;
    result.rows = a.rows;
    result.cols = a.cols;

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix result;
    result.rows = a.cols;
    result.cols = a.rows;

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}

double det_matrix(Matrix a)
{
    if (a.cols != a.rows)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    if (a.rows == 1)
    {
        return a.data[0][0];
    }

    double result = 0.0;

    for (int j = 0; j < a.cols; j++)
    {
        Matrix diguimatrix; // 创建新矩阵储存递归用的矩阵 rows与cols都是原矩阵减1
        diguimatrix.rows = a.rows - 1;
        diguimatrix.cols = a.cols - 1;

        for (int i = 1; i < a.rows; i++) // 跳过第j列 生成Laplace定理所需的矩阵
        {
            int k = 0;
            for (int t = 0; t < a.cols; t++)
            {
                if (t != j)
                {
                    diguimatrix.data[i - 1][k++] = a.data[i][t];
                }
            }
        }

        double diguimatrix_det = det_matrix(diguimatrix);
        result += (j % 2 == 0 ? 1 : -1) * a.data[0][j] * diguimatrix_det;
    }

    return result;
}

Matrix inv_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    if (det_matrix(a) == 0.0)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }

    Matrix bansui_Matric;
    bansui_Matric.rows = a.rows;
    bansui_Matric.cols = a.rows;

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.rows; j++)
        {
            Matrix newmatrix; // 用于计算代数余子式
            newmatrix.rows = a.rows - 1;
            newmatrix.cols = a.rows - 1;

            int new_i = 0;
            int new_j = 0;

            for (int k = 0; k < a.rows; k++)
            {
                if (k != i)
                {
                    new_j = 0;
                    for (int m = 0; m < a.rows; m++)
                    {
                        if (m != j)
                        {
                            newmatrix.data[new_i][new_j] = a.data[k][m];
                            new_j++;
                        }
                    }
                    new_i++;
                }
            }

            double num = (((i + j) % 2 == 0 ? 1 : -1) * det_matrix(newmatrix));

            bansui_Matric.data[j][i] = num;
        }
    }
    Matrix result;
    result.rows = a.rows;
    result.cols = a.cols;

    double det_result = det_matrix(a);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.rows; j++)
        {
            result.data[i][j] = bansui_Matric.data[i][j] / det_result;
        }
    }
    return result;
}

int rank_matrix(Matrix a)
{
    int rank = 0;
    for (int i = 0; i < a.rows; i++)
    {
        int flag = -1; // 用于判定当前列是否全为0
        for (int j = i; j < a.rows; j++)
        {
            if (a.data[j][i] != 0.0)
            {
                flag = j;
                break;
            }
        }

        if (flag == -1)
        {
            // 若当前列全为零，直接跳过
            continue;
        }

        if (flag != i) // 若下一行此列为零且之后某行此列不为零，则将此行与下一行调换位置
        {
            for (int j = 0; j < a.cols; j++)
            {
                double temp = a.data[i][j];
                a.data[i][j] = a.data[flag][j];
                a.data[flag][j] = temp;
            }
        }
        for (int j = i + 1; j < a.rows; j++) // 将该列往后的数字都减为0
        {
            double beishu = a.data[j][i] / a.data[i][i];
            for (int k = i; k < a.cols; k++)
            {
                a.data[j][k] -= beishu * a.data[i][k];
            }
        }
        rank++;
    }
    return rank;
}

double trace_matrix(Matrix a)
{
    if (a.cols != a.rows)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    double trace = 0.0;
    for (int i = 0; i < a.rows; i++)
    {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}