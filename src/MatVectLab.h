
#ifndef _MATVECTLAB_H_
#define _MATVECTLAB_H_
#endif

#include <math.h>

typedef unsigned int SIZE;






class VectorOperations
{
public:
    // this function prints the vector
    template <typename dtype, SIZE N>
    void print(dtype (&vector1)[N])
    {
        Serial.println(String(N) + "D Vector");
        Serial.print("[");
        for (SIZE c = 0; c < N; c += 1)
        {
            if (c == N - 1)
                Serial.print(String(vector1[c]));
            else
                Serial.print(String(vector1[c]) + ", ");
        }
        Serial.println("]");
        Serial.println();
    }

    template <typename dtype>
    void print3D(dtype (&vector1)[3])
    {
        Serial.println(String(vector1[0]) + "i + " + String(vector1[1]) + "j + " + String(vector1[2]) + "k");
    }

    template <typename dtype>
    void print2D(dtype (&vector1)[2])
    {
        Serial.println(String(vector1[0]) + "i + " + String(vector1[1]) + "j");
    }

    // this function copies a vector to another vector of the same size
    template <typename dtype, typename dtype1, SIZE N>
    void copy(dtype (&destinationVector)[N], dtype1 (&sourceVector)[N])
    {
        for (SIZE c = 0; c < N; c += 1)
        {
            destinationVector[c] = (dtype)sourceVector[c];
        }
    }

    // this function converts a matrix to a vector
    template <typename dtype, typename dtype1, SIZE N>
    void mat2vect(dtype (&destinationVector)[N], dtype1 (&sourceMatrix)[1][N])
    {
        for (SIZE n = 0; n < N; n += 1)
        {
            destinationVector[n] = (dtype)sourceMatrix[0][n];
        }
    }

    template <typename dtype, typename dtype1, SIZE N>
    void mat2vect(dtype (&destinationVector)[N], dtype1 (&sourceMatrix)[N][1])
    {
        for (SIZE n = 0; n < N; n += 1)
        {
            destinationVector[n] = (dtype)sourceMatrix[n][0];
        }
    }

    // this function copies a vector to a matrix to create a column vector matrix
    template <typename dtype, typename dtype1, SIZE N>
    void vect2mat_C(dtype (&destinationMatrix)[N][1], dtype1 (&sourceVector)[N])
    {
        for (SIZE n = 0; n < N; n += 1)
        {
            destinationMatrix[n][1] = (dtype)sourceVector[n];
        }
    }

    // this function clears a vector
    template <typename dtype, SIZE N>
    void clear(dtype (&vector1)[N])
    {
        for (SIZE c = 0; c < N; c += 1)
        {
            vector1[c] = 0;
        }
    }

    // this function adds two vectors together
    template <typename dtype, typename dtype1, typename dtype2, SIZE N>
    void add(dtype (&result)[N], dtype1 (&vector1)[N], dtype2 (&vector2)[N])
    {
        double buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = (double)vector1[c] + (double)vector2[c];
        }
        clear(result);
        copy(result, buffer);
    }

    // this function subtracts a vector from another
    template <typename dtype, typename dtype1, typename dtype2, SIZE N>
    void subtract(dtype (&result)[N], dtype1 (&vector1)[N], dtype2 (&vector2)[N])
    {
        double buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = (double)vector1[c] - (double)vector2[c];
        }
        clear(result);
        copy(result, buffer);
    }

    // this function multiply and divide two vectors together element by element
    template <typename dtype, typename dtype1, typename dtype2, SIZE N>
    void multiply(dtype (&result)[N], dtype1 (&vector1)[N], dtype2 (&vector2)[N])
    {
        double buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = (double)vector1[c] * (double)vector2[c];
        }
        clear(result);
        copy(result, buffer);
    }

    template <typename dtype, typename dtype1, typename dtype2, SIZE N>
    void divide(dtype (&result)[N], dtype1 (&vector1)[N], dtype2 (&vector2)[N])
    {
        double buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = (double)vector1[c] / (double)vector2[c];
        }
        clear(result);
        copy(result, buffer);
    }


    // this function scale a vector
    template <typename dtype, typename dtype1, typename dtype2, SIZE N>
    void scale(dtype (&result)[N], dtype1 (&vector1)[N], dtype2 num)
    {
        double buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = (double)vector1[c] * (double)num;
        }
        clear(result);
        copy(result, buffer);
    }

    template <typename dtype, typename dtype1, typename dtype2, SIZE N>
    void scaleDiv(dtype (&result)[N], dtype1 (&vector1)[N], dtype2 num)
    {
        double buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = (double)vector1[c] / (double)num;
        }
        clear(result);
        copy(result, buffer);
    }

    // this function performs dot product on two vector
    template <typename dtype, typename dtype1, typename dtype2, SIZE N>
    dtype dot(dtype1 (&vector1)[N], dtype2 (&vector2)[N])
    {
        double dot_prod;
        for (int c = 0; c < N; c += 1)
        {
            dot_prod += (double)vector1[c] * (double)vector2[c];
        }
        return dot_prod;
    }

    // this function transforms a vector using a transformation matrix
    template <typename dtype, typename dtype1, typename dtype2, SIZE N>
    void transform(dtype (&result_vector)[N], dtype1 (&transformationMatrix)[N][N], dtype2 (&vector)[N])
    {
        double buffer[N];
        double val;

        for (SIZE row = 0; row < N; row += 1)
        {
            for (SIZE col = 0; col < 1; col += 1)
            {
                for (SIZE count = 0; count < N; count += 1)
                {
                    val += (double)transformationMatrix[row][count] * (double)vector[count];
                }
                buffer[row] = val;
                val = 0;
            }
        }

        clear(result_vector);
        copy(result_vector, buffer);
    }

    // this function performs cross product on two vector
    template <typename dtype, typename dtype1, typename dtype2>
    void cross(dtype (&result)[3], dtype1 (&vector1)[2], dtype2 (&vector2)[2])
    {
        double buffer[3];

        buffer[0] = ((double)vector1[1] * 0) - (0 * (double)vector2[1]);
        buffer[1] = -1 * (((double)vector1[0] * 0) - (0 * (double)vector2[0]));
        buffer[2] = ((double)vector1[0] * (double)vector2[1]) - ((double)vector1[1] * (double)vector2[0]);

        clear(result);
        copy(result, buffer);
    }

    template <typename dtype, typename dtype1, typename dtype2>
    void cross(dtype (&result)[3], dtype1 (&vector1)[3], dtype2 (&vector2)[3])
    {
        double buffer[3];

        buffer[0] = ((double)vector1[1] * (double)vector2[2]) - ((double)vector1[2] * (double)vector2[1]);
        buffer[1] = -1 * (((double)vector1[0] * (double)vector2[2]) - ((double)vector1[2] * (double)vector2[0]));
        buffer[2] = ((double)vector1[0] * (double)vector2[1]) - ((double)vector1[1] * (double)vector2[0]);

        clear(result);
        copy(result, buffer);
    }

    // this function calculates the magnitude of a vector
    template <typename dtype, typename dtype1, SIZE N>
    dtype magnitude(dtype1 (&vector1)[N])
    {
        double mag;
        for (int c = 0; c < N; c += 1)
        {
            mag += pow((double)vector1[c], 2.0);
        }
        mag = sqrt(mag);

        return (dtype)mag;
    }

    // this function normalize a vector (it generates unit vector for a 3D or 2D vector)
    template <typename dtype, typename dtype1, SIZE N>
    void normalize(dtype (&result)[N], dtype1 (&vector1)[N])
    {
        double buffer[N];
        double mag = magnitude<double>(vector1);

        for (int c = 0; c < N; c += 1)
        {
            buffer[N] = (double)vector1[N] / mag;
        }

        clear(result);
        copy(result, buffer);
    }

    // this function calculates the cosine of angle between two 3Dvectors
    template <typename dtype, typename dtype1, typename dtype2>
    dtype cosineOfAngleBtw(dtype1 (&vector1)[3], dtype2 (&vector2)[3])
    {
        double unit_vector1[3], unit_vector2[3];

        normalize(unit_vector1, vector1);
        normalize(unit_vector2, vector2);

        double cosine_of_angle = dot<double>(unit_vector1, unit_vector2);
        return (dtype)cosine_of_angle;
    }

    // this function calculates the angle between two 3Dvectors
    template <typename dtype, typename dtype1, typename dtype2>
    dtype angleBtwRad(dtype1 (&vector1)[3], dtype2 (&vector2)[3])
    {
        double angle = acos(cosineOfAngleBtw<double>(vector1, vector2));
        return (dtype)angle;
    }

    // this function calculates the angle between two 3Dvectors
    template <typename dtype, typename dtype1, typename dtype2>
    dtype angleBtwDeg(dtype1 (&vector1)[3], dtype2 (&vector2)[3])
    {
        double angle = acos(cosineOfAngleBtw<double>(vector1, vector2)) * 180/M_PI;
        return (dtype)angle;
    }

    // this function rounds a vector to the nearest decimal place
    template <SIZE N>
    double round_dp(double (&rounded_vector)[N], double (&vector1)[N], int dp)
    {
        double buffer[N];
        for (int c = 0; c < N; c += 1)
        {
            buffer[c] = rounded(vector1[c], dp);
        }
        clear(rounded_vector);
        copy(rounded_vector, buffer);
    }

private:
    double rounded(double val, int dp)
    {
        if (val >= 0)
        {
            double rounded_val = (int)((val * pow(10, dp)) + 0.5);
            return (double)rounded_val / pow(10, dp);
        }
        else
        {
            double rounded_val = (int)((val * pow(10, dp)) - 0.5);
            return (double)rounded_val / pow(10, dp);
        }
    }
};

VectorOperations vectOp;









class MatrixOperations
{
public:
    // this function prints the matrix
    template <typename dtype, SIZE R, SIZE C>
    void print(dtype (&matrix)[R][C])
    {
        Serial.println(String(R) + " by " + String(C) + " Matrix");
        for (int r = 0; r < R; r += 1)
        {
            for (int c = 0; c < C; c += 1)
            {
                Serial.print(String(matrix[r][c]) + "\t");
            }
            Serial.println();
        }
        Serial.println();
    }

    // this function copies a matrix to another matrix of the same size
    template <typename dtype, typename dtype1, SIZE R, SIZE C>
    void copy(dtype (&destinationMatrix)[R][C], dtype1 (&sourceMatrix)[R][C])
    {
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                destinationMatrix[r][c] = (dtype)sourceMatrix[r][c];
            }
        }
    }


    // this function copies a vector to a matrix to create a row vector matrix
    template <typename dtype, typename dtype1, SIZE N>
    void vect2mat_R(dtype (&destinationMatrix)[1][N], dtype1 (&sourceVector)[N])
    {
        for (SIZE n = 0; n < N; n += 1)
        {
            destinationMatrix[0][n] = (dtype)sourceVector[n];
        }
    }

    // this function copies a vector to a matrix to create a column vector matrix
    template <typename dtype, typename dtype1, SIZE N>
    void vect2mat_C(dtype (&destinationMatrix)[N][1], dtype1 (&sourceVector)[N])
    {
        for (SIZE n = 0; n < N; n += 1)
        {
            destinationMatrix[n][0] = (dtype)sourceVector[n];
        }
    }

    // this function clears a matrix
    template <typename dtype, SIZE R, SIZE C>
    void clear(dtype (&matrix)[R][C])
    {
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                matrix[r][c] = 0;
            }
        }
    }


    // this function adds two matrices or vectors together
    template <typename dtype, typename dtype1, typename dtype2, SIZE R, SIZE C>
    void add(dtype (&result)[R][C], dtype1 (&matrix1)[R][C], dtype2 (&matrix2)[R][C])
    {
        double buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (double)matrix1[r][c] + (double)matrix2[r][c];
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function subtracts a matrix or vector from another
    template <typename dtype, typename dtype1, typename dtype2, SIZE R, SIZE C>
    void subtract(dtype (&result)[R][C], dtype1 (&matrix1)[R][C], dtype2 (&matrix2)[R][C])
    {
        double buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (double)matrix1[r][c] - (double)matrix2[r][c];
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function multiplies matrices together element by element
    template <typename dtype, typename dtype1, typename dtype2, SIZE R, SIZE C>
    void multiply(dtype (&result)[R][C], dtype1 (&matrix1)[R][C], dtype2 (&matrix2)[R][C])
    {
        double buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (double)matrix1[r][c] * (double)matrix2[r][c];
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function divide matrices together element by element
    template <typename dtype, typename dtype1, typename dtype2, SIZE R, SIZE C>
    void divide(dtype (&result)[R][C], dtype1 (&matrix1)[R][C], dtype2 (&matrix2)[R][C])
    {
        double buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (double)matrix1[r][c] / (double)matrix2[r][c];
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function transpose a matrix or vector
    template <typename dtype, typename dtype1, SIZE R, SIZE C>
    void transpose(dtype (&T_matrix)[C][R], dtype1 (&matrix)[R][C])
    {
        dtype buffer[C][R];
        for (SIZE r = 0; r < C; r += 1)
        {
            for (SIZE c = 0; c < R; c += 1)
            {
                buffer[r][c] = (dtype)matrix[c][r];
            }
        }
        clear(T_matrix);
        copy(T_matrix, buffer);
    }

    // this function multiplies two matrices together
    template <typename dtype, typename dtype1, typename dtype2, SIZE N, SIZE R, SIZE C>
    void dot(dtype (&result)[R][C], dtype1 (&matrix1)[R][N], dtype2 (&matrix2)[N][C])
    {
        double buffer[R][C];
        double val;
        dtype2 T_matrix2[C][N];
        transpose(T_matrix2, matrix2); // transpose matrix2 and store ans in T_matrix2

        for (SIZE row = 0; row < R; row += 1)
        {
            for (SIZE col = 0; col < C; col += 1)
            {
                for (SIZE count = 0; count < N; count += 1)
                {
                    val += (double)matrix1[row][count] * (double)T_matrix2[col][count];
                }
                buffer[row][col] = val;
                val = 0;
            }
        }

        clear(result);
        copy(result, buffer);
    }


    // this function multiplies a matrix or vector by a scaler
    template <typename dtype, typename dtype1, typename dtype2, SIZE R, SIZE C>
    void scale(dtype (&result)[R][C], dtype1 (&matrix1)[R][C], dtype2 num)
    {
        double buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (double)matrix1[r][c] * (double)num;
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function divides a matrix or vector by a scaler
    template <typename dtype, typename dtype1, typename dtype2, SIZE R, SIZE C>
    void scaleDiv(dtype (&result)[R][C], dtype1 (&matrix1)[R][C], dtype2 num)
    {
        double buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (double)matrix1[r][c] / (double)num;
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this generates inverse of a given matrix
    template <typename dtype, typename dtype1, SIZE N>
    void inverse(dtype (&inv_matrix)[N][N], dtype1 (&matrix)[N][N])
    {
        double buffer[N][N];
        double d = det(matrix);
   
        minorMatrix(buffer, matrix);
        cofactorMatrix(buffer, buffer);
        adjointMatrix(buffer, buffer);

        scaleDiv(inv_matrix, buffer, d);
    }

    // this rounds a given matrix to the nearest decimal place , dp=0 gives the nearest whole number
    template <SIZE N>
    void round_dp(double (&rounded_matrix)[N][N], double (&matrix)[N][N], int dp)
    {
        double buffer[N][N];
        for (int r = 0; r < N; r += 1)
        {
            for (int c = 0; c < N; c += 1)
            {
                buffer[r][c] = rounded(matrix[r][c], dp);
            }
        }
        clear(rounded_matrix);
        copy(rounded_matrix, buffer);
    }

private:
    double rounded(double val, int dp)
    {
        if (val >= 0)
        {
            double rounded_val = (int)((val * pow(10, dp)) + 0.5);
            return (double)rounded_val / pow(10, dp);
        }
        else
        {
            double rounded_val = (int)((val * pow(10, dp)) - 0.5);
            return (double)rounded_val / pow(10, dp);
        }
    }

    // this function calculates the determinant of a given matrix
    template <typename dtype>
    double det(dtype (&matrix)[1][1])
    {
        return (double)matrix[0][0];
    }

    template <typename dtype, SIZE N>
    double det(dtype (&matrix)[N][N])
    {
        double d = 0;
        double minor_matrix[N - 1][N - 1];

        for (SIZE c = 0; c < N; c += 1)
        {
            minor(minor_matrix, matrix, 0, c);
            if (c == 0 || (c % 2) == 0)
                d += ((double)matrix[0][c] * det(minor_matrix));
            else
                d -= ((double)matrix[0][c] * det(minor_matrix));
        }
        return d;
    }

    // this functon gets the minor matrix element for a giving matrix position
    template <typename dtype, typename dtype1, SIZE N>
    void minor(dtype (&minor_matrix)[N - 1][N - 1], dtype1 (&major_matrix)[N][N], int pos_r, int pos_c)
    {
        SIZE minor_r = 0, minor_c = 0;

        dtype buffer[N - 1][N - 1];

        if (N >= 2)
        {
            for (SIZE r = 0; r < N; r += 1)
            {

                for (SIZE c = 0; c < N; c += 1)
                {
                    if ((r == pos_r) || (c == pos_c))
                    {
                        continue;
                    }
                    else
                    {
                        buffer[minor_r][minor_c] = (dtype)major_matrix[r][c];
                        minor_c += 1;
                        if (minor_c >= N - 1)
                            minor_c = 0;
                    }
                }

                if (r == pos_r)
                {
                    continue;
                }
                else
                {
                    minor_r += 1;
                    if (minor_r >= N - 1)
                        minor_r = 0;
                }
            }
        }

        clear(minor_matrix);
        copy(minor_matrix, buffer);
    }

    // this generates a new matrices of the det of minors of a given matrix
    template <typename dtype, typename dtype1, SIZE N>
    void minorMatrix(dtype (&minor_matrix)[N][N], dtype1 (&major_matrix)[N][N])
    {
        double buffer[N][N];
        dtype1 Minor[N - 1][N - 1];

        for (SIZE r = 0; r < N; r += 1)
        {
            for (SIZE c = 0; c < N; c += 1)
            {
                minor(Minor, major_matrix, r, c);
                buffer[r][c] = det(Minor);
            }
        }

        clear(minor_matrix);
        copy(minor_matrix, buffer);
    }

    // this generates cofactor of a given matrix
    template <typename dtype, typename dtype1, SIZE N>
    void cofactorMatrix(dtype (&co_matrix)[N][N], dtype1 (&minor_matrix)[N][N])
    {
        dtype buffer[N][N];

        for (SIZE r = 0; r < N; r += 1)
        {
            for (SIZE c = 0; c < N; c += 1)
            {
                if ((r + c) == 0 || ((r + c) % 2) == 0)
                    buffer[r][c] = (dtype)minor_matrix[r][c];
                else
                    buffer[r][c] = -1 * (dtype)minor_matrix[r][c];
            }
        }

        clear(co_matrix);
        copy(co_matrix, buffer);
    }

    // this generates adjoint of a given cofactor matrix
    template <typename dtype, typename dtype1, SIZE N>
    void adjointMatrix(dtype (&ad_matrix)[N][N], dtype1 (&co_matrix)[N][N])
    {
        transpose(ad_matrix, co_matrix);
    }
};

MatrixOperations matOp;
