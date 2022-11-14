
#ifndef _MATRIXLAB_H_
#define _MATRIXLAB_H_
#endif

#include <math.h>

typedef unsigned int SIZE;

/* the rounding to decimal place algorithm was gotten from geeks for geeks and edited by me*/
// float round_dp(float val, int dp)
// {
//     if (val >= 0)
//     {
//         float rounded_val = (int)((val * pow(10, dp)) + 0.5);
//         return (float)rounded_val / pow(10, dp);
//     }
//     else
//     {
//         float rounded_val = (int)((val * pow(10, dp)) - 0.5);
//         return (float)rounded_val / pow(10, dp);
//     }
// }



class VectorOperations
{
public:
    // this function prints the vector
    template <typename dataType, SIZE N>
    void print(dataType (&vector1)[N])
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

    template <typename dataType>
    void print3D(dataType (&vector1)[3])
    {
        Serial.println(String(vector1[0]) + "i + " + String(vector1[1]) + "j + " + String(vector1[2]) + "k");
    }

    template <typename dataType>
    void print2D(dataType (&vector1)[2])
    {
        Serial.println(String(vector1[0]) + "i + " + String(vector1[1]) + "j");
    }

    // this function copies a vector to another vector of the same size
    template <typename dataType, SIZE N>
    void copy(dataType (&destinationVector)[N], dataType (&sourceVector)[N])
    {
        for (SIZE c = 0; c < N; c += 1)
        {
            destinationVector[c] = sourceVector[c];
        }
    }

    // this function clears a vector
    template <typename dataType, SIZE N>
    void clear(dataType (&vector1)[N])
    {
        for (SIZE c = 0; c < N; c += 1)
        {
            vector1[c] = 0;
        }
    }

    // this function adds two vectors together
    template <typename dataType, SIZE N>
    void add(dataType (&result)[N], dataType (&vector1)[N], dataType (&vector2)[N])
    {
        dataType buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = vector1[c] + vector2[c];
        }
        clear(result);
        copy(result, buffer);
    }

    // this function subtracts a vector from another
    template <typename dataType, SIZE N>
    void subtract(dataType (&result)[N], dataType (&vector1)[N], dataType (&vector2)[N])
    {
        dataType buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = vector1[c] - vector2[c];
        }
        clear(result);
        copy(result, buffer);
    }

    // this function multiply two vectors together element by element
    template <typename dataType, SIZE N>
    void multiply(dataType (&result)[N], dataType (&vector1)[N], dataType (&vector2)[N])
    {
        dataType buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = vector1[c] * vector2[c];
        }
        clear(result);
        copy(result, buffer);
    }

    // this function scale a vector
    template <typename dataType, SIZE N>
    void scale(dataType (&result)[N], dataType (&vector1)[N], dataType num)
    {
        dataType buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = vector1[c] * num;
        }
        clear(result);
        copy(result, buffer);
    }

    template <typename dataType, SIZE N>
    void scaleDiv(dataType (&result)[N], dataType (&vector1)[N], dataType num)
    {
        dataType buffer[N];
        for (SIZE c = 0; c < N; c += 1)
        {
            buffer[c] = vector1[c] / num;
        }
        clear(result);
        copy(result, buffer);
    }

    // this function performs dot product on two vector
    template <typename dataType, SIZE N>
    float dot(dataType (&vector1)[N], dataType (&vector2)[N])
    {
        dataType dot_prod;
        for (int c = 0; c < N; c += 1)
        {
            dot_prod += (vector1[c] * vector2[c]);
        }
        return dot_prod;
    }

    // this function performs cross product on two vector
    template <typename dataType>
    void cross(dataType (&result)[3], dataType (&vector1)[2], dataType (&vector2)[2])
    {
        dataType buffer[3];

        buffer[0] = (vector1[1] * 0) - (0 * vector2[1]);
        buffer[1] = -1 * ((vector1[0] * 0) - (0 * vector2[0]));
        buffer[2] = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);

        clear(result);
        copy(result, buffer);
    }

    template <typename dataType>
    void cross(dataType (&result)[3], dataType (&vector1)[3], dataType (&vector2)[3])
    {
        dataType buffer[3];

        buffer[0] = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
        buffer[1] = -1 * ((vector1[0] * vector2[2]) - (vector1[2] * vector2[0]));
        buffer[2] = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);

        clear(result);
        copy(result, buffer);
    }

    // this function calculates the magnitude of a vector
    template <typename dataType, SIZE N>
    float magnitude(dataType (&vector1)[N])
    {
        dataType mag;
        for (int c = 0; c < N; c += 1)
        {
            mag += pow(vector1[c], 2);
        }
        mag = sqrt(mag);

        return mag;
    }

    // this function normalize a vector (it generates unit vector for a 3D or 2D vector)
    template <typename dataType, SIZE N>
    void normalize(dataType (&result)[N], dataType (&vector1)[N])
    {
        dataType buffer[N];
        dataType mag = magnitude(vector1);

        for (int c = 0; c < N; c += 1)
        {
            buffer[N] = vector1[N] / mag;
        }

        clear(result);
        copy(result, buffer);
    }

    // this function calculates the Direction Cosine in radians
    void directionCosinesRad(float (&result)[3], float (&vector1)[3])
    {
        float buffer[3];
        float unit_vector[3];
        normalize(unit_vector, vector1);

        buffer[0] = acos(unit_vector[0]);
        buffer[1] = acos(unit_vector[1]);
        buffer[2] = acos(unit_vector[2]);

        clear(result);
        copy(result, buffer);
    }

    // this function calculates the Direction Cosine in radians
    void directionCosinesDeg(float (&result)[3], float (&vector1)[3])
    {
        float buffer[3];
        float unit_vector[3];
        normalize(unit_vector, vector1);

        buffer[0] = acos(unit_vector[0]) * 180 / M_PI;
        buffer[1] = acos(unit_vector[1]) * 180 / M_PI;
        buffer[2] = acos(unit_vector[2]) * 180 / M_PI;

        clear(result);
        copy(result, buffer);
    }

    // this function calculates the angle between two vectors
    float angleBtwRad(float (&vector1)[3], float (&vector2)[3])
    {
        float unit_vector1[3], unit_vector2[3];

        normalize(unit_vector1, vector1);
        normalize(unit_vector2, vector2);

        float angle = acos((unit_vector1[0] * unit_vector2[0]) + (unit_vector1[1] * unit_vector2[1]) + (unit_vector1[2] * unit_vector2[2]));
        return angle;
    }

    // this function calculates the angle between two vectors
    float angleBtwDeg(float (&vector1)[3], float (&vector2)[3])
    {
        float unit_vector1[3], unit_vector2[3];

        normalize(unit_vector1, vector1);
        normalize(unit_vector2, vector2);

        float angle = acos((unit_vector1[0] * unit_vector2[0]) + (unit_vector1[1] * unit_vector2[1]) + (unit_vector1[2] * unit_vector2[2])) * 180 / M_PI;
        return angle;
    }

    // this function rounds a vector to the nearest decimal place
    template <SIZE N>
    void round_dp(float (&rounded_vector)[N], float (&vector1)[N], int dp)
    {
        float buffer[N];
        for (int c = 0; c < N; c += 1)
        {
            buffer[c] = rounded(vector1[c], dp);
        }
        clear(rounded_vector);
        copy(rounded_vector, buffer);
    }

private:
    float rounded(float val, int dp)
    {
        if (val >= 0)
        {
            float rounded_val = (int)((val * pow(10, dp)) + 0.5);
            return (float)rounded_val / pow(10, dp);
        }
        else
        {
            float rounded_val = (int)((val * pow(10, dp)) - 0.5);
            return (float)rounded_val / pow(10, dp);
        }
    }
};

VectorOperations vectOp;







class MatrixOperations
{
public:
    // this function prints the matrix
    template <typename dataType, SIZE R, SIZE C>
    void print(dataType (&matrix)[R][C])
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
    template <typename dataType, SIZE R, SIZE C>
    void copy(dataType (&destinationMatrix)[R][C], dataType (&sourceMatrix)[R][C])
    {
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                destinationMatrix[r][c] = sourceMatrix[r][c];
            }
        }
    }

    // this function clears a matrix
    template <typename dataType, SIZE R, SIZE C>
    void clear(dataType (&matrix)[R][C])
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
    template <typename dataType, SIZE R, SIZE C>
    void add(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType (&matrix2)[R][C])
    {
        dataType buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = matrix1[r][c] + matrix2[r][c];
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function subtracts a matrix or vector from another
    template <typename dataType, SIZE R, SIZE C>
    void subtract(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType (&matrix2)[R][C])
    {
        dataType buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = matrix1[r][c] - matrix2[r][c];
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function transpose a matrix or vector
    template <typename dataType, SIZE R, SIZE C>
    void transpose(dataType (&T_matrix)[C][R], dataType (&matrix)[R][C])
    {
        dataType buffer[C][R];
        for (SIZE r = 0; r < C; r += 1)
        {
            for (SIZE c = 0; c < R; c += 1)
            {
                buffer[r][c] = matrix[c][r];
            }
        }
        clear(T_matrix);
        copy(T_matrix, buffer);
    }

    // this function multiplies two matrices together
    template <typename dataType, SIZE R, SIZE C, SIZE N>
    void dot(dataType (&result)[R][C], dataType (&matrix1)[R][N], dataType (&matrix2)[N][C])
    {
        dataType buffer[R][C];
        dataType val;
        dataType T_matrix2[C][N];
        transpose(T_matrix2, matrix2); // transpose matrix2 and store ans in T_matrix2

        for (SIZE row = 0; row < R; row += 1)
        {
            for (SIZE col = 0; col < C; col += 1)
            {
                for (SIZE count = 0; count < N; count += 1)
                {
                    val += matrix1[row][count] * T_matrix2[col][count];
                }
                buffer[row][col] = val;
                val = 0;
            }
        }

        clear(result);
        copy(result, buffer);
    }

    // this function multiplies a matrix and a vector, stores the ans in a vector
    template <typename dataType, SIZE N>
    void dot(dataType (&result)[N], dataType (&matrix)[N][N], dataType (&vector1)[N])
    {
        dataType buffer[N];
        dataType val;

        for (SIZE row = 0; row < N; row += 1)
        {
            for (SIZE col = 0; col < 1; col += 1)
            {
                for (SIZE count = 0; count < N; count += 1)
                {
                    val += matrix[row][count] * vector1[count];
                }
                buffer[row] = val;
                val = 0;
            }
        }

        vectOp.clear(result);
        vectOp.copy(result, buffer);
    }

    // this function multiplies a matrix or vector by a scaler
    template <typename dataType, SIZE R, SIZE C>
    void scale(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType num)
    {
        dataType buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = matrix1[r][c] * num;
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function divides a matrix or vector by a scaler
    template <typename dataType, SIZE R, SIZE C>
    void scaleDiv(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType num)
    {
        dataType buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = matrix1[r][c] / num;
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function multiplies matrices together element by element
    template <typename dataType, SIZE R, SIZE C>
    void multiply(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType (&matrix2)[R][C])
    {
        dataType buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = matrix1[r][c] * matrix2[r][c];
            }
        }
        clear(result);
        copy(result, buffer);
    }

    // this function calculates the determinant of a given matrix
    template <typename dataType>
    dataType det(dataType (&matrix)[1][1])
    {
        return matrix[0][0];
    }

    template <typename dataType, SIZE N>
    dataType det(dataType (&matrix)[N][N])
    {
        dataType d = 0;
        dataType minor_matrix[N - 1][N - 1];

        for (SIZE c = 0; c < N; c += 1)
        {
            minor(minor_matrix, matrix, 0, c);
            if (c == 0 || (c % 2) == 0)
                d += (matrix[0][c] * det(minor_matrix));
            else
                d -= (matrix[0][c] * det(minor_matrix));
        }
        return d;
    }

    // this generates inverse of a given matrix
    template <typename dataType, SIZE N>
    void inverse(double (&inv_matrix)[N][N], dataType (&matrix)[N][N])
    {
        double _matrix[N][N];
        toDouble(_matrix, matrix);

        double buffer[N][N];

        double d = det(_matrix);
   
        minorMatrix(buffer, _matrix);
        cofactorMatrix(buffer, buffer);
        adjointMatrix(buffer, buffer);

        scaleDiv(inv_matrix, buffer, d);
    }

    // this rounds a given matrix to the nearest decimal place , dp=0 gives the nearest whole number
    template <SIZE N>
    void round_dp(float (&rounded_matrix)[N][N], float (&matrix)[N][N], int dp)
    {
        float buffer[N][N];
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

    // this function copies a matrix to another matrix of the same size
    template <typename dataType, SIZE R, SIZE C>
    void toFloat(float (&destinationMatrix)[R][C], dataType (&sourceMatrix)[R][C])
    {
        float buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (float)sourceMatrix[r][c];
            }
        }

        clear(destinationMatrix);
        copy(destinationMatrix, buffer);
    }


    template <typename dataType, SIZE R, SIZE C>
    void toDouble(double (&destinationMatrix)[R][C], dataType (&sourceMatrix)[R][C])
    {
        double buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (double)sourceMatrix[r][c];
            }
        }

        clear(destinationMatrix);
        copy(destinationMatrix, buffer);
    }


    // this function copies a matrix to another matrix of the same size
    template <typename dataType, SIZE R, SIZE C>
    void toInt(int (&destinationMatrix)[R][C], dataType (&sourceMatrix)[R][C])
    {
        int buffer[R][C];
        for (SIZE r = 0; r < R; r += 1)
        {
            for (SIZE c = 0; c < C; c += 1)
            {
                buffer[r][c] = (int)sourceMatrix[r][c];
            }
        }

        clear(destinationMatrix);
        copy(destinationMatrix, buffer);
    }


    

private:
    float rounded(float val, int dp)
    {
        if (val >= 0)
        {
            float rounded_val = (int)((val * pow(10, dp)) + 0.5);
            return (float)rounded_val / pow(10, dp);
        }
        else
        {
            float rounded_val = (int)((val * pow(10, dp)) - 0.5);
            return (float)rounded_val / pow(10, dp);
        }
    }

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

    // this functon gets the minor matrix element for a giving matrix position
    template <typename dataType, SIZE N>
    void minor(dataType (&minor_matrix)[N - 1][N - 1], dataType (&major_matrix)[N][N], int pos_r, int pos_c)
    {
        SIZE minor_r = 0, minor_c = 0;

        dataType buffer[N - 1][N - 1];

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
                        buffer[minor_r][minor_c] = major_matrix[r][c];
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
    template <typename dataType, SIZE N>
    void minorMatrix(dataType (&minor_matrix)[N][N], dataType (&major_matrix)[N][N])
    {
        dataType buffer[N][N];
        dataType Minor[N - 1][N - 1];

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
    template <typename dataType, SIZE N>
    void cofactorMatrix(dataType (&co_matrix)[N][N], dataType (&minor_matrix)[N][N])
    {
        dataType buffer[N][N];

        for (SIZE r = 0; r < N; r += 1)
        {
            for (SIZE c = 0; c < N; c += 1)
            {
                if ((r + c) == 0 || ((r + c) % 2) == 0)
                    buffer[r][c] = minor_matrix[r][c];
                else
                    buffer[r][c] = -1 * minor_matrix[r][c];
            }
        }

        clear(co_matrix);
        copy(co_matrix, buffer);
    }

    // this generates adjoint of a given cofactor matrix
    template <typename dataType, SIZE N>
    void adjointMatrix(dataType (&ad_matrix)[N][N], dataType (&co_matrix)[N][N])
    {
        transpose(ad_matrix, co_matrix);
    }
};

MatrixOperations matOp;
