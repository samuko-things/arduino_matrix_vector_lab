#include <MatVectLab.h>


void setup() {
  // put your setup code here, to run once:
Serial.begin(9600);

  /*BASIC ELEMENT WISE MATRIX OPERATIONS*/
  //////////////////////////////////////////////////////////////
  int mat_A[2][4] = {
    { 6, 5, 4, 1 },
    { 2, 3, -7, 8 },
  };
    
  int mat_B[2][4] = {
     {1, 4, 2, 3},
     {6, -1, 0, 5},
  };

  int mat_C[3][3] = {
    {8, 3, 6},
    {5, 2, 7},
    {1, 0, 4},
  };

  int mat_D[3][3] = {
    {1, 2, 3},
    {4, 5, 6},
    {7, 8, 9},
  };

  double mat_E[3][3];


/* element wise operations */
matOp.add(mat_A, mat_A, mat_B); // mat_A = mat_A + mat_B
Serial.println("mat_A: ");
matOp.print(mat_A);


Serial.println("mat_B: ");
matOp.print(mat_B);

matOp.multiply(mat_B, mat_B, mat_A); // mat_B = mat_B * mat_A (element wise)
Serial.println("mat_A * mat_B (element_wise): ");
matOp.print(mat_B);

matOp.subtract(mat_D, mat_C, mat_D); // mat_D = mat_C - mat_D 
Serial.println("mat_D: ");
matOp.print(mat_D);

matOp.clear(mat_C);
matOp.copy(mat_C, mat_D); // copy contents of mat_D into mat_C
Serial.println("mat_C: ");

/////////////////////////////////////////////////////////////////




/*SCALING MATRICES*/
//////////////////////////////////////////////////////////////////////

// /* scalar operations */
// matOp.copy(mat_E, mat_C);
// Serial.println("mat_E: ");
// matOp.print(mat_E);

// matOp.scale(mat_E, mat_E, 2.5);
// Serial.println("mat_E*2.5: ");
// matOp.print(mat_E);

// matOp.scaleDiv(mat_E, mat_E, 2.5);
// Serial.println("mat_E/2.5: ");
// matOp.print(mat_E);

//////////////////////////////////////////////////////////////////////




/* DOT PRODUCT AND TRANSPOSE */
/////////////////////////////////////////////////////////////////////

// int A[3][2] = {
//   {5, 2},
//   {7, 4},
//   {3, 1},   
// };

// int B[2][3] = {
//   {9, 2, 4},
//   {-2, 3, 6},    
// };

// int AB[3][3];
// int BA[2][2];

// Serial.println("matrix A: ");
// matOp.print(A);
// Serial.println("matrix B: ");
// matOp.print(B);

// // matOp.dot(AB, B, A); // AB != B.A (2x3 dot 3x2 != 3x3 matrix) -> ERROR
// matOp.dot(AB, A, B); // AB = A.B (3x2 dot 2x3 = 3x3 matrix) -> CORRECT
// Serial.println("A.B: ");
// matOp.print(AB);

// // matOp.dot(BA, A, B); // BA != A.B (3x2 dot 2x3 != 2x2 matrix) -> ERROR
// matOp.dot(BA, B, A); // BA = B.A (2x3 dot 3x2 = 2x2 matrix) -> CORRECT
// Serial.println("B.A: ");
// matOp.print(BA);

// // transpose A(3x2) and put values in B (2x3)
// matOp.transpose(B,A); // B (2x3) = transpose of A(3x2)
// Serial.println("B (A_T): ");
// matOp.print(B);

//////////////////////////////////////////////////////////////////////////



/* INVERSE OF MATRIX */
/////////////////////////////////////////////////////////////////////////

// double A[3][3] = {
// {2, 7, 4},
// {3, 1, 6},
// {5, 0, 8},  
// };

// double A_inv[3][3];

// // inverse of matrix
// matOp.inverse(A_inv, A);

// Serial.println("A: ");
// matOp.print(A);

// Serial.println("A_inv: ");
// matOp.print(A_inv);
/////////////////////////////////////////////////////////////////////////

}

void loop() {
  // put your main code here, to run repeatedly:
}
