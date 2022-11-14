#include <MatVectLab.h>


void setup() {
  // put your setup code here, to run once:
Serial.begin(9600);

// /*BASIC ELEMENT WISE MATRIX OPERATIONS*/
//////////////////////////////////////////////////////////////
int vect_A[4] = { 6, 5, 4, 1 };
  
int vect_B[4] = {1,4,2,3};

int vect_C[3] = {8, 3, 6};

double vect_D[3] = {1, 2, 3};

double vect_E[3];


Serial.println("vect_A: ");
vectOp.print(vect_A);

vectOp.add(vect_A, vect_A, vect_B); // vect_A = vect_A + vect_B
Serial.println("vect_A: ");
vectOp.print(vect_A);


Serial.println("vect_B: ");
vectOp.print(vect_B);

vectOp.multiply(vect_B, vect_B, vect_A); // vect_B = vect_B * vect_A (element wise)
Serial.println("vect_A * vect_B (element_wise): ");
vectOp.print(vect_B);

vectOp.subtract(vect_D, vect_C, vect_D); // vect_D = vect_C - vect_D 
Serial.println("vect_D: ");
vectOp.print(vect_D);

vectOp.clear(vect_C);
vectOp.copy(vect_C, vect_D); // copy contents of vect_D into vect_C
Serial.println("vect_C: ");

/////////////////////////////////////////////////////////////////





// /*SCALING MATRICES*/
//////////////////////////////////////////////////////////////////////

// /* scalar operations */
// vectOp.copy(vect_E, vect_C);
// Serial.println("vect_E: ");
// vectOp.print(mat_E);

// vectOp.scale(vect_E, vect_E, 2.5);
// Serial.println("vect_E*2.5: ");
// vectOp.print(vect_E);

// vectOp.scaleDiv(vect_E, vect_E, 2.5);
// Serial.println("vect_E/2.5: ");
// vectOp.print(vect_E);

//////////////////////////////////////////////////////////////////////





// /* DOT PRODUCT AND MAGNITUDE AND NORMALIZING */
/////////////////////////////////////////////////////////////////////

// double A[4] = { 5.45, 7.89, -3.05, 8.20};

// int B[4] = { 9, 2, 4, 8};

// double normalized_A[4];


// int dotAB = vectOp.dot<int>(A,B);
// double dotBA = vectOp.dot<double>(B,A);

// double magA = vectOp.magnitude<double>(A);
// int magA = vectOp.magnitude<int>(B);


// Serial.println("vect A: ");
// vectOp.print(A);
// Serial.println("vect B: ");
// vectOp.print(B);

// Serial.print("A.B: ");
// Serial.println(dotAB);

// Serial.print("B.A: ");
// Serial.println(dotBA);

// Serial.print("A magnitude: ");
// Serial.println(magA);

// Serial.print("B magnitude: ");
// Serial.println(magB);

// Serial.print("A normalized: ");
// vectOp.normalize(normalized_A, A);
// matOp.print(normalized_A);


// // convert vector to a column vector matrix
// double A_c [4][1];
// vectOp.vect2mat_C(A_c, A);
// matOp.print(A_c);

// // convert vector to a row vector matrix
// int B_r [1][4];
// vectOp.vect2mat_R(B_r, B);
// matOp.print(B_r);

//////////////////////////////////////////////////////////////////////////






// /* VECTOR TRANSFORMATION */
/////////////////////////////////////////////////////////////////////

// int A[4] = { 5, 7, -3, 8};
// int At[4][4] = { // transformation matrix
//   {1, -1, 2, 1},
//   {5, 8, 4, 2},
//   {6, 8, -3, 5},
//   {2, 1, -1, 0},  
// };

// Serial.println("vect A: ");
// vectOp.print(A);

// vectOp.transorm(A, At, A); // transform vectorA copy the ans to vector A

// Serial.println("vect A transformed: ");
// vectOp.print(A);


// double B[4] = { 9, 2, 4};
// double Bt[3][3] = { // transormation matrix
//   {0.5, -2.5, -1},
//   {1, -0.75, 3.45},
//   {0.875, 1.05, -1},      
// };


// Serial.println("vect B: ");
// vectOp.print(B);

// vectOp.transorm(B, Bt, B); // transform vectorB copy the ans to vector B

// Serial.println("vect B transformed: ");
// vectOp.print(B);

//////////////////////////////////////////////////////////////////////////







// /* 2D and 3D vectors SPECIFIC FUNCTION OPERATIONS */
/////////////////////////////////////////////////////////////////////

// int vect_A[3] = {8, 3, 6};
// double vect_B[3] = {1, 2, 3};
// float result3D[3];

// int vect_C[2] = {8, 3};
// double vect_D[2] = {2, 3};


// vectOp.print2D(vect_D);
// vectOp.print3D(vect_B);

// vectOp.cross(result3D, vect_A, vect_B); // result = A cross B
// vectOp.print3D(result3D);
// vectOp.cross(result3D, vect_B, vect_A); // result = B cross A
// vectOp.print3D(result3D);


// vectOp.cross(result3D, vect_C, vect_D); // result = C cross D
// vectOp.print3D(result3D);
// vectOp.cross(result3D, vect_D, vect_C); // result = D cross C
// vectOp.print3D(result3D);


// // compute unit vector direction
// double unit_vect_B[3];
// vectOp.normalize(unit_vect_B, B);

// // angle between 3D vectors
//  double cos_of_angle_btw_AB = cosineOfAngleBtw<double>(A, B);
// Serial.print("cos_of_angle_btw_AB: ");
// Serial.println(cos_of_angle_btw_AB);  

//   double angle_btw_AB_rad = angleBtwRad<double>(A, B); 
//   Serial.print("angle_btw_AB_rad: ");
// Serial.println(angle_btw_AB_rad); 

//   double angle_btw_AB_deg = angleBtwDeg<double>(A, B);
//   Serial.print("angle_btw_AB_deg: ");
// Serial.println(angle_btw_AB_deg);




//////////////////////////////////////////////////////////////////////////


}

void loop() {
  // put your main code here, to run repeatedly:
}