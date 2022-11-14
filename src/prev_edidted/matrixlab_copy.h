#ifndef _MATRIXLAB_H_
#define _MATRIXLAB_H_
#endif

namespace mat{

	// this function copies a matrix to another matrix of the same size
	template<typename dataType, size_t R, size_t C>
	void copy(dataType (&destinationMatrix)[R][C], dataType (&sourceMatrix)[R][C]) {
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      destinationMatrix[r][c] = sourceMatrix[r][c];
	    }
	  }
	}



	// this function clears a matrix
	template<typename dataType, size_t R, size_t C>
	void clear(dataType (&matrix)[R][C]) {
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      matrix[r][c] = 0;
	    }
	  }
	}

		// this function adds two matrices or vectors together
	template<typename dataType, size_t R, size_t C>
	void add(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType (&matrix2)[R][C]) {
	  dataType Buffer[R][C];
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      Buffer[r][c] = matrix1[r][c] + matrix2[r][c];
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function adds a number to all elements of a matrix or vector
	template<typename dataType, size_t R, size_t C>
	void add(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType num) {
	  dataType Buffer[R][C];
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      Buffer[r][c] = matrix1[r][c] + num;
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function adds a number to all elements of a matrix or vector
	template<typename dataType, size_t R, size_t C>
	void add(dataType (&result)[R][C], dataType num, dataType (&matrix1)[R][C]) {
	  dataType Buffer[R][C];
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      Buffer[r][c] = num + matrix1[r][c];
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function subtracts a matrix or vector from another
	template<typename dataType, size_t R, size_t C>
	void sub(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType (&matrix2)[R][C]) {
	  dataType Buffer[R][C];
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      Buffer[r][c] = matrix1[r][c] - matrix2[r][c];
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function subtracts a number from all elements of a matrix
	template<typename dataType, size_t R, size_t C>
	void sub(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType num) {
	  dataType Buffer[R][C];
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      Buffer[r][c] = matrix1[r][c] - num;
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function subtracts a number from all elements of a matrix
	template<typename dataType, size_t R, size_t C>
	void sub(dataType (&result)[R][C], dataType num, dataType (&matrix1)[R][C]) {
	  dataType Buffer[R][C];
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      Buffer[r][c] = num - matrix1[r][c];
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function transpose a matrix or vector
	template<typename dataType, size_t R, size_t C>
	void transpose(dataType (&T_matrix)[C][R], dataType (&matrix)[R][C]) {
	  dataType Buffer[C][R];
	  for (int r = 0; r < C; r += 1) {
	    for (int c = 0; c < R; c += 1) {
	      Buffer[r][c] = matrix[c][r];
	    }
	  }

	  clear(T_matrix);
	  copy(T_matrix, Buffer);
	}



	// this function multiplies two matrices together
	template<typename dataType, size_t R, size_t C, size_t N>
	void prod(dataType (&result)[R][C], dataType (&matrix1)[R][N], dataType (&matrix2)[N][C]) {
	  dataType Buffer[R][C];
	  dataType val;

	  dataType T_matrix2[C][N];
	  transpose(T_matrix2, matrix2); // transpose matrix2 and store ans in T_matrix2

	  for (int row = 0; row < R; row += 1) {
	    for (int col = 0; col < C; col += 1) {
	      for (int count = 0; count < N; count += 1) {
	        val += matrix1[row][count] * T_matrix2[col][count];
	      }
	      Buffer[row][col] = val;
	      val = 0;
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function multiplies a matrix or vector by a scaler
	template<typename dataType, size_t R, size_t C>
	void scale(dataType (&result)[R][C], dataType num, dataType (&matrix1)[R][C]) {
	  dataType Buffer[R][C];
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      Buffer[r][c] = num * matrix1[r][c];
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function multiplies matrices together element by element
	template<typename dataType, size_t R, size_t C>
	void eProd(dataType (&result)[R][C], dataType (&matrix1)[R][C], dataType (&matrix2)[R][C]) {
	  dataType Buffer[R][C];
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      result[r][c] = matrix1[r][c] * matrix2[r][c];
	    }
	  }

	  clear(result);
	  copy(result, Buffer);
	}



	// this function calculates the determinant of a 2*2 or 3*3 matrix
	template<typename dataType, size_t N>
	dataType det(dataType (&matrix)[N][N]) {
	  if (N == 2) {
	    dataType d = (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
	    return d;
	  }

	  else if (N == 3) {
	    dataType d = ( (matrix[0][0] * matrix[1][1] * matrix[2][2]) + (matrix[0][1] * matrix[1][2] * matrix[2][0]) + (matrix[0][2] * matrix[1][0] * matrix[2][1]) ) - ( (matrix[0][2] * matrix[1][1] * matrix[2][0]) + (matrix[0][0] * matrix[1][2] * matrix[2][1]) + (matrix[0][1] * matrix[1][0] * matrix[2][2]) );
	    return d;
	  }

	  else return 0;
	}



	// this function calculates the inverse of a 1*1, 2*2 or 3*3 matrix
	template<typename dataType, size_t N>
	void inverse(dataType (&inv_matrix)[N][N], dataType (&matrix)[N][N]) {
	  dataType Buffer[N][N];
	  if (N == 1) { // for 1 by 1 matrix
	    Buffer[0][0] = 1 / matrix[0][0];
	    clear(inv_matrix);
	    copy(inv_matrix, Buffer);
	  }

	  else if (N == 2) { // for 2 by 2 matrix
	    dataType d = det(matrix);

	    Buffer[0][0] = matrix[1][1] / d;
	    Buffer[0][1] = (-1 * matrix[0][1]) / d;
	    Buffer[1][0] = (-1 * matrix[1][0]) / d;
	    Buffer[1][1] = matrix[0][0] / d;

	    clear(inv_matrix);
	    copy(inv_matrix, Buffer);
	  }

	  else if (N == 3) { // for 3 by 3 matrix
	    dataType d = det(matrix);

	    Buffer[0][0] = ( (matrix[1][1] * matrix[2][2]) - (matrix[1][2] * matrix[2][1]) ) / d;
	    Buffer[0][1] = -1 * ( (matrix[0][1] * matrix[2][2]) - (matrix[0][2] * matrix[2][1]) ) / d;
	    Buffer[0][2] = ( (matrix[0][1] * matrix[1][2]) - (matrix[0][2] * matrix[1][1]) ) / d;

	    Buffer[1][0] = -1 * ( (matrix[1][0] * matrix[2][2]) - (matrix[1][2] * matrix[2][0]) ) / d;
	    Buffer[1][1] = ( (matrix[0][0] * matrix[2][2]) - (matrix[0][2] * matrix[2][0]) ) / d;
	    Buffer[1][2] = -1 * ( (matrix[0][0] * matrix[1][2]) - (matrix[0][2] * matrix[1][0]) ) / d;

	    Buffer[2][0] = ( (matrix[1][0] * matrix[2][1]) - (matrix[1][1] * matrix[2][0]) ) / d;
	    Buffer[2][1] = -1 * ( (matrix[0][0] * matrix[2][1]) - (matrix[0][1] * matrix[2][0]) ) / d;
	    Buffer[2][2] = ( (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]) ) / d;

	    clear(inv_matrix);
	    copy(inv_matrix, Buffer);
	  }

	}



	// this function prints the matrix
	template<typename dataType, size_t R, size_t C>
	void printMatrix(dataType (&matrix)[R][C]) {
	  Serial.println(String(R) + " by " + String(C) + " Matrix");
	  for (int r = 0; r < R; r += 1) {
	    for (int c = 0; c < C; c += 1) {
	      Serial.print(String(matrix[r][c]) + "\t");
	    }
	    Serial.println();
	  }
	  Serial.println();
	}



	








	// this function performs dot product on two vector
	template<typename dataType>
	dataType dot(dataType (&vector1)[3][1], dataType (&vector2)[3][1]) {
	  dataType dot_prod;

	  dot_prod = (vector1[0][0] * vector2[0][0]) + (vector1[1][0] * vector2[1][0]) + (vector1[2][0] * vector2[2][0]);
	  //  for (int count = 0; count < 3; count += 1) {
	  //    dotprod += (vector1[count][0] * vector2[count][0]);
	  //  }
	  return dot_prod;
	}


	// this function performs cross product on two vector
	template<typename dataType>
	void cross(dataType (&result)[3][1], dataType (&vector1)[3][1], dataType (&vector2)[3][1]) {
	  dataType Buffer[3][1];

	  Buffer[0][0] = (vector1[1][0] * vector2[2][0]) - (vector1[2][0] * vector2[1][0]);
	  Buffer[1][0] = -1 * ( (vector1[0][0] * vector2[2][0]) - (vector1[2][0] * vector2[0][0]) );
	  Buffer[2][0] = (vector1[0][0] * vector2[1][0]) - (vector1[1][0] * vector2[0][0]);

	  clear(result);
	  copy(result, Buffer);
	}


	// this function calculates the mag of a vector
	template<typename dataType>
	dataType mag(dataType (&vector)[3][1]) {
	  dataType magnitude = sqrt( (vector[0][0] * vector[0][0]) + (vector[1][0] * vector[1][0]) + (vector[2][0] * vector[2][0]) );
	  return magnitude;
	}


	// this function calculates the unit vector direction
	template<typename dataType>
	void DCM(dataType (&result)[3][1], dataType (&vector)[3][1]) {
	  dataType Buffer[3][1];
	  dataType magnitude = mag(vector);

	  Buffer[0][0] = vector[0][0] / magnitude;
	  Buffer[1][0] = vector[1][0] / magnitude;
	  Buffer[2][0] = vector[2][0] / magnitude;

	  clear(result);
	  copy(result, Buffer);
	}

	// this function calculates the Direction Cosine in deg
	template<typename dataType>
	void DCM_deg(dataType (&result)[3][1], dataType (&vector)[3][1]) {
	  dataType Buffer[3][1];

	  dataType unit_vector[3][1];
	  DCM(unit_vector, vector);

	  Buffer[0][0] = acos(unit_vector[0][0]) * 180 / PI;
	  Buffer[1][0] = acos(unit_vector[1][0]) * 180 / PI;
	  Buffer[2][0] = acos(unit_vector[2][0]) * 180 / PI;

	  clear(result);
	  copy(result, Buffer);
	}


	// this function calculates the Direction Cosine in radians
	template<typename dataType>
	void DCM_rad(dataType (&result)[3][1], dataType (&vector)[3][1]) {
	  dataType Buffer[3][1];

	  dataType unit_vector[3][1];
	  DCM(unit_vector, vector);

	  Buffer[0][0] = acos(unit_vector[0][0]);
	  Buffer[1][0] = acos(unit_vector[1][0]);
	  Buffer[2][0] = acos(unit_vector[2][0]);

	  clear(result);
	  copy(result, Buffer);
	}


	// this function calculates the angle between two vectors
	template<typename dataType>
	dataType vectorAngle_deg(dataType (&vector1)[3][1], dataType (&vector2)[3][1]) {
	  dataType unit_vector1[3][1], unit_vector2[3][1];

	  DCM(unit_vector1, vector1);
	  DCM(unit_vector2, vector2);

	  dataType angle = acos( (unit_vector1[0][0] * unit_vector2[0][0]) + (unit_vector1[1][0] * unit_vector2[1][0]) + (unit_vector1[2][0] * unit_vector2[2][0]) ) * 180 / PI;
	  return angle;
	}


	// this function calculates the angle between two vectors
	template<typename dataType>
	dataType vectorAngle_rad(dataType (&vector1)[3][1], dataType (&vector2)[3][1]) {
	  dataType unit_vector1[3][1], unit_vector2[3][1];

	  DCM(unit_vector1, vector1);
	  DCM(unit_vector2, vector2);

	  dataType angle = acos( (unit_vector1[0][0] * unit_vector2[0][0]) + (unit_vector1[1][0] * unit_vector2[1][0]) + (unit_vector1[2][0] * unit_vector2[2][0]) );
	  return angle;
	}


	// this function prints the vector
	template<typename dataType>
	void printVector(dataType (&vector)[3][1]) {
	  Serial.println( String(vector[0][0]) + "i + " + String(vector[1][0]) + "j + " + String(vector[2][0]) + "k" + "\n");
	}

}
