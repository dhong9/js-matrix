/** Class representing a matrix. */
class Matrix {
  
  /**
   * Creates a matrix of dimensions rows x cols. 
   * If the number of columns is not specified, then
   * the matrix will be a square matrix of dimensions
   * rows x rows.
   * @param {Number} rows - The number of rows in the matrix.
   * @param {Number} cols - The number of cols in the matrix.
   */
  constructor(rows, cols) {
    this.rows = rows;
    this.cols = cols ? cols : rows;
    
    this.matrix = []; // Matrices are represented with arrays
    // Initialize all entries of the matrix to zeros
    for (let i = 0; i < this.rows; i++) {
      let tempArr = [];
      for (let j = 0; j < this.cols; j++) {
        tempArr.push(0);
      }
      this.matrix.push(tempArr);
    }
  }
  
  /**
   * Sets all entries of the matrix to 0.
   */
  zeros() {
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        if (this.matrix[i][j] !== 0) this.matrix[i][j] = 0;
      }
    }
  }
  
  /**
   * Static method that makes new matrix with specified number 
   * of rows of zeros and specified number of columns of zeros. 
   * If the number of columns is not specified, then the matrix 
   * will be a square matrix of dimensions rows x rows.
   * @param {Number} rows - The number of rows of zeros.
   * @param {Number} cols - The number of columns of zeros.
   * @return {Matrix} The new matrix with the specified rows of zeros and specified columns of zeros.
   */
  static zeros(rows, cols) {
    return new Matrix(rows, cols ? cols : rows);
  }
  
  /**
   * Sets all entries of the matrix to 1.
   */
  ones() {
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        if (this.matrix[i][j] !== 1) this.matrix[i][j] = 1;
      }
    }
  }
  
  /**
   * Static method to create a new matrix with specified number 
   * of rows of ones ad specified number of columns of ones.
   * If the number of columns is not specified, then the matrix 
   * will be a square matrix of dimensions rows x rows.
   * @param {Number} rows = The number of rows of ones.
   * @param {Number} cols - The number of columns of ones.
   * @return {Matrix} The new matrix with the specified rows of ones and specified columns of ones.
   */
  static ones(rows, cols) {
    var matrix = new Matrix(rows, cols ? cols : rows);
    matrix.ones();
    return matrix;
  }
  
  /**
   * Turns the current matrix into an identity matrix.
   */
  identity() {
    this.zeros(); // Set all elements of the matrix to 0
    
    // Set the diagonal of zero-matrix to 1
    for (let i = 0; i < this.rows; i++) {
      this.set(i, i, 1);
    }
  }
  
  /**
   * Creates an identity matrix of specified size
   * @param {Number} size - The size of the identity matrix.
   * @return {Matrix} The identity matrix of given size.
   */
  static identity(size) {
    var I = new Matrix(size); // Create new zero matrix
    
    // Set the diagonal of new matrix to 1
    for (let i = 0; i < size; i++) {
      I.set(i, i, 1);
    }
    
    return I;
  }
  
  /**
   * Transposes the current matrix.
   */
  transpose() {
    var trans = []; // Array to store transpose values in
    
    // Loop to apply the definition of transpose
    // a[j, i] = a[i, j]
    for (let i = 0; i < this.cols; i++) {
      let temp = [];
      for (let j = 0; j < this.rows; j++) {
        temp.push(this.get(j, i));
      }
      trans.push(temp);
    }
    this.matrix = trans;
    [this.rows, this.cols] = [this.cols, this.rows];
  }
  
  /**
   * Static method to transpose a given matrix.
   * @param {Matrix} matrix - The matrix to take the transpose of.
   * @return {Matrix} The resulting matrix when the given matrix is transposed.
   */
  static transpose(matrix) {
    var transMatrix = new Matrix(matrix.cols, matrix.rows); // Matrix that is the transpose of input matrix
    var trans = []; // Array to store transpose values in
    
    // Loop to apply the definition of transpose 
    // a[j, i] = a[i, j]
    for (let i = 0; i < matrix.cols; i++) {
      let temp = [];
      for (let j = 0; j < matrix.rows; j++) {
        temp.push(matrix.get(j, i));
      }
      trans.push(temp);
    }
    transMatrix.matrix = trans;
    return transMatrix;
  }
  
  /**
   * Adds a given matrix to the current matrix.
   * @param {Matrix} matrix - The matrix to add to the current matrix.
   */
  add(matrix) {
    // Loop to add the matrices component-wise
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        this.matrix[i][j] += matrix.get(i, j);
      }
    }  
  }
  
  /**
   * Static method to create a new matrix that is the sum 
   * of the two given matrices. The order that the matrices 
   * does not matter because matrix addition is commutative.
   * @param {Matrix} matrix1 - One of the matrices to add.
   * @param {Matrix} matrix2 - Another matrix to add.
   * @return {Matrix} The sum of the two given matrices.
   */
  static add(matrix1, matrix2) {
    var rows = matrix1.rows, cols = matrix1.cols; // Get number of  rows and columns
    var sum = new Matrix(rows, cols); // New matrix to store the sum in
    
    // Loop to add the matrices component-wise
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        sum.set(i, j, matrix1.get(i, j) + matrix2.get(i, j));
      }
    }
    
    return sum;
  }
  
  /**
   * Subtracts the given matrix from the current matrix
   * @param {Matrix} matrix - The matrix to subtract the current matrix from.
   */
  subtract(matrix) {
    // Loop to subtract the matrices component-wise
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        this.matrix[i][j] -= matrix[i][j];
      }
    }  
  }
  
  /**
   * Static method for subtracting one matrix from another. 
   * Note that the order that the matrices are given in 
   * matters because subtraction is not commutative.
   * @param {Matrix} matrix1 - The matrix to be subtracted from.
   * @param {Matrix} matrix2 - The matrix to subtract the first matrix from.
   * @return {Matrix} The result of the second matrix subtracted from the first matrix.
   */
  static subtract(matrix1, matrix2) {
    var rows = matrix1.rows, cols = matrix1.cols; // Get number of  rows and columns
    var diff = new Matrix(rows, cols); // New matrix to store the difference in
    
    // Loop to subtract the matrices component-wise
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        diff.set(i, j, matrix1.get(i, j) - matrix2.get(i, j));
      }
    }
    return matrix1;
  }
  
  /**
   * Multiplies the current matrix by another matrix 
   * row by column.
   * @param {Matrix} matrix - The matrix to multiply the current matrix by.
   */
  multiply(matrix) {
    var product = []; // Array to store product elements in
    
    // Loop to multiply matrices row by column
    for (let i = 0; i < this.rows; i++) {
      let arr = [];
      for (let j = 0; j < matrix.cols; j++) {
        let prodSum = 0;
        for (let k = 0; k < this.cols; k++) {
          prodSum += this.get(i, k) * matrix.get(k, j);
        }
        arr.push(prodSum);
      }
      product.push(arr);
    }
    
    this.matrix = product; // Update current matrix with new array representation
    this.cols = matrix.cols; // The number of columns has changed
  }
  
  /**
   * Static method to multiply 2 matrices. 
   * Note that the order that the matrices are given in matters 
   * because matrix multiplication is not commuative.
   * @param {Matrix} matrix1 - The left matrix.
   * @param {Matrix} matrix2 - The right matrix.
   * @return {Matrix} The result of the left matrix multiplied by the right matrix.
   */
  static multiply(matrix1, matrix2) {
    var prodMatrix = new Matrix(matrix1.rows, matrix2.cols); // Initialize resulting matrix
    var product = []; // Array to store product in
    
    // Loop to multiply matrices row by column
    for (let i = 0; i < matrix1.rows; i++) {
      let arr = [];
      for (let j = 0; j < matrix2.cols; j++) {
        let prodSum = 0;
        for (let k = 0; k < matrix1.cols; k++) {
          prodSum += matrix1.get(i, k) * matrix2.get(k, j);
        }
        arr.push(prodSum);
      }
      product.push(arr);
    }
    
    prodMatrix.matrix = product; // Update resulting matrix's array
    return prodMatrix;
  }
  
  /**
   * Scales a matrix by a given scalar.
   * @param {Number} k - What to scale the matrix by.
   */
  scale(k) {
    this.matrix = this.matrix.map(x => x.map(y => k * y));
  }
  
  /**
   * Finds the matrix minor at a given row and column of the matrix. 
   * The current matrix will be modified to the calculated matrix minor.
   * @param {Number} row - The row to remove.
   * @param {Number} col - The column to remove.
   */
  matrixMinor(row, col) {
    this.matrix.splice(row, 1); // Removes input row
    this.rows--;
    
    // Removes input column
    for (let i = 0; i < this.rows; i++) {
      this.matrix[i].splice(col, 1);
    }
    this.cols--;
  }
  
  /**
   * Static method for finding the matrix minor for a 
   * given matrix at its specified row and column.
   * @param {Matrix} matrix - The matrix to find the matrix minor of.
   * @param {Number} row - The row to remove.
   * @param {Number} col - The column to remove.
   * @return {Matrix} The minor of the input matrix at the given row and column.
   */
  static matrixMinor(matrix, row, col) {
    // Copy input matrix into another matrix that will be modified
    var minor = new Matrix(matrix.rows, matrix.cols);
    minor.matrix = JSON.parse(JSON.stringify(matrix.matrix));
    
    // Find matrix minor of copied matrix
    minor.matrixMinor(row, col);
    return minor;
  }
  
  /**
   * Finds the determinant of the current matrix.
   * @return {Number} The determinant of the matrix.
   */
  determinant() {
    // If the matrix is 1x1...
    if (this.rows === 1) return this.get(0, 0);
    
    // If the matrix is 2x2...
    if (this.rows === 2) 
      return this.get(0, 0) * this.get(1, 1) - this.get(0, 1) * this.get(1, 0);
    
    // General case
    var det = 0;
    for (let i = 0; i < this.rows; i++) {
      det += Math.pow(-1, i) * this.get(0, i) * Matrix.matrixMinor(this, 0, i).determinant();
    }
    return det;
  }
  
  /**
   * Finds the determinant of a specified matrix.
   * @param {Matrix} matrix - The matrix to find the determinant of.
   * @return {Number} The determinant of the input matrix.
   */
  static determinant(matrix) {
    return matrix.determinant();
  }
  
  /**
   * Inverts the current matrix
   */
  inverse() {
    var det = this.determinant();
    
    // Special case for 1x1 matrix
    if (this.rows === 1) {
      this.set(0, 0, 1 / get(0, 0));
    }
    
    // Special case for 2x2 matrix
    else if (this.rows === 2) {
      [this.matrix[0][0], this.matrix[1][1]] = [this.matrix[1][1], this.matrix[0][0]];
      this.matrix[0][1] *= -1;
      this.matrix[1][0] *= -1;
      this.scale(1 / det);
    }
    
    // General case
    else {
      let cofactors = [];
      
      for (let r = 0; r < this.rows; r++) {
        let cofactorRow = [];
        for (let c = 0; c < this.rows; c++) {
          let minor = Matrix.matrixMinor(this, r, c);
          cofactorRow.push(Math.pow(-1, r + c) * minor.determinant());
        }
        cofactors.push(cofactorRow);
      }
      
      this.matrix = cofactors;
      this.transpose();
      this.scale(1 / det);
    }
  }
  
  /**
   * Static method for computing the inverse of a given matrix.
   * @param {Matrix} matrix - The matrix to find the inverse of.
   * @return {Matrix} matrix - The inverse of the input matrix.
   */
  static inverse(matrix) {
    // Copy input matrix into another matrix that will be modified
    var inv = new Matrix(matrix.rows, matrix.cols);
    inv.matrix = JSON.parse(JSON.stringify(matrix.matrix));
    
    // Find inverse of copied matrix
    inv.inverse();
    return inv;
  }
  
  /**
   * Gets the element of the matrix at the given row and column
   * @param {Number} row - The row of the matrix location.
   * @param {Number} col - The column of the matrix location.
   * @return {Number} The matrix element at the specified location.
   */
  get(row, col) {
    return this.matrix[row][col];
  }
  
  /**
   * Sets the matrix at the given row and column to the given value.
   * @param {Number} row - The row of the matrix location.
   * @param {Number} col - The column of the matrix location.
   * @param {Number} val - The value to set the matrix at the specified location to.
   */
  set(row, col, val) {
    this.matrix[row][col] = val;
  }
  
  /**
   * Converts the matrix object to a readable string format.
   * @return {string} The string representation of the matrix.
   */
  toString() {
    var str = "[";
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        str += `${this.matrix[i][j]} `;
      }
      str = str.slice(0, -1);
      str += "\n ";
    }
    str = str.slice(0, -2);
    str += "]";
    return str;
  }
}