/**
 * Custom exception that gets thrown when the matrix cannot 
 * be created with the given dimensions.
 * @extends Error
 */
class InvalidDimension extends Error {
  
  /**
   * Creates exception that will display what dimension is not 
   * valid and what the reason is.
   * @param {string} argName - The name of the argument.
   * @param {*} value - The value of the argument.
   * @param {string} reason - Why the argument value is invalid.
   */
  constructor(argName, value, reason) {
    var message = `Cannot create matrix with ${value} ${argName} because ${reason}.`;
    super(message);
    this.name = this.constructor.name;
    if (typeof Error.captureStackTrace === 'function') {
      Error.captureStackTrace(this, this.constructor);
    } else { 
      this.stack = (new Error(message)).stack; 
    }
  }
  
}

/**
 * Custom exception that gets thrown when the matrix 
 * has to be square.
 * @extends Error
 */
class NonsquareMatrix extends Error {
  
  /**
   * Creates exception that will display what cannot be done 
   * without a nonsquare matrix.
   * @param {Matrix} matrix - The nonsquare matrix.
   * @param {string} inability - What could not be done (due to nonsquare matrix).
   */
  constructor(matrix, inability) {
    var message = `${inability} from nonsquare matrix. The number of rows and the number of columns must be the same. Your matrix has ${matrix.rows} rows and ${matrix.cols} columns.`;
    super(message);
    this.name = this.constructor.name;
    if (typeof Error.captureStackTrace === 'function') {
      Error.captureStackTrace(this, this.constructor);
    } else { 
      this.stack = (new Error(message)).stack; 
    }
  }
  
}

/**
 * Custom exception that gets thrown when matrices being 
 * added or subtracted have dimensions.
 * @extends Error
 */
class DifferentDimensions extends Error {
  
  /**
   * Creates exception that will display dimensions of the 
   * matrices with disagreeing dimensions.
   * @param {Matrix} matrix1 - The left matrix with dimensions that disagree with the right matrix.
   * @param {Matrix} matrix2 - The right matrix with dimensions that disagree with the left matrix.
   */
  constructor(matrix1, matrix2, op) {
    var message = `Cannot ${op} matrices with different dimensions. The left matrix has ${matrix1.rows} row(s) and ${matrix1.cols} column(s) whereas the right matrix has ${matrix2.rows} row(s) and ${matrix2.cols} column(s)`;
    super(message);
    this.name = this.constructor.name;
    if (typeof Error.captureStackTrace === 'function') {
      Error.captureStackTrace(this, this.constructor);
    } else { 
      this.stack = (new Error(message)).stack; 
    }
  }
  
}

/**
 * Custom exception that gets thrown when matrices being 
 * multiplied together have disagreeing inner dimensions. 
 * In other words, the columns of the left matrix is not the 
 * same as the rows of the right matrix.
 * @extends Error
 */
class InnerDimensions extends Error {
  
  /**
   * Creates exception that will display where the inner 
   * dimensions disagree of two matrices that are being 
   * multiplied together.
   * @param {Matrix} matrix1 - The left matrix with a different 
   * number of columns from he right matrix's rows.
   * @param {Matrix} matrix2 - The right matrix with a different 
   * number of rows as the left matrix's columns.
   */
  constructor(matrix1, matrix2) {
    var message = `Cannot multiply matrices with different inner dimensions. The left matrix has ${matrix1.cols} column(s) whereas the right matrix has ${matrix2.rows} row(s).`;
    super(message);
    this.name = this.constructor.name;
    if (typeof Error.captureStackTrace === 'function') {
      Error.captureStackTrace(this, this.constructor);
    } else { 
      this.stack = (new Error(message)).stack; 
    }
  }
  
}

/**
 * Custom exception for matrices with determinants of 0.
 * @extends Error
 */
class ZeroDeterminant extends Error {
  
  /**
   * Creates exception that will throw for non-invertible matrices.
   */
  constructor() {
    var message = "Cannot invert matrix with determinant of 0!";
    super(message);
    this.name = this.constructor.name;
    if (typeof Error.captureStackTrace === 'function') {
      Error.captureStackTrace(this, this.constructor);
    } else { 
      this.stack = (new Error(message)).stack; 
    }
  }
  
}

/** Class representing a matrix. */
class Matrix {
  
  /**
   * Creates a matrix of dimensions rows x cols. 
   * If the number of columns is not specified, then
   * the matrix will be a square matrix of dimensions
   * rows x rows. 
   * All entries are automatically initialized to 0.
   * @param {Number} rows - The number of rows in the matrix.
   * @param {Number} cols - The number of cols in the matrix.
   */
  constructor(rows, cols) {
    this.rows = rows;
    this.cols = cols ? cols : rows;
    
    // Error handling when creating matrix:
    
    // First, make sure that the dimensions are numbers
    if (isNaN(this.rows))
      throw new InvalidDimension("rows", this.rows, `'${this.rows}' is not a number`);
    if (isNaN(this.cols))
      throw new InvalidDimension("columns", this.cols, `'${this.cols}' is not a number`);
    
    // Since we know that the arguments are numbers,
    // convert the numbers into a readable format
    // (in case of strings)
    this.rows = +`${this.rows}`;
    this.cols = +`${this.cols}`;
    
    // Second, make sure that the dimensions are integers
    if (this.rows !== parseInt(this.rows, 10))
      throw new InvalidDimension("rows", this.rows, `'${this.rows}' is not an integer`);
    if (this.cols !== parseInt(this.cols, 10))
      throw new InvalidDimension("columns", this.cols, `'${this.cols}' is not an integer`);
    
    // Finally, make sure that the dimensions are positive
    if (this.rows < 1)
      throw new InvalidDimension("rows", this.rows, `'${this.rows}' is not positive`);
    if (this.cols < 1)
      throw new InvalidDimension("columns", this.cols, `'${this.cols}' is not positive`);
    
    
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
   * @example <caption>Example of setting all entries of a 2x2
   * matrix of ones to 0.</caption>
   * var m = Matrix.ones(2); // Matrix m: [1 1; 1 1]
   * m.zeros();              // Matrix m: [0 0; 0 0]
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
   * @example <caption>Example of making a zero matrix using 
   * both parameters.</caption>
   * // Matrix: [0 0]
   * Matrix.zeros(1, 2);
   * @example <caption>Example of making a square matrix using 
   * only the first parameter.</caption>
   * // Matrix: [0 0; 0 0]
   * Matrix.zeros(2);
   * @param {Number} rows - The number of rows of zeros.
   * @param {Number} cols - The number of columns of zeros.
   * @return {Matrix} The new matrix with the specified rows of 
   * zeros and specified columns of zeros.
   */
  static zeros(rows, cols) {
    return new Matrix(rows, cols ? cols : rows);
  }
  
  /**
   * Sets all entries of the matrix to 1.
   * @example <caption>Example of setting all entries of a 2x3 
   * zero matrix to 1.</caption>
   * var m = new Matrix(2, 3); // Matrix m: [0 0 0; 0 0 0]
   * m.ones();                 // Matrix m: [1 1 1; 1 1 1]
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
   * @example <caption>Example of making a matrix of ones using 
   * both parameters.</caption>
   * // Matrix: [1 1; 1 1; 1 1]
   * Matrix.ones(3, 2);
   * @example <caption>Example of making a square matrix using 
   * only the first parameter.</caption>
   * // Matrix: [1 1; 1 1]
   * Matrix.ones(2);
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
   * @throws {NonsquareMatrix} Identity matrices can only be created from square matrices.
   */
  identity() {
    
    // Make sure that the matrix is square
    if (this.rows !== this.cols)
      throw new NonsquareMatrix(this, "Unable to create identity matrix");
    
    this.zeros(); // Set all elements of the matrix to 0
    
    // Set the diagonal of zero-matrix to 1
    for (let i = 0; i < this.rows; i++) {
      this.set(i, i, 1);
    }
  }
  
  /**
   * Creates an identity matrix of specified size
   * @example <caption>Example of making an 2x2 identity 
   * matrix.</caption>
   * 
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
   * @throws {DifferentDimensions} Matrices being added 
   * together must have the same dimensions as each other.
   */
  add(matrix) {
    // Make sure that the matrices being added together have
    // the same dimensions
    if (this.rows !== matrix.rows || this.cols !== matrix.cols)
      throw new DifferentDimensions(this, matrix, "add");
    
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
   * @throws {DifferentDimensions} Matrices being added 
   * together must have the same dimensions as each other.
   */
  static add(matrix1, matrix2) {
    // Make sure that the matrices being added together have
    // the same dimensions
    if (matrix1.rows !== matrix2.rows || matrix1.cols !== matrix2.cols)
      throw new DifferentDimensions(matrix1, matrix2, "add");
    
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
   * @throws {DifferentDimensions} Matrices being subtracted 
   * from each other must have the same dimensions as each other.
   */
  subtract(matrix) {
    // Make sure that the matrices being subtracted from 
    // each other the same dimensions
    if (this.rows !== matrix.rows || this.cols !== matrix.cols)
      throw new DifferentDimensions(this, matrix, "subtract");
      
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
   * @throws {DifferentDimensions} Matrices being subtracted 
   * from each other must have the same dimensions as each other.
   */
  static subtract(matrix1, matrix2) {
    // Make sure that the matrices being subtracted from 
    // each other the same dimensions
    if (matrix1.rows !== matrix2.rows || matrix1.cols !== matrix2.cols)
      throw new DifferentDimensions(matrix1, matrix2, "subtract");
      
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
   * @throws {InnerDimensions} The left matrix must have the 
   * same number of columns as the right matrix's rows.
   */
  multiply(matrix) {
    // Make sure that inner dimensions agree
    if (this.cols !== matrix.rows)
      throw new InnerDimensions(this, matrix);
    
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
    // Make sure that inner dimensions agree
    if (matrix1.cols !== matrix2.rows)
      throw new InnerDimensions(matrix1, matrix2);
    
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
   * Scales the current matrix by a given scalar.
   * @param {Number} k - What to scale the matrix by.
   */
  scale(k) {
    this.matrix = this.matrix.map(x => x.map(y => k * y));
  }
  
  /**
   * Scales a given matrix by a given scalar.
   * @param {Matrix} matrix - The matrix to scale.
   * @param {Number} k - What the scale the input matrix by.
   * @return {Matrix} The resulting matrix after scaling the input matrix.
   */
  static scale(matrix, k) {
    // Save a copy of the input matrix
    var scaled = new Matrix(matrix.rows, matrix.cols);
    scaled.matrix = JSON.parse(JSON.stringify(matrix.matrix));
    
    // Scale the copied matrix and return it
    scaled.matrix = scaled.matrix.map(x => x.map(y => k * y));
    return scaled;
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
   * @throws {NonsquareMatrix} Cannot find determinant of nonsquare matrices.
   */
  determinant() {
    // We can only find determinants of square matrices
    if (this.rows !== this.cols)
      throw new NonsquareMatrix(this, "Cannot find determinant");
    
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
   * @throws {ZeroDeterminant} Matrices with determinants of 0 
   * cannot be inverted.
   */
  inverse() {
    var det = this.determinant();
    if (det === 0)
      throw new ZeroDeterminant();
    
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

var m = Matrix.ones(2);
console.log(m);