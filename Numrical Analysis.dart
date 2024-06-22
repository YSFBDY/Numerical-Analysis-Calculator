import 'dart:math';

//<<<<<< CHAPTER 1 >>>>>>//

//__Bisection Method__//

double fB(double x) {
  return 4 * pow(x, 3) - 6 * pow(x, 2) + 7 * x - 2.3; //user input
}

double bisect(double xl, double xu) {
  double eps = 1.0; //user input
  double xr = 0;
  double xrOld = 0;
  double error = 0;
  int iter = 0;

  if (fB(xl) * fB(xu) < 0) {
    do {
      xrOld = xr;
      xr = (xl + xu) / 2;
      error = (xr - xrOld).abs() / xr * 100;

      if (iter == 0) {
        print(
            "iteration: $iter | xl= ${xl.toStringAsFixed(3)} | f(xl)= ${fB(xl).toStringAsFixed(3)} | xu= ${xu.toStringAsFixed(3)} | f(xu)= ${fB(xu).toStringAsFixed(3)} | xr= ${xr.toStringAsFixed(3)} | f(xr)= ${fB(xr).toStringAsFixed(3)}");
      } else {
        print(
            "iteration: $iter | xl= ${xl.toStringAsFixed(3)} | f(xl)= ${fB(xl).toStringAsFixed(3)} | xu = ${xu.toStringAsFixed(3)} | f(xu) = ${fB(xu).toStringAsFixed(3)} | xr = ${xr.toStringAsFixed(3)} | f(xr) = ${fB(xr).toStringAsFixed(3)} | error = ${error.toStringAsFixed(3)} %");
      }

      double m = fB(xl) * fB(xr);
      if (m > 0) {
        xl = xr;
      } else if (m == 0) {
        return xr;
      } else {
        xu = xr;
      }
      iter++;
    } while (error > eps);

    return xr;
  } else {
    print("Not correct xl and xu");
    return 0;
  }
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//__False Position Method__//

double fF(double x) {
  return -26 +
      82.3 * x -
      88 * pow(x, 2) +
      45.4 * pow(x, 3) -
      9 * pow(x, 4) +
      0.65 * pow(x, 5); //user input
}

double falseP(double xl, double xu) {
  double eps = 0.2; // user input
  double xr = 0;
  double xrOld = 0;
  double error = 0;
  int iter = 0;

  if (fF(xl) * fF(xu) < 0) {
    do {
      xrOld = xr;
      xr = xu - (fF(xu) * (xl - xu)) / (fF(xl) - fF(xu));
      error = (xr - xrOld).abs() / xr * 100;

      if (iter == 0) {
        print(
            "iteration: $iter | xl= ${xl.toStringAsFixed(3)} | f(xl)= ${fF(xl).toStringAsFixed(3)} | xu= ${xu.toStringAsFixed(3)} | f(xu)= ${fF(xu).toStringAsFixed(3)} | xr= ${xr.toStringAsFixed(3)} | f(xr)= ${fF(xr).toStringAsFixed(3)}");
      } else {
        print(
            "iteration: $iter | xl= ${xl.toStringAsFixed(3)} | f(xl)= ${fF(xl).toStringAsFixed(3)} | xu = ${xu.toStringAsFixed(3)} | f(xu) = ${fF(xu).toStringAsFixed(3)} | xr = ${xr.toStringAsFixed(3)} | f(xr) = ${fF(xr).toStringAsFixed(3)} | error = ${error.toStringAsFixed(3)} %");
      }

      double m = fF(xl) * fF(xr);
      if (m > 0) {
        xl = xr;
      } else if (m == 0) {
        return xr;
      } else {
        xu = xr;
      }
      iter++;
    } while (error > eps);

    return xr;
  } else {
    print("Not correct xl and xu");
    return 0;
  }
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//__Simple Fixedpoint Method__//

double fs(double x) {
  return sqrt(((1.7 * x) + 2.5) / 0.9); //user input
}

double fixedpoint(double x) {
  double eps = 0.7; // user input
  int iter = 0;
  double xiPlus1 = 0;
  double xi = x;
  double error = 0;

  do {
    xiPlus1 = fs(xi);
    error = (100 * (xiPlus1 - xi) / xiPlus1).abs();

    print(
        "i=$iter | Xi= ${xi.toStringAsFixed(3)} | Xi+1= ${xiPlus1.toStringAsFixed(3)} | E= ${error.toStringAsFixed(3)}%");

    xi = xiPlus1; // Update xi for the next iteration

    iter++;
  } while (error > eps);

  return xi;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//__Newton Method__//

double fN(double x) {
  return -0.9 * pow(x, 2) + 1.7 * x + 2.5; //user input
}

double numericalDerivative(double Function(double) f, double x, double h) {
  return (f(x + h) - f(x)) / h;
}

double newton(double xo, double h) {
  double eps = 0.7; //user input
  double error = 0;
  int iter = 0;
  double xiPlus1 = 0;
  double xi = xo;

  do {
    print(
        "i=$iter | Xi= ${xi.toStringAsFixed(3)} | F(Xi)= ${fN(xi).toStringAsFixed(3)} | F'(Xi)= ${numericalDerivative(fN, xi, h).toStringAsFixed(3)} | E= ${error.toStringAsFixed(3)}%");

    xiPlus1 = xi - (fN(xi) / numericalDerivative(fN, xi, h));

    error = (100 * (xiPlus1 - xi) / xiPlus1).abs();

    xi = xiPlus1;
    iter++;
  } while (error > eps);
  print(
      "i=$iter | Xi= ${xi.toStringAsFixed(3)} | F(Xi)= ${fN(xi).toStringAsFixed(3)} | F'(Xi)= ${numericalDerivative(fN, xi, h).toStringAsFixed(3)} | E= ${error.toStringAsFixed(3)}%");

  return xi;
}


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//__Secant Method__//

double fSC(double x) {
  return 0.95 * pow(x, 3) - 5.9 * pow(x, 2) + 10.9 * x - 6; //user input
}

double secant(double xu, double xl) {
  double eps = 0.5; //user input
  double error = 0;
  int iter = 0;
  double Xi = xu;
  double xiMin1 = xl;
  double xiPlus1 = 0;

  do {
    print(
        "i=$iter | Xi-1=${xiMin1.toStringAsFixed(3)} | f(Xi-1)= ${fSC(xiMin1).toStringAsFixed(3)} | Xi = ${Xi.toStringAsFixed(3)} | f(Xi) ${fSC(Xi).toStringAsFixed(3)} | Error  ${error.toStringAsFixed(3)}%");

    xiPlus1 = Xi - (fSC(Xi) * (xiMin1 - Xi)) / (fSC(xiMin1) - fSC(Xi));
    error = (100 * (xiPlus1 - Xi) / xiPlus1).abs();
    iter++;
    xiMin1 = Xi;
    Xi = xiPlus1;
  } while (error > eps);
  print(
      "i=$iter | Xi-1=${xiMin1.toStringAsFixed(3)} | f(Xi-1)= ${fSC(xiMin1).toStringAsFixed(3)} | Xi = ${Xi.toStringAsFixed(3)} | f(Xi) ${fSC(Xi).toStringAsFixed(3)} | Error  ${error.toStringAsFixed(3)}%");

  return Xi;
}

//_________________________________________________________________________________________________________//

//<<<<<< CHAPTER 2 >>>>>>//

double m21 = 0;
double m31 = 0;
double m32 = 0;

void display(List<List<double>> matrix) {
  for (int i = 0; i < matrix.length; i++) {
    String row = "[ ";
    for (int j = 0; j < matrix[i].length; j++) {
      row += "${matrix[i][j].toStringAsFixed(0)} ";
    }
    row += "]";
    print(row);
  }
  print("");
}

void copyMatrix(List<List<double>> src, List<List<double>> dest) {
  for (int i = 0; i < src.length; i++) {
    for (int j = 0; j < src[i].length; j++) {
      dest[i][j] = src[i][j];
    }
  }
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//__Gaussian Elimination__//

void gje(List<List<double>> matrix, List<double> mValues) {
  display(matrix);
  m21 = matrix[1][0] / matrix[0][0];
  m31 = matrix[2][0] / matrix[0][0];

  for (int j = 0; j < matrix[0].length; j++) {
    double e2 = matrix[1][j];
    double e1 = m21 * matrix[0][j];
    matrix[1][j] = e2 - e1;
  }
  display(matrix);

  for (int j = 0; j < matrix[0].length; j++) {
    double e3 = matrix[2][j];
    double e1 = m31 * matrix[0][j];
    matrix[2][j] = e3 - e1;
  }

  m32 = matrix[2][1] / matrix[1][1];

  for (int j = 0; j < matrix[0].length; j++) {
    double e3 = matrix[2][j];
    double e1 = m32 * matrix[1][j];
    matrix[2][j] = e3 - e1;
  }

  display(matrix);

  double x3 = matrix[2][3] / matrix[2][2];
  double x2 = (matrix[1][3] - (matrix[1][2] * x3)) / matrix[1][1];
  double x1 = (matrix[0][3] - ((matrix[0][1] * x2) + (matrix[0][2] * x3))) /
      matrix[0][0];

  print("Gauss Jordan Result\nX1 = $x1\nX2 = $x2\nX3 = $x3");
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//__LU Decomposition __//

void luDecomposition(List<List<double>> matrix) {
  List<List<double>> u = List.generate(3, (_) => List.generate(3, (_) => 0.0));
  List<double> b = List.generate(3, (_) => 0.0);
  List<double> mValues = [0, 0, 0];
  for (int i = 0; i < 3; i++) {
    b[i] = matrix[i][3];
  }
  gje(matrix, mValues);

  List<List<double>> l = List.generate(3, (_) => List.generate(3, (_) => 0.0));
  for (int i = 0; i < 3; i++) {
    l[i][i] = 1; // Set diagonal elements of L to 1
  }
  l[1][0] =
      m21; // Set the elements below the diagonal to the multipliers used during elimination
  l[2][0] = m31;
  l[2][1] = m32;

  print("\nU matrix");
  display(matrix); // U matrix is the result of Gaussian Jordan Elimination

  print("L matrix");
  display(l);

  // Solve Ly = b for y
  List<double> y = [0, 0, 0];
  y[0] = b[0] / l[0][0];
  y[1] = (b[1] - l[1][0] * y[0]) / l[1][1];
  y[2] = (b[2] - l[2][0] * y[0] - l[2][1] * y[1]) / l[2][2];

  // Solve Ux = y for x
  List<double> x = [0, 0, 0];
  x[2] = y[2] / matrix[2][2];
  x[1] = (y[1] - matrix[1][2] * x[2]) / matrix[1][1];
  x[0] = (y[0] - matrix[0][1] * x[1] - matrix[0][2] * x[2]) / matrix[0][0];

  print("LU decomposition\nX1 = ${x[0]}\nX2 = ${x[1]}\nX3 = ${x[2]}");
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//__Cramer’s rule __//

void cramersRule(List<List<double>> matrix) {
  List<List<double>> tempMatrix =
      List.generate(3, (_) => List.generate(3, (_) => 0.0));
  List<double> detA = [0, 0, 0];
  double det = 0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      tempMatrix[i][j] = matrix[i][j];
    }
  }
  det = (matrix[0][0] *
          (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])) -
      (matrix[0][1] *
          (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])) +
      (matrix[0][2] *
          (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]));
  print("Determinant of A = $det");

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix[j][i] = matrix[j][3];
    }

    detA[i] = ((matrix[0][0] *
                (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])) -
            (matrix[0][1] *
                (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])) +
            (matrix[0][2] *
                (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))) /
        det;
    copyMatrix(tempMatrix, matrix);
    print("A[${i + 1}] = ${detA[i]}");
  }

  for (int i = 0; i < 3; i++) {
    print("x${i + 1} = ${detA[i]}");
  }
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void main() {
  //<<<<<< CHAPTER 1 >>>>>>//

  double xl = 3.5; //user input      (bisection & false & secant)
  double xu = 2.5; //user input      (bisection & false & secant)
  double xo = 5; //user input      (simple fixed point & Newton)

  //bisection print
  //print("Root = ${bisect(xl, xu).toStringAsFixed(3)}");

  //false postion print
  //print("Root = ${falseP(xl, xu).toStringAsFixed(3)}");

  //simple fixed point print
 // print("Root = ${fixedpoint(xo).toStringAsFixed(3)}");

  //Newton print
  //print("Root = ${newton(xo,0.0001).toStringAsFixed(3)}");

  //Secant Method
  //print("Root = ${secant(xl, xu).toStringAsFixed(3)}");

///////////////////////////////////////////////////////////////

//<<<<<< CHAPTER 2 >>>>>>//

  List<List<double>> matrix = [
    [2, 1, 1, 8],
    [4, 1, 0, 11],
    [-2, 2, 1, 3]
  ]; //user input

  List<double> mValues = [0, 0, 0];

  //Gaussian Elimination print
  //gje(matrix, mValues);

  //LU Decomposition print
  //luDecomposition(matrix);

  //Cramer’s rule print
  //cramersRule(matrix);
}
