#include <iostream>
#include "windows.h"
#include <vector>
#include "string"
#include <cmath>

using namespace std;

class Matrix {
private:
	vector<vector<double>> matrix;
	vector<int> width_form;
	int m = 0, n = 0;

	void Format() {
		width_form.clear();
		int width, buf_width, width_null;
		string str_width;
		width = 1;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				str_width = to_string(matrix[j][i]);
				for (int i = str_width.size() - 1; i >= 0; i--) {
					if (str_width[i] == '0' || str_width[i] == ',') {
						str_width.erase(i);
					}
					else break;
				}
				buf_width = str_width.size();
				if (width < buf_width) width = buf_width;
			}
			width_form.push_back(width);
		}
	}
	void Size() {
		m = matrix.size();
		n = matrix[0].size();
		for (int i = 1; i < m; i++) {
			if (matrix[i].size() != n) {
				n = 0;
				break;
			}
		}
	}
public:
	void Input_Matrix(int m, int n)
	{
		double number;
		vector<double> str_matrix;
		cout << "Enter the [" << m << "," << n << "] matrix:" << endl;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cin >> number;
				str_matrix.push_back(number);
			}
			matrix.push_back(str_matrix);
			str_matrix.clear();
		}
	}
	void Output_Matrix() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	void Output_SLAU() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
				if (j == n - 2) cout << "| ";
			}
			cout << endl;
		}
		cout << endl;
	}
	vector<vector<double>> get_Matrix() {
		return matrix;
	}
	void set_SLAU(vector<vector<double>> slau) {
		matrix = slau;
	}
};

class Determinant {
private:
	vector<vector<double>> matrix_det;
	double det = 0;
	Matrix M;

	double Accur(double n) {
		double buf = abs(n);
		if (buf != 0) {
			if (fmod(buf, 0.001) < 1) {
				n = round(n * 1000) / 1000;
			}
		}
		return n;
	}
	double Det_2x(vector<vector<double>> matrix_2x) {
		double right = matrix_2x[0][0] * matrix_2x[1][1];
		double left = matrix_2x[0][1] * matrix_2x[1][0];
		return Accur(right - left);
	}
	double Det_3x(vector<vector<double>> matrix_3x) {
		double right = 0;
		double buf = 1, buf1 = 1;
		for (int i = 0; i < 3; i++) {
			buf *= matrix_3x[i][i];
		}
		right += buf;
		buf = 1;
		for (int i = 0; i < 3; i++) {
			if (i + 1 < 3) {
				buf *= matrix_3x[i][i + 1];
				buf1 *= matrix_3x[i + 1][i];
			}
			else {
				buf *= matrix_3x[i][0];
				buf1 *= matrix_3x[0][i];
			}
		}
		right += buf;
		right += buf1;
		double left = 0;
		buf = matrix_3x[0][2] * matrix_3x[2][0] * matrix_3x[1][1];
		left += buf;
		buf = matrix_3x[0][1] * matrix_3x[1][0] * matrix_3x[2][2];
		left += buf;
		buf = matrix_3x[1][2] * matrix_3x[2][1] * matrix_3x[0][0];
		left += buf;
		return Accur(right - left);
	}
	double Det_mx(vector<vector<double>> matrix_mx) {
		if (matrix_mx.size() == 2) return Det_2x(matrix_mx);
		else if (matrix_mx.size() == 3) return Det_3x(matrix_mx);
		else {
			double det_mx = 0;
			for (int i = 0; i < matrix_mx.size(); i++) {
				vector<vector<double>> buf_m = matrix_mx;
				buf_m.erase(buf_m.begin());
				int buf_size = buf_m.size();
				for (int buf_index = 0; buf_index < buf_size; buf_index++) {
					buf_m[buf_index].erase(buf_m[buf_index].begin() + i);
				}
				if (i % 2 == 0) {
					det_mx += matrix_mx[0][i] * Det_mx(buf_m);
				}
				else det_mx -= matrix_mx[0][i] * Det_mx(buf_m);
			}
			return Accur(det_mx);
		}
	}
public:
	void Input_Det_Matrix(vector<vector<double>> matrix) {
		matrix_det = matrix;
	}
	void Decision() {
		det = Det_mx(matrix_det);
		cout << "det(A) = " << det << endl;
	}
	double get_det() {
		return det;
	}
};

class SLAU {
private:
	vector<vector<double>> matrix_slau, A, B, T, T_, D;
	vector<double> Y, X;
	Matrix M;
	int m;
	void Accur_str(int str_index) {
		double buf;
		vector<double> str_slau = matrix_slau[str_index];
		for (int i = 0; i < str_slau.size(); i++) {
			buf = abs(str_slau[i]);
			if (buf != 0) {
				if (fmod(buf, 0.001) < 1) {
					str_slau[i] = round(str_slau[i] * 1000) / 1000;
				}
			}
		}
		matrix_slau[str_index] = str_slau;
	}
	double Accur(double n) {
		double buf = abs(n);
		if (buf != 0) {
			if (fmod(buf, 0.001) < 1) {
				n = round(n * 1000) / 1000;
			}
		}
		return n;
	}
	void Output_SLAU() {
		M.set_SLAU(matrix_slau);
		M.Output_SLAU();
	}
	void X_Record(int number) {
		double buf_x = 0;
		for (int i = 0; i < matrix_slau.size(); i++) {
			if (i != number) {
				if (abs(matrix_slau[number][i]) >= 0.001) {
					buf_x += (matrix_slau[number][i] * X[matrix_slau.size() - i - 1]);
				}
			}
		}
		buf_x -= matrix_slau[number][matrix_slau.size()];
		buf_x /= -(matrix_slau[number][number]);
		buf_x = Accur(buf_x);
		X.push_back(buf_x);
	}
	void Y_Record(int number) {
		double buf_x = 0;
		for (int i = 0; i < matrix_slau.size(); i++) {
			if (i != number) {
				if (abs(matrix_slau[number][i]) >= 0.001) {
					buf_x += (matrix_slau[number][i] * Y[i]);
				}
			}
		}
		buf_x -= matrix_slau[number][matrix_slau.size()];
		buf_x /= -(matrix_slau[number][number]);
		buf_x = Accur(buf_x);
		Y.push_back(buf_x);
	}
	double sign(double x) {
		if (x > 0) return 1;
		if (x < 0) return -1;
		else return 0;
	}
	double dii(int i_index) {
		if (i_index <= 0 || i_index >= m) return 0;
		double sum = 0;
		for (int i = 0; i < i_index; i++) {
			sum += (T[i][i_index] * T[i][i_index] * T[i][i]);
		}
		double buf = A[i_index][i_index] - sum;
		return Accur(sign(buf));
	}
	double t11() {
		if (A[0][0] < 0) return Accur(sqrt(-A[0][0]));
		return Accur(sqrt(A[0][0]));
	}
	double t1j(int j_index) {
		if (j_index > 0 && j_index < m) return Accur(A[0][j_index] / T[0][0]);
		else return 0;
	}
	double tii(int i_index) {
		if (i_index <= 0 || i_index >= m) return 0;
		double sum = 0;
		for (int i = 0; i < i_index; i++) {
			sum += (T[i][i_index] * T[i][i_index] * D[i][i]);
		}
		double buf = A[i_index][i_index] - sum;
		if(buf < 0) return Accur((sqrt(-buf)));
		return Accur(sqrt(buf));
	}
	double tij(int i_index, int j_index) {
		if (i_index <= 0 || j_index >= m || j_index <= i_index) return 0;
		double sum = 0;
		double d = dii(i_index);
		for (int i = 0; i < i_index; i++) {
			sum += (T[i][i_index] * T[i][j_index] * D[i][i]);
		}
		double buf = A[i_index][j_index] - sum;
		return Accur(buf / (T[i_index][i_index] * d));
	}
	void Create_T() {
		T[0][0] = t11();
		D[0][0] = sign(A[0][0]);
		for (int i = 1; i < m; i++) {
			T[0][i] = t1j(i);
		}
		for (int i = 1; i < m; i++) {
			D[i][i] = dii(i);
			T[i][i] = tii(i);
			for (int j = m - 1; j > i; j--) {
				T[i][j] = tij(i, j);
			}
		}
	}
	void Create_T_() {
		T_ = T;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				T_[i][j] = T[j][i];
			}
		}
	}
	bool Test() {
		bool test = true;
		vector<vector<double>> Buf = A;
		while (true) {
			Buf.erase(Buf.end() - 1);
			for (int i = 0; i < Buf.size(); i++) {
				Buf[i].erase(Buf[i].end() - 1);
			}
			if (Buf.size() == 1) 
				break;
			Determinant Det;
			Det.Input_Det_Matrix(Buf);
			Det.Decision();
			if (Det.get_det() <= 0) {
				test = false;
				break;
			}
		}
		return test;
	}
public:
	bool Input_SLAU(vector<vector<double>> matrix_A, vector<vector<double>> matrix_B) {
		m = matrix_A.size();
		A = matrix_A, B = matrix_B;
		if (m != matrix_B.size() || matrix_B[0].size() != 1) {
			return false;
		}
		vector<double> str_matrix;
		for (int i = 0; i < m; i++) {
			str_matrix.push_back(0);
		}
		for (int i = 0; i < m; i++) {
			T.push_back(str_matrix);
		}
		D = T;
		return Test();
	}
	void Decision() {
		Create_T();
		Create_T_();
		matrix_slau = T_;
		for (int i = 0; i < m; i++) {
			matrix_slau[i].push_back(B[i][0]);
		}
		Output_SLAU();
		for (int i = 0; i < matrix_slau.size(); i++) {
			Y_Record(i);
		}
		matrix_slau = T;
		for (int i = 0; i < m; i++) {
			matrix_slau[i].push_back(Y[i]);
		}
		Output_SLAU();
		for (int i = matrix_slau.size() - 1; i >= 0; i--) {
			X_Record(i);
		}
		reverse(X.begin(), X.end());
	}
	void Rezult() {
		cout << "Answer: ";
		for (int i = 0; i < X.size(); i++) {
			cout << "x" << i + 1 << " = " << X[i] << ", ";
		}
		cout << endl;
	}
	
};

class Task {
private:
	Matrix M;
	Matrix M1;
	SLAU Slau;
	Determinant Det;
	int m, n;
	void Create_SLAU() {
		M.Input_Matrix(m, m);
		Det.Input_Det_Matrix(M.get_Matrix());
		n = 1;
		M1.Input_Matrix(m, n);
	}
public:
	void SLAU_Task() { //почему не работают отрицательные числа???
		cout << "Input the size of the matrix." << endl << "m = ";
		cin >> m;
		if (m < 2 || m > 49) {
			cout << "Incorrect size of the matrix." << endl;
			return;
		}
		Create_SLAU();
		Det.Decision();
		if (Det.get_det() == 0) {
			cout << "Error! The matrix det = 0." << endl;
			return;
		}
		M.Output_Matrix();
		M1.Output_Matrix();
		bool flag = Slau.Input_SLAU(M.get_Matrix(), M1.get_Matrix());
		if (!flag) {
			cout << "Matrix isn't positiv." << endl;
			return;
		}
		Slau.Decision();
		Slau.Rezult();
	}
};

int main() {
	setlocale(LC_ALL, "Rus");
	SetConsoleCP(1251);
	int ans, exit = 1;
	while (exit == 1) {
		Task T;
		cout << "1.SLAU" << endl << "2.Exit" << endl << "Choose a way:" << endl;
		cin >> ans;
		switch (ans)
		{
		case 1:
			T.SLAU_Task();
			break;
		case 2:
			exit = 0;
			break;
		default:
			cout << "This task does not exist" << endl;
			break;
		}
	}
	system("pause");
	return 0;
}

/*
0.32 -0.42 0.85
0.63 -1.43 -0.58
0.84 -2.23 -0.52
1.32
-0.44
0.64

1 2 4
2 13 23
4 23 77
10
50
150

2 1 4
1 1 3
4 3 14
16
12
52

3.78 1.08 -1.35
1.08 -2.28 0.37
-1.35 0.37 2.86
0.35
1.27
0.47

3.23 1.62 0.65  
1.62 2.33 -0.43
0.65 -0.43 2.16
1.28
0.87
-2.87


¬ариант 11
2.74 -1.18 0.23
-1.18 2.71 -0.52
0.23 -0.52 1.62
0.16
1.81
-1.25

*/