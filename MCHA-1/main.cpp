#include <iostream>
#include <cmath>
#include <vector>

double epsilon = 0.001;

double f(double x, double y) {
    return -y * cos(x) + sin(x) * cos(x);
}

double fs(double x, double y) {
    return 2 * f(x, y) * cos(x) + y * pow(cos(x), 2) - sin(x) * pow(cos(x), 2);
}

double F(double x1, double y0, double y1, double h) {
    return y1 - y0 - h * f(x1, y1);
}

double Fs(double x1, double y1, double h) {
    return 1 + h * f(x1, y1) * cos(x1);
}

std::vector <double> methodRungeKutta(double a, double b, double h, double y0) {
    std::vector <double> y(11);
    double fi_0 = 0.0;
    double fi_1 = 0.0;
    y[0] = y0;
    for (size_t j = 0; j < y.size(); ++j) {
        fi_0 = h * f(a + h * j, y[j]);
        fi_1 = h * f(a + h * (j + 0.5), y[j] + 0.5 * fi_0);
        y[j + 1] = y[j] + fi_1;
    }
    
    std::cout << "Метод Рунге-Кутта:\n";
    
    for (int i = 0; i < y.size(); ++i) {
        std::cout << "y(" << a + h * i << ") = " << y[i] << std::endl;
    }
    
    std::cout << "\n\n";
    return y;
}

std::vector <double> methodPPPT(double a, double b, double h, double y0) {
    std::vector <double> y(11);
    y[0] = y0;
    for (size_t j = 0; j < y.size(); ++j) {
        y[j + 1] = y[j] + h * f(a + h * j, y[j]);
        y[j + 1] = y[j] + h * 0.5 * (f(a + h * j, y[j]) + f(a + h * (j + 1), y[j + 1]));
    }
    std::cout << " Метод ПППТ-2:\n";
    
    for (int i = 0; i < y.size(); ++i) {
        std::cout << " y(" << a + h * i << ") = " << y[i] << std::endl;
    }
    
    std::cout << "\n\n";
    return y;
}

std::vector <double> methodEiler(double a, double b, double h, double y0) {
    std::vector <double> y(11);
    y[0] = y0;
    for (size_t j = 0; j < y.size(); ++j) {
        double y_prev = 0.0;
        y[j + 1] = y[j];
        do {
            y_prev = y[j];
            y[j + 1] = y_prev - F(a + (j + 1) * h, y_prev, y[j + 1], h) / Fs(a + (j + 1) * h, y[j + 1], h);
        } while (fabs(y[j] - y_prev) >= epsilon);
    }
    
    std::cout << " Метод Эйлера:\n";
    
    for (int i = 0; i < y.size(); ++i) {
        std::cout << " y(" << a + h * i << ") = " << y[i] << std::endl;
    }
    
    std::cout << "\n\n";
    
    return y;
}

std::vector <double> methodAdams(double a, double b, double h, double y0) {
    std::vector <double> y(11);
    y[0] = y0;
    y[1] = y[0] + h * f(a, y[0]);
    y[2] = y[1] + h * f(a + 0.1, y[1]);
    for (size_t j = 2; j < y.size(); ++j) {
        y[j + 1] = y[j] + (h / 12) * (23 * f(a + h * j, y[j]) - 16 * f(a + h * (j - 1), y[j - 1]) + 5 * f(a + h * (j - 2), y[j - 2]));
    }
    std::cout << " Метод Адамса 3-его порядка:\n";
    
    for (int i = 0; i < y.size(); ++i) {
        std::cout << " y(" << a + h * i << ") = " << y[i] << std::endl;
    }
    
    std::cout << "\n\n";
    return y;
}

void countAccuracy (std::vector <double> y1, std::vector <double> y2) {
    std::cout << " Погрешность:\n";
    for (size_t i = 0; i < y1.size(); ++i) {
        std::cout << " " << y1[i] - y2[i] << "\n";
    }
    std::cout << "\n";
}

int main() {
    double a = 0;          // Задано по условию
    double b = 1;
    double y0 = -1;
    double n = 10;
    double h = (b - a) / n;
    std::vector <double> yr = methodRungeKutta(a, b, h, y0);
    std::vector <double> yp = methodPPPT(a, b, h, y0);
    std::vector <double> ya = methodAdams(a, b, h, y0);
    std::vector <double> ye = methodEiler(a, b, h, y0);
    countAccuracy(ya, yr);
    countAccuracy(ya, yp);
    countAccuracy(ya, ye);
    return 0;
}
