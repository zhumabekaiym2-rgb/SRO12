

package kz.kaznu.aiym.first;

import java.util.Locale;

public class MainApp {

    public static void main(String[] args) {

        Locale.setDefault(Locale.US);

        System.out.println("======================================");
        System.out.println("       СРО № 12");
        System.out.println(" Вариант 4: ВАХ диода");
        System.out.println("======================================");

        // Демонстрационные узлы интерполяции
        double[] x = {-1.0, -0.5, 0.0, 0.5, 1.0};
        double[] y = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            y[i] = Math.exp(-x[i] * x[i]);
        }

        double point = 0.25;

        // Вычисление интерполяции
        double lagrangeResult = lagrange(x, y, point);

        double[] coefficients = dividedDifferences(x, y);
        double newtonResult = newton(x, coefficients, point);

        System.out.println("\nИНТЕРПОЛЯЦИЯ ЛАГРАНЖА И НЬЮТОНА");
        System.out.printf("Точка x = %.2f%n", point);
        System.out.printf("Метод Лагранжа: %.8f%n", lagrangeResult);
        System.out.printf("Метод Ньютона:  %.8f%n", newtonResult);
    }

    // Полином Лагранжа
    public static double lagrange(double[] x, double[] y,
                                  double point) {

        checkData(x, y);

        double result = 0.0;
        int n = x.length;

        for (int i = 0; i < n; i++) {

            double term = y[i];

            for (int j = 0; j < n; j++) {

                if (i != j) {
                    term *= (point - x[j]) / (x[i] - x[j]);
                }
            }

            result += term;
        }

        return result;
    }

    // Таблица разделенных разностей Ньютона
    public static double[] dividedDifferences(double[] x,
                                              double[] y) {

        checkData(x, y);

        int n = x.length;
        double[] coefficients = y.clone();

        for (int j = 1; j < n; j++) {

            for (int i = n - 1; i >= j; i--) {

                coefficients[i] =
                        (coefficients[i] - coefficients[i - 1])
                                / (x[i] - x[i - j]);
            }
        }

        return coefficients;
    }

    // Вычисление полинома Ньютона
    public static double newton(double[] x,
                                double[] coefficients,
                                double point) {

        if (x == null || coefficients == null
                || x.length == 0
                || x.length != coefficients.length) {
            throw new IllegalArgumentException(
                    "Некорректные данные для метода Ньютона");
        }

        double result = coefficients[x.length - 1];

        for (int i = x.length - 2; i >= 0; i--) {

            result = result * (point - x[i])
                    + coefficients[i];
        }

        return result;
    }

    // Проверка исходных данных
    private static void checkData(double[] x, double[] y) {

        if (x == null || y == null
                || x.length == 0
                || x.length != y.length) {
            throw new IllegalArgumentException(
                    "Массивы x и y должны иметь одинаковую длину");
        }

        for (int i = 0; i < x.length; i++) {

            for (int j = i + 1; j < x.length; j++) {

                if (x[i] == x[j]) {
                    throw new IllegalArgumentException(
                            "Узлы интерполяции должны быть различными");
                }
            }
        }
    }
}
