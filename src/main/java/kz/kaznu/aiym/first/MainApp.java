



package kz.kaznu.aiym.first;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;

public class MainApp {

    private static final int ERROR_GRID_POINTS = 10001;

    public static void main(String[] args) {

        Locale.setDefault(Locale.US);

        System.out.println("======================================");
        System.out.println("       СРО № 12");
        System.out.println(" Вариант 4: ВАХ диода");
        System.out.println("======================================");

        // Демонстрационные данные для интерполяции
        double[] x = {-1.0, -0.5, 0.0, 0.5, 1.0};
        double[] y = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            y[i] = Math.exp(-x[i] * x[i]);
        }

        double point = 0.25;

        double lagrangeResult = lagrange(x, y, point);

        double[] coefficients = dividedDifferences(x, y);
        double newtonResult = newton(x, coefficients, point);

        double splineResult = cubicSpline(x, y, point);

        System.out.println("\nРЕЗУЛЬТАТЫ ИНТЕРПОЛЯЦИИ");

        System.out.printf("Точка x = %.2f%n", point);
        System.out.printf("Метод Лагранжа: %.8f%n", lagrangeResult);
        System.out.printf("Метод Ньютона: %.8f%n", newtonResult);
        System.out.printf("Кубический сплайн: %.8f%n", splineResult);

        // Исследование эффекта Рунге
        rungeStudy();
    }

    // Полином Лагранжа
    public static double lagrange(double[] x, double[] y,
                                  double point) {

        checkData(x, y);

        return lagrangeUnchecked(x, y, point);
    }

    private static double lagrangeUnchecked(double[] x,
                                            double[] y,
                                            double point) {

        double result = 0.0;

        for (int i = 0; i < x.length; i++) {

            double term = y[i];

            for (int j = 0; j < x.length; j++) {

                if (i != j) {
                    term *= (point - x[j]) / (x[i] - x[j]);
                }
            }

            result += term;
        }

        return result;
    }

    // Разделенные разности Ньютона
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

    // Полином Ньютона
    public static double newton(double[] x,
                                double[] coefficients,
                                double point) {

        if (x == null || coefficients == null
                || x.length == 0
                || x.length != coefficients.length) {

            throw new IllegalArgumentException(
                    "Некорректные данные метода Ньютона");
        }

        double result = coefficients[x.length - 1];

        for (int i = x.length - 2; i >= 0; i--) {

            result = result * (point - x[i])
                    + coefficients[i];
        }

        return result;
    }

    // Естественный кубический сплайн
    public static double cubicSpline(double[] x,
                                     double[] y,
                                     double point) {

        checkData(x, y);

        int n = x.length;

        if (n < 2) {
            throw new IllegalArgumentException(
                    "Для сплайна нужны минимум две точки");
        }

        for (int i = 0; i < n - 1; i++) {

            if (x[i] >= x[i + 1]) {
                throw new IllegalArgumentException(
                        "Узлы должны идти по возрастанию");
            }
        }

        if (point < x[0] || point > x[n - 1]) {
            throw new IllegalArgumentException(
                    "Точка вне интервала интерполяции");
        }

        double[] h = new double[n - 1];

        for (int i = 0; i < n - 1; i++) {
            h[i] = x[i + 1] - x[i];
        }

        double[] m = new double[n];

        int size = n - 2;

        if (size > 0) {

            double[] diagonal = new double[size];
            double[] upper = new double[size];
            double[] rhs = new double[size];

            for (int i = 1; i < n - 1; i++) {

                int k = i - 1;

                diagonal[k] =
                        2.0 * (h[i - 1] + h[i]);

                rhs[k] = 6.0 * (
                        (y[i + 1] - y[i]) / h[i]
                                - (y[i] - y[i - 1]) / h[i - 1]
                );

                if (k < size - 1) {
                    upper[k] = h[i];
                }
            }

            // Прямой ход метода прогонки
            for (int i = 1; i < size; i++) {

                double factor = h[i] / diagonal[i - 1];

                diagonal[i] -= factor * upper[i - 1];
                rhs[i] -= factor * rhs[i - 1];
            }

            // Обратный ход
            m[n - 2] =
                    rhs[size - 1] / diagonal[size - 1];

            for (int i = size - 2; i >= 0; i--) {

                m[i + 1] =
                        (rhs[i] - upper[i] * m[i + 2])
                                / diagonal[i];
            }
        }

        int interval = n - 2;

        for (int i = 0; i < n - 1; i++) {

            if (point <= x[i + 1]) {
                interval = i;
                break;
            }
        }

        double hi = h[interval];

        double left = x[interval + 1] - point;
        double right = point - x[interval];

        return m[interval] * Math.pow(left, 3) / (6.0 * hi)
                + m[interval + 1] * Math.pow(right, 3) / (6.0 * hi)
                + (y[interval] - m[interval] * hi * hi / 6.0)
                * left / hi
                + (y[interval + 1]
                - m[interval + 1] * hi * hi / 6.0)
                * right / hi;
    }

    // Исследование эффекта Рунге
    private static void rungeStudy() {

        System.out.println("\nИССЛЕДОВАНИЕ ЭФФЕКТА РУНГЕ");

        System.out.println(
                "f(x) = 1 / (1 + 25*x^2), интервал [-1; 1]");

        System.out.println(
                "n — степень полинома, число узлов = n + 1");

        System.out.println("-----------------------------------------------");

        System.out.printf("%8s %12s %20s %20s%n",
                "n", "Узлы",
                "Равномерная сетка",
                "Сетка Чебышёва");

        Path outputFile = Path.of("runge_errors.csv");

        try (BufferedWriter writer = Files.newBufferedWriter(
                outputFile, StandardCharsets.UTF_8)) {

            writer.write(
                    "degree,node_count,max_error_equidistant,max_error_chebyshev");

            writer.newLine();

            for (int degree = 4; degree <= 20; degree += 2) {

                double[] uniformNodes =
                        equidistantNodes(degree);

                double[] chebyshevNodes =
                        chebyshevLobattoNodes(degree);

                double[] uniformValues =
                        functionValues(uniformNodes);

                double[] chebyshevValues =
                        functionValues(chebyshevNodes);

                double uniformError =
                        maxInterpolationError(
                                uniformNodes, uniformValues);

                double chebyshevError =
                        maxInterpolationError(
                                chebyshevNodes, chebyshevValues);

                System.out.printf(Locale.US,
                        "%8d %12d %20.8e %20.8e%n",
                        degree,
                        degree + 1,
                        uniformError,
                        chebyshevError);

                writer.write(String.format(Locale.US,
                        "%d,%d,%.12e,%.12e",
                        degree,
                        degree + 1,
                        uniformError,
                        chebyshevError));

                writer.newLine();
            }

            System.out.println("\nТаблица ошибок сохранена:");
            System.out.println(outputFile.toAbsolutePath());

            System.out.println(
                    "Для графика используйте логарифмическую ось Y.");

        } catch (IOException e) {

            System.out.println(
                    "Ошибка сохранения файла: " + e.getMessage());
        }
    }

    // Функция Рунге
    private static double rungeFunction(double x) {

        return 1.0 / (1.0 + 25.0 * x * x);
    }

    // Равномерная сетка
    private static double[] equidistantNodes(int degree) {

        double[] nodes = new double[degree + 1];

        for (int i = 0; i <= degree; i++) {

            nodes[i] = -1.0 + 2.0 * i / degree;
        }

        return nodes;
    }

    // Узлы Чебышёва — Лобатто
    private static double[] chebyshevLobattoNodes(int degree) {

        double[] nodes = new double[degree + 1];

        for (int i = 0; i <= degree; i++) {

            nodes[i] = -Math.cos(Math.PI * i / degree);
        }

        return nodes;
    }

    // Значения функции в узлах
    private static double[] functionValues(double[] nodes) {

        double[] values = new double[nodes.length];

        for (int i = 0; i < nodes.length; i++) {

            values[i] = rungeFunction(nodes[i]);
        }

        return values;
    }

    // Максимальная абсолютная погрешность
    private static double maxInterpolationError(
            double[] nodes, double[] values) {

        double maxError = 0.0;

        for (int i = 0; i < ERROR_GRID_POINTS; i++) {

            double point = -1.0
                    + 2.0 * i / (ERROR_GRID_POINTS - 1);

            double exact = rungeFunction(point);

            double approximate =
                    lagrangeUnchecked(nodes, values, point);

            double error = Math.abs(exact - approximate);

            if (error > maxError) {
                maxError = error;
            }
        }

        return maxError;
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
                            "Узлы должны быть различными");
                }
            }
        }
    }
}
