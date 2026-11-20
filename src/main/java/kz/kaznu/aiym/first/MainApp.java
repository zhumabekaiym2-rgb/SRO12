
package kz.kaznu.aiym.first;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;

public class MainApp {

    private static final int ERROR_GRID_POINTS = 10001;

    private static final double K_B = 1.380649e-23;
    private static final double Q = 1.602176634e-19;
    private static final double TEMPERATURE = 300.0;

    public static void main(String[] args) {

        Locale.setDefault(Locale.US);

        System.out.println("======================================");
        System.out.println("              СРО № 12");
        System.out.println("       Основы вычислительной физики");
        System.out.println("          Вариант 4: ВАХ диода");
        System.out.println("======================================");

        // Демонстрационная функция для интерполяции
        double[] x = {-1.0, -0.5, 0.0, 0.5, 1.0};
        double[] y = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            y[i] = Math.exp(-x[i] * x[i]);
        }

        double point = 0.25;

        // Метод Лагранжа
        double lagrangeResult = lagrange(x, y, point);

        // Метод Ньютона
        double[] coefficients = dividedDifferences(x, y);
        double newtonResult = newton(x, coefficients, point);

        // Кубический сплайн
        double splineResult = cubicSpline(x, y, point);

        System.out.println("\nРЕЗУЛЬТАТЫ ИНТЕРПОЛЯЦИИ");
        System.out.printf("Точка x = %.2f%n", point);
        System.out.printf("Лагранж: %.8f%n", lagrangeResult);
        System.out.printf("Ньютон: %.8f%n", newtonResult);
        System.out.printf("Кубический сплайн: %.8f%n", splineResult);

        // Исследование эффекта Рунге
        rungeStudy();

        // Обработка данных диода методом МНК
        diodeLeastSquaresStudy();
    }

    // ==========================================
    // 1. ИНТЕРПОЛЯЦИЯ ЛАГРАНЖА
    // ==========================================

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

    // ==========================================
    // 2. ИНТЕРПОЛЯЦИЯ НЬЮТОНА
    // ==========================================

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

    // ==========================================
    // 3. ЕСТЕСТВЕННЫЙ КУБИЧЕСКИЙ СПЛАЙН
    // ==========================================

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

    // ==========================================
    // 4. ИССЛЕДОВАНИЕ ЭФФЕКТА РУНГЕ
    // ==========================================

    private static double rungeFunction(double x) {

        return 1.0 / (1.0 + 25.0 * x * x);
    }

    private static double[] equidistantNodes(int degree) {

        double[] nodes = new double[degree + 1];

        for (int i = 0; i <= degree; i++) {
            nodes[i] = -1.0 + 2.0 * i / degree;
        }

        return nodes;
    }

    private static double[] chebyshevLobattoNodes(int degree) {

        double[] nodes = new double[degree + 1];

        for (int i = 0; i <= degree; i++) {
            nodes[i] = -Math.cos(Math.PI * i / degree);
        }

        return nodes;
    }

    private static double[] functionValues(double[] nodes) {

        double[] values = new double[nodes.length];

        for (int i = 0; i < nodes.length; i++) {
            values[i] = rungeFunction(nodes[i]);
        }

        return values;
    }

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

    private static void rungeStudy() {

        System.out.println("\nИССЛЕДОВАНИЕ ЭФФЕКТА РУНГЕ");

        System.out.printf("%8s %12s %20s %20s%n",
                "n", "Узлы",
                "Равномерная",
                "Чебышёва");

        try (BufferedWriter writer = Files.newBufferedWriter(
                Path.of("runge_errors.csv"),
                StandardCharsets.UTF_8)) {

            writer.write(
                    "degree,node_count,max_error_equidistant,max_error_chebyshev");

            writer.newLine();

            for (int degree = 4; degree <= 20; degree += 2) {

                double[] uniformNodes =
                        equidistantNodes(degree);

                double[] chebyshevNodes =
                        chebyshevLobattoNodes(degree);

                double uniformError =
                        maxInterpolationError(
                                uniformNodes,
                                functionValues(uniformNodes));

                double chebyshevError =
                        maxInterpolationError(
                                chebyshevNodes,
                                functionValues(chebyshevNodes));

                System.out.printf(Locale.US,
                        "%8d %12d %20.8e %20.8e%n",
                        degree, degree + 1,
                        uniformError, chebyshevError);

                writer.write(String.format(Locale.US,
                        "%d,%d,%.12e,%.12e",
                        degree, degree + 1,
                        uniformError, chebyshevError));

                writer.newLine();
            }

            System.out.println(
                    "Таблица ошибок сохранена в runge_errors.csv");

        } catch (IOException e) {

            System.out.println(
                    "Ошибка записи файла Рунге: " + e.getMessage());
        }
    }

    // ==========================================
    // 5. МЕТОД НАИМЕНЬШИХ КВАДРАТОВ
    // ==========================================

    private static final class DiodeFit {

        double intercept;
        double slope;

        double i0;
        double ideality;

        double rSquared;
        double correlation;

        double seIntercept;
        double seSlope;

        double seI0;
        double seIdeality;

        double sse;

        DiodeFit(double intercept, double slope,
                 double i0, double ideality,
                 double rSquared, double correlation,
                 double seIntercept, double seSlope,
                 double seI0, double seIdeality,
                 double sse) {

            this.intercept = intercept;
            this.slope = slope;
            this.i0 = i0;
            this.ideality = ideality;
            this.rSquared = rSquared;
            this.correlation = correlation;
            this.seIntercept = seIntercept;
            this.seSlope = seSlope;
            this.seI0 = seI0;
            this.seIdeality = seIdeality;
            this.sse = sse;
        }
    }

    // Линейная регрессия ln(I) = A + B*U
    private static DiodeFit linearLeastSquares(
            double[] voltage,
            double[] current,
            double thermalVoltage) {

        int n = voltage.length;

        if (n <= 2 || current.length != n) {
            throw new IllegalArgumentException(
                    "Для МНК нужны минимум три согласованные точки");
        }

        double[] logCurrent = new double[n];

        double meanU = 0.0;
        double meanLogI = 0.0;

        for (int i = 0; i < n; i++) {

            if (current[i] <= 0.0) {
                throw new IllegalArgumentException(
                        "Ток должен быть положительным");
            }

            logCurrent[i] = Math.log(current[i]);

            meanU += voltage[i];
            meanLogI += logCurrent[i];
        }

        meanU /= n;
        meanLogI /= n;

        double sxx = 0.0;
        double sxy = 0.0;
        double syy = 0.0;

        for (int i = 0; i < n; i++) {

            double dx = voltage[i] - meanU;
            double dy = logCurrent[i] - meanLogI;

            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }

        if (sxx == 0.0 || syy == 0.0) {
            throw new IllegalArgumentException(
                    "Недостаточный разброс исходных данных");
        }

        // Коэффициенты прямой
        double slope = sxy / sxx;
        double intercept = meanLogI - slope * meanU;

        if (slope <= 0.0) {
            throw new IllegalArgumentException(
                    "Наклон регрессии должен быть положительным");
        }

        // Остаточная сумма квадратов
        double sse = 0.0;

        for (int i = 0; i < n; i++) {

            double predicted =
                    intercept + slope * voltage[i];

            double residual = logCurrent[i] - predicted;

            sse += residual * residual;
        }

        // Коэффициент детерминации
        double rSquared = 1.0 - sse / syy;

        // Коэффициент корреляции
        double correlation =
                sxy / Math.sqrt(sxx * syy);

        // Оценка дисперсии остатков
        double variance = sse / (n - 2);

        // Стандартные ошибки коэффициентов
        double seSlope = Math.sqrt(variance / sxx);

        double seIntercept = Math.sqrt(
                variance * (1.0 / n + meanU * meanU / sxx)
        );

        // Параметры диода
        double i0 = Math.exp(intercept);
        double ideality = 1.0 / (slope * thermalVoltage);

        // Приближенные стандартные ошибки параметров
        double seI0 = i0 * seIntercept;

        double seIdeality =
                seSlope / (slope * slope * thermalVoltage);

        return new DiodeFit(
                intercept, slope,
                i0, ideality,
                rSquared, correlation,
                seIntercept, seSlope,
                seI0, seIdeality,
                sse
        );
    }

    // Модель диода по найденным параметрам
    private static double diodeCurrent(
            double voltage,
            double i0,
            double ideality,
            double thermalVoltage) {

        return i0 * Math.expm1(
                voltage / (ideality * thermalVoltage)
        );
    }

    // Генерация и обработка смоделированных данных
    private static void diodeLeastSquaresStudy() {

        System.out.println("\nМНК: ВОЛЬТ-АМПЕРНАЯ ХАРАКТЕРИСТИКА ДИОДА");

        double thermalVoltage =
                K_B * TEMPERATURE / Q;

        // Истинные параметры для моделирования
        double trueI0 = 5.0e-9;
        double trueN = 1.8;

        // Напряжение, В
        double[] voltage = {
                0.30, 0.32, 0.34, 0.36, 0.38,
                0.40, 0.42, 0.44, 0.46, 0.48, 0.50
        };

        // Факторы небольшого детерминированного шума
        double[] noise = {
                1.012, 0.991, 1.018, 0.986,
                1.006, 0.995, 1.014, 0.989,
                1.009, 0.994, 1.011
        };

        // Формирование синтетического набора данных
        double[] current = new double[voltage.length];

        for (int i = 0; i < voltage.length; i++) {

            current[i] =
                    trueI0
                            * Math.expm1(
                            voltage[i] / (trueN * thermalVoltage))
                            * noise[i];
        }

        // Аппроксимация методом наименьших квадратов
        DiodeFit fit = linearLeastSquares(
                voltage, current, thermalVoltage);

        System.out.println("Тип данных: смоделированные с шумом");
        System.out.printf(Locale.US,
                "Температура: %.1f K%n", TEMPERATURE);
        System.out.printf(Locale.US,
                "Тепловое напряжение: %.8f V%n", thermalVoltage);

        System.out.println("\nТАБЛИЦА РЕЗУЛЬТАТОВ");

        System.out.printf("%10s %18s %18s %18s%n",
                "U, V", "I измер., A",
                "ln(I)", "I расч., A");

        for (int i = 0; i < voltage.length; i++) {

            double calculated = diodeCurrent(
                    voltage[i],
                    fit.i0,
                    fit.ideality,
                    thermalVoltage
            );

            System.out.printf(Locale.US,
                    "%10.3f %18.8e %18.8f %18.8e%n",
                    voltage[i],
                    current[i],
                    Math.log(current[i]),
                    calculated
            );
        }

        System.out.println("\nНАЙДЕННЫЕ ПАРАМЕТРЫ ДИОДА");

        System.out.printf(Locale.US,
                "I0 = %.8e A%n", fit.i0);

        System.out.printf(Locale.US,
                "n  = %.8f%n", fit.ideality);

        System.out.printf(Locale.US,
                "Стандартная ошибка I0 = %.8e A%n",
                fit.seI0);

        System.out.printf(Locale.US,
                "Стандартная ошибка n = %.8f%n",
                fit.seIdeality);

        System.out.printf(Locale.US,
                "Коэффициент корреляции r = %.8f%n",
                fit.correlation);

        System.out.printf(Locale.US,
                "Коэффициент детерминации R^2 = %.8f%n",
                fit.rSquared);

        System.out.printf(Locale.US,
                "Сумма квадратов остатков = %.8e%n",
                fit.sse);

        System.out.println("\nКонтрольные параметры модели:");
        System.out.printf(Locale.US,
                "Заданный I0 = %.8e A%n", trueI0);
        System.out.printf(Locale.US,
                "Заданный n = %.4f%n", trueN);

        saveDiodeResults(voltage, current, fit, thermalVoltage);
    }

    // Сохранение данных для Excel и построения графика
    private static void saveDiodeResults(
            double[] voltage,
            double[] current,
            DiodeFit fit,
            double thermalVoltage) {

        Path outputFile = Path.of("diode_fit.csv");

        try (BufferedWriter writer = Files.newBufferedWriter(
                outputFile, StandardCharsets.UTF_8)) {

            writer.write(
                    "U_V,I_measured_A,ln_I,I_calculated_A,residual_A");

            writer.newLine();

            for (int i = 0; i < voltage.length; i++) {

                double calculated = diodeCurrent(
                        voltage[i],
                        fit.i0,
                        fit.ideality,
                        thermalVoltage
                );

                double residual = current[i] - calculated;

                writer.write(String.format(Locale.US,
                        "%.6f,%.12e,%.10f,%.12e,%.12e",
                        voltage[i],
                        current[i],
                        Math.log(current[i]),
                        calculated,
                        residual
                ));

                writer.newLine();
            }

            System.out.println("\nФайл diode_fit.csv сохранён.");

        } catch (IOException e) {

            System.out.println(
                    "Ошибка сохранения данных: " + e.getMessage());
        }
    }

    // ==========================================
    // ПРОВЕРКА ИСХОДНЫХ ДАННЫХ
    // ==========================================

    private static void checkData(double[] x, double[] y) {

        if (x == null || y == null
                || x.length == 0
                || x.length != y.length) {

            throw new IllegalArgumentException(
                    "Массивы x и y должны иметь одинаковую длину");
        }

        for (int i = 0; i < x.length; i++) {

            if (!Double.isFinite(x[i])
                    || !Double.isFinite(y[i])) {

                throw new IllegalArgumentException(
                        "Данные должны быть конечными числами");
            }

            for (int j = i + 1; j < x.length; j++) {

                if (x[i] == x[j]) {
                    throw new IllegalArgumentException(
                            "Узлы интерполяции должны быть различными");
                }
            }
        }
    }
}
