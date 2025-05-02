#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

void drawGraph(const std::vector<double>& x, const std::vector<double>& y, int width = 60, int height = 20) {
    if (x.size() != y.size() || x.empty()) {
        std::cerr << "Ошибка: векторы x и y должны быть одинаковой длины и не пустыми.\n";
        return;
    }

    // Находим диапазоны значений
    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());
    double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());

    // Инициализация пустого поля
    std::vector<std::string> grid(height, std::string(width, ' '));

    for (size_t i = 0; i < x.size(); ++i) {
        int col = static_cast<int>((x[i] - x_min) / (x_max - x_min) * (width - 1));
        int row = static_cast<int>((y[i] - y_min) / (y_max - y_min) * (height - 1));
        row = height - 1 - row; // инвертировать по y для корректного отображения
        if (row >= 0 && row < height && col >= 0 && col < width)
            grid[row][col] = '*';
    }

    // Вывод графика
    for (const auto& line : grid)
        std::cout << "|" << line << "|\n";
    
    // Ось x
    std::cout << "+" << std::string(width, '-') << "+\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "x: [" << x_min << ", " << x_max << "]  ";
    std::cout << "y: [" << y_min << ", " << y_max << "]\n";
}
