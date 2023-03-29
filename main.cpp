
#include <functional>
#include <random>
#include <ctime>
#include <memory>
#include <omp.h>

#include "solver.hpp"

constexpr int N = 1000000000;

int main()
{

    std::cout << "Задаем вектор и нормаль нашей плоскости" << std::endl;
    Vector V0(1, -3, 4), normal(2, 1, 3); //Векторы задающие плоскость
    std::cout << "Vector = " << V0 << "; Normal = " << normal << std::endl;

    auto plane = Plane::FromPointNormal(V0, normal);    //Создаем плоскость
    std::cout << "Plane = " << plane << std::endl;

    std::cout << "Задаем вектор и нормаль плоскости верхней грани призмы" << std::endl;
    Vector V0_top(2, 4, -7), normal_t(-20, 1, -3); //Векторы задающие плоскость верхней грани призмы
    std::cout << "Vector = " << V0_top << "; Normal = " << normal_t << std::endl;
    auto plane_prizm_top = Plane::FromPointNormal(V0_top, normal_t); //Создаем плоскость верхней грани
    std::cout << "Plane = " << plane_prizm_top << std::endl;

    //Задаем многогранник на этой плоскости из шести граней
    std::cout << "Задаем многогранник на этой плоскости из шести граней" << std::endl;
    std::time_t now = std::time(0);
    std::mt19937 gen{static_cast<std::uint32_t>(now)};
    std::uniform_real_distribution<double> dist{-20, 20};
//    std::vector<Vector> V0_top_prizm;
    auto V0_top_prizm = std::make_unique<std::vector<Vector>>();
    for (int i = 0; i < 6; i++) {
        const double tmp = dist(gen);
        const double xx  = tmp;
        const double yy  = tmp + 2*i;
        V0_top_prizm->push_back({xx, yy, plane_prizm_top.get_z(xx, yy)});
    }

    std::cout << "Координаты точек многогранника: " << std::endl;
    for(const auto& x : *V0_top_prizm) {
        std::cout << x << "; ";
    }
    std::cout << std::endl;
    std::cout << "Вычисляем центроид многоугольника: " << std::endl;
    auto barycentr = Vector::Barycentric(*V0_top_prizm); //Вычисляем центроид многоугольника.
    std::cout << "Барицентр: " << barycentr << std::endl;

    std::cout << "Вычисляем переменную плоскость нижней грани, паралельную плоскости верхней грани: " << std::endl;

    auto plane_prizm_bot = [=](double x, double y, double z) {//Создаем плоскость нижней грани
        Vector tmp(x, y, z);
        const double length_prizm_ = barycentr.length(tmp);
        const double D2 = length_prizm_*std::sqrt(plane_prizm_top[0]*plane_prizm_top[0] +
                                                  plane_prizm_top[1]*plane_prizm_top[1] +
                                                  plane_prizm_top[2]*plane_prizm_top[2]) + plane_prizm_top[3];
        return Plane(normal_t, D2);
    };
    std::cout << "Задаем конечный вектор плоскости нижней грани и его нормаль: " << std::endl;
    Vector V0_b_(3, 5, 7);
    std::cout << "Vector = " << V0_b_ << "; Normal = " << normal_t << std::endl;
    std::cout << "плоскость нижней грани = " << plane_prizm_bot(3, 5, 7) << std::endl;

    std::cout << "Задаем функцию растояния от координат центра элементарного объема до плоскости: " << std::endl;
    std::function<double(double, double, double)> func = [=](double x, double y, double z) {
        const auto plane_tmp = plane_prizm_bot(x, y, z);
        double arg = std::sqrt(plane[0]*plane[0] +
                               plane[1]*plane[1] +
                               plane[2]*plane[2]);
        return std::abs(plane[0]*x + plane[1]*y + plane[2]*plane_tmp.get_z(x, y) + plane[3])/arg;
    };

    std::cout << "Границы интегрирования нижняя: " << barycentr << std::endl;
    const double lower_bound_x = barycentr[0];
    const double lower_bound_y = barycentr[1];
    const double lower_bound_z = barycentr[2];

    std::cout << "Границы интегрирования верхняя: " << V0_b_ << std::endl;
    const double upper_bound_x = V0_b_[0];
    const double upper_bound_y = V0_b_[1];
    const double upper_bound_z = V0_b_[2];

    std::cout << "Интегрируем методом Монте Карло: " << std::endl;
    std::random_device rd; // Генератор случайных чисел
    double sum = 0.0;
    #pragma omp parallel
    {
        int seed;
        #pragma omp critical
        {
            seed = rd();
        }

        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist_x(lower_bound_x, upper_bound_x);
        std::uniform_real_distribution<double> dist_y(lower_bound_y, upper_bound_y);
        std::uniform_real_distribution<double> dist_z(lower_bound_z, upper_bound_z);

        #pragma omp for reduction(+ : sum)
        for (int i = 0; i < N; i++) {
            double x = lower_bound_x + (upper_bound_x - lower_bound_x) * dist_x(gen);
            double y = lower_bound_y + (upper_bound_y - lower_bound_y) * dist_y(gen);
            double z = lower_bound_z + (upper_bound_z - lower_bound_z) * dist_z(gen);
            sum += func(x, y, z);
        }
    }

    double result = sum/**(upper_bound_x - lower_bound_x)*(upper_bound_y - lower_bound_y)*(upper_bound_z - lower_bound_z)*//N;
    std::cout << "Результат интегрирования: " << result << std::endl;
    return 0;
}


