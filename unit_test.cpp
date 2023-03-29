#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_plane

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <random>
#include <memory>
#include "solver.hpp"


BOOST_AUTO_TEST_CASE(test_plane) {
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

    std::cout << "Есть ли пересечения плоскости и верхней грани: " << std::endl;
    BOOST_CHECK_EQUAL(Plane::Intersection(plane, plane_prizm_top), true);

    std::cout << "Задаем многогранник на этой плоскости из шести граней" << std::endl;
    std::time_t now = std::time(0);
    std::mt19937 gen{static_cast<std::uint32_t>(now)};
    std::uniform_real_distribution<double> dist{-20, 20};

    auto V0_top_prizm = std::make_unique<std::vector<Vector>>();
    for (int i = 0; i < 6; i++) {
        const double tmp = dist(gen);
        const double xx  = tmp;
        const double yy  = tmp + 2*i;
        V0_top_prizm->push_back({xx, yy, plane_prizm_top.get_z(xx, yy)});
    }
    std::cout << "Вычисляем центроид многоугольника: " << std::endl;
    auto barycentr = Vector::Barycentric(*V0_top_prizm); //Вычисляем центроид многоугольника.

    std::cout << "Проверка барицентра, значения должно быть равно нулю или около нуля " << plane_prizm_top.SignedDistance(barycentr) << std::endl;
    BOOST_CHECK_SMALL(plane_prizm_top.SignedDistance(barycentr), 1.0e-10);
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

    std::cout << "Проверка пересечения плоскостей верхней и нижней грани: " << std::endl;
    BOOST_CHECK_EQUAL(Plane::Intersection(plane_prizm_top, plane_prizm_bot(3, 5, 7)), false);
    std::cout << "Проверка на исключение " << std::endl;
    BOOST_CHECK_THROW(plane_prizm_top[-1], std::runtime_error);
}
