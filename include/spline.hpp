#ifndef SPLINE_HPP
#define SPLINE_HPP

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <algorithm>
#include <limits>

class Spline{
private:
    std::vector<double> x_ = {};
    std::vector<double> y_ = {};

public:
    Spline(std::vector<double>&& gx, std::vector<double>&& gy):
        x_(gx), y_(gy) {}
    explicit Spline(std::istream& input);
    double operator()(double gx) const;
    std::pair<size_t, size_t> bin_search_range(double g_x) const;
};



#endif //SPLINE_HPP
