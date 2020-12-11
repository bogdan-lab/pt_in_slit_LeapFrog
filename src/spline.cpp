#include "spline.hpp"


Spline::Spline(std::istream& input){
    double tmp_x, tmp_y;
    input >> std::ws;
    char c = static_cast<char>(input.peek());
    if(c=='#') input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while(input>>tmp_x >> tmp_y){
        x_.push_back(tmp_x);
        y_.push_back(tmp_y);
    }
    if(!std::is_sorted(x_.cbegin(), x_.cend())){
        fprintf(stderr, "Expect spline to be sorted!");
        exit(EXIT_FAILURE);
    }
}

double Spline::operator()(double g_x) const {
    auto range = bin_search_range(g_x);
    size_t l = range.first;
    size_t r = range.second;
    double k = (y_[r]-y_[l])/(x_[r] - x_[l]);
    double b = y_[r] - k*x_[r];
    return k*g_x + b;
}

std::pair<size_t, size_t> Spline::bin_search_range(double g_x) const {
    if(g_x<x_.front() || g_x>x_.back()){
        fprintf(stderr, "Given value %.10e is out of spline bounds\n", g_x);
        exit(EXIT_FAILURE);
    }
    size_t left = 0u;
    size_t right = x_.size()-1;
    while(right-left>1){
        size_t mid = (right+left)/2;
        if(x_[mid]>g_x) right = mid;
        else left = mid;
    }
    return std::make_pair(left, right);
}

