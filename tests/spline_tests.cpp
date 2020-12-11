#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "spline.hpp"


TEST_CASE("testing binary search", "[spline]"){
    Spline tst({1,2,3,4,5}, {2,3,4,5,6});
    auto range1 = tst.bin_search_range(1);
    CHECK(range1.first==0);
    CHECK(range1.second==1);
    auto range2 = tst.bin_search_range(5);
    CHECK(range2.first==3);
    CHECK(range2.second==4);
    auto range3 = tst.bin_search_range(3.5);
    CHECK(range3.first==2);
    CHECK(range3.second==3);
}


TEST_CASE("testing spline calls", "[spline]"){
    Spline tst({1,2,3,4,5}, {2,3,4,5,6});
    CHECK(tst(1)==2);
    CHECK(tst(5)==6);
    CHECK(tst(3)==4);
    CHECK(tst(2.5)==3.5);
}

