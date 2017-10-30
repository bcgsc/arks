#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "Common/StatUtil.h"
#include <ctgmath>
#include <iostream>

using namespace std;

template <class T>
bool approxEqual(const T& a, const T& b, const T& epsilon) {
	return fabs(a - b) < epsilon;
}

TEST_CASE("quantile calculation", "[StatUtil]")
{
	array<int, 4> data = { 1, 2, 3, 4 };
	REQUIRE(approxEqual(2.5, quantile(data.begin(), data.end(), 0.5), 0.001));
}
