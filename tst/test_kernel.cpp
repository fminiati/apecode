#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cassert>
#include <type_traits>
#define USE_APED
#include "Util.h"

using namespace fm::aped;

struct LinearBroadening
{
    static constexpr inline Real fwhm(const Real, const Real, const Real) { return 1; }
    static constexpr spacing_t spacing() { return Spacing::linear_uniform; }
};

template <typename Profile>
void go_through_k(const Real a_len, const Real a_tol)
{
    fm::aped::Kernel k{build_kernel<Profile>(a_len, a_tol)};

    std::cout << "    left wing size =" << k.left_wing_size() << ", right wing size=" << k.right_wing_size() << '\n';
    std::cout << "    weights: \n";
    for (int w = -k.left_wing_size(); w <= k.right_wing_size(); ++w)
        std::cout << "        w[" << w << "]=" << k.weight(w) << '\n';
}

int main(int argc, char* argv[])
{
    const Real line = 1.5;
    const Real g_tolerance = 1.e-6;
    const Real l_tolerance = 1.e-3;

    std::cout << "Gaussian - linear-uniform: \n";
    go_through_k<LineProfile<Gaussian,LinearBroadening>>(line,g_tolerance);

    std::cout << "Lorentzian - linear-uniform: \n";
    go_through_k<LineProfile<Lorentzian, LinearBroadening>>(line, l_tolerance);

    std::cout << "Gaussian - log-uniform: \n";
    go_through_k<LineProfile<Gaussian, ThermalBroadening>>(line,g_tolerance);

    std::cout << "Lorentzian - log-uniform: \n";
    go_through_k<LineProfile<Lorentzian, ThermalBroadening>>(line, l_tolerance);
}
