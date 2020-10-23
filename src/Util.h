//
// Copyright (C) 2020 Francesco Miniati <francesco.miniati@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#ifdef USE_APED

#ifndef APED_UTIL_H
#define APED_UTIL_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <cassert>

#ifdef USE_TIMER
#include "Timer.h" // available at https://github.com/fminiati/mthread-timer
using namespace fm::profiling;
#else
#ifndef TIMER_H
#define TIMER_H
template <unsigned T=0> struct Timer_t {
    Timer_t(std::string &&){};
};
#endif
#endif

namespace fm::aped
{
    using Real = double;

    constexpr Real MIN_W_INCR_TO_KERNEL_TOL = 0.1;

    constexpr Real zero = 0.0e0;
    constexpr Real half = 0.5e0;
    constexpr Real one = 1.0e0;
    constexpr Real two = 2.0e0;
    constexpr Real ten = 1.0e1;
    constexpr Real pi = 3.1415926535897932e0;
    constexpr Real sqrt_two = 1.4142135623730951e0;
    constexpr Real sqrt_ln2 = 0.8325546111576977e0; //std::sqrt(std::log(2));

    // Boltzmann's constant
    constexpr Real kB_cgs = 1.380648528e-16; // Xspec values: 1.3806511609069063e-16 
    // Speed of light in cm s^-1
    constexpr Real c_light_cgs = 2.99792458e10;
    // conversion from keV to Angstrom
    constexpr Real keVToAngstrom = 12.39841974e0;
    // conversion from keV to Kelvin
    constexpr Real keVToKelvin = 1.1604505e7;

    enum class LineShape : char
    {
        delta = 0,
        gaussian = 1,
        lorentzian = 2,
        pseudovoigt = 3
    };

    enum class LineBroadening : char
    {
        none = 0,
        thermal = 1
        //, turbulent = 2
        //, collisional = 3
    };

    struct Spacing
    {
        enum type : char
        {
            undetermined = 0,
            irregular = 1,
            linear_uniform = 2,
            log_uniform = 3
        };

        static constexpr bool is_uniform(const Spacing::type a_type)
        {
            return a_type == linear_uniform || a_type == log_uniform;
        }

        template <typename T>
        static auto grid_spacing(const std::vector<T> &a_x)
        {
            Timer_t<3> t("grid_spacing");
            constexpr T eps = std::numeric_limits<T>::epsilon();

            if (a_x.size() < 3)
                return Spacing::undetermined;

            auto all_adjacent_of = [](const std::vector<T> x, const auto test) {
                size_t i = 0;
                while (++i < x.size() && test(x[i - 1], x[i]))
                {}
                return i == x.size();
            };

            // check linear
            const T dx = a_x[1] - a_x[0];
            if (all_adjacent_of(a_x, [dx](const T a, const T b) { return std::fabs(b - a - dx) < eps * dx; }))
                return Spacing::linear_uniform;

            // check logarithmic
            const T dlx = a_x[1] / a_x[0];
            if (all_adjacent_of(a_x, [dlx](const T a, const T b) { return std::fabs(a / b * dlx - one) < eps * dlx; }))
                return Spacing::log_uniform;

            return Spacing::irregular;
        }
    };
    using spacing_t = Spacing::type;

    // broadening mechanism determining FWHM of emission lines and its corresponding type of spacing
    struct NoBroadening
    {
        static constexpr inline Real fwhm(const Real, const Real, const Real) { return 0; }
        static constexpr spacing_t spacing() { return Spacing::undetermined; }
    };
    struct ThermalBroadening
    {
        static constexpr Real sqrt_8_ln2_to_c_cgs = two * sqrt_two * sqrt_ln2 / c_light_cgs;
        static inline Real fwhm(const Real a_temp, const Real a_atomic_mass, const Real a_ph_energy)
        {
            return sqrt_8_ln2_to_c_cgs * std::sqrt(kB_cgs * a_temp / a_atomic_mass) * a_ph_energy;
        }
        static constexpr spacing_t spacing() { return Spacing::log_uniform; }
    };

    // define voigt broadening two mechanism through templates
    template <typename GaussianBroadening, typename LorentzianBroadening>
    struct PseudoVoigtBroadening
    {
        static inline Real fwhm(const Real a_t, const Real a_m, const Real a_e)
        {
            const auto [f,eta] = fwhm_and_eta(a_t, a_m, a_e);
            return f;
        }
        static inline auto fwhm_and_eta(const Real a_t, const Real a_m, const Real a_e)
        {
            const Real fG = GaussianBroadening::fwhm(a_t, a_m, a_e);
            const Real fL = LorentzianBroadening::fwhm(a_t, a_m, a_e);
            const Real fL2 = fL * fL, fG2 = fG * fG;
            const Real fL3 = fL2 * fL, fG3 = fG * fG2;
            const Real fL4 = fL2 * fL2, fG4 = fG2 * fG2;
            const Real fL5 = fL3 * fL2, fG5 = fG3 * fG2;
            const Real f = std::pow(fG5 + 2.69269 * fG4 * fL + 2.42843 * fG3 * fL2 + 4.47163 * fG2 * fL3 + 0.07842 * fG * fL4 + fL5, 0.2);
            const Real x = fL / f;
            const Real eta = x * (1.36603 + x * (-0.47719 + x * 0.11116));

            return std::pair(f, eta);
        }
        static constexpr spacing_t spacing()
        {
            return GaussianBroadening::spacing() == LorentzianBroadening::spacing() ? GaussianBroadening::spacing() : Spacing::irregular;
        }
    };

    // line shapes
    struct Delta
    {
        static inline constexpr Real tail_integral(const Real) { return 1; }
    };
   
    struct Gaussian
    {
        static inline Real tail_integral(const Real a_x)  { return half * std::erf(sqrt_ln2 * a_x); }
    };

    struct Lorentzian
    {
        static constexpr Real one_over_pi = one / pi;
        static inline Real tail_integral(const Real a_x)  { return one_over_pi * std::atan(a_x); }
    };

    struct PseudoVoigt
    {
        static inline Real tail_integral(const Real a_x) 
        {
            return (one - m_eta) * Gaussian::tail_integral(a_x) + m_eta * Lorentzian::tail_integral(a_x);
        }
        static void set_eta(const Real a_eta) { m_eta = a_eta; }
    protected:
        static thread_local Real m_eta;
    };
    inline thread_local Real PseudoVoigt::m_eta{};


    // Broadening describe the mechanism determing the full-width-at-half-maximum
    template <typename Shape, typename Broadening, bool PseudoContBrd=false>
    struct LineProfile
    {
        static inline Real tail_integral(const Real a_x) { return Shape::tail_integral(a_x); }
        static inline Real fwhm(const Real a_t, const Real a_m, const Real a_e) { return Broadening::fwhm(a_t, a_m, a_e); }
    };

    // specialisation for Voigt
    template <typename Voigt, typename G, typename L, template<typename...> typename Broadening, bool PseudoContBrd>
    struct LineProfile<Voigt, Broadening<G,L>, PseudoContBrd>
    {
        static inline Real tail_integral(const Real a_x) { return Voigt::tail_integral(a_x); }
        static inline Real fwhm(const Real a_t, const Real a_m, const Real a_e)
        {
            const auto [f, eta] = Broadening<G,L>::fwhm_and_eta(a_t, a_m, a_e);
            Voigt::set_eta(eta);
            return f;
        }
    };

    // type traits for line profile
    template <typename T> 
    struct line_profile_type_traits {};
    template <typename S, typename B, bool Z>
    struct line_profile_type_traits<LineProfile<S, B, Z>>
    {
        static_assert(std::is_same_v<S, Delta> == std::is_same_v<B, NoBroadening>,
                      "Inconsistent Line Profile: choose Delta profile if and only if NoBroadening is the broadening mechanism.");
        static_assert(std::disjunction_v<std::negation<std::bool_constant<Z>>, std::negation<std::is_same<S, Delta>>>,
                      "Inconsistent line profile: apply pseudo continuum broadening if and only if line shape is not a delta.");

        using shape_t = S;
        using broadening_t = B;
        static constexpr bool pseudo_cont_broadening_v = Z;
    };
    template <typename T>
    using line_shape_t = typename line_profile_type_traits<T>::shape_t;
    template <typename T>
    using line_broadening_t = typename line_profile_type_traits<T>::broadening_t;
    template <typename T>
    inline constexpr bool apply_pseudo_cont_broadening_v = line_profile_type_traits<T>::pseudo_cont_broadening_v;
    template <typename T>
    inline constexpr bool apply_line_broadening_v = std::is_same_v<line_shape_t<T>, Delta>;
    template <typename T>
    inline constexpr spacing_t broadening_spacing_v = std::integral_constant<spacing_t, line_broadening_t<T>::spacing()>();

    struct Kernel
    {
        ~Kernel() {}

        Kernel(const int a_left_size, const int a_right_size, const std::vector<Real> a_w)
            : m_left_wing_size{a_left_size}, m_right_wing_size{a_right_size}, m_w{a_w}
        {
        }

        virtual inline int left_wing_size() const
        {
            return m_left_wing_size;
        }

        virtual inline int right_wing_size() const
        {
            return m_right_wing_size;
        }

        // peak is centered at m_wing_size
        virtual inline Real weight(const int a_index) const
        {
            return m_w[a_index + m_left_wing_size];
        }

    protected:
        int m_left_wing_size, m_right_wing_size;
        std::vector<Real> m_w;
    };

    // convolution
    struct Convolution
    {
        Convolution() {}
        ~Convolution() {}

         // convolve spectrum with this kernel
        static inline void convolve(std::vector<Real> &a_f, Kernel &&a_kernel)
        {
            Timer_t<5> t("conv_spect_kern");

            std::vector<Real> convolution(a_f.size());
            std::memset(&convolution[0], 0, sizeof(Real) * convolution.size());
            for (size_t i = 0; i < a_f.size(); ++i)
            {
                convolve(convolution, a_f[i], i, a_kernel);
            }
            a_f.swap(convolution);
        }

        // convolve single line emission with kernel
        static inline void convolve(std::vector<Real> &a_c,
                                    const Real a_I0,
                                    const size_t a_bin,
                                    const Kernel &a_kernel)
        {
            Timer_t<5> t("conv_line_kern");

            const size_t lo = std::max(a_bin - a_kernel.left_wing_size(), (size_t)0);
            const size_t hi = std::min(a_bin + a_kernel.right_wing_size(), a_c.size() - 1);
            for (size_t pos = lo; pos <= hi; ++pos)
            {
                a_c[pos] += a_I0 * a_kernel.weight(pos - a_bin);
            }
        }

        // use this and the following function to directly compute the convolution
        // when building the weights is not advantegeous
        template <typename Shape>
        static void convolve(std::vector<Real> &a_f,
                             const std::vector<Real> &a_x,
                             const std::vector<Real> &a_fwhm,
                             const Real a_tolerance)
        {
            Timer_t<5> t("conv_spectrum");

            std::vector<Real> convolution(a_f.size());
            std::memset(&convolution[0], 0, sizeof(Real) * convolution.size());
            for (size_t i = 0; i < a_f.size(); ++i)
            {
                convolve<Shape>(convolution, a_f[i], half * (a_x[i] + a_x[i + 1]), a_fwhm[i], i, a_x, a_tolerance);
            }
            a_f.swap(convolution);
        }

        template <typename Shape>
        static inline void convolve(std::vector<Real> &a_c,
                                    const Real a_I0,
                                    const Real a_centre,
                                    const Real a_fwhm,
                                    const size_t a_bin,
                                    const std::vector<Real> &a_x,
                                    const Real a_tolerance)
        {
            Timer_t<5> t("conv_line");
            assert(a_bin + 1 < a_x.size() && a_centre > a_x[a_bin] && a_centre < a_x[a_bin + 1]);

            const Real delta_tol = MIN_W_INCR_TO_KERNEL_TOL * a_tolerance;

            const Real norm = two / a_fwhm;
            { // left wing
                const Real asymptote = -Shape::tail_integral(-std::numeric_limits<Real>::infinity());
                Real w = -Shape::tail_integral(norm * (a_x[a_bin] - a_centre));
                a_c[a_bin] += w * a_I0;

                Real wm{};
                Real err = asymptote - w;
                for (int i = a_bin - 1; err > a_tolerance && i >= 0 && w >= delta_tol; --i)
                {
                    wm += w;
                    w = -Shape::tail_integral(norm * (a_x[i] - a_centre)) - wm;
                    a_c[i] += w * a_I0;
                    err -= w;
                }
#ifndef NDEBUG
                if (err >  a_tolerance)
                    std::clog << "Convolution::convolve: left-wing weight calculation did not converge... "
                              << "the residual error is " << err << '\n';
#endif
            }
            { // right wing
                const Real asymptote = Shape::tail_integral(std::numeric_limits<Real>::infinity());
                Real w = Shape::tail_integral(norm * (a_x[a_bin + 1] - a_centre));
                a_c[a_bin] += w * a_I0;

                Real wm{};
                Real err = asymptote - w;
                for (size_t i = a_bin + 1; err > a_tolerance && i < a_c.size() && w >= delta_tol; ++i)
                {
                    wm += w;
                    w = Shape::tail_integral(norm * (a_x[i + 1] - a_centre)) - wm;
                    a_c[i] += w * a_I0;
                    err -= w;
                }
#ifndef NDEBUG
                if (err >  a_tolerance)
                    std::clog << "Convolution::convolve: right-wing weight calculation did not converge... "
                              << "the residual error is " << err << '\n';
#endif
            }
        }
    };

    template <typename Profile>
    auto build_kernel(const Real a_length, const Real a_kernel_tol)
    {
        static_assert(Spacing::is_uniform(broadening_spacing_v<Profile>), "build_kernel requires uniform Spacing type");
        Timer_t<3> t("build_kernel");

        // compute integral of shape from centre to a mesh nodes along a wing. a_next_node is a function 
        // taking a node and returning the next in line
        auto kernel_weights = [a_kernel_tol](const Real a_centre, const Real a_node, const auto a_next_node) {
            const Real s = two * (a_node > a_centre ? one : -one);
            const Real asymptote = s * line_shape_t<Profile>::tail_integral(s * std::numeric_limits<Real>::infinity());
            const Real delta_tol = MIN_W_INCR_TO_KERNEL_TOL * a_kernel_tol;

            std::vector<Real> w(1, s * line_shape_t<Profile>::tail_integral(a_node - a_centre));
            Real node = a_node;
            Real dw = one;
            while (asymptote - w.back() > a_kernel_tol && dw >= delta_tol)
            {
                dw = -w.back();
                node = a_next_node(node);
                w.emplace_back(s * line_shape_t<Profile>::tail_integral(node - a_centre));
                dw += w.back();
            }
#ifndef NDEBUG
            if (const Real err = asymptote - w.back(); err > a_kernel_tol)
                std::clog << "build_kernel: " << (s < zero ? "left" : "right")
                          << " wing weight calculation did not converge... the residual error is " << err << '\n';
#endif
            return w;
        };

        std::vector<Real> wm, wp, w;
        {
            const Real lo = one, mid = one + half * a_length, hi = one + a_length;
            if constexpr (broadening_spacing_v<Profile> == Spacing::log_uniform)
            {
                wm = kernel_weights(mid, lo, [f = one / a_length](const Real x) { return f * x; });
                wp = kernel_weights(mid, hi, [f = a_length](const Real x) { return f * x; });
            }
            else if constexpr (broadening_spacing_v<Profile> == Spacing::linear_uniform)
            {
                wm = kernel_weights(mid, lo, [dx = a_length](const Real x) { return x - dx; });
                wp = kernel_weights(mid, hi, [dx = a_length](const Real x) { return x + dx; });
            }
        }
        // now go backward, and get the proper weights by successive subtraction
        // first left wing
        for (size_t i = wm.size() - 1; i > 0; --i)
            w.emplace_back(wm[i] - wm[i - 1]);
        // centre
        w.emplace_back(wm[0] + wp[0]);
        // right wing
        for (size_t i = 1; i < wp.size(); ++i)
            w.emplace_back(wp[i] - wp[i - 1]);

#ifndef NDEBUG
        Real sum = zero;
        for (auto x : w)
            sum += x;
        if (const Real err = abs(one - sum); err > two * a_kernel_tol)
            std::clog << "build_kernel: kernel not normalized, error is :" << std::setprecision(10) << err << '\n';
#endif

        return Kernel(wm.size() - 1, wp.size() - 1, w);
    }
} // namespace fm::aped

#endif // APED_UTIL_H
#endif // USE_APED
