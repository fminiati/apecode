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

    // conversion from keV to Angstrom
    constexpr double keVToAngstrom = 12.39854;
    // conversion from keV to Kelvin
    constexpr double keVToKelvin = 1.1604505e7;

    // set tolerance for convolution ops to float precision as atomdb is float
    constexpr double KERNEL_TOL = std::numeric_limits<float>::epsilon();

    constexpr double zero = 0.0e0;
    constexpr double half = 0.5e0;
    constexpr double one = 1.0e0;
    constexpr double two = 2.0e0;
    constexpr double ten = 1.0e1;
    constexpr double pi = 3.1415926535897932e0;
    constexpr double sqrt_two = 1.4142135623730951e0;
    constexpr double sqrt_ln2 = 0.8325546111576977e0; //std::sqrt(std::log(2));

    // Boltzmann's constant
    constexpr double kB_cgs = 1.3806488e-16;
    // Speed of light in cm s^-1
    constexpr double c_light_cgs = 2.9979246e10;

    enum class LineShape
    {
        delta = 0,
        gaussian = 1,
        lorentzian = 2,
        pseudovoigt = 3
    };

    enum class LineBroadening
    {
        none = 0,
        thermal = 1,
        turbulent = 2
        //, collisional = 3
    };

    struct Spacing
    {
        enum type
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
        static constexpr auto grid_spacing(const std::vector<T> &a_x)
        {
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
            if (all_adjacent_of(a_x, [dx, eps](const T a, const T b) { return std::fabs(b - a - dx) < eps * dx; }))
                return Spacing::linear_uniform;

            // check logarithmic
            const T dlx = a_x[1] / a_x[0];
            if (all_adjacent_of(a_x, [dlx, eps](const T a, const T b) { return std::fabs(a / b * dlx - one) < eps * dlx; }))
                return Spacing::log_uniform;

            return Spacing::irregular;
        }
    };
    using spacing_t = Spacing::type;

    // broadening mechanism determining FWHM of emission lines and its corresponding type of spacing
    struct ThermalBroadening
    {
        static constexpr Real sqrt_8_ln2_to_c_cgs = two * sqrt_two * sqrt_ln2 / c_light_cgs;
        static inline Real fwhm(const Real a_temp, const Real a_atomic_mass, const Real a_ph_energy)
        {
            return sqrt_8_ln2_to_c_cgs * std::sqrt(kB_cgs * a_temp / a_atomic_mass);
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
        static inline Real area(const Real) { return 1; }
        static inline Real fwhm(const Real, const Real, const Real) { return 0; }
    };
   
    struct Gaussian
    {
        static inline Real area(const Real a_x)  { return half * std::erf(sqrt_ln2 * a_x); }
    };

    struct Lorentzian
    {
        static constexpr Real one_over_pi = one / pi;
        static inline Real area(const Real a_x)  { return one_over_pi * std::atan(a_x); }
    };

    struct PseudoVoigt
    {
        static inline Real area(const Real a_x) 
        {
            return (one - m_eta) * Gaussian::area(a_x) + m_eta * Lorentzian::area(a_x);
        }
        static void set_eta(const Real a_eta) { m_eta = a_eta; }
    protected:
        static thread_local Real m_eta;
    };
    inline thread_local Real PseudoVoigt::m_eta{};


    // Broadening describe the mechanism determing the full-width-at-half-maximum
    template <typename Shape, typename Broadening>
    struct LineProfile
    {
        static inline Real area(const Real a_x) { return Shape::area(a_x); }
        static inline Real fwhm(const Real a_t, const Real a_m, const Real a_e) { return Broadening::fwhm(a_t, a_m, a_e); }
    };

    // specialisation for Voigt
    template <typename Voigt, typename G, typename L, template<typename...> typename Broadening>
    struct LineProfile<Voigt, Broadening<G,L>>
    {
        static inline Real area(const Real a_x) { return Voigt::area(a_x); }
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
    template <typename S, typename B>
    struct line_profile_type_traits<LineProfile<S, B>>
    {
        using shape_t = S;
        using broadening_t = B;
    };
    template <typename T>
    using line_shape_t = typename line_profile_type_traits<T>::shape_t;
    template <typename T>
    using line_broadening_t = typename line_profile_type_traits<T>::broadening_t;
    template <typename T>
    inline constexpr spacing_t broadening_spacing_v = std::integral_constant<spacing_t, line_broadening_t<T>::spacing()>();

    // kernel object to compute convolutio in case of uniformly spaced grids such that the weights are const.
    template <typename Profile>
    struct LineKernel
    {
        ~LineKernel() {}

        LineKernel(const Real a_mesh_size)
        {
            static_assert(Spacing::is_uniform(broadening_spacing_v<Profile>), "LineKernel requires const Spacing type");
            Timer_t<4> t("Aped::LineKernel");

            // compute integral of line shape from line centre to a mesh node along a wing. it is expected to asymptote to
            // half since the line shape is normalised to one. a_next_mesh_node is a function taking a node
            // and returning the next thus encapsulating information about the grid structure
            auto _weights = [](const Real a_centre, const Real a_next_node, const auto a_next_mesh_node) {
                const Real s = a_next_node > a_centre ? one : -one;
                std::vector<Real> w(1, line_shape_t<Profile>::area(s * (a_centre - a_next_node)));
                Real mesh_node = a_next_node;
                while (half - w.back() > KERNEL_TOL)
                {
                    mesh_node = a_next_mesh_node(mesh_node);
                    w.emplace_back(line_shape_t<Profile>::area(s * (a_centre - mesh_node)));
                }
                return w;
            };

            std::vector<Real> wm, wp;
            {
                const Real lo = one, mid = one + half * a_mesh_size, hi = one + a_mesh_size;
                if constexpr(broadening_spacing_v<Profile> == Spacing::log_uniform)
                {
                    wm = _weights(mid, lo, [f = one / a_mesh_size](const Real x) { return f * x; });
                    wp = _weights(mid, hi, [f = a_mesh_size](const Real x) { return f * x; });
                }
                else if constexpr (broadening_spacing_v<Profile> == Spacing::linear_uniform)
                {
                    wm = _weights(mid, lo, [dx = a_mesh_size](const Real x) { return x - dx; });
                    wp = _weights(mid, hi, [dx = a_mesh_size](const Real x) { return x + dx; });
                }
            }

            // now go backward, and get the proper weights by subtracting successive erfs
            // first left wing
            for (size_t i = wm.size() - 1; i > 0; --i)
                m_w.emplace_back(wm[i] - wm[i - 1]);
            // centre
            m_w.emplace_back(wm[0] + wp[0]);
            // right wing
            for (size_t i = 1; i < wp.size(); ++i)
                m_w.emplace_back(wp[i] - wp[i - 1]);

            m_left_wing_size = wm.size() - 1;
            m_right_wing_size = wp.size() - 1;
#ifndef NDEBUG
            Real sum = zero;
            for (auto w : wm)
                sum += w; 
            for (auto w : wp)
                sum += w; 
            assert(abs(one - sum) < KERNEL_TOL);
#endif
        }

        virtual inline size_t left_wing_size() const 
        {
            return m_left_wing_size;
        }

        virtual inline size_t right_wing_size() const 
        {
            return m_right_wing_size;
        }

        // peak is centered at m_wing_size
        virtual inline Real weight(const size_t a_index) const 
        {
            return (Real)m_w[a_index + m_left_wing_size];
        }

    protected:
        size_t m_left_wing_size, m_right_wing_size;
        std::vector<Real> m_w;
    };

    // convolution
    struct Convolution
    {
        Convolution() {}
        ~Convolution() {}

         // convolve spectrum with this kernel
        template <typename Kernel>
        static inline void convolve(std::vector<Real> &a_f, const Kernel &a_kernel)
        {
            Timer_t<4> t("LineKernel::convolve_spectrum_kernel");

            std::vector<Real> convolution(a_f.size());
            std::memset(&convolution[0], 0, sizeof(Real) * convolution.size());
            for (size_t i = 0; i < a_f.size(); ++i)
            {
                convolve(convolution, a_f[i], i, a_kernel);
            }
            a_f.swap(convolution);
        }

        // convolve single line emission with kernel
        template <typename Kernel>
        static inline void convolve(std::vector<Real> &a_c,
                                    const Real a_I0,
                                    const size_t a_bin,
                                    const Kernel &a_kernel)
        {
            Timer_t<4> t("LineKernel::convolve_line_kernel");

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
                             const std::vector<Real> &a_fwhm)
        {
            Timer_t<4> t("Aped::GaussianKernel::convolve_spectrum_profile");

            std::vector<Real> convolution(a_f.size());
            std::memset(&convolution[0], 0, sizeof(Real) * convolution.size());
            for (size_t i = 0; i < a_f.size(); ++i)
            {
                convolve<Shape>(convolution, a_f[i], half * (a_x[i] + a_x[i + 1]), a_fwhm[i], i, a_x);
            }
            a_f.swap(convolution);
        }

        template <typename Shape>
        static void convolve(std::vector<Real> &a_c,
                             const Real a_I0,
                             const Real a_centre,
                             const Real a_fwhm,
                             const size_t a_bin,
                             const std::vector<Real> &a_x)
        {
            Timer_t<4> t("Aped::LineKernel::convolve_once_line_profile");
            assert(a_bin + 1 < a_x.size() && a_centre > a_x[a_bin] && a_centre < a_x[a_bin + 1]);

            const Real norm = two / a_fwhm;
            { // left wing
                double w = Shape::area(norm * (a_centre - a_x[a_bin]));
                a_c[a_bin] += w * a_I0;

                double wm = 0;
                double err = half - w;
                for (int i = a_bin - 1; err > KERNEL_TOL && i >= 0; --i)
                {
                    wm += w;
                    w = Shape::area(norm * (a_centre - a_x[i])) - wm;
                    a_c[i] += w * a_I0;
                    err -= w;
                }
            }
            { // right wing
                double w = Shape::area(norm * (a_x[a_bin + 1] - a_centre));
                a_c[a_bin] += w * a_I0;

                double wm = 0;
                double err = half - w;
                for (size_t i = a_bin + 1; err > KERNEL_TOL && i < a_c.size(); ++i)
                {
                    wm += w;
                    w = Shape::area(norm * (a_x[i + 1] - a_centre)) - wm;
                    a_c[i] += w * a_I0;
                    err -= w;
                }
            }
        }
    };
} // namespace fm::aped

#endif // APED_UTIL_H
#endif // USE_APED
