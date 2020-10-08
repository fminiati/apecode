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

#ifndef APEC_H
#define APEC_H

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "Aped.h"
#include "SutherlandDopita.h"

namespace fm::aped
{

    // return interpolated function (f(x)) value at input position xp
    template <class T, class R>
    void linear_interp(std::vector<T> &fp, const std::vector<T> &xp,
                       const std::vector<R> &f, const std::vector<R> &x)
    {
        Timer_t<4> t("Apec::linear_interp");

        std::vector<T> xh(f.size());
        for (size_t i = 0; i<f.size(); ++i)
            xh[i] = half * (x[i] + *x[i + 1]);

        size_t ip = 1;
        for (size_t i = 0; i < fp.size(); ++i)
        {
            const T hp = half * (xp[i] + xp[i + 1]);
            //if (hp>x.front() &&  hp<x.back()) {
            if (hp > xh.front() && hp < xh.back())
            {
                while (hp > xh[ip])
                    ++ip;

                const T w = (xh[ip] - hp) / (xh[ip] - xh[ip - 1]);
                fp[i] += (one - w) * f[ip] + w * f[ip - 1];
            }
        }
    }

    // compute emissivity spectra using lookup table based on APED
    struct Apec
    {
        // null constructor
        Apec() {}

        // overloaded version of emission spectrum
        Apec(const std::string a_aped_path,
             const std::string a_aped_version,
             const Real a_energy_min,
             const Real a_energy_max,
             const int a_num_energy_bins,
             const std::string a_energy_bin_spacing,
             const std::string a_abundances_model = "AndersGrevesse",
             const std::string a_line_broadening = "convolution",
             const bool a_line_emission = true,
             const bool a_cont_emission = true,
             const int a_verbosity = 0)
            : m_spectrum_size(a_num_energy_bins),
              m_abundances_model(a_abundances_model),
              m_line_broadening(a_line_broadening),
              m_line_emission(a_line_emission),
              m_cont_emission(a_cont_emission),
              m_verbosity(a_verbosity)
        {
            assert( a_line_broadening == "none" ||
                   (a_line_broadening == "convolution" && a_energy_bin_spacing == "log") ||
                   (a_line_broadening == "linebyline" && a_energy_bin_spacing == "const"));

            // energy log spacing
            if (a_energy_bin_spacing == "log")
            {
                m_d_en = log10(a_energy_max / a_energy_min) / a_num_energy_bins;
                const unsigned num_buf_lo = (unsigned)ceil(-log10(one - MAX_DOPPLER_SHIFT) / m_d_en);
                const unsigned num_buf_hi = (unsigned)ceil(log10(one + MAX_DOPPLER_SHIFT) / m_d_en);
                // build energy
                for (int i = 0; i <= a_num_energy_bins; ++i)
                {
                    m_energy.emplace_back(a_energy_min * pow(ten, (i * m_d_en)));
                }
                // and buffer energy
                const Real min_buf_en = a_energy_min * pow(ten, -(num_buf_lo * m_d_en));

                for (size_t i = 0; i <= (num_buf_lo + a_num_energy_bins + num_buf_hi); ++i)
                {
                    m_buf_energy.emplace_back(min_buf_en * pow(ten, (i * m_d_en)));
                }
            }
            // energy const spacing
            else if (a_energy_bin_spacing == "const")
            {
                m_d_en = (a_energy_max - a_energy_min) / a_num_energy_bins;
                const unsigned num_buf_lo = (unsigned)ceil(a_energy_min * MAX_DOPPLER_SHIFT / m_d_en);
                const unsigned num_buf_hi = (unsigned)ceil(a_energy_max * MAX_DOPPLER_SHIFT / m_d_en);
                // build energy range and buffers
                for (int i = 0; i <= a_num_energy_bins; ++i)
                {
                    m_energy.emplace_back(a_energy_min + i * m_d_en);
                }
                // and buffer energy
                const Real min_buf_en = a_energy_min - num_buf_lo * m_d_en;
                for (size_t i = 0; i <= (num_buf_lo + a_num_energy_bins + num_buf_hi); ++i)
                {
                    m_buf_energy.emplace_back(min_buf_en + i * m_d_en);
                }
            }
            assert(m_spectrum_size == m_energy.size() - 1);

            // build aped
            Aped aped(a_aped_path, a_aped_version);
            // set temperature table
            m_temperature = aped.temperatures();
            m_dlog_temp = aped.temp_log_interv();

            // temperature range
            m_temp_min = (Real)*std::min_element(m_temperature.begin(), m_temperature.end());
            m_temp_max = (Real)*std::max_element(m_temperature.begin(), m_temperature.end());

            // separate BBN elements ...
            std::vector<unsigned> bbn_el(APED_atomic_numbers, APED_atomic_numbers + NBBNELEMENTS);
            build_emissivity_table(m_j_bbn_el, bbn_el, aped);

            // ... from metals
            std::vector<unsigned> metals(APED_atomic_numbers + NBBNELEMENTS, APED_atomic_numbers + NUM_APEC_ATOMS);
            build_emissivity_table(m_j_metals, metals, aped);

            if (m_verbosity > 0)
            {
                std::cout << "\nApec::Apec built successfully " << '\n';
            }
        }

        // copy constructor
        Apec(const Apec &) = default;

        // spectral energy
        std::vector<Real> spectral_energy() const
        {
            return m_energy;
        }

        // ph cm^3 s^-1
        void emission_spectrum(std::vector<Real> &a_spectrum,
                               const Real a_temperature,
                               const Real a_metallicity,
                               const Real a_doppler_shift) const
        {
            Timer_t<> t("APED::Apec::emission_spectrum");

            // no emissivity below temperature floor
            if (a_temperature > m_temp_min && a_temperature < m_temp_max)
            {
                // identify temperature bin
                const int it_lo = (int)floor(log10(a_temperature / m_temp_min) / m_dlog_temp);
                const int it_hi = std::min((size_t)it_lo + 1, m_temperature.size() - 1);
                assert(a_temperature >= (Real)m_temperature[it_lo] && a_temperature <= (Real)m_temperature[it_hi]);

                // temperature interpolation coefficient
                const Real f = std::abs(log10(a_temperature / (Real)m_temperature[it_lo])) / m_dlog_temp;

                // multiply by ne so that the spectrum remains simply normalized to nH^2
                const Real ne = 1.0;//m_sd.n_e(a_temperature, a_metallicity);

                // rename
                const Real z = a_metallicity;

                const unsigned buf_spectrum_size = m_buf_energy.size() - 1;
                std::vector<Real> js(buf_spectrum_size);
                {
                    Timer_t<> t("APED::Apec::emission_spectrum:temp_interpolation");
                    for (size_t i = 0; i < buf_spectrum_size; ++i)
                    {
                        js[i] = ne * ((one - f) * (m_j_bbn_el[it_lo][i] + z * m_j_metals[it_lo][i]) + f * (m_j_bbn_el[it_hi][i] + z * m_j_metals[it_hi][i]));
                    }
                }

                // Interpolate from rest frame to lab frame spectrum
                {
                    Timer_t<> t("APED::Apec::emission_spectrum:ph_en_interpolation");
                    std::vector<Real> shifted_energy(m_buf_energy.size());

                    const Real df = one + a_doppler_shift;
                    for (size_t i = 0; i < m_buf_energy.size(); ++i)
                    {
                        shifted_energy[i] = df * m_buf_energy[i];
                    }

                    a_spectrum.resize(m_spectrum_size, zero);
                    linear_interp(a_spectrum, m_energy, js, shifted_energy);
                }
            }
            else
            {
                Timer_t<> t("APED::Apec::emission_spectrum:out_of_temp_range");
                a_spectrum.resize(m_spectrum_size, zero);
            }
        }

        // ph cm^3 s^-1
        Real emissivity(const Real a_temperature, const Real a_metallicity, const Real a_doppler_shift)
        {
            Timer_t<> t("APED::Apec::emissivity");

            // call aped
            std::vector<Real> spectrum;
            emission_spectrum(spectrum, a_temperature, a_metallicity, a_doppler_shift);

            Real sum = zero;
            for (const auto &s : spectrum)
                sum += s;

            return sum;
        }

    protected:
        //
        void build_emissivity_table(std::vector<std::vector<Real>> &a_jt,
                                    const std::vector<unsigned> &a_elements,
                                    const Aped &a_aped)
        {
            // here use these dummy values and rescaled properly later
            Real dummy_metallicity = 1.0;
            Real no_doppler_shift = 0.0;

            // resize table
            a_jt.clear();

            // compute relative abundances
            std::vector<ElementAbundance> rel_abundances;
            AbundanceUtil::relative_abundances(rel_abundances, a_elements, dummy_metallicity, m_abundances_model);

            // temperature dimensions
            a_jt.resize(m_temperature.size());

            // build call aped
            auto it = m_temperature.begin();
            for (auto &jt : a_jt)
            {
                std::vector<Real> j;
                a_aped.emission_spectrum(j,
                                         m_buf_energy,
                                         rel_abundances,
                                         *it,
                                         no_doppler_shift,
                                         m_line_broadening,
                                         m_line_emission,
                                         m_cont_emission);
                // move emissivity to table
                jt = std::move(j);
                // push temperature iterator
                ++it;
            }
        }

    protected:
        // tables
        std::vector<std::vector<Real>> m_j_bbn_el;
        std::vector<std::vector<Real>> m_j_metals;
        // SD object to get ne(T,Z)
        fm::cooling_tables::SutherlandDopita m_sd;
        // energy, buffered energy and temperature grid
        std::vector<Real> m_energy;
        std::vector<Real> m_buf_energy;
        std::vector<double> m_temperature;

    protected:
        // size of spectrum
        unsigned m_spectrum_size;
        // energy and (log of) temperature intervals
        Real m_d_en;
        double m_dlog_temp;
        // temp bounds
        double m_temp_min, m_temp_max;
        // emissivity features
        std::string m_abundances_model;
        std::string m_line_broadening;
        bool m_line_emission;
        bool m_cont_emission;
        int m_verbosity;
    };

} // namespace fm::aped

#endif // APED_H
#endif // USE_APED
