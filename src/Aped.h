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

#ifndef APED_H
#define APED_H

#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cstdlib>
#include <cmath>
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

    constexpr double zero = 0.0e0;
    constexpr double half = 0.5e0;
    constexpr double one = 1.0e0;
    constexpr double two = 2.0e0;
    constexpr double ten = 1.0e1;
    constexpr double sqrt_two = std::sqrt(two);

    constexpr double TOLERANCE = (1.e-5);

    // aped database contains data for 28 species
    constexpr unsigned NUM_APEC_ATOMS = 30;

    // Boltzmann's constant
    constexpr double k_B = 1.3806488e-16;
    // conversion from keV to Angstrom
    constexpr double keVToAngstrom = 12.39854;
    // conversion from Angstrom to keV
    constexpr double AngstromTokeV = 0.080654657725829;
    // conversion from keV to Kelvin
    constexpr double keVToKelvin = 1.1604505e7;
    // kT in erg corresponding to temperature of 1 keV
    constexpr double keVToerg = 1.60217646e-9;

    // Speed of light in cm s^-1
    constexpr double C_cm_s = 2.9979246e10;
    // Unified atomic mass constant in g
    constexpr double AMU_g = 1.660538e-24;
    // km in cm
    constexpr double KM_cm = 1.0e5;
    // equivalent to v/c with v=15000 km/sec
    constexpr double MAX_DOPPLER_SHIFT = (0.05);

    // virtual base class for kernel
    struct LineKernel
    {
        LineKernel() {}
        ~LineKernel() {}

        virtual int wing_size() const = 0;
        virtual Real val(const int index) const = 0;
    };

    // gaussian kernel
    struct GaussianKernel : public LineKernel
    {
        GaussianKernel(const Real e_line, const Real e_lo, const Real e_hi)
        {
            Timer_t<4> t("Aped::GaussianKernel");

            // builds kernel with size sufficient to have negligible residuals in the wings
            const double de = e_hi - e_lo;
            const double de_hi = e_hi - e_line;
            const double de_lo = e_line - e_lo;

            std::vector<double> wm, wp;

            double err = one;
            unsigned i = 0;
            while (err > TOLERANCE)
            {
                wm.emplace_back(half * erf(de_lo + i * de));
                wp.emplace_back(half * erf(de_hi + i * de));
                err = abs(one - (wm[i] + wp[i]));
                ++i;
            }
            const unsigned N = i;
            // now go backward, and get the proper weights by subtracting successive erfs
            // first left wing
            while (--i > 0)
            {
                m_w.emplace_back(wm[i] - wm[i - 1]);
            }
            assert(i == 0);
            // centre
            m_w.emplace_back(wm[i] + wp[i]);
            // right wing
            while (++i < N)
            {
                m_w.emplace_back(wp[i] - wp[i - 1]);
            }
            // sanity checks
            assert(m_w.size() % 2 == 1);
            m_wing_size = (m_w.size() - 1) / 2;

            // shift index from N -> N-1
            double sum = m_w[--i];
            // loop over wings
            while (i-- > 0)
            {
                sum += m_w[i] + m_w[N + i];
            }
            assert(abs(one - sum) < TOLERANCE);
        }

        virtual ~GaussianKernel() {}

        virtual int wing_size() const
        {
            return m_wing_size;
        }

        // peak is centered at m_wing_size
        virtual Real val(const int a_index) const
        {
            return (Real)m_w.at(a_index + m_wing_size);
        }

    protected:
        size_t m_wing_size;
        std::vector<double> m_w;
    };

    //
    // data
    static const unsigned APED_atomic_numbers[NUM_APEC_ATOMS] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    static const char *APED_atom_names[NUM_APEC_ATOMS] = {
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
        "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"};
    static const float APED_atomic_masses[NUM_APEC_ATOMS] = {
        1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
        22.989770, 24.3050, 26.981538, 28.0855, 30.973761, 32.065, 35.4527, 39.948, 39.0983, 40.078,
        44.955910, 47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.933200, 58.6934, 63.456, 65.39};

    static const float AndersGrevesseAbundances[NUM_APEC_ATOMS] = {
        12.0, 10.99, 1.16, 1.15, 2.60, 8.56, 8.05, 8.93, 4.56, 8.09, 6.33, 7.58, 6.47, 7.55, 5.45,
        7.21, 5.50, 6.56, 5.12, 6.36, 3.10, 4.99, 4.00, 5.67, 5.39, 7.67, 4.92, 6.25, 4.21, 4.60};

    static const float LoddersAbundances[NUM_APEC_ATOMS] = {
        12.0, 10.984, 3.35, 1.48, 2.85, 8.46, 7.90, 8.76, 4.53, 7.95, 6.37, 7.62, 6.54, 7.61, 5.54,
        7.26, 5.33, 6.62, 5.18, 6.41, 3.15, 5.00, 4.07, 5.72, 5.58, 7.54, 4.98, 6.29, 4.34, 4.70};

    // BBN elements
    constexpr unsigned NBBNELEMENTS =4;
    // static const char *BBNElements[NBBNELEMENTS] = {"H", "He", "Li", "Be"};

    // single element abundance
    struct ElementAbundance
    {
        // null
        ElementAbundance() {}
        // full constructor
        ElementAbundance(const unsigned a_n, const Real a_a)
            : atomic_number(a_n), abundance(a_a)
        {}
        // copy
        ElementAbundance(const ElementAbundance &) = default;

        ~ElementAbundance() {}

        unsigned atomic_number;
        Real abundance;
    };

    // Compute abundances relative to Anders and Grevesse which is assumed in Aped emissivities
    struct AbundanceUtil
    {
        // constructor
        AbundanceUtil() {}
        // destructor
        ~AbundanceUtil() {}

        // compute abundances relative to AndersGrevesse
        static void relative_abundances(std::vector<ElementAbundance> &a_relative_abundances,
                                        const std::vector<unsigned> &a_elements,
                                        const Real a_metallicity,
                                        const std::string a_abundances_model)
        {
            // make sure elements are within range
            for (const auto a : a_elements)
            {
                if (a > NUM_APEC_ATOMS)
                    throw std::runtime_error("AbundanceUtil:: atomic number is out of range ! ");
            }

            // if elements.size()==0, use all elements
            std::vector<unsigned> elements;
            if (a_elements.size() > 0)
                elements.assign(a_elements.begin(), a_elements.end());
            else
                elements.assign(APED_atomic_numbers, APED_atomic_numbers + NUM_APEC_ATOMS);

            // reference abundance values
            std::vector<float> ref_ab(AndersGrevesseAbundances, AndersGrevesseAbundances + NUM_APEC_ATOMS);

            // model abundance values
            std::vector<float> mod_ab;
            if (a_abundances_model == "AndersGrevesse")
                mod_ab.assign(AndersGrevesseAbundances, AndersGrevesseAbundances + NUM_APEC_ATOMS);
            else if (a_abundances_model == "Lodders")
                mod_ab.assign(LoddersAbundances, LoddersAbundances + NUM_APEC_ATOMS);
            else
                throw std::runtime_error("abundance model " + a_abundances_model + " is not recongnised ! ");

            a_relative_abundances.clear();
            for (const auto &a : elements)
            {
                const Real abundance = pow(ten, mod_ab[a - 1] - ref_ab[a - 1]) * (a <= NBBNELEMENTS ? 1 : a_metallicity);
                a_relative_abundances.emplace_back(ElementAbundance(a, abundance));
            }
        }
    }; // namespace APED

    // Ion Entry
    struct Ion
    {
        // null constructor
        Ion() {}
        // copy constructor
        Ion(const Ion &) = default;
        // destructor
        ~Ion() {}
        // partial constructor (set ion species)
        Ion(const unsigned a_ion)
            : ion(a_ion)
        {}
        // assignemt operator
        Ion &operator=(const Ion &) = default;

        unsigned ion;
        std::vector<float> cont_energy;
        std::vector<float> continuum;
        std::vector<float> continuum_err;
        std::vector<float> pseudo_cont_energy;
        std::vector<float> pseudo_cont;
        std::vector<float> pseudo_cont_err;
        std::vector<float> line_energy;
        std::vector<float> line_energy_err;
        std::vector<float> line_emissivity;
        std::vector<float> line_emissivity_err;
        std::vector<int> elem_driver;
        std::vector<int> ion_driver;
        std::vector<int> lower_level;
        std::vector<int> upper_level;
    };

    //
    struct Element
    {
        // constructor
        Element() {}
        // destructor
        ~Element() {}
        // copy constructors
        Element(const Element &) = default;
        // partial constructor
        Element(const std::string a_name, const unsigned a_A, const float a_M)
            : name(a_name), atomic_number(a_A), atomic_mass(a_M)
        {}
        // assignment operator
        Element &operator=(const Element &) = default;

        void check_ion(const unsigned a_ion)
        {
            if (ions.find(a_ion) == ions.end())
            {
                ions.insert(std::pair<unsigned, Ion>(a_ion, Ion(a_ion)));
            }
        }

        unsigned num_ions() const
        {
            return ions.size();
        }

        //
        std::string name;
        unsigned atomic_number;
        float atomic_mass;

        // map ionization state to ion
        std::map<unsigned, Ion> ions;
    };

    // aped database
    struct TemperatureRecord
    {
        // null constructor
        TemperatureRecord() {}
        // copy
        TemperatureRecord(const TemperatureRecord &) = default;
        // destructor
        ~TemperatureRecord() {}
        // assignment operator
        TemperatureRecord &operator=(const TemperatureRecord &) = default;

        // add new element
        void check_element(const unsigned a_atomic_number)
        {
            if (elements.find(a_atomic_number) == elements.end())
            {
                assert(a_atomic_number <= NUM_APEC_ATOMS);

                Element E(APED_atom_names[a_atomic_number - 1],
                          a_atomic_number,
                          APED_atomic_masses[a_atomic_number - 1]);

                elements.insert(std::pair<unsigned, Element>(a_atomic_number, E));
            }
        }

        // data
        Real temperature;
        // map atomic number to element
        std::map<unsigned, Element> elements;
    };

    //
    struct Aped
    {
        // null constructor
        Aped()
            : m_verbosity(0)
        {}
        // full constructor
        Aped(const std::string a_aped_path, const std::string a_version, const int a_verbosity = 0);
        // copy constructor
        Aped(const Aped &) = default;
        // destructor
        ~Aped() {}

        // number of tabulated temperatures
        int num_temperatures() const
        {
            return m_temperatures.size();
        }
        // tabulated temperatures
        std::vector<Real> temperatures() const
        {
            return m_temperatures;
        }
        // temperature table step
        double temp_log_interv() const
        {
            return m_dlog_temp;
        }

        //
        int num_elements() const
        {
            return m_num_elements_line.size();
        }

        // photon emission spectrum in ph cm^3 s^-1
        void emission_spectrum(std::vector<Real> &a_spectrum,
                               const std::vector<Real> &a_energy,
                               const std::vector<unsigned> &a_elements, // can be empty to mean all
                               const std::string a_abundances_model,
                               const Real a_metallicity,
                               const Real a_temperature,
                               const Real a_doppler_shift,
                               const std::string a_line_broadening,
                               const bool a_line_emission,
                               const bool a_cont_emission) const
        {
            // utility to compute abundances re
            std::vector<ElementAbundance> rel_ab;
            AbundanceUtil::relative_abundances(rel_ab, a_elements, a_metallicity, a_abundances_model);

            // compute spectrum
            this->emission_spectrum(a_spectrum, a_energy, rel_ab, a_temperature,
                                    a_doppler_shift, a_line_broadening,
                                    a_line_emission, a_cont_emission);
        }

        // overloaded version of emission spectrum in ph cm^3 s^-1
        void emission_spectrum(std::vector<Real> &a_spectrum,
                               const std::vector<Real> &a_energy,
                               const std::vector<ElementAbundance> &a_atom_abundances,
                               const Real a_temperature,
                               const Real a_doppler_shift,
                               const std::string a_line_broadening,
                               const bool a_line_emission,
                               const bool a_cont_emission) const;

    protected:
        // single ion emission line spectrum
        void ion_line_emission(std::vector<Real> &a_spectrum,
                               const std::vector<Real> &a_energy,
                               const int a_atomic_number,
                               const int a_rmJ,
                               const int a_temp_idx,
                               const Real a_doppler_shift,
                               const std::string a_line_broadening) const;
        // continuum
        void ion_continuum_emission(std::vector<Real> &a_spectrum,
                                    const std::vector<Real> &a_energy,
                                    const int a_atomic_number,
                                    const int a_rmJ,
                                    const int a_temp_idx,
                                    const Real a_doppler_shift,
                                    const bool a_cont_emission,
                                    const bool a_pseudo_cont_emission) const;

    protected:
        // do simple convolution with gaussian weights
        void simple_convolution(std::vector<Real> &a_spectrum,
                                const LineKernel *a_kernel) const
        {
            Timer_t<4> t("Aped::simple_convolution");

            const int kwing = a_kernel->wing_size();

            std::vector<Real> convolution(a_spectrum.size(), 0);

            for (size_t i = 0; i < a_spectrum.size(); ++i)
            {
                // don't convolve for nothing
                if (a_spectrum[i] > 0)
                {
                    for (int j = -kwing; j <= kwing; ++j)
                    {
                        if (i + j >= 0 && i + j < a_spectrum.size())
                            convolution[i + j] += a_spectrum[i] * a_kernel->val(j);
                    }
                }
            }
            // replace spectrum with convolution
            a_spectrum.assign(convolution.begin(), convolution.end());
        }

        // add gaussian profile of line emission
        void line_profile(std::vector<Real> &a_spectrum,
                          std::vector<Real> &a_energy,
                          const Real a_line_energy,
                          const Real a_line_emissivity,
                          const LineKernel *a_kernel) const
        {
            Timer_t<4> t("Aped::line_profile");

            if (a_line_energy >= a_energy.front() && a_line_energy < a_energy.back())
            {
                // find line energy bin
                int ie = 0;
                while (a_energy[ie] < a_line_energy)
                    ++ie;
                ie = std::max(0, ie - 1);

                // use gaussian kernel
                const size_t kwing = a_kernel->wing_size();

                for (size_t j = -kwing; j <= kwing; ++j)
                {
                    if (ie + j >= 0 && ie + j < a_spectrum.size())
                        a_spectrum[ie + j] += a_line_emissivity * a_kernel->val(j);
                }
            }
        }

    protected:
        // temperature data
        std::vector<TemperatureRecord> m_aped_data;
        // log step of temp database
        double m_dlog_temp;
        // verbosity, zero by default
        int m_verbosity;

        std::vector<double> m_temperatures;
        std::vector<double> m_density;
        std::vector<int> m_num_lines;
        std::vector<int> m_num_elements_line;
        std::vector<int> m_num_elements_coco;
        std::vector<int> m_num_continuum;
        std::vector<int> m_num_pseudo;
    };

} // namespace fm::aped

#endif // _APED_
#endif // USE_APED
