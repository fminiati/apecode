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

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <map>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "Util.H"

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
    // Unified atomic mass constant in g
    constexpr Real AMU_cgs = 1.660539040e-24; 

    // aped database contains data for 28 species
    constexpr unsigned NUM_APED_ATOMS = 30;
    // data
    static const unsigned atomic_numbers[NUM_APED_ATOMS] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    static const char *atom_names[NUM_APED_ATOMS] = {
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
        "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"};
    static const Real atomic_masses[NUM_APED_ATOMS] = {
        1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
        22.989770, 24.3050, 26.981538, 28.0855, 30.973761, 32.065, 35.4527, 39.948, 39.0983, 40.078,
        44.955910, 47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.933200, 58.6934, 63.456, 65.39};
    static const Real AndersGrevesseAbundances[NUM_APED_ATOMS] = {
        12.0, 10.99, 1.16, 1.15, 2.60, 8.56, 8.05, 8.93, 4.56, 8.09, 6.33, 7.58, 6.47, 7.55, 5.45,
        7.21, 5.50, 6.56, 5.12, 6.36, 3.10, 4.99, 4.00, 5.67, 5.39, 7.67, 4.92, 6.25, 4.21, 4.60};
    static const Real LoddersAbundances[NUM_APED_ATOMS] = {
        12.0, 10.984, 3.35, 1.48, 2.85, 8.46, 7.90, 8.76, 4.53, 7.95, 6.37, 7.62, 6.54, 7.61, 5.54,
        7.26, 5.33, 6.62, 5.18, 6.41, 3.15, 5.00, 4.07, 5.72, 5.58, 7.54, 4.98, 6.29, 4.34, 4.70};
    // BBN elements
    constexpr unsigned NBBNELEMENTS =4;
    // static const char *BBNElements[NBBNELEMENTS] = {"H", "He", "Li", "Be"};

    // single element abundance
    struct ElementAbundance
    {
        ElementAbundance() {}
        ~ElementAbundance() {}
        ElementAbundance(const unsigned a_n, const Real a_a)
            : m_atomic_number(a_n), m_abundance(a_a) {}
        ElementAbundance(const ElementAbundance &) = default;

        unsigned m_atomic_number;
        Real m_abundance;
    };

    // Compute abundances relative to Anders and Grevesse which is assumed in Aped emissivities
    struct AbundanceUtil
    {
        // compute abundances relative to AndersGrevesse
        static void relative_abundances(std::vector<ElementAbundance> &a_relative_abundances,
                                        const std::vector<unsigned> &a_elements,
                                        const Real a_metallicity,
                                        const std::string a_abundances_model)
        {
            // make sure elements are within range
            for (const auto a : a_elements)
            {
                if (a > NUM_APED_ATOMS)
                    throw std::runtime_error("AbundanceUtil:: atomic number is out of range ! ");
            }

            // if elements.size()==0, use all elements
            std::vector<unsigned> elements;
            if (a_elements.size() > 0)
                elements.assign(a_elements.begin(), a_elements.end());
            else
                elements.assign(atomic_numbers, atomic_numbers + NUM_APED_ATOMS);

            // reference abundance values
            std::vector<Real> ref_ab(AndersGrevesseAbundances, AndersGrevesseAbundances + NUM_APED_ATOMS);

            // model abundance values
            std::vector<Real> mod_ab;
            if (a_abundances_model == "AndersGrevesse")
                mod_ab.assign(AndersGrevesseAbundances, AndersGrevesseAbundances + NUM_APED_ATOMS);
            else if (a_abundances_model == "Lodders")
                mod_ab.assign(LoddersAbundances, LoddersAbundances + NUM_APED_ATOMS);
            else
                throw std::runtime_error("abundance model " + a_abundances_model + " is not recongnised ! ");

            a_relative_abundances.clear();
            for (const auto &a : elements)
            {
                const Real abundance = pow(ten, mod_ab[a - 1] - ref_ab[a - 1]) * (a <= NBBNELEMENTS ? 1 : a_metallicity);
                a_relative_abundances.emplace_back(ElementAbundance(a, abundance));
            }
        }
    };

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
            : m_ion(a_ion)
        {}
        // assignemt operator
        Ion &operator=(const Ion &) = default;

        unsigned m_ion;
        std::vector<Real> m_cont_energy;
        std::vector<Real> m_continuum;
        std::vector<Real> m_continuum_err;
        std::vector<Real> m_pseudo_cont_energy;
        std::vector<Real> m_pseudo_cont;
        std::vector<Real> m_pseudo_cont_err;
        std::vector<Real> m_line_energy;
        std::vector<Real> m_line_energy_err;
        std::vector<Real> m_line_emissivity;
        std::vector<Real> m_line_emissivity_err;
        std::vector<int> m_elem_driver;
        std::vector<int> m_ion_driver;
        std::vector<int> m_lower_level;
        std::vector<int> m_upper_level;
    };

    struct Element
    {
        // constructor
        Element() {}
        // destructor
        ~Element() {}
        // copy constructors
        Element(const Element &) = default;
        // partial constructor
        Element(const std::string a_name, const unsigned a_A, const Real a_M)
            : m_name(a_name), m_atomic_number(a_A), m_atomic_mass(a_M)
        {}
        // assignment operator
        Element &operator=(const Element &) = default;

        void check_ion(const unsigned a_ion)
        {
            if (m_ions.find(a_ion) == m_ions.end())
            {
                m_ions.insert(std::pair<unsigned, Ion>(a_ion, Ion(a_ion)));
            }
        }

        unsigned num_ions() const
        {
            return m_ions.size();
        }

        std::string m_name;
        unsigned m_atomic_number;
        Real m_atomic_mass;

        // map ionization state to ion
        std::map<unsigned, Ion> m_ions;
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
            if (m_elements.find(a_atomic_number) == m_elements.end())
            {
                assert(a_atomic_number <= NUM_APED_ATOMS);

                Element E(atom_names[a_atomic_number - 1],
                          a_atomic_number,
                          atomic_masses[a_atomic_number - 1]);

                m_elements.insert(std::pair<unsigned, Element>(a_atomic_number, E));
            }
        }

        // data
        Real m_temperature;
        // map atomic number to element
        std::map<unsigned, Element> m_elements;
    };

    struct Aped
    {
        // null constructor
        Aped()
            : m_verbosity{0}, m_energy_spacing{Spacing::undetermined}
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
        Real temp_log_interv() const
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
                               const bool a_cont_emission,
                               const bool a_line_emission,
                               const LineShape a_line_profile,
                               const LineBroadening a_line_broadening,
                               const Real a_kernel_tolerance = 1.e-6) const
        {
            // utility to compute abundances re
            std::vector<ElementAbundance> rel_ab;
            AbundanceUtil::relative_abundances(rel_ab, a_elements, a_metallicity, a_abundances_model);

            // compute spectrum
            switch (a_line_profile)
            {
                case LineShape::delta :
                    emission_spectrum<LineProfile<Delta, NoBroadening>>(a_spectrum, a_energy, rel_ab,
                                                                        a_temperature, a_doppler_shift,
                                                                        a_cont_emission, a_line_emission,
                                                                        a_kernel_tolerance);
                    break;
                case LineShape::gaussian:
                    emission_spectrum<LineProfile<Gaussian, ThermalBroadening>>(a_spectrum, a_energy, rel_ab,
                                                                                a_temperature, a_doppler_shift,
                                                                                a_cont_emission, a_line_emission,
                                                                                a_kernel_tolerance);
                    break;
                case LineShape::lorentzian:
                    emission_spectrum<LineProfile<Lorentzian, ThermalBroadening>>(a_spectrum, a_energy, rel_ab,
                                                                                  a_temperature, a_doppler_shift,
                                                                                  a_cont_emission, a_line_emission,
                                                                                  a_kernel_tolerance);
                    break;
                case LineShape::pseudovoigt:
                    using LP = LineProfile<PseudoVoigt,PseudoVoigtBroadening<ThermalBroadening,ThermalBroadening>>;
                    //using PVB = PseudoVoigtBroadening;
                    emission_spectrum<LP>(a_spectrum, a_energy, rel_ab,
                                          a_temperature, a_doppler_shift,
                                          a_cont_emission, a_line_emission, a_kernel_tolerance);
                    break;
                default:
                    std::cerr << "\n LineShape " << static_cast<int>(a_line_profile) << "was not recognised, aped will return!\n\n";
            }
        }


        // photon emission spectrum in ph cm^3 s^-1
        template <typename Profile>
        void emission_spectrum(std::vector<Real> &a_spectrum,
                               const std::vector<Real> &a_energy,
                               const std::vector<unsigned> &a_elements, // can be empty to mean all
                               const std::string a_abundances_model,
                               const Real a_metallicity,
                               const Real a_temperature,
                               const Real a_doppler_shift,
                               const bool a_cont_emission,
                               const bool a_line_emission,
                               const Real a_kernel_tolerance) const
        {
            // utility to compute abundances re
            std::vector<ElementAbundance> rel_ab;
            AbundanceUtil::relative_abundances(rel_ab, a_elements, a_metallicity, a_abundances_model);

            // compute spectrum
            emission_spectrum<Profile>(a_spectrum, a_energy, rel_ab,
                                       a_temperature, a_doppler_shift,
                                       a_cont_emission, a_line_emission,
                                       a_kernel_tolerance);
        }

        // overloaded version of emission spectrum in ph cm^3 s^-1
        template <typename Profile>
        void emission_spectrum(std::vector<Real> &a_spectrum,
                               const std::vector<Real> &a_energy,
                               const std::vector<ElementAbundance> &a_atom_abundances,
                               const Real a_temperature,
                               const Real a_doppler_shift,
                               const bool a_cont_emission,
                               const bool a_line_emission,
                               const Real a_kernel_tolerance) const;

    protected:
        // single ion emission line spectrum
        template <typename Profile>
        void ion_line_emission(std::vector<Real> &a_spectrum,
                               const std::vector<Real> &a_energy,
                               const int a_atomic_number,
                               const int a_rmJ,
                               const int a_temp_bin,
                               const Real a_doppler_shift,
                               const Real a_kernel_tolerance) const;
        // continuum
        template <typename Profile>
        void ion_continuum_emission(std::vector<Real> &a_spectrum,
                                    const std::vector<Real> &a_energy,
                                    const int a_atomic_number,
                                    const int a_rmJ,
                                    const int a_temp_bin,
                                    const Real a_doppler_shift,
                                    const bool a_cont_emission,
                                    const bool a_pseudo_cont_emission,
                                    const Real a_kernel_tolerance) const;

    protected:
        // temperature data
        std::vector<TemperatureRecord> m_aped_data;
        // log step of temp database
        Real m_dlog_temp;
        // verbosity, zero by default
        int m_verbosity;
        // spacing of spectral energy nodes
        mutable spacing_t m_energy_spacing;

        std::vector<Real> m_temperatures;
        std::vector<Real> m_density;
        std::vector<int> m_num_lines;
        std::vector<int> m_num_elements_line;
        std::vector<int> m_num_elements_coco;
        std::vector<int> m_num_continuum;
        std::vector<int> m_num_pseudo;
    };

    //
    // Implementation of function templates
    //

    // continuum spectrum including pseudo continuum
    template <typename Profile>
    void Aped::emission_spectrum(std::vector<Real> &a_spectrum,
                                 const std::vector<Real> &a_energy,
                                 const std::vector<ElementAbundance> &a_atom_abundances,
                                 const Real a_temperature,
                                 const Real a_doppler_shift,
                                 const bool a_cont_emission,
                                 const bool a_line_emission,
                                 const Real a_kernel_tolerance) const
    {
        static_assert(std::conjunction_v<std::is_same<line_shape_t<Profile>, Delta>,
                                         std::is_same<line_broadening_t<Profile>, NoBroadening>> ||
                          std::negation_v<std::disjunction<std::is_same<line_shape_t<Profile>, Delta>,
                                                           std::is_same<line_broadening_t<Profile>, NoBroadening>>>,
                      "Inconsistent Line Profile: choose Delta profile if and only if NoBroadening is the broadening mechanism.");
        Timer_t<> t("Aped::emission_spectrum");

        // initialize spectrum
        a_spectrum.resize(a_energy.size() - 1);
        std::memset(&a_spectrum[0], 0, sizeof(Real) * a_spectrum.size());

        // check temperature range
        const Real T_min = (Real)*std::min_element(m_temperatures.begin(), m_temperatures.end());
        const Real T_max = (Real)*std::max_element(m_temperatures.begin(), m_temperatures.end());

        // set spacing of energy grid
        m_energy_spacing = Spacing::grid_spacing(a_energy);

        // if T is within allowed range
        if (a_temperature >= T_min && a_temperature <= T_max)
        {
            // identify right bin
            const int temp_bin_lo = (int)floor(log10(a_temperature / T_min) / m_dlog_temp);
            const int temp_bin_hi = std::min((size_t)temp_bin_lo + 1, m_temperatures.size() - 1);
            assert(a_temperature >= (Real)m_temperatures[temp_bin_lo] && a_temperature <= (Real)m_temperatures[temp_bin_hi]);

            // loop over all elements
            for (const auto &A : a_atom_abundances)
            {
                // loop over temp bins
                for (int temp_bin = temp_bin_lo; temp_bin <= temp_bin_hi; ++temp_bin)
                {
                    // Linear interpolation coefficient as adopted by original Aped code
                    const Real f = one - std::abs(a_temperature - (Real)m_temperatures[temp_bin]) / (m_temperatures[temp_bin_lo] * (pow(ten, m_dlog_temp) - one));

                    // add ion emission to spectrum taking into account io abundance and ionization fraction
                    auto add_emission_to_spectrum = [x = f * A.m_abundance](std::vector<Real> &j, const std::vector<Real> &i) {
                        for (size_t k = 0; k < j.size(); ++k)
                            j[k] += x * i[k];
                    };

                    if (a_line_emission)
                    {
                        Timer_t<2> t("Aped::emission_spectrum:line_emission");

                        std::vector<Real> line_emission;
                        ion_line_emission<Profile>(line_emission,
                                                   a_energy,
                                                   A.m_atomic_number,
                                                   0,
                                                   temp_bin,
                                                   a_doppler_shift,
                                                   a_kernel_tolerance);

                        // add ion line emission
                        add_emission_to_spectrum(a_spectrum, line_emission);
                    }

                    if (a_cont_emission || a_line_emission)
                    {
                        Timer_t<2> t("Aped::emission_spectrum:cont_emission");

                        // continuum
                        std::vector<Real> emission;
                        ion_continuum_emission<Profile>(emission,
                                                        a_energy,
                                                        A.m_atomic_number,
                                                        0, temp_bin,
                                                        a_doppler_shift,
                                                        a_cont_emission,
                                                        a_line_emission,
                                                        a_kernel_tolerance);

                        // add true and/or psd cont emission
                        add_emission_to_spectrum(a_spectrum, emission);
                    }
                }
            }
        }
    }

    template <typename Profile>
    void Aped::ion_line_emission(std::vector<Real> &a_spectrum,
                                 const std::vector<Real> &a_energy,
                                 const int a_atomic_number,
                                 const int a_rmJ,
                                 const int a_temp_bin,
                                 const Real a_doppler_shift,
                                 const Real a_kernel_tolerance) const
    {
        // resize spectrum vector
        a_spectrum.resize(a_energy.size() - 1);
        std::memset(&a_spectrum[0], 0, sizeof(Real) * a_spectrum.size());

        // count num lines
        size_t num_lines = 0;

        // find the element at the input temperature bin
        if (const auto el_iterator = m_aped_data[a_temp_bin].m_elements.find(a_atomic_number);
            el_iterator != m_aped_data[a_temp_bin].m_elements.end())
        {
            const Element &atom = el_iterator->second;

            // if ionization state (rmJ) == 0 then add up all ions, else select ionization state according to input
            const auto [beg, end] = (a_rmJ == 0 ? std::pair(atom.m_ions.cbegin(), atom.m_ions.cend()) : std::pair(atom.m_ions.find(a_rmJ), std::next(atom.m_ions.find(a_rmJ))));
            for (auto ion_iterator = beg; ion_iterator != end; ++ion_iterator)
            {
                // set ion pointer
                const Ion &ion = ion_iterator->second;
                assert(ion.m_line_emissivity.size() == ion.m_line_energy.size());

                // loop thorugh emission lines
                int energy_bin = 0;
                for (size_t i_line = 0; i_line < ion.m_line_emissivity.size(); ++i_line)
                {
                    // doppler shifted photon energy
                    const Real line_energy = ion.m_line_energy[i_line] * (one + a_doppler_shift);

                    if (line_energy >= a_energy.front() && line_energy < a_energy.back())
                    {
                        ++num_lines;

                        while (a_energy[energy_bin] < line_energy)
                            ++energy_bin;
                        energy_bin = (energy_bin > 0 ? energy_bin - 1 : energy_bin);

                        if constexpr (std::is_same_v<line_shape_t<Profile>,Delta>)
                        {
                            a_spectrum[energy_bin] += ion.m_line_emissivity[i_line];
                        }
                        else
                        {
                            Timer_t<4> t("Aped::ion_line_emission:convolution");
                            const Real line_fwhm = Profile::fwhm(m_temperatures[a_temp_bin],
                                                                 AMU_cgs * atom.m_atomic_mass, line_energy);
                            Convolution::convolve<Profile>(a_spectrum, (Real)ion.m_line_emissivity[i_line],
                                                           line_energy, line_fwhm, energy_bin,
                                                           a_energy, a_kernel_tolerance);
                        }
                    }
                }
            }
        }
        else
        {
            std::cerr << " strangely element " << a_atomic_number << " was not found " << std::endl;
        }

        if (m_verbosity > 0)
        {
            std::cout << " added " << num_lines << " emission lines to the spectrum " << std::endl;
        }
    }

    // continuum
    template <typename Profile>
    void Aped::ion_continuum_emission(std::vector<Real> &a_continuum,
                                      const std::vector<Real> &a_energy,
                                      const int a_atomic_number,
                                      const int a_rmJ,
                                      const int a_temp_bin,
                                      const Real a_doppler_shift,
                                      const bool a_cont_emission,
                                      const bool a_psd_cont_emission,
                                      const Real a_kernel_tolerance) const
    {
        auto add_cont_emission_to_spectrum = [&j = a_continuum, &e = a_energy](const std::vector<Real> &js,
                                                                               const std::vector<Real> &es) {
            Timer_t<4> t("Aped::ion_continuum_emission::add_to_spectrum");

            if (es.back() < e.front() && es.front() > e.back())
                return;

            size_t k = 0;
            Real ef = e[0], jf = js[0]; // foot point values
            for (size_t i = 0; i < j.size(); ++i)
            {
                if (es.back() > e[i])
                {
                    while (es[k] < e[i])
                        ++k;

                    // reset foot emissivity if need be
                    if (i == 0 && k > 0)
                        jf = js[k - 1] + (js[k] - js[k - 1]) / (es[k] - es[k - 1]) * (ef - es[k - 1]);

                    // loop through source contributions within this e-bin
                    while (es[k] <= e[i + 1] && k < js.size())
                    {
                        j[i] += half * (jf + js[k]) * (es[k] - ef);
                        jf = js[k];
                        ef = es[k];
                        ++k;
                    }
                    // add final contribution and reset foot values for next e-bin
                    if (k < js.size()) // --> es[k] > e[i + 1]
                    {
                        const Real jh = jf + (js[k] - jf) / (es[k] - ef) * (e[i + 1] - ef);
                        j[i] += half * (jf + jh) * (e[i + 1] - ef);
                        jf = jh;
                        ef = e[i + 1];
                    }
                }
            }
        };

        // initialise vectors
        a_continuum.resize(a_energy.size() - 1);
        std::memset(&a_continuum[0], 0, sizeof(Real) * a_continuum.size());

        // find the element at the input temperature bin
        if (const auto el_iterator = m_aped_data[a_temp_bin].m_elements.find(a_atomic_number); el_iterator != m_aped_data[a_temp_bin].m_elements.end())
        {
            const Element &atom = el_iterator->second;

            // if ionization state (rmJ) == 0 then add up all ions, else select ionization state according to input
            if (const auto &ion_iterator = atom.m_ions.find(a_rmJ); ion_iterator != atom.m_ions.end())
            {
                const Ion &ion = ion_iterator->second;

                if (a_psd_cont_emission)
                {
                    std::vector<Real> pseudo_cont_energy(ion.m_pseudo_cont_energy);
                    if (a_doppler_shift != 0.e0)
                    {
                        for (auto &e : pseudo_cont_energy)
                            e *= (one + a_doppler_shift);
                    }
                    add_cont_emission_to_spectrum(ion.m_pseudo_cont, pseudo_cont_energy);

                    if constexpr (!std::is_same_v<line_shape_t<Profile>,Delta>)
                    {
                        Timer_t<3> t("Aped::ion_continuum_emission:convolution");

                        if (Spacing::is_uniform(m_energy_spacing) &&  broadening_spacing_v<Profile> == m_energy_spacing)
                        {
                            const Real x_fwhm = Profile::fwhm(m_temperatures[a_temp_bin], AMU_cgs * atom.m_atomic_mass, one);
                            const Real mesh_size = m_energy_spacing == Spacing::log_uniform ? a_energy[1] / a_energy[0] : a_energy[1] - a_energy[0];
                            Convolution::convolve(a_continuum, build_kernel<Profile>(mesh_size / x_fwhm, a_kernel_tolerance));
                        }
                        else
                        {
                            // need to calculate fwhm for eacg energy to properly set eta in case of voigt profile
                            const Real x_fwhm = Profile::fwhm(m_temperatures[a_temp_bin], AMU_cgs * atom.m_atomic_mass, one);
                            std::vector<Real> fwhm(a_continuum.size());
                            for (size_t i = 0; i < fwhm.size(); ++i)
                                fwhm[i] =  x_fwhm* half * (a_energy[i] + a_energy[i + 1]);
                            Convolution::convolve<line_shape_t<Profile>>(a_continuum, a_energy, fwhm, a_kernel_tolerance);
                        }
                    }
                }

                if (a_cont_emission)
                {
                    std::vector<Real> cont_energy(ion.m_cont_energy);
                    if (a_doppler_shift != 0.e0)
                    {
                        for (auto &e : cont_energy)
                            e *= (one + a_doppler_shift);
                    }
                    add_cont_emission_to_spectrum(ion.m_continuum, cont_energy);
                }
            }
            else
            {
                std::cerr << " strangely ion " << a_rmJ << " was not found " << std::endl;
            }
        }
        else
        {
            std::cerr << " strangely element " << a_atomic_number << " was not found " << std::endl;
        }
    }

} // namespace fm::aped

#endif // APED_H
#endif // USE_APED
