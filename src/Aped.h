#ifdef USE_APED

#ifndef APED_H
#define APED_H

#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "Timer.h"
#include "SutherlandDopita.h"

namespace fm::aped
{
    using namespace fm::profiling;
    using Real = double;

    constexpr double zero = 0.0e0;
    constexpr double half = 0.5e0;
    constexpr double one = 1.0e0;
    constexpr double two = 2.0e0;
    constexpr double ten = 1.0e1;

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

    // sum - add second vector to first
    template <class T>
    void add(std::vector<T> &va, const std::vector<T> &vb, const T f)
    {
        assert(va.size() == vb.size());

        unsigned i = 0;
        unsigned n = va.size();
        T *a = &va[0];
        const T *b = &vb[0];
        while (i++ < n)
            *a++ += f * *b++;
    }

    // return interpolated function (f(x)) value at input position xp
    template <class T>
    T linear_interp(const T xp, const std::vector<T> &x, const std::vector<T> &f)
    {
        //
        if (xp < x.front() || xp >= x.back())
            return T(0);

        size_t i = 0;
        while (xp > x[i])
            ++i;

        return ((f[i] - f[i - 1]) / (x[i] - x[i - 1]) * (xp - x[i]) + f[i]);
    }

    // return interpolated function (f(x)) value at input position xp
    template <class T, class R>
    void linear_interp(std::vector<T> &fp, const std::vector<T> &xp,
                       const std::vector<R> &f, const std::vector<R> &x)
    {
        //
        std::vector<T> xh(f.size());
        {
            const R *xl = &x[0];
            const R *xr = &x[1];
            Real *x = &xh[0];
            unsigned i = xh.size();
            while (i-- > 0)
            {
                *x++ = half * (*xl++ + *xr++);
            }
        }

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

    // return interpolated function (f(x)) value at input position xp
    template <class T, class R>
    void linear_integrate(std::vector<T> &fp, const std::vector<T> &xp,
                          const std::vector<R> &f, const std::vector<R> &x)
    {
        //
        size_t ip = 0;
        for (size_t i = 0; i < fp.size(); ++i)
        {

            const T hp = half * (xp[i] + xp[i + 1]);
            if (hp >= x.front() && hp <= x.back())
            {
                while (hp > x[ip])
                    ++ip;
                fp[i] += ((f[ip] - f[ip - 1]) / (x[ip] - x[ip - 1]) * (hp - x[ip]) + f[ip]) * (xp[i + 1] - xp[i]);
            }
        }
    }

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
        // simple constructor
        ElementAbundance(const unsigned a_n, const Real a_a)
            : atomic_number(a_n), abundance(a_a)
        {
        }
        // null
        ElementAbundance(const ElementAbundance &ea)
            : atomic_number(ea.atomic_number), abundance(ea.abundance)
        {
        }
        //
        ~ElementAbundance() {}

        unsigned atomic_number;
        Real abundance;
    };

    // utility to compute abundances to be used in aped API
    struct AbundanceUtil
    {
        // constructor
        AbundanceUtil() {}
        // destructor
        ~AbundanceUtil() {}

        // compute abundances relative to AndersGrevesse
        void relative_abundances(std::list<ElementAbundance> &a_relative_abundances,
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
        // full constructor
        Ion(const Ion &a_ion)
        {
            *this = a_ion;
        }
        // destructor
        ~Ion() {}
        // partial constructor (set ion species)
        Ion(const unsigned a_ion)
            : ion(a_ion)
        {
        }

        // copy constructor
        Ion &operator=(const Ion &a_rhs)
        {
            this->ion = a_rhs.ion;
            this->cont_energy = a_rhs.cont_energy;
            this->continuum = a_rhs.continuum;
            this->continuum_err = a_rhs.continuum_err;
            this->pseudo_cont_energy = a_rhs.pseudo_cont_energy;
            this->pseudo_cont = a_rhs.pseudo_cont;
            this->pseudo_cont_err = a_rhs.pseudo_cont_err;
            this->line_energy = a_rhs.line_energy;
            this->line_energy_err = a_rhs.line_energy_err;
            this->line_emissivity = a_rhs.line_emissivity;
            this->line_emissivity_err = a_rhs.line_emissivity_err;
            this->elem_driver = a_rhs.elem_driver;
            this->ion_driver = a_rhs.ion_driver;
            this->lower_level = a_rhs.lower_level;
            this->upper_level = a_rhs.upper_level;
            return *this;
        }

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
        // partial constructors
        Element(const Element &A)
            : name(A.name), atomic_number(A.atomic_number), atomic_mass(A.atomic_mass), ions(A.ions)
        {
        }
        //
        Element(const std::string a_name, const unsigned a_A, const float a_M)
            : name(a_name), atomic_number(a_A), atomic_mass(a_M)
        {
        }

        // copy constructor
        Element &operator=(const Element &a_rhs)
        {
            this->name = a_rhs.name;
            this->atomic_number = a_rhs.atomic_number;
            this->atomic_mass = a_rhs.atomic_mass;
            this->ions.insert(a_rhs.ions.begin(), a_rhs.ions.end());
            return *this;
        }

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
        //
        TemperatureRecord(const TemperatureRecord &a_tr)
            : temperature(a_tr.temperature), elements(a_tr.elements)
        {
        }

        // destructor
        ~TemperatureRecord() {}
        // copy constructor
        TemperatureRecord &operator=(const TemperatureRecord &a_rhs)
        {
            this->temperature = a_rhs.temperature;
            this->elements.insert(a_rhs.elements.begin(), a_rhs.elements.end());
            return *this;
        }

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
        // constructor
        Aped()
            : m_verbosity(0)
        {
        }

        Aped(const std::string a_aped_path, const std::string a_version, const int a_verbosity = 0)
            : m_verbosity(a_verbosity)
        {
            define(a_aped_path, a_version);
        }

        Aped(const Aped &a_aped)
        {
            // temperature data
            this->m_aped_data = a_aped.m_aped_data;
            // log step of temp database
            this->m_dlog_temp = a_aped.m_dlog_temp;
            // verbosity, zero by default
            this->m_verbosity = a_aped.m_verbosity;
            //
            this->m_temperatures = a_aped.m_temperatures;
            this->m_density = a_aped.m_density;
            this->m_num_lines = a_aped.m_num_lines;
            this->m_num_elements_line = a_aped.m_num_elements_line;
            this->m_num_elements_coco = a_aped.m_num_elements_coco;
            this->m_num_continuum = a_aped.m_num_continuum;
            this->m_num_pseudo = a_aped.m_num_pseudo;
        }

        // full construtor
        void define(const std::string a_aped_path, const std::string a_version);

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
            AbundanceUtil ab_ut;
            std::list<ElementAbundance> rel_ab;
            ab_ut.relative_abundances(rel_ab, a_elements, a_metallicity, a_abundances_model);

            // compute spectrum
            this->emission_spectrum(a_spectrum, a_energy, rel_ab, a_temperature,
                                    a_doppler_shift, a_line_broadening,
                                    a_line_emission, a_cont_emission);
        }

        // overloaded version of emission spectrum in ph cm^3 s^-1
        void emission_spectrum(std::vector<Real> &a_spectrum,
                               const std::vector<Real> &a_energy,
                               const std::list<ElementAbundance> &a_atom_abundances,
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
                                    const Real a_doppler_shift) const;

    protected:
        // do simple convolution with gaussian weights
        void simple_convolution(std::vector<Real> &a_spectrum,
                                const LineKernel *a_kernel) const
        {
            //
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
        //
        std::vector<double> m_temperatures;
        std::vector<double> m_density;
        std::vector<int> m_num_lines;
        std::vector<int> m_num_elements_line;
        std::vector<int> m_num_elements_coco;
        std::vector<int> m_num_continuum;
        std::vector<int> m_num_pseudo;
    };

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
             const bool a_cont_emission = true)
            : m_spectrum_size(a_num_energy_bins),
              m_abundances_model(a_abundances_model),
              m_line_broadening(a_line_broadening),
              m_line_emission(a_line_emission),
              m_cont_emission(a_cont_emission)
        {
            assert((a_line_broadening == "convolution" && a_energy_bin_spacing == "log") ||
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

            std::cout << "  " << '\n';
            std::cout << "Apec: Apec built successfully " << '\n';
            std::cout << "  " << '\n';
        }

        // copy constructor
        Apec(const Apec &a_apec)
            : m_j_bbn_el(a_apec.m_j_bbn_el),
              m_j_metals(a_apec.m_j_metals),
              m_energy(a_apec.m_energy),
              m_buf_energy(a_apec.m_buf_energy),
              m_temperature(a_apec.m_temperature),
              m_spectrum_size(a_apec.m_spectrum_size),
              m_d_en(a_apec.m_d_en),
              m_dlog_temp(a_apec.m_dlog_temp),
              m_temp_min(a_apec.m_temp_min),
              m_temp_max(a_apec.m_temp_max),
              m_abundances_model(a_apec.m_abundances_model),
              m_line_broadening(a_apec.m_line_broadening),
              m_line_emission(a_apec.m_line_emission),
              m_cont_emission(a_apec.m_cont_emission)
        {
        }

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
            // timer
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
                const Real ne = m_sd.n_e(a_temperature, a_metallicity);

                // rename
                const Real z = a_metallicity;

                const unsigned buf_spectrum_size = m_buf_energy.size() - 1;
                std::vector<Real> js(buf_spectrum_size);
                {
                    Timer_t<> t("APED::Apec::emission_spectrum:temp_interpolation");

                    unsigned i = 0;
                    const Real *jblo = &m_j_bbn_el[it_lo][i]; // bbn low  T
                    const Real *jbhi = &m_j_bbn_el[it_hi][i]; // bbn high T
                    const Real *jzlo = &m_j_metals[it_lo][i]; // met low  T
                    const Real *jzhi = &m_j_metals[it_hi][i]; // met high T
                    Real *j = &js[i];
                    while (i++ < buf_spectrum_size)
                    {
                        *j++ = ne * ((one - f) * (*jblo++ + z * *jzlo++) + f * (*jbhi++ + z * *jzhi++));
                    }
                }

                // Interpolate from rest frame to lab frame spectrum
                {
                    Timer_t<> t("APED::Apec::emission_spectrum:ph_en_interpolation");
                    std::vector<Real> shifted_energy(m_buf_energy.size());

                    const Real df = one + a_doppler_shift;
                    const Real *be = &m_buf_energy[0];
                    Real *se = &shifted_energy[0];
                    unsigned i = 0;
                    while (i++ < m_buf_energy.size())
                    {
                        *se++ = df * *be++;
                    }
                    //

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
            // timer
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

            // utility to compute abundances re
            AbundanceUtil ab_ut;
            //
            std::list<ElementAbundance> rel_abundances;
            ab_ut.relative_abundances(rel_abundances, a_elements, dummy_metallicity, m_abundances_model);

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
    };

} // namespace fm::aped

#endif // APED_H
#endif // USE_APED
