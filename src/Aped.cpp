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

#define MAXSTRLEN 1024
#define READONLY 0 // options when opening a file

#include <algorithm>
#include <cassert>
#include "Aped.h"
#include "FitsUtil.h"

namespace fm::aped
{
    using namespace fm::fits_util;

    // full construtor
    Aped::Aped(const std::string a_aped_path, const std::string a_version, const int a_verbosity)
        : m_verbosity(a_verbosity)
    {
        Timer_t<> t("Aped::define");

        const std::string apec_line_file = a_aped_path + "/apec_v" + a_version + "_line.fits";
        const std::string apec_coco_file = a_aped_path + "/apec_v" + a_version + "_coco.fits";

        // open line file
        char comment[MAXSTRLEN];
        char value[MAXSTRLEN];
        int status = 0;

        // file pointer
        fitsfile *f_line_ptr = NULL;
        fitsfile *f_coco_ptr = NULL;

        if (fits_open_file(&f_line_ptr, apec_line_file.c_str(), READONLY, &status))
        {
            print_fits_error(status);
        }
        if (fits_open_file(&f_coco_ptr, apec_coco_file.c_str(), READONLY, &status))
        {
            print_fits_error(status);
        }

        // number of hdus
        int num_hdus = -1, num_chdus = -1;
        {
            Timer_t<> t("Aped::define:fits_read_in_params");

            fits_get_num_hdus(f_line_ptr, &num_hdus, &status);
            fits_get_num_hdus(f_coco_ptr, &num_chdus, &status);
            assert(num_hdus == num_chdus);

            // get preliminary parameters
            int num_ltemps = -1, num_ctemps;
            fits_read_key(f_line_ptr, TINT, "INUM_TEMP", &num_ltemps, comment, &status);
            fits_read_key(f_coco_ptr, TINT, "INUM_TEMP", &num_ctemps, comment, &status);
            assert(num_ltemps == num_ctemps);
            assert(status == 0);
            double ltemp_start, ltemp_step, ldens_start, ldens_step;
            fits_read_key(f_line_ptr, TDOUBLE, "DTEMP_START", &ltemp_start, comment, &status);
            fits_read_key(f_line_ptr, TDOUBLE, "DTEMP_STEP", &ltemp_step, comment, &status);
            fits_read_key(f_line_ptr, TDOUBLE, "DDENSITY_START", &ldens_start, comment, &status);
            fits_read_key(f_line_ptr, TDOUBLE, "DDENSITY_STEP", &ldens_step, comment, &status);
            double ctemp_start, ctemp_step, cdens_start, cdens_step;
            fits_read_key(f_coco_ptr, TDOUBLE, "DTEMP_START", &ctemp_start, comment, &status);
            fits_read_key(f_coco_ptr, TDOUBLE, "DTEMP_STEP", &ctemp_step, comment, &status);
            fits_read_key(f_coco_ptr, TDOUBLE, "DDENSITY_START", &cdens_start, comment, &status);
            fits_read_key(f_line_ptr, TDOUBLE, "DDENSITY_STEP", &cdens_step, comment, &status);
            assert(ltemp_start == ctemp_start && status == 0);
            assert(ltemp_step == ctemp_step && status == 0);
            assert(ldens_start == cdens_start && status == 0);
            assert(ldens_step == cdens_step && status == 0);

            // set log temperature step
            m_dlog_temp = (double)ltemp_step;

            char PARAMETERS[] = "PARAMETERS";
            // now move to parameters HDU in line & coco file
            fits_movnam_hdu(f_line_ptr, ANY_HDU, PARAMETERS, 0, &status);
            fits_movnam_hdu(f_coco_ptr, ANY_HDU, PARAMETERS, 0, &status);

            // read parameters from input fits files:
            int num_dims = 0;
            fits_read_key(f_line_ptr, TINT, "NAXIS", &num_dims, comment, &status);
            assert(num_dims == 2);
            int num_cols = -1;
            fits_read_key(f_line_ptr, TINT, "NAXIS1", &num_cols, comment, &status);
            assert(status == 0);
            int num_rows = -1;
            fits_read_key(f_line_ptr, TINT, "NAXIS2", &num_rows, comment, &status);
            assert(status == 0);

            // read parameters data: store temp(K) and dens read from HDUs below
            std::vector<float> lkT;
            read_fits_column(f_line_ptr, lkT, "kT", num_rows);
            //read_fits_column(f_line_ptr,m_density ,"EDensity",num_rows);
            read_fits_column(f_line_ptr, m_num_elements_line, "NElement", num_rows);
            read_fits_column(f_line_ptr, m_num_lines, "Nline", num_rows);
            //
            read_fits_column(f_coco_ptr, m_num_elements_coco, "NElement", num_rows);
            read_fits_column(f_coco_ptr, m_density, "EDensity", num_rows);
            read_fits_column(f_coco_ptr, m_num_continuum, "NCont", num_rows);
            read_fits_column(f_coco_ptr, m_num_pseudo, "NPseudo", num_rows);

            { // sanity check: same temperatures for line and cont emission
                std::vector<float> ckT;
                read_fits_column(f_coco_ptr, ckT, "kT", num_rows);
                for (size_t i = 0; i < ckT.size(); ++i)
                {
                    assert(ckT[i] == lkT[i]);
                }
            }
        }

        {
            Timer_t<> t("Aped::define:fits_read_in_data");

            // loop over fits hdus starting from second
            for (int hdu = 3; hdu <= num_hdus; ++hdu)
            {
                // advance hdu of line file
                int hdu_type;
                fits_movabs_hdu(f_line_ptr, hdu, &hdu_type, &status);
                print_fits_error(status);

                fits_read_keyword(f_line_ptr, "EXTNAME", value, comment, &status);
                // read in numbers from header
                double temperature, density, time;
                int num_lines;
                fits_read_key(f_line_ptr, TDOUBLE, "HIERARCH TEMPERATURE", &temperature, comment, &status);
                fits_read_key(f_line_ptr, TDOUBLE, "DENSITY", &density, comment, &status);
                fits_read_key(f_line_ptr, TDOUBLE, "TIME", &time, comment, &status);
                fits_read_key(f_line_ptr, TINT, "N_LINES", &num_lines, comment, &status);
                // sanity checks
                // std::cout << " num-L = " << num_lines << ", num L_hdu= " << m_num_lines.at(hdu-3) << endl;
                // assert(num_lines==m_num_lines.at(hdu-3));
                // store temperature and density for this HDU
                m_temperatures.emplace_back(temperature);
                m_density.emplace_back(density);

                std::vector<int> element, ion, upper_lev, lower_lev;
                std::vector<float> lambda, lambda_err, epsilon, epsilon_err;
                read_fits_column(f_line_ptr, lambda, "Lambda", num_lines);
                read_fits_column(f_line_ptr, lambda_err, "Lambda_Err", num_lines);
                read_fits_column(f_line_ptr, epsilon, "Epsilon", num_lines);
                read_fits_column(f_line_ptr, epsilon_err, "Epsilon_Err", num_lines);
                read_fits_column(f_line_ptr, lower_lev, "LowerLev", num_lines);
                read_fits_column(f_line_ptr, upper_lev, "UpperLev", num_lines);
                read_fits_column(f_line_ptr, element, "Element", num_lines);
                read_fits_column(f_line_ptr, ion, "Ion", num_lines);

                // now move forward hdu of coco file
                fits_movabs_hdu(f_coco_ptr, hdu, &hdu_type, &status);
                print_fits_error(status);

                // read from header
                fits_read_keyword(f_coco_ptr, "EXTNAME", value, comment, &status);
                // read in numbers from header
                double ctemperature; //,density,time;
                int num_rows;
                fits_read_key(f_coco_ptr, TDOUBLE, "HIERARCH TEMPERATURE", &ctemperature, comment, &status);
                // fits_read_key(f_coco_ptr,TDOUBLE, "DENSITY"             , &density    , comment, &status);
                // fits_read_key(f_coco_ptr,TDOUBLE, "TIME"                , &time       , comment, &status);
                fits_read_key(f_coco_ptr, TINT, "NAXIS2", &num_rows, comment, &status);

                //
                std::vector<int> Z, rmJ, num_cont, num_pseudo;
                read_fits_column(f_coco_ptr, Z, "Z", num_rows);
                read_fits_column(f_coco_ptr, rmJ, "rmJ", num_rows);
                read_fits_column(f_coco_ptr, num_cont, "N_Cont", num_rows);
                read_fits_column(f_coco_ptr, num_pseudo, "N_Pseudo", num_rows);
                // std::cout << " hdu " << hdu << " num-C = " << num_cont.size() << ", num C_hdu= " << m_num_continuum.at(hdu-3) << endl;
                // sanity check
                //assert(num_cont.size()==(size_t)m_num_continuum.at(hdu-3));
                //assert(num_pseudo.size()==(size_t)m_num_pseudo.at(hdu-3));
                assert(ctemperature == m_temperatures.at(hdu - 3));

                //
                std::vector<std::vector<float>> enrg_cont(num_rows), continuum(num_rows), enrg_pseudo(num_rows), pseudo(num_rows);
                for (int r = 0; r < num_rows; ++r)
                {
                    enrg_cont[r].resize(num_cont[r]);
                    continuum[r].resize(num_cont[r]);
                    enrg_pseudo[r].resize(num_pseudo[r]);
                    pseudo[r].resize(num_pseudo[r]);

                    read_fits_column(f_coco_ptr, enrg_cont[r], "E_Cont", num_cont[r], r + 1);
                    read_fits_column(f_coco_ptr, continuum[r], "Continuum", num_cont[r], r + 1);
                    read_fits_column(f_coco_ptr, enrg_pseudo[r], "E_Pseudo", num_pseudo[r], r + 1);
                    read_fits_column(f_coco_ptr, pseudo[r], "Pseudo", num_pseudo[r], r + 1);
                }

                // loop over elements
                TemperatureRecord temp_record;
                temp_record.temperature = temperature;

                // parse emission line data
                for (int l = 0; l < num_lines; ++l)
                {
                    // make sure record contains relevant element and ion
                    temp_record.check_element(element[l]);
                    temp_record.elements[element[l]].check_ion(ion[l]);
                    //if (element[l]==3) std::cout << " A(3) " << temp_record.elements[element[l]].atomic_number << std::endl;

                    temp_record.elements[element[l]].ions[ion[l]].line_energy.emplace_back(keVToAngstrom / lambda[l]);
                    temp_record.elements[element[l]].ions[ion[l]].line_emissivity.emplace_back(epsilon[l]);
                    temp_record.elements[element[l]].ions[ion[l]].lower_level.emplace_back(lower_lev[l]);
                    temp_record.elements[element[l]].ions[ion[l]].upper_level.emplace_back(upper_lev[l]);
                    // temp_record.elements[element[l]].ions[ion[l]].line_energy_err    .emplace_back( keVToAngstrom/lambda_err[l] );
                    // temp_record.elements[element[l]].ions[ion[l]].line_emissivity_err.emplace_back( epsilon_err[l] );
                }

                // parse continuum emission
                for (int r = 0; r < num_rows; ++r)
                {
                    // make sure record contains relevant element and ion
                    temp_record.check_element(Z[r]);
                    temp_record.elements[Z[r]].check_ion(rmJ[r]);

                    temp_record.elements[Z[r]].ions[rmJ[r]].cont_energy = enrg_cont[r];
                    temp_record.elements[Z[r]].ions[rmJ[r]].continuum = continuum[r];
                    temp_record.elements[Z[r]].ions[rmJ[r]].pseudo_cont_energy = enrg_pseudo[r];
                    temp_record.elements[Z[r]].ions[rmJ[r]].pseudo_cont = pseudo[r];
                    // std::cout << " e_c size= " << enrg_cont.size() << ", n rows= " << num_cont[r] <<
                    // " tr-e_c size= " << temp_record.elements[Z[r]].ions[rmJ[r]].cont_energy.size()  << std::endl;
                }
                assert(temp_record.elements.size() == NUM_APED_ATOMS);

                // finally add data to table
                m_aped_data.emplace_back(temp_record);

            } // loop over hdus
        } // timer

        if (m_verbosity > 0)
        {
            std::cout << "  " << std::endl;
            std::cout << " Aped Code was built successfully using fits files in " << a_aped_path << std::endl;
            std::cout << "  " << std::endl;
        }
    }

    /// continuum spectrum including pseudo continuum
    void Aped::emission_spectrum(std::vector<Real> &a_spectrum,
                                 const std::vector<Real> &a_energy,
                                 const std::vector<ElementAbundance> &a_atom_abundances,
                                 const Real a_temperature,
                                 const Real a_doppler_shift,
                                 const std::string a_line_broadening,
                                 const bool a_line_emission,
                                 const bool a_cont_emission) const
    {
        Timer_t<> t("Aped::emission_spectrum");
        // sanity check, catch miscommunications
        assert(a_line_broadening == "none" || a_line_broadening == "convolution" || a_line_broadening == "linebyline");

        // initialize spectrum
        a_spectrum.resize(a_energy.size() - 1, 0);

        // check temperature range
        const Real T_min = (Real)*std::min_element(m_temperatures.begin(), m_temperatures.end());
        const Real T_max = (Real)*std::max_element(m_temperatures.begin(), m_temperatures.end());

        // if T is within allowed range
        if (a_temperature >= T_min && a_temperature <= T_max)
        {
            // identify right bin
            const int it_lo = (int)floor(log10(a_temperature / T_min) / m_dlog_temp);
            const int it_hi = std::min((size_t)it_lo + 1, m_temperatures.size() - 1);
            assert(a_temperature >= (Real)m_temperatures[it_lo] && a_temperature <= (Real)m_temperatures[it_hi]);

            // loop over all elements
            for (const auto &A : a_atom_abundances)
            {
                // loop over temp bins
                for (int it = it_lo; it <= it_hi; ++it)
                {
                    // temperature interpolation coefficient: interp. in log space, following log spacing of temperature tabulation
                    // const Real f=one - std::abs(log10(a_temperature/(Real)m_temperatures[it]))/m_dlog_temp;
                    // Linear interpolation adopted by original Aped code
                    const Real f = one - std::abs(a_temperature - (Real)m_temperatures[it]) / (m_temperatures[it_lo] * (pow(ten, m_dlog_temp) - one));

                    // add ion emission to spectrum taking into accout io abundance and ionization fraction
                    auto add_emission_to_spectrum = [&j = a_spectrum, x = f * A.abundance](const std::vector<Real> &i) {
                        for (size_t k = 0; k < j.size(); ++k)
                            j[k] += x * i[k];
                    };

                    if (a_line_emission)
                    {
                        Timer_t<2> t("Aped::emission_spectrum:line_emission");

                        std::vector<Real> line_emission(a_spectrum.size(), 0);
                        ion_line_emission(line_emission, a_energy, A.atomic_number, 0, it, a_doppler_shift, a_line_broadening);

                        // add ion line emission
                        add_emission_to_spectrum(line_emission);
                    }
                    //
                    if (a_cont_emission || a_line_emission)
                    {
                        Timer_t<2> t("Aped::emission_spectrum:cont_emission");

                        // continuum
                        std::vector<Real> emission(a_spectrum.size(), 0);
                        ion_continuum_emission(emission,
                                               a_energy,
                                               A.atomic_number,
                                               0, it,
                                               a_doppler_shift,
                                               a_cont_emission,
                                               a_line_emission);

                        // add psd-cont ion emission
                        add_emission_to_spectrum(emission);
                    }
                }
            }
        }
    }

    //
    void Aped::ion_line_emission(std::vector<Real> &a_spectrum,
                                 const std::vector<Real> &a_energy,
                                 const int a_atomic_number,
                                 const int a_rmJ,
                                 const int a_temp_idx,
                                 const Real a_doppler_shift,
                                 const std::string a_line_broadening) const
    {
        // resize spectrum vector
        a_spectrum.resize(a_energy.size() - 1, 0);

        // count num lines
        size_t num_elines = 0;

        // find the element at the input temperature bin
        if (const auto ei = m_aped_data[a_temp_idx].elements.find(a_atomic_number);
            ei != m_aped_data[a_temp_idx].elements.end())
        {
            const Element &atom = ei->second;

            // thermal velocity (kT/m)^1/2
            const Real x_thermal = sqrt_two * std::sqrt(kB_cgs * m_temperatures[a_temp_idx] / (AMU_cgs * atom.atomic_mass)) / c_cgs;

            // if ionization state (rmJ) == 0 then add up all ions, else select ionization state according to input
            const auto beg = a_rmJ == 0 ? atom.ions.begin() : atom.ions.find(a_rmJ);
            const auto end = a_rmJ == 0 ? atom.ions.end() : (beg != atom.ions.end() ? std::next(beg) : beg);
            for (auto ii=beg; ii != end; ++ii)
            {
                // set ion pointer
                const Ion &ion = ii->second;
                assert(ion.line_emissivity.size() == ion.line_energy.size());

                // loop thorugh emission lines; ie is the energy bin of line
                int ie = 0;
                for (size_t il = 0; il < ion.line_emissivity.size(); ++il)
                {
                    // doppler shifted photon energy
                    const Real e_line = ion.line_energy[il] * (one + a_doppler_shift);

                    if (e_line >= a_energy.front() && e_line < a_energy.back())
                    {
                        ++num_elines;

                        // find e-bin
                        while (a_energy[ie] < e_line)
                            ++ie;
                        ie = (ie > 0 ? ie - 1 : ie);

                        // add thermal broadening
                        if (a_line_broadening == "linebyline")
                        {
                            Timer_t<4> t("Aped::ion_line_emission:linebyline_broadening");
                            // use gaussian kernel: thermal broadening c_th/c * e_line
                            const Real sqr2_sigma = x_thermal * e_line;

                            GaussianKernel kernel(e_line/sqr2_sigma, a_energy[ie]/sqr2_sigma, a_energy[ie+1]/sqr2_sigma);
                            kernel.convolve(a_spectrum, ion.line_emissivity[il], ie);
                        }
                        else
                        {
                            // add to spectrum
                            a_spectrum[ie] += ion.line_emissivity[il]; 
                        }
                    }
                }
                // if a_rmJ==0 we are done
                if (a_rmJ != 0)
                    break;
            }
            // thermal broadening
            if (a_line_broadening == "convolution")
            {
                Timer_t<3> t("Aped::emission_spectrum:line_convolution");
                const Real Elo = one / x_thermal;
                const Real Ehi = Elo * a_energy[1] / a_energy[0];
                const Real Eline = half * (Elo + Ehi);

                GaussianKernel kernel(Eline, Elo, Ehi);
                kernel.convolve(a_spectrum);
            }
        }
        else
        {
            std::cerr << " strangely element " << a_atomic_number << " was not found " << std::endl;
        }

        if (m_verbosity > 0)
        {
            std::cout << " added " << num_elines << " to spectrum " << std::endl;
        }
    }

    // continuum
    void Aped::ion_continuum_emission(std::vector<Real> &a_spectrum,
                                      const std::vector<Real> &a_energy,
                                      const int a_atomic_number,
                                      const int a_rmJ,
                                      const int a_temp_idx,
                                      const Real a_doppler_shift,
                                      const bool a_cont_emission,
                                      const bool a_pseudo_cont_emission) const
    {
        // resize spectrum vector
        a_spectrum.resize(a_energy.size() - 1, 0);

        auto add_cont_emission_to_spectrum = [&j = a_spectrum, &e = a_energy](const std::vector<float> &js,
                                                                              const std::vector<float> &es) {
            Timer_t<4> t("Aped::ion_continuum_emission::add_to_spectrum");

            if (es.back() < e.front() && es.front() > e.back())
                return;

            size_t k = 0;
            float ef=e[0], jf = js[0]; // foot point values
            for (size_t i = 0; i < j.size(); ++i)
            {
                if (es.back() > e[i])
                {
                    while (es[k] < e[i])
                        ++k;

                    // reset foot emissivity if need be
                    if (i==0 && k>0)
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
                        const float jh = jf + (js[k] - jf) / (es[k] - ef) * (e[i + 1] - ef);
                        j[i] += half * (jf + jh) * (e[i + 1] - ef);
                        jf = jh;
                        ef = e[i + 1];
                    }
                }
            }
        };

        // find the element at the input temperature bin
        if (const auto ei = m_aped_data[a_temp_idx].elements.find(a_atomic_number); ei != m_aped_data[a_temp_idx].elements.end())
        {
            const Element &atom = ei->second;

            // if ionization state (rmJ) == 0 then add up all ions, else select ionization state according to input
            if (const auto &it = atom.ions.find(a_rmJ); it != atom.ions.end())
            {
                const Ion &ion = it->second;

                if (a_cont_emission)
                {
                    std::vector<float> cont_energy(ion.cont_energy);
                    if (a_doppler_shift != 0.e0)
                    {
                        for (auto &e : cont_energy)
                            e *= (one + a_doppler_shift);
                    }
                    add_cont_emission_to_spectrum(ion.continuum, cont_energy);
                }

                if (a_pseudo_cont_emission)
                {
                    std::vector<float> pseudo_cont_energy(ion.pseudo_cont_energy);
                    if (a_doppler_shift != 0.e0)
                    {
                        for (auto &e : pseudo_cont_energy)
                            e *= (one + a_doppler_shift);
                    }
                    add_cont_emission_to_spectrum(ion.pseudo_cont, pseudo_cont_energy);
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
#endif // USE_APED
