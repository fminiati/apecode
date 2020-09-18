#ifdef USE_APED

#define MAXSTRLEN 1024
#define READONLY 0 // options when opening a file

#include <algorithm>
#include <cassert>
#include "Aped.h"
#include "Util.h"
#include "FitsUtil.h"

using namespace fm::fits_util;
using namespace fm::aped;

// full construtor
void Aped::define(const std::string a_aped_path, const std::string a_version)
{
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
    fits_get_num_hdus(f_line_ptr, &num_hdus, &status);
    fits_get_num_hdus(f_coco_ptr, &num_chdus, &status);
    assert(num_hdus == num_chdus);

    // get preliminary parameters
    //
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
        m_temperatures.push_back(temperature);
        m_density.push_back(density);

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

            temp_record.elements[element[l]].ions[ion[l]].line_energy.push_back(keVToAngstrom / lambda[l]);
            temp_record.elements[element[l]].ions[ion[l]].line_emissivity.push_back(epsilon[l]);
            temp_record.elements[element[l]].ions[ion[l]].lower_level.push_back(lower_lev[l]);
            temp_record.elements[element[l]].ions[ion[l]].upper_level.push_back(upper_lev[l]);
            // temp_record.elements[element[l]].ions[ion[l]].line_energy_err    .push_back( keVToAngstrom/lambda_err[l] );
            // temp_record.elements[element[l]].ions[ion[l]].line_emissivity_err.push_back( epsilon_err[l] );
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
        assert(temp_record.elements.size() == NUM_APEC_ATOMS);

        // finally add data to table
        m_aped_data.push_back(temp_record);

    } // loop over hdus

    std::cout << "  " << std::endl;
    std::cout << " Aped Code was built successfully using fits files in " << a_aped_path << std::endl;
    std::cout << "  " << std::endl;
}

/// continuum spectrum including pseudo continuum
void Aped::emission_spectrum(std::vector<Real> &a_spectrum,
                             const std::vector<Real> &a_energy,
                             const std::list<ElementAbundance> &a_atom_abundances,
                             const Real a_temperature,
                             const Real a_doppler_shift,
                             const std::string a_line_broadening,
                             const bool a_line_emission,
                             const bool a_cont_emission) const
{
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

                if (a_line_emission)
                {
                    // lines
                    std::vector<Real> spectrum(a_spectrum.size(), 0);
                    ion_line_emission(spectrum, a_energy, A.atomic_number, 0, it, a_doppler_shift, a_line_broadening);

                    // thermal broadening
                    if (a_line_broadening == "convolution")
                    {
                        const auto &el = m_aped_data[it].elements.find(A.atomic_number);
                        //             const Real atomic_mass= AMU_g*el->second.atomic_mass; //m_aped_data[it].elements[A->atomic_number].atomic_mass;
                        //             const Real sigma=sqrt(k_B*a_temperature/atomic_mass) / C_cm_s;
                        //             const Real ksize=(a_energy[1]-a_energy[0])/a_energy[0] / (sqrt(two)*sigma);
                        const Real sqrt2_sigma = sqrt(2 * k_B * a_temperature / (AMU_g * el->second.atomic_mass)) / C_cm_s;
                        const Real Elo = one / (sqrt2_sigma);
                        const Real Ehi = a_energy[1] / a_energy[0] / (sqrt2_sigma);
                        const Real Eline = half * (Elo + Ehi);

                        GaussianKernel *kernel = new GaussianKernel(Eline, Elo, Ehi);
                        simple_convolution(spectrum, kernel);
                        delete kernel;
                    }

                    // include abundance factor
                    add(a_spectrum, spectrum, f * A.abundance);
                }
                //
                if (a_cont_emission)
                {
                    // continuum
                    std::vector<Real> spectrum(a_spectrum.size(), 0);
                    ion_continuum_emission(spectrum,
                                           a_energy,
                                           A.atomic_number,
                                           0, it,
                                           a_doppler_shift);

                    // include abundance factor
                    add(a_spectrum, spectrum, f * A.abundance);
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
    long num_elines = 0;

    // find the element at the input temperature bin
    if (const auto ei = m_aped_data[a_temp_idx].elements.find(a_atomic_number);
        ei != m_aped_data[a_temp_idx].elements.end())
    {
        const Element &atom = ei->second;

        // if ionization state (rmJ) == 0 then add up all ions, else select ionization state according to input
        auto ii = a_rmJ == 0 ? atom.ions.begin() : atom.ions.find(a_rmJ);
        for (; ii != atom.ions.end(); ++ii)
        {
            // set ion pointer
            const Ion &ion = ii->second;
            assert(ion.line_emissivity.size() == ion.line_energy.size());

            // loop through energy of emission line
            //const auto j=ion.line_emissivity.begin();
            for (auto hn = ion.line_energy.begin(), j = ion.line_emissivity.begin(); hn != ion.line_energy.end(); ++hn, ++j)
            {
                // photon energy shifted by doppler effect
                const Real e_line = (Real)(*hn) * (one + a_doppler_shift);

                if (e_line >= a_energy.front() && e_line < a_energy.back())
                {
                    ++num_elines;

                    // find e-bin
                    int ie = 0;
                    while (a_energy[ie] < e_line)
                        ++ie;
                    ie = (ie > 0 ? ie - 1 : ie);
                    assert(ie < (int)a_spectrum.size());

                    // add thermal broadening
                    if (a_line_broadening == "linebyline")
                    {
                        // use gaussian kernel: thermal broadening (kT/mc2)^1/2 * e_line
                        const Real sqrt2_sigma = sqrt(2 * k_B * m_temperatures[a_temp_idx] / (AMU_g * atom.atomic_mass)) / C_cm_s * e_line;
                        const Real Elo = a_energy[ie] / (sqrt2_sigma);
                        const Real Ehi = a_energy[ie + 1] / (sqrt2_sigma);
                        const Real Eline = e_line / (sqrt2_sigma);
                        // energy
                        GaussianKernel kernel(Eline, Elo, Ehi);
                        const int kwing = kernel.wing_size();
                        // line profile
                        for (int k = -kwing; k <= kwing; ++k)
                        {
                            if (ie + k >= 0 && ie + k < (int)a_spectrum.size())
                                a_spectrum[ie + k] += (Real)(*j) * kernel.val(k);
                        }
                    }
                    else
                    {
                        // add to spectrum
                        a_spectrum.at(ie) += (Real)(*j);
                    }
                }
            }
            // if a_rmJ==0 we are done
            if (a_rmJ != 0)
                break;
        }
    }
    else
    {
        std::cerr << " strangely element " << a_atomic_number << " was not found " << std::endl;
    }
    //
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
                                  const Real a_doppler_shift) const
{
    // resize spectrum vector
    a_spectrum.resize(a_energy.size() - 1, 0);

    // find the element at the input temperature bin
    if (const auto ei = m_aped_data[a_temp_idx].elements.find(a_atomic_number); ei != m_aped_data[a_temp_idx].elements.end())
    {
        // element
        const Element &atom = ei->second;

        // if ionization state (rmJ) == 0 then add up all ions, else select ionization state according to input
        if (const auto &it = atom.ions.find(a_rmJ); it != atom.ions.end())
        {
            // ion
            const Ion &ion = it->second;
            // shift spectra
            std::vector<float> cnt_enrg_shifted(ion.cont_energy), psd_enrg_shifted(ion.pseudo_cont_energy);
            if (a_doppler_shift != 0.e0)
            {
                for (auto &vi : cnt_enrg_shifted)
                    vi *= (one + a_doppler_shift);
            }
            for (auto &vi : psd_enrg_shifted)
                vi *= (one + a_doppler_shift);

            linear_integrate(a_spectrum, a_energy, ion.continuum, cnt_enrg_shifted);
            linear_integrate(a_spectrum, a_energy, ion.pseudo_cont, psd_enrg_shifted);
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

#endif // USE_APED
