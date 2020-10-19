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

#include <cassert>
#include "Aped.h"
#include "FitsUtil.h"

namespace fm::aped
{
    using namespace fm::fits_util;

    // full construtor
    Aped::Aped(const std::string a_aped_path, const std::string a_version, const int a_verbosity)
        : m_verbosity{a_verbosity}, m_energy_spacing{Spacing::undetermined}
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
            m_dlog_temp = ltemp_step;

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
            std::vector<Real> lkT;
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
                std::vector<Real> ckT;
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
                std::vector<Real> lambda, lambda_err, epsilon, epsilon_err;
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

                std::vector<std::vector<Real>> enrg_cont(num_rows), continuum(num_rows), enrg_pseudo(num_rows), pseudo(num_rows);
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
                temp_record.m_temperature = temperature;

                // parse emission line data
                for (int l = 0; l < num_lines; ++l)
                {
                    // make sure record contains relevant element and ion
                    temp_record.check_element(element[l]);
                    temp_record.m_elements[element[l]].check_ion(ion[l]);
                    //if (element[l]==3) std::cout << " A(3) " << temp_record.elements[element[l]].m_atomic_number << std::endl;

                    temp_record.m_elements[element[l]].m_ions[ion[l]].m_line_energy.emplace_back(keVToAngstrom / lambda[l]);
                    temp_record.m_elements[element[l]].m_ions[ion[l]].m_line_emissivity.emplace_back(epsilon[l]);
                    temp_record.m_elements[element[l]].m_ions[ion[l]].m_lower_level.emplace_back(lower_lev[l]);
                    temp_record.m_elements[element[l]].m_ions[ion[l]].m_upper_level.emplace_back(upper_lev[l]);
                    // temp_record.m_elements[element[l]].m_ions[ion[l]].line_energy_err    .emplace_back( keVToAngstrom/lambda_err[l] );
                    // temp_record.m_elements[element[l]].m_ions[ion[l]].line_emissivity_err.emplace_back( epsilon_err[l] );
                }

                // parse continuum emission
                for (int r = 0; r < num_rows; ++r)
                {
                    // make sure record contains relevant element and ion
                    temp_record.check_element(Z[r]);
                    temp_record.m_elements[Z[r]].check_ion(rmJ[r]);

                    temp_record.m_elements[Z[r]].m_ions[rmJ[r]].m_cont_energy = enrg_cont[r];
                    temp_record.m_elements[Z[r]].m_ions[rmJ[r]].m_continuum = continuum[r];
                    temp_record.m_elements[Z[r]].m_ions[rmJ[r]].m_pseudo_cont_energy = enrg_pseudo[r];
                    temp_record.m_elements[Z[r]].m_ions[rmJ[r]].m_pseudo_cont = pseudo[r];
                    // std::cout << " e_c size= " << enrg_cont.size() << ", n rows= " << num_cont[r] <<
                    // " tr-e_c size= " << temp_record.m_elements[Z[r]].m_ions[rmJ[r]].m_cont_energy.size()  << std::endl;
                }
                assert(temp_record.m_elements.size() == NUM_APED_ATOMS);

                // finally add data to table
                m_aped_data.emplace_back(temp_record);

            } // loop over hdus
        }     // timer

        if (m_verbosity > 0)
        {
            std::cout << "  " << std::endl;
            std::cout << " Aped Code was built successfully using fits files in " << a_aped_path << std::endl;
            std::cout << "  " << std::endl;
        }
    }
} // namespace fm::aped
#endif // USE_APED
