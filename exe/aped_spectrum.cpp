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
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

// this variable would be usually set at compilation, however this is safe
// as this code wouldn't make sense without such definition
#ifndef USE_APED
#define USE_APED
#endif
#include "Aped.h"
#include "FileParser.h"

int main(int argc, char *argv[])
{
    using namespace fm::aped;

    // get input file
    std::string input_file = "unknown";

    for (int i = 1; i < argc; ++i)
    {
        if (strncmp(argv[i], "-f", 3) == 0)
            input_file = argv[i + 1];
    }
    std::cout << "input file name is " << input_file << '\n';
    if (input_file == "unknown")
    {
        const std::string prog(argv[0]);
        std::cout << "\nrun this prog with: " + prog + " -f input_filename\n\n";
        return 0;
    }

    // read input file
    fm::FileParser input(input_file);

    // read in path and version
    std::string aped_file_path;
    input.get_item(aped_file_path, "aped.aped_file_path");

    std::string aped_version;
    input.get_item(aped_version, "aped.aped_version");

    // define Aped object
    Aped aped(aped_file_path, aped_version);

    // min and max ph energy
    double emin = 0.0;
    double emax = -1.0;
    input.get_item(emin, "aped.min_ph_energy[keV]");
    input.get_item(emax, "aped.max_ph_energy[keV]");

    // spectral resolution
    double de = 99999;
    input.get_item(de, "aped.spectral_energy_resolution");

    // using log scale for ph energy
    bool log_scale = true;
    input.get_item(log_scale, "aped.use_log_energy_scale");

    if (log_scale)
        de = pow(ten, de);

    // plasma temperature in keV
    double temperature = 1.0;
    input.get_item(temperature, "aped.plasma_temperature[keV]");
    // convert from keV to K
    temperature *= keVToKelvin;

    // plasma velocity in cm/sec
    double velocity = 0.0;
    input.get_item(velocity, "aped.plasma_velocity[cm/s]");

    // Doppler shift
    const double doppler_shift = velocity / c_light_cgs;

    // metallicity
    double metallicity = 0.3;
    input.get_item(metallicity, "aped.plasma_metallicity");

    // abundance model
    int int_abundance_model = 0;
    input.get_item(int_abundance_model, "aped.plasma_abundance_model");
    fm::aped::AbundanceModel abundance_model=static_cast<fm::aped::AbundanceModel>(int_abundance_model);

    // line emission
    bool line_emission = false;
    input.get_item(line_emission, "aped.use_line_emission");

    // continuum emission
    bool cont_emission = false;
    input.get_item(cont_emission, "aped.use_cont_emission");

    // output file
    std::string output_file;
    input.get_item(output_file, "aped.output_file_name");

    // default is delta function
    int int_line_shape = 0;
    input.get_item(int_line_shape, "aped.line_shape");
    fm::aped::LineShape line_shape = static_cast<fm::aped::LineShape>(int_line_shape);

    // number of elements
    int num_elements;
    input.get_item(num_elements, "aped.num_elements");

    // empty vector means all elements
    std::vector<unsigned> elements(num_elements);
    input.get_items(elements, "aped.elements");

    double kernel_tolerance;
    input.get_item(kernel_tolerance, "aped.kernel_tolerance");

    // spectral energy and emission
    std::vector<double> ph_energy(1, emin), spectrum;

    double e = emin;
    while (e < emax)
    {
        e = (log_scale ? e * de : e + de);
        ph_energy.push_back(e);
    }
    // adjust max energy
    emax = ph_energy.back();
    std::cout << " Max spectral energy[keV]=" << emax << '\n';

    // compute spectrum from database
    aped.emission_spectrum(spectrum,
                           ph_energy,
                           elements,
                           abundance_model,
                           metallicity,
                           temperature,
                           doppler_shift,
                           cont_emission,
                           line_emission,
                           line_shape,
                           kernel_tolerance);

    std::cout << " size of emission spectrum is " << spectrum.size() << '\n';

    { // open file, pring header, momentum grid and initial function
        std::ofstream file(output_file, std::ios_base::trunc);

        file << "#  " << '\n'
             << "# Plasma Emissivity:" << '\n'
             << "#  " << '\n'
             << "#     Aped Database - v" << aped_version << '\n'
             << "#  " << '\n'
             << "# Parameters: " << '\n'
             << "#  n[cm^-3] =1.000 " << '\n'
             << "#  Temp[keV]=" << std::setprecision(4) << temperature << '\n'
             << "#  Emin[keV]= " << std::setprecision(4) << emin << '\n'
             << "#  Emax[keV]= " << std::setprecision(4) << emax << '\n'
             << "#  Doppler Shft = " << std::setprecision(4) << doppler_shift << '\n'
             << "#  Metallicity  = " << std::setprecision(4) << metallicity << '\n'
             << "#  Abund. Model = " << (int_abundance_model==0?"AndersGravesse":"Lodders") << '\n'
             << "#  Line Emission= " << (line_emission ? "Yes" : "No") << '\n'
             << "#  Cont Emission= " << (cont_emission ? "Yes" : "No") << '\n'
             << "#  Spectral Res.= " << std::setprecision(4) << de << '\n'
             << "#  En. log-scale= " << (log_scale ? "Yes" : "No") << '\n'
             << "# " << '\n';

        file << "#    " << '\n'
             << "#    E                  j  " << '\n'
             << "#  [keV]          [ph cm^3 s^-1] " << '\n'
             << "# "
             << '\n';

        for (size_t i = 0; i < spectrum.size(); ++i)
        {
            file << std::setw(10) << std::setprecision(4) << std::scientific << ph_energy.at(i)
                 << "          "
                 << std::setw(10) << std::setprecision(4) << std::scientific << spectrum.at(i)
                 << '\n';
        }
    }

    return 0;
}
