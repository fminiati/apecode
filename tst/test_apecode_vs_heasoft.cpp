#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <xsFortran.h>
#include <xsTypes.h>
#include <FunctionUtility.h>
#include "Aped.h" // this is from $HEADAS/include

#include "../src/Aped.h" // header files with same name but from ../src
#include "FitsUtil.h"
#include "FileParser.h"

int main(int argc, char *argv[])
{
    // get input file
    bool verbose = false;
    std::string input_file = "unknown";

    for (int i = 1; i < argc; ++i)
    {
        if (strncmp(argv[i], "-v", 3) == 0)
            verbose = true;
        if (strncmp(argv[i], "-f", 3) == 0)
            input_file = argv[i + 1];
    }
    if (input_file == "unknown")
    {
        const std::string prog(argv[0]);
        std::cout << "\nrun this prog with: " + prog + " -f input_filename\n\n";
        return 1;
    }

    // read input file
    fm::FileParser input(input_file);

    // read in path and version
    std::string aped_file_path;
    input.get_item(aped_file_path, "aped.aped_file_path");

    std::string aped_version;
    input.get_item(aped_version, "aped.aped_version");

    const std::string adb_path = "/Users/francesco.miniati/Work/clab/atomdb_v3.0.9";
    const std::string cocofile = aped_file_path + "/apec_v"+aped_version+"_coco.fits";
    const std::string linefile = aped_file_path + "/apec_v"+aped_version+"_line.fits";

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
        de = std::pow(fm::aped::ten, de);

    // plasma temperature in keV
    double temperature = 1.0;
    input.get_item(temperature, "aped.plasma_temperature[keV]");
    // convert from keV to K
    temperature *= fm::aped::keVToKelvin;

    // plasma velocity in cm/sec
    double velocity = 0.0;
    input.get_item(velocity, "aped.plasma_velocity[cm/s]");

    // Doppler shift
    const double doppler_shift = velocity / fm::aped::C_cm_s;

    // metallicity
    double metallicity = 0.3;
    input.get_item(metallicity, "aped.plasma_metallicity");

    // abundance model
    std::string abundance_model = "Lodders";
    input.get_item(abundance_model, "aped.plasma_abundance_model");

    // line emission
    bool line_emission = false;
    input.get_item(line_emission, "aped.use_line_emission");

    // line emission
    bool cont_emission = false;
    input.get_item(cont_emission, "aped.use_cont_emission");

    // no line broadening
    std::string line_broadening;
    input.get_item(line_broadening, "aped.line_broadening");

    // number of elements
    int num_elements;
    input.get_item(num_elements, "aped.num_elements");

    // empty vector means all elements
    std::vector<unsigned> elements(num_elements);
    input.get_items(elements, "aped.elements");

    // spectral energy and emission
    std::vector<double> ph_energy(1, emin);
    {
        double e = emin;
        while (e < emax)
        {
            e = (log_scale ? e * de : e + de);
            ph_energy.push_back(e);
        }
    }
    // adjust max energy
    emax = ph_energy.back();
    std::cout << "\nmax spectral energy in keV reset to=" << emax << '\n';

    std::cout << "\nbuilding fm::aped::aped... ";
    fm::aped::Aped fm_aped(aped_file_path, aped_version);
    std::cout << " done! \n";

    // Initializes data directory locations needed by the models.
    std::cout << "\ninitialising Xspec's data directory... \n";
    FNINIT();
    std::cout << " done! \n";

    // std::cout << "\nbuilding Xspec's aped... ";
    // Aped std_aped;
    // const int status = std_aped.Read(cocofile, linefile);
    // if (status != 0)
    //     return status;
    // std::cout << " done! \n";

    // std_aped.SetThermalBroadening(line_broadening != "none");
    // std_aped.SetVelocityBroadening(velocity);
    // std_aped.SetNoLines(!line_emission);
    // std_aped.SetNoResonanceLines(!line_emission);
    // std_aped.SetMinimumLinefluxForBroadening(zero);
    // std_aped.SetBroadenPseudoContinuum(line_broadening != "none");
    //std_aped.SetLogTempInterpolation(log_scale);

    for (const auto A : elements)
    {
        std::cout << "\nTemp[K]= " << std::setw(10) << std::setprecision(4) << std::scientific << temperature
                  << ",  atomic number = " << A << std::endl;
        // std::cout << " xspec abund = " << FunctionUtility::getAbundance(A)
        //   << " xspec A&G abund = " << FunctionUtility::getAbundance("angr", A) << '\n';

        std::vector<fm::aped::ElementAbundance> el_abundance;
        fm::aped::AbundanceUtil::relative_abundances(el_abundance, std::vector<unsigned>(1, A), metallicity, abundance_model);

        std::vector<double> fm_aped_spectrum;
        if (verbose)
            std::cout << " computing emission_spectrum with fm::aped... ";
        fm_aped.emission_spectrum(fm_aped_spectrum,
                                  ph_energy,
                                  el_abundance,
                                  temperature,
                                  doppler_shift,
                                  line_broadening,
                                  line_emission,
                                  cont_emission);
        if (verbose)
            std::cout << " done! \n";

        const Real temperature_kev = temperature / fm::aped::keVToKelvin;
        //const Real temperature_broadening = 0.000861739; //0.e0; //temperature / fm::aped::keVtoKelvin;
        const Real emission_measure = 1.e-14;
        // IntegerArray is  a std::vector<int>
        const IntegerArray ia_element(1,el_abundance[0].atomic_number);
        // RealArray is a std::valarray<Real>
        const RealArray ra_el_abundance( el_abundance[0].abundance, 1);
        const RealArray ra_ph_energy(&ph_energy[0], ph_energy.size());
        RealArray ra_xspec_spectrum(zero, ph_energy.size()-1), ra_spectrum_err(zero, ph_energy.size()-1);

        if (verbose)
            std::cout << " computing emission_spectrum with xspec's aped... ";
        // std_aped.SumEqSpectra(ra_ph_energy, ia_element, ra_el_abundance, doppler_shift,
        //                       temperature_kev, temperature_broadening, emission_measure,
        //                       ra_xspec_spectrum, ra_spectrum_err);
        const int status = calcCEISpectrum(ra_ph_energy, ia_element, ra_el_abundance, doppler_shift,
                                           temperature_kev, emission_measure, false,
                                           doppler_shift, true, ra_xspec_spectrum, ra_spectrum_err);
        if (status) std::cerr << " Aped status = " << status << '\n';
        if (verbose)
            std::cout << " done! \n";

        auto rel_diff = [](const auto a, const auto b) {
            const auto c = 0.5 * (a + b);
            return (c > 0 ? std::abs(a - b) / b : zero);
        };

        double aped_vs_xspec_max_rel_diff = 0.0;
        for (size_t i = 0; i < fm_aped_spectrum.size(); ++i)
        {
            aped_vs_xspec_max_rel_diff = std::max(aped_vs_xspec_max_rel_diff,rel_diff(fm_aped_spectrum[i],ra_xspec_spectrum[i]));
        }
        std::cout << " Max Relative Diff fm::aped vs xspec's aped= " << std::setw(10) << std::setprecision(4) << std::scientific << aped_vs_xspec_max_rel_diff
                  << '\n';

        if (verbose)
        {
            for (size_t i = 0; i < fm_aped_spectrum.size(); ++i)
            {
                std::cout << " E[keV]=" << std::setw(10) << std::setprecision(4) << std::scientific << ph_energy[i]
                          << " -->   aped= " << std::setw(10) << std::setprecision(4) << std::scientific << fm_aped_spectrum[i]
                          << " -->   xspec= " << std::setw(10) << std::setprecision(4) << std::scientific << ra_xspec_spectrum[i]
                          //      << "  --->   eps[ph cm^3 s^-1]= " << std::setw(10) << std::setprecision(4) << std::scientific << spc_va2[i]
                          << " --->  Daped-xspec(\%)= " << std::setw(10) << std::setprecision(4) << std::scientific << rel_diff(fm_aped_spectrum[i], ra_xspec_spectrum[i])
                          //   << "  --->   Deps/eps= " << std::setw(10) << std::setprecision(4) << std::scientific << (spc_va[i]-spc_va2[i])/spc_va[i]
                          << std::endl;
            }
        }
    }
    // read parameters from input fits files:
    // read temperature data

    return 0;
}
