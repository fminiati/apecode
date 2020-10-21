#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include <xsFortran.h>
#include <xsTypes.h>
#include <FunctionUtility.h>
#include "Aped.h" // this is from $HEADAS/include

#include "../src/Aped.h" // header files with same name but from ../src
//#include "../src/XspecAped.h" // this is from $HEADAS/include
#include "FileParser.h"

int main(int argc, char *argv[])
{
    using Clock = std::chrono::steady_clock;

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

usa    auto rel_diff = [](const auto a, const auto b) {
        const auto c = 0.5 * (a + b);
        return (c > 0 ? std::abs(a - b) / b : fm::aped::zero);
    };
    auto max_rel_diff = [rel_diff](const auto &a, const auto &b) {
        double max_diff = 0;
        for (size_t i = 0; i < a.size(); ++i)
        {
            max_diff = std::max(max_diff, rel_diff(a[i], b[i]));
        }
        return max_diff;
    };

    // read input file
    fm::FileParser input(input_file);

    // read in path and version
    std::string aped_file_path;
    input.get_item(aped_file_path, "aped.aped_file_path");

    std::string aped_version;
    input.get_item(aped_version, "aped.aped_version");

    const std::string cocofile = aped_file_path + "/apec_v"+aped_version+"_coco.fits";
    const std::string linefile = aped_file_path + "/apec_v"+aped_version+"_line.fits";

    int verbosity;
    input.get_item(verbosity, "aped.verbosity");

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
    const double doppler_shift = velocity / fm::aped::c_light_cgs;

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
    int int_line_broadening = 0;
    input.get_item(int_line_broadening, "aped.line_broadening");
    fm::aped::LineBroadening line_broadening=static_cast<fm::aped::LineBroadening>(int_line_broadening);

    // default is delta function
    fm::aped::LineShape line_shape = static_cast<fm::aped::LineShape>(0);
    if (int_line_broadening)
    {
        int int_line_shape = 0;
        input.get_item(int_line_shape, "aped.line_shape");
        line_shape = static_cast<fm::aped::LineShape>(int_line_shape);
    }

    // number of elements
    int num_elements;
    input.get_item(num_elements, "aped.num_elements");

    // empty vector means all elements
    std::vector<unsigned> elements(num_elements);
    input.get_items(elements, "aped.elements");

    double kernel_tolerance;
    input.get_item(kernel_tolerance, "aped.kernel_tolerance");

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
    fm::aped::Aped fm_aped(aped_file_path, aped_version, verbosity);
    std::cout << " done! \n";

    // Initializes data directory locations needed by the models.
    std::cout << "\ninitialising Xspec's data directory... \n";
    FNINIT();
    std::cout << " done! \n";

    std::cout << "\nbuilding Xspec's aped... ";
    Aped std_aped;
    if (const int status = std_aped.Read(cocofile, linefile); status != 0)
        return status;
    std::cout << " done! \n";

    std_aped.SetThermalBroadening(int_line_broadening==1);
    std_aped.SetVelocityBroadening(velocity);
    std_aped.SetNoLines(!line_emission);
    //std_aped.SetNoResonanceLines(!line_emission);
    std_aped.SetMinimumLinefluxForBroadening(fm::aped::zero);
    //std_aped.SetBroadenPseudoContinuum(line_broadening != "none");
    std_aped.SetLogTempInterpolation(false);

    // Xspec data
    const Real temperature_kev = temperature / fm::aped::keVToKelvin;
    const Real emission_measure = 1.e-14;
    const RealArray ra_ph_energy(&ph_energy[0], ph_energy.size());

    std::vector<fm::aped::ElementAbundance> el_abundance;
    fm::aped::AbundanceUtil::relative_abundances(el_abundance,
                                                 elements,
                                                 metallicity,
                                                 abundance_model);

    // performance
    std::chrono::duration<double, std::micro> tot_fm_dur{}, tot_fm_dur1{};
    std::chrono::duration<double, std::micro> tot_xspec_dur{}, tot_xspec_dur1{};
    for (const auto A : elements)
    {
        std::cout << "\nTemp[K]= " << std::setw(10) << std::setprecision(4) << std::scientific << temperature
                  << ",  atomic number = " << A << std::endl;

        std::vector<unsigned> one_element(1, A);
        std::vector<double> fm_aped_spectrum;
        if (verbose)
            std::cout << " computing emission_spectrum with fm::aped... ";
        std::chrono::duration<double, std::micro> fm_aped_dur;
        {
#ifdef USE_TIMER
            Timer_t<> t("Aped.Ind");
#endif
            const auto t_i{Clock::now()};
            fm_aped.emission_spectrum(fm_aped_spectrum,
                                      ph_energy,
                                      one_element,
                                      abundance_model,
                                      metallicity,
                                      temperature,
                                      doppler_shift,
                                      cont_emission,
                                      line_emission,
                                      line_shape,
                                      line_broadening,
                                      kernel_tolerance);
            const auto t_e{Clock::now()};
            fm_aped_dur = t_e - t_i;
            tot_fm_dur += fm_aped_dur;
            if (A != elements[0]) tot_fm_dur1 += fm_aped_dur;
        }
        if (verbose)
            std::cout << " done! \n";

        // IntegerArray is  a std::vector<int>
        const IntegerArray ia_element(1, el_abundance[A-1].m_atomic_number);
        // RealArray is a std::valarray<Real>
        const RealArray ra_el_abundance(el_abundance[A-1].m_abundance, 1);
        RealArray ra_xspec_spectrum(0.0, ph_energy.size() - 1), ra_spectrum_err(0.0, ph_energy.size() - 1);

        if (verbose)
            std::cout << " computing emission_spectrum with xspec's aped... ";

        std::chrono::duration<double, std::micro> xspec_dur;
        {
            const auto t_i{Clock::now()};
            std_aped.SumEqSpectra(ra_ph_energy, ia_element, ra_el_abundance,
                                  doppler_shift, temperature_kev, emission_measure,
                                  ra_xspec_spectrum, ra_spectrum_err);
            const auto t_e{Clock::now()};
            xspec_dur = t_e - t_i;
            tot_xspec_dur += xspec_dur;
            if (A != elements[0]) tot_xspec_dur1 += xspec_dur;
        }
        // const int status = xspec::calcCEISpectrum(ra_ph_energy, ia_element, ra_el_abundance, doppler_shift,
        //                                    temperature_kev, emission_measure, int_line_broadening,
        //                                    doppler_shift, !line_emission, ra_xspec_spectrum, ra_spectrum_err);
        // if (status)
        //     std::cerr << " Aped status = " << status << '\n';
        if (verbose)
            std::cout << " done! \n";

        std::cout
            << " Timing: fm_aped: " << fm_aped_dur.count() << "us, xspec_aped: " << xspec_dur.count() << " us\n";
        std::cout << " Max Relative Diff fm::aped vs xspec's aped= "
                  << std::setw(10) << std::setprecision(4) << std::scientific << max_rel_diff(fm_aped_spectrum,ra_xspec_spectrum)
                  << '\n';

        if (verbose)
        {
            for (size_t i = 0; i < fm_aped_spectrum.size(); ++i)
            {
                std::cout << " E[keV]=" << std::setw(10) << std::setprecision(4) << std::scientific << ph_energy[i]
                          << " -->   aped= " << std::setw(10) << std::setprecision(4) << std::scientific << fm_aped_spectrum[i]
                          << " -->   xspec= " << std::setw(10) << std::setprecision(4) << std::scientific << ra_xspec_spectrum[i]
                          << " --->  Daped-xspec(\%)= " << std::setw(10) << std::setprecision(4) << std::scientific 
                          << rel_diff(fm_aped_spectrum[i], ra_xspec_spectrum[i])
                          << std::endl;
            }
        }
    }
    std::cout << "\nFinal timing: \n";
    std::cout << "    total.......... : fm_aped: " << tot_fm_dur.count() << "us, xspec_aped: " << tot_xspec_dur.count() << " us\n";
    std::cout << "    total w/o first : fm_aped: " << tot_fm_dur1.count() << "us, xspec_aped: " << tot_xspec_dur1.count() << " us\n";
    {
        std::vector<double> fm_aped_spectrum;
        RealArray ra_xspec_spectrum(0.0, ph_energy.size() - 1), ra_spectrum_err(0.0, ph_energy.size() - 1);
        std::chrono::duration<double, std::micro> full_xspec_dur;
        {
            // RealArray is a std::valarray<Real>
            RealArray ra_el_abundance(num_elements);
            IntegerArray ia_element(num_elements);
            for (const auto& A : el_abundance)
            {
                ra_el_abundance[A.m_atomic_number-1] = A.m_abundance;
                ia_element[A.m_atomic_number-1] = A.m_atomic_number;
            }
            const auto t_i{Clock::now()};
            std_aped.SumEqSpectra(ra_ph_energy, ia_element, ra_el_abundance,
                                  doppler_shift, temperature_kev, emission_measure,
                                  ra_xspec_spectrum, ra_spectrum_err);
            const auto t_e{Clock::now()};
            full_xspec_dur = t_e - t_i;
        }
        std::chrono::duration<double, std::micro> fm_full_dur;
        {
#ifdef USE_TIMER
            Timer_t<> t("Aped.All");
#endif
            const auto t_i{Clock::now()};
            fm_aped.emission_spectrum(fm_aped_spectrum,
                                      ph_energy,
                                      elements,
                                      abundance_model,
                                      metallicity,
                                      temperature,
                                      doppler_shift,
                                      cont_emission,
                                      line_emission,
                                      line_shape,
                                      line_broadening,
                                      kernel_tolerance);
            const auto t_e{Clock::now()};
            fm_full_dur = t_e - t_i;
        }
        std::cout << "\nFull calculation using " << elements.size() << "elements\n";
        std::cout << "    Timing..... : fm_aped: " << fm_full_dur.count() << "us, xspec_aped : " << full_xspec_dur.count() << " us\n";
        std::cout << "    Max Relative Diff fm::aped vs xspec's aped= "
                  << std::setw(10) << std::setprecision(4) << std::scientific << max_rel_diff(fm_aped_spectrum,ra_xspec_spectrum)
                  << '\n';
    }

#ifdef USE_TIMER
    Timer_t<>::print_record(std::cout);
#endif

    return 0;
}
