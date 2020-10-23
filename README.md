Apecode is a code to compute astrophysical plasma spectral emission based on the AtomDB database.

The main code is contained in the class object Aped (Aped.h, Aped.cpp) while additional utility 
functions are in Util.h and FitsUtil.h. The code structure reflects the AtomDB database as 
represented in the XSPEC libraries. It is meant to be a standalone version to be used and linked 
with personal analysis or simulation codes without the dependency on large libraries. The main 
code and some support classes/structs maintain similar names to the original implementations
for easy of recognition but are wrapped in a namespace (fm) to avoid name conflict. Thus in 
addition to the support classes there is an Aped code which computes the emission spectra
straight from the database, by summing up the contribution of each atomic species according
to its abundance and ionization state as determined by the emitting plasma's temperature in 
accord with collisional ionization equilibrium. The case of non-collisional equilibrium is 
left for future implementation.

Aped uses a user defined (in Util.h) datatype "Real" for both the database data, input parameter 
and returned spectrum.

In addition to C++17, compilation requires cfitsio libraries given the FITS binary file format 
of the AtomDB database (I have used so far clang and gnu compilers). Execution requires obviously
the above database files. The code works seamlessly with either the traditional format using 51 
temperature bins or the upgraded version with 201 temperature bins, between 10^4 and 10^9 K. I have 

The code tst/test_apecode_vs_heasoft.cpp carries out an accuracy and performance comparison with the
original Xspsec's Aped.h code both element by element as well as for the full 28 element spectral 
calculation, showing a max relative difference in accuracy of order 10^-5 (partly affected by 
Xspsec's optimizations of the erf function). (tst/test_kernel.cpp is a very simple code looking 
at the  behaviour of the build_kernel function for Gaussian and Lorentzian shapes, to get an idea 
of the slow convergence of the Lorentzian function, while exe/aped_spectrum.cpp is a simple use 
example of Aped.h).

There are three API's to Aped.h: the one shown here is an ordinary member function which 
specifies the spectrum calculation through a set of input parameters:

    // photon emission spectrum in ph cm^3 s^-1
    void emission_spectrum(std::vector<Real> &a_spectrum,           // output emission spectrum
                           const std::vector<Real> &a_energy,       // input energy bins. a_spectrum[i] is the
                                                                    // emissivity between a_energy[i] and a_energy[i+1]
                           const std::vector<unsigned> &a_elements, // atomic numbers of emitting elements; empty vector
                                                                    // means all 28 aped elements
                           const AbundanceModel a_abundances_model, // select abundance model from enum type AbundanceModel
                           const Real a_metallicity,                // metallicity
                           const Real a_temperature,                // temperature in keV
                           const Real a_doppler_shift,              // v_plasma/c_light
                           const bool a_cont_emission,              // whether to include continuum emission
                           const bool a_line_emission,              // whether to include line (and pseudo cont) emission
                           const LineShape a_line_profile,          // select line profile from enum type LineShape
                           const LineBroadening a_line_broadening,  // select broadening mechanism from enum type LineBroadening
                           const Real a_kernel_tolerance = 1.e-6) const; // tolerated error for convolution of line emission on the 
                                                                    // spectrum's grid when a_line_shape is not delta function

In particular the a_abundances_model parameter selects the abundance model from the following enum types:

    enum class AbundanceModel : char  { AndersGrevesse = 0, Lodders = 1 };

This is quite limited but as shown further below there is the possibility to specify any input abundance model.
Likewise, a_line_profile and a_line_broadening select, respectively, the line shape and broadening 
mechanism from a limited set of choices defined by the following enum types:

    enum class LineShape : char { delta = 0,  gaussian = 1, lorentzian = 2, pseudovoigt = 3 };

    enum class LineBroadening : char { none = 0, thermal = 1 /*, turbulent = 2, collisional = 3 */ };

LineShape determines the shape of the emission line and LineBroadening the function calculating 
the full-width-at-half-maximum. This polymorphic behavior is expressed through the implementation 
of a function template whose template parameters are object classes with static function providing
the required functionality. This takes us to the second main API for computing the plasma emission
spectrum shown below. Notice that objects classes for the default cases enlisted in the above enum types
(LineShape and LineBroadening) are already provided (in Util.h). However, this second API enables 
the user to specify arbitrary line shape and broadening mechanism through template parameters of 
their own definition. The second API is as follows (function parameters with the same name have 
the same meaning as above):

    template <typename LineProfile<typename Shape, typename Broadening, bool PseudoContBrd>>
    void emission_spectrum(std::vector<Real> &a_spectrum,
                           const std::vector<Real> &a_energy,
                           const std::vector<unsigned> &a_elements, 
                           const AbundanceModel a_abundances_model,
                           const Real a_metallicity,
                           const Real a_temperature,
                           const Real a_doppler_shift,
                           const bool a_cont_emission,
                           const bool a_line_emission,
                           const Real a_kernel_tolerance) const;

Here Shape is a template parameter which contains (at least) a function "tail_integral(const Real x)"
returning the integral of the line shape from the line centre to x, the signed distance therefrom 
normalised to half the FWHM. For example, for a Gaussian shape we have defined the following object:

    struct Gaussian
    {
        static inline Real tail_integral(const Real a_x)  { return half * std::erf(sqrt_ln2 * a_x); }
    };

It is possible to have skewed distributions, the only requirement is that the shape has normalised
total area, i.e.: 

    tail_integral(+inf)-tail_integral(-inf) = 1.

Likewise the Broadening template parameter contains (at least) a function returning the 
the line's full-width-at-half-maximum. For thermal broadening we use:

    struct ThermalBroadening
    {
        static constexpr Real sqrt_8_ln2_to_c_cgs = two * sqrt_two * sqrt_ln2 / c_light_cgs;
        static inline Real fwhm(const Real a_temp, const Real a_atomic_mass, const Real a_ph_energy)
        {
            return sqrt_8_ln2_to_c_cgs * std::sqrt(kB_cgs * a_temp / a_atomic_mass) * a_ph_energy;
        }
        static constexpr spacing_t spacing() { return Spacing::log_uniform; }
    };

Here the spacing function returns the type of spacing associated to this mechanism and is used for 
optimising the convolution in case the energy grid shares the same spacing type. In this case
it's log_uniform because it is proportional to the energy itself (and would be uniform in log space).

The PseudoContBrd is a boolean type template which specifies whether or not the broadening 
should also be applied to pseudo continuum. Its default value is set to false as in Xspec (see below).

Finally LineProfile is the following class template:

    template <typename Shape, typename Broadening, bool PseudoContBrd=false> struct LineProfile
    {
        static inline Real tail_integral(const Real a_x) { return Shape::tail_integral(a_x); }
        static inline Real fwhm(const Real a_t, const Real a_m, const Real a_e) { return Broadening::fwhm(a_t, a_m, a_e); }
    };

In the case of a Voigt profile we use a slightly more complex object because now there are 
two broadening mechanisms at work. So for a pseudo Voigt type of profile given by a linear 
combination of a Gaussian and Lorentzian shape with linear parameter eta, we can write:

    template <typename Voigt, typename G, typename L, template<typename...> typename Broadening, bool PseudoContBrd>
    struct LineProfile<Voigt, Broadening<G,L>, PseudoContBrd>
    {
        static inline Real tail_integral(const Real a_x) { return Voigt::tail_integral(a_x); }
        static inline Real fwhm(const Real a_t, const Real a_m, const Real a_e)
        {
            const auto [f, eta] = Broadening<G,L>::fwhm_and_eta(a_t, a_m, a_e);
            Voigt::set_eta(eta);
            return f;
        }
    };

and provide a 

    template <typename GaussianBroadening, typename LorentzianBroadening>
    struct PseudoVoigtBroadening
    {
        static inline Real fwhm(const Real a_t, const Real a_m, const Real a_e)
        {
            const auto [f,eta] = fwhm_and_eta(a_t, a_m, a_e);
            return f;
        }
        // define fwhm and eta
        static inline auto fwhm_and_eta(const Real a_t, const Real a_m, const Real a_e)
        {
            const Real fwhmG = GaussianBroadening::fwhm(a_t, a_m, a_e);
            const Real fwhmL = LorentzianBroadening::fwhm(a_t, a_m, a_e);
            const Real fwhm = fwhm_USER_DEF_FUNCTION(fwhmG,fwhmL);
            const Real eta = eta_USER_DEF_FUNCTION(fwhmG,fwhmL);
            return std::pair(f, eta);
        }
        static constexpr spacing_t spacing()
        {
            // sometimes Lorentizian are used with thermal broadening...
            return GaussianBroadening::spacing() == LorentzianBroadening::spacing() ? GaussianBroadening::spacing() : Spacing::irregular;
        }
    };

A few additional examples are available in Util.H. However, at this point it should be straightforward for the user
to define their own shape and broadening objects.

There is a third, final and most general API to Aped which is still a function template as the previous one
but offers the possibility to specify the abundance of each element:

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

Thus instead of a_abundance_model and a std::vector<unsigned> of atomic numbers here the API
takes as input a std::vector<ElementAbundance>. ElementAbundance is the following struct type:

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

containint an unsigned member data for the atomic number and a Real member data for its abundance.
Notice that the abundances must be expressed as relative values with respect to the Anders & Gravesse
values, which is the assumed model for the emissivities of the AtomDB database. So if H_ag is the 
Anders & Gravesse value for the hydrogen abundance and you want to specify H_foo, you have to define 
Element(1, H_foo/H_ag).
