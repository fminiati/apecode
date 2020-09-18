Apecode is a code to compute astrophysical plasma spectral emission base on the AtomDB database.

The code structure reflects the AtomDB database as represented in the XSPEC libraries.
It is meant to be a standalone version to be used and linked with personal analysis or 
simulation codes without the dependency on large libraries.
The main code and support classes/structs maintain similar names to the original implementations
for easy of recognition but are wrapped in a namespace (fm) to avoid name conflict.
Thus in addition to the support classes there is an Aped code which computes the emission spectra
straight from the database, by summing up the contribution of each atomic species according
to its abundance and ionization state as determined by the emitting plasma temperature and
density parameters. And there is an Apec class which computes the emission spectrum using a 
temperature and density dependent lookup table, previously created with the Aped code itself.
This is more efficient if one needs to calculate the emission spectrum many times. The plasma 
is supposed to be in thermal equilibrium.

Both Aped and Apec use native float data for database objects but extract spectra in a user 
defined datatype Real.

Compilation requires cfitsio librieries given the fits file format of the AtomDB database.
Execution currently requires the apec_3.0.9_51 files version of the database, namely the 
one using the standard binning with 51 temeprature between 10^4 and 10^9 K.
