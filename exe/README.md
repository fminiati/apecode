# Enviromental variables
Set the following variables to the appropriate values in order to run xspec code
(change accordingly if not using bash shell type):

        export HEADAS=absolute_path_to_heasoft/platform
        export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH":$HEADAS/lib"

For example:

        export HEADAS=/usr/local/lib/heasoft/heasoft-6.28/x86_64-apple-darwin19.6.0/
        export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH":$HEADAS/lib"
