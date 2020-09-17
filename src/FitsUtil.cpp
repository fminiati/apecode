
#ifdef USE_APED
#include "Fits_util.H"

void FitsUtil::print_fits_error(const int status, const std::string s)
{
  if (status)
    {
      // print error report
      fits_report_error(stderr, status);
      throw std::runtime_error("A FITS error occurred"+s);
    }
}
#endif
