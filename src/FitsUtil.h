
#ifdef USE_APED

#ifndef FITSTOOLS_H
#define FITSTOOLS_H

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstring>
#include <typeinfo>
#include <cassert>

#ifdef zero
#undef zero
#endif
#include "fitsio.h"
#ifndef zero
#define zero (0.0e0)
#endif

#define MAXDIM 7

namespace fm::fits_util {

  void print_fits_error(const int status, const std::string s={});
  
  template <class T> int fits_data_type(T t)
  {
    if      (typeid(T).name() == typeid(char) .name()) return TSBYTE;
    else if (typeid(T).name() == typeid(short).name()) return TSHORT;
    else if (typeid(T).name() == typeid(bool) .name()) return TLOGICAL;
    else if (typeid(T).name() == typeid(int)  .name()) return TINT;
    else if (typeid(T).name() == typeid(long) .name()) return TLONG;
    else if (typeid(T).name() == typeid(unsigned char).name()) return TBYTE;
    else if (typeid(T).name() == typeid(unsigned short).name()) return TUSHORT;
    else if (typeid(T).name() == typeid(unsigned).name()) return TUINT;
    else if (typeid(T).name() == typeid(unsigned long).name()) return TULONG;
    else if (typeid(T).name() == typeid(long long).name()) return TLONGLONG;
    else if (typeid(T).name() == typeid(float).name()) return TFLOAT;
    else if (typeid(T).name() == typeid(double).name()) return TDOUBLE;
    //else if (typeid(T).name() == typeid(complex).name()) return TCOMPLEX;
    //else if (typeid(T).name() == typeid(double complex).name()) return TDBLCOMPLEX;
    else
    {
      std::cout << " fits_data_type " << typeid(T).name() << std::endl;
      std::cerr << "cannot use this datatype with this function" << std::endl;
    }
    // just to avoid the compiler's complains
    return TINT;
  }

  //
  template <class T> void read_fits_column(fitsfile* a_file_ptr,
                                           std::vector<T>& a_v,
                                           const std::string a_col_name,
                                           const int a_num_elements,
                                           const int a_first_row=1)
  {
    int column=-1, status=0;
    {
      char col_name[100];//=(char*)a_col_name.data();
      strcpy(col_name, a_col_name.c_str());
      fits_get_colnum(a_file_ptr, CASEINSEN, col_name, &column, &status);
      print_fits_error(status," reading column - "+a_col_name+" - in fits_get_colnum");
    }

    int anynul;
    a_v.resize(a_num_elements);
    
    fits_read_col(a_file_ptr,
                  fits_data_type(T(0)),
                  column,
                  (long)a_first_row,
                  1,
                  a_num_elements,
                  0,
                  a_v.data(),
                  &anynul,
                  &status);

    print_fits_error(status," reading data of "+a_col_name+" in fits_read_col");
  }
}

#endif // FITSTOOLS_H
#endif // USE_APED
