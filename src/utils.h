#ifndef JSGEODA_UTILS
#define JSGEODA_UTILS

#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>

// c++14
inline bool iequals(const std::string& a, const std::string& b)
{
    return boost::iequals(a, b);
    //return std::equal(a.begin(), a.end(),
    //                  b.begin(), b.end(),
    //                  [](char a, char b) {
    //                      return tolower(a) == tolower(b);
    //                  });
}
#endif