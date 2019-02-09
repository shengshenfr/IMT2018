#include <iostream>
#include <ql/quantlib.hpp>  // Ce fichier seul suffit pour utiliser QuantLib

using namespace QuantLib ;

int main() {
  SimpleQuote q ;
  q.setValue (10) ;

  std::cout << q.value() << std::endl ;

 return 0 ;
}
