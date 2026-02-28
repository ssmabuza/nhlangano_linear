// ==============================================================================
//                Hlangana: Copyright S. Mabuza Enterprises LLC
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#include <iostream>

#include <Teuchos_ParameterList.hpp>

int main(int argc, char **argv)
{

  Teuchos::ParameterList params;
  params.set("Parameter1", 42);
  params.print(std::cout);

  std::cout << "Hello World!" << std::endl;
  return 0;
}