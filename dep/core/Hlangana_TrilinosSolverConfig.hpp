// @HEADER
// ****************************************************************************
//                Hlangana: Copyright S. Mabuza Enterprises LLC
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hlangana_TrilinosSolverConfig_HPP__
#define __Hlangana_TrilinosSolverConfig_HPP__

#include "Hlangana_OptionHandler.hpp"

namespace hlangana
{

class TrilinosSolverConfig
  : public Optionable
{
public:

  TrilinosSolverConfig(const std::shared_ptr<OptionHandler> &optionHandler)
      : Optionable(optionHandler)
  {
    mOptionHandler = optionHandler;
    optionHandler->parse();
  }

  ~TrilinosSolverConfig() {}

  void initializeConfig(int argc, char **argv, MPI_Comm rawMpiComm);

  Teuchos::RCP<const Teuchos::MpiComm<int>> getComm() const {
    TEUCHOS_ASSERT(mInitialized);
    return mComm;
  }

private:

  virtual void addOptionsCallback(po::options_description &config) override
  {
    config.add_options()("trilinos-parameters", 
                          po::value<std::string>(&mTrilinosParams)->is_required(), 
                          "Get the XML file for all Trilinos Options");
  }

  std::shared_ptr<OptionHandler> mOptionHandler;
  std::string mTrilinosParams;
  Teuchos::RCP<const Teuchos::MpiComm<int>> mComm;
  bool mInitialized;
};

}
// namespace hlangana

#endif /** __Hlangana_TrilinosSolverConfig_HPP__ */