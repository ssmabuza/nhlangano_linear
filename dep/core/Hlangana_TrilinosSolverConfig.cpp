// @HEADER
// ****************************************************************************
//                Hlangana: Copyright S. Mabuza Enterprises LLC
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_StackedTimer.hpp>

#include "Hlangana_TrilinosSolverConfig.hpp"

namespace hlangana
{

  void TrilinosSolverConfig::initializeConfig(int argc, char **argv, MPI_Comm rawMpiComm)
  {
    mComm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper(rawMpiComm));
    Teuchos::RCP<Teuchos::ParameterList> inputParams 
      = Teuchos::rcp(new Teuchos::ParameterList("Hlangana Trilinos Parameters"));
    {
      const std::string checkXML = ".xml";
      const auto seachXML = mTrilinosParams.find(checkXML);
      if (seachXML != std::string::npos)
      {
        Teuchos::updateParametersFromXmlFileAndBroadcast(mTrilinosParams, inputParams.ptr(), *mComm);
      }
      else
      {

        const std::string checkYaml = ".yaml";
        const auto searchYaml = mTrilinosParams.find(checkYaml);
        if (searchYaml != std::string::npos)
        {
          Teuchos::updateParametersFromYamlFileAndBroadcast(mTrilinosParams, inputParams.ptr(), *mComm);
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                                     "ERROR: Input file named: "
                                     << mTrilinosParams
                                     << "\nrequires a suffix of type \".xml\" "
                                     <<"or \".yaml\" to determine parser type!\n");
        }
      }
    }
  }

} // namespace hlangana
