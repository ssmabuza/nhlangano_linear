// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>

#include <Kokkos_Core.hpp>

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
// #include <Kokkos_DefaultNode.hpp>
#include <Phalanx_MDField.hpp>

#include "Flujo_DriverFactory.hpp"

using namespace flujo;

// Bring in the version file, so we can print out the string
// note that escape sequences in this file: see
// https://stackoverflow.com/questions/410980/include-a-text-file-in-a-c-program-as-a-char
const char * version_string = {
#include "flujo_version_file.txt"
};

int main(int argc, char *argv[])
{
  int status = 0;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  if (mpiSession.getNProc() > 1)
  {
    out->setShowProcRank(true);
    out->setOutputToRootOnly(0);
    pout->setShowProcRank(true);
  }

  const auto stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Flujo Stacked Timer"));
  Teuchos::TimeMonitor::setStackedTimer(stackedTimer);

  // print out valiant versioning information
  (*out) << std::endl;
  (*out) << "                                Welcome to Flujo!                          " << std::endl;
  (*out) << std::endl;
  (*out) << "***************************************************************************" << std::endl;
  out->pushTab(3);
  (*out) << version_string << "\n" ;
  out->popTab();
  (*out) << "***************************************************************************\n" << std::endl;
  
  bool output_timings = true;
  std::string timingsFile;
  
  try {
    
    Teuchos::RCP<Teuchos::Time> total_time = Teuchos::TimeMonitor::getNewTimer("Flujo: Total Time");
    Teuchos::TimeMonitor timer(*total_time); 
    
    // *********************
    // Parse the command line arguments
    // *********************
    std::string input_file_name = "flujo.xml";
    int num_threads=-1;
    
    // Hard coded to MPI - we do not support a serial build!
    Teuchos::RCP<const Teuchos::MpiComm<int>> comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper(Teuchos::as<MPI_Comm>(MPI_COMM_WORLD)));
    {
      Teuchos::CommandLineProcessor clp;
      clp.setOption("i", &input_file_name, "Flujo input xml/yaml filename");

      clp.setOption("kokkos-threads", &num_threads);
      clp.setOption("timings-file",&timingsFile, "File containing the YAML output for the timing data, default don't produce a file");
      
      Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = 
          clp.parse(argc,argv,&std::cerr);

      TEUCHOS_TEST_FOR_EXCEPTION(parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL, 
                                 std::runtime_error, "Failed to parse command line!");
    }
    
    Kokkos::initialize(argc, argv);
    {
      // *********************
      // Parse the xml input file and broadcast to other processes
      // *********************
      Teuchos::RCP<Teuchos::ParameterList> input_params = Teuchos::rcp(new Teuchos::ParameterList("Flujo Parameters"));
      {
        const std::string check_xml = ".xml";
        const auto seach_xml = input_file_name.find(check_xml);
        if (seach_xml != std::string::npos)
        {
          Teuchos::updateParametersFromXmlFileAndBroadcast(input_file_name, input_params.ptr(), *comm);
        } else {
          
          const std::string check_yaml = ".yaml";
          const auto search_yaml = input_file_name.find(check_yaml);
          if (search_yaml != std::string::npos) {
            Teuchos::updateParametersFromYamlFileAndBroadcast(input_file_name, input_params.ptr(), *comm);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                      "ERROR: Input file named: " 
                                      << input_file_name 
                                      << "\nrequires a suffix of type \".xml\" or \".yaml\" to determine parser type!\n");
          }
        }
      }

      // Create the driver and solve the problem
      Teuchos::RCP<Driver> driver = DriverFactory::build(input_params,comm);
      // set up the driver using the given parameters
      driver->setup();
      // solve the problem and save data
      driver->solve();
    }
    Kokkos::finalize();
    
  }
  catch (std::exception& e) {
    *pout << "*********** Caught Exception std::exception: Begin Error Report ***********" << std::endl;
    *pout << e.what() << std::endl;
    *pout << "************ Caught Exception std::exception: End Error Report ************" << std::endl;
    status = -1;
  }
  catch (std::string& msg) {
    *pout << "*********** Caught Exception std::string: Begin Error Report ***********" << std::endl;
    *pout << msg << std::endl;
    *pout << "************ Caught Exception std::string: End Error Report ************" << std::endl;
    status = -1;
  }
  catch (...) {
    *pout << "*********** Caught Exception: Begin Error Report ***********" << std::endl;
    *pout << "Caught UNKOWN exception" << std::endl;
    *pout << "************ Caught Exception: End Error Report ************" << std::endl;
    status = -1;
  }
  
  if (output_timings)
    Teuchos::TimeMonitor::summarize(*out,false,true,false,Teuchos::Union);
  // *********************
  // Output time monitor information
  // *********************
  // use YAML b/c its easier to parse in Python
  if (!timingsFile.empty())
  {
    std::ofstream fout(timingsFile.c_str());
    Teuchos::RCP<Teuchos::ParameterList> reportParams = parameterList(* (Teuchos::TimeMonitor::getValidReportParameters()));
    reportParams->set("Report format", "YAML");
    reportParams->set("YAML style", "spacious");
    Teuchos::TimeMonitor::report(fout,reportParams);
  }
  
  stackedTimer->stop("Flujo Stacked Timer");
  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = true;
  options.output_minmax = true;
  options.output_histogram = true;
  options.num_histogram = 5;
  stackedTimer->report(std::cout, Teuchos::DefaultComm<int>::getComm(), options);

  if (0 == status)
    *out << "Flujo run completed." << std::endl;

  return status;
}
