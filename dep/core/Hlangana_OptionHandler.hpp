// @HEADER
// ****************************************************************************
//                Hlangana: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hlangana_OptionHandler_HPP__
#define __Hlangana_OptionHandler_HPP__

#include "Hlangana.hpp"

namespace hlangana
{

  namespace po = boost::program_options;

  class Optionable;

  /**
   * @brief A class to handle options and objects which need to be initialised using these options.
   */
  class OptionHandler
  {
  public:
    /**
     * @brief Constructs a OptionHandler.
     *
     * @param argc Number of command line arguments.
     * @param argv The command line arguments.
     * @param generic The generic options.
     * @param hidden The hidden options.
     * @param configfilename The name of the configuration file.
     */
    OptionHandler(int argc,
                  char **argv,
                  po::options_description generic,
                  po::options_description hidden,
                  std::string configfilename)
        : argc(argc), argv(argv), generic(generic), hidden(hidden), configfilename(configfilename)
    {
      adddouble = false;
    }

    /**
     * @brief Parses the options.
     *
     * @param checkunused If true, it checks whether there are options which are specified by the user but not used by the program.
     */
    void parse(bool checkunused = false);

    /** @brief Retrieve the double value connect to the option name with the givendefaultvalue*/
    double getdoublevalue(std::string name, double defaultvalue)
    {
      adddouble = true;
      doublename = name;
      doubledefaultvalue = defaultvalue;
      parse();
      adddouble = false;
      return doublevalue;
    }

  private:
    bool adddouble;
    std::string doublename;
    double doubledefaultvalue;
    double doublevalue;

  public:
    /**
     * @brief Adds the Optionable Object to control of this OptionHandler.
     *
     * @param callback The Optionable object to add.
     */
    void addAddOptionsCallback(Optionable &callback)
    {
      callbacks.push_back(&callback);
    }

    /**
     * @brief Removes the Optionable Object from the control of this OptionHandler.
     *
     * @param callback The Optionable object to remove.
     */
    void removeAddOptionsCallback(Optionable &callback)
    {
      for (std::size_t i = 0; i < callbacks.size(); ++i)
      {
        if (&callback == callbacks[i])
        {
          callbacks[i] = nullptr;
          break;
        }
      }
    }

    virtual ~OptionHandler() {}

  private:
    std::vector<Optionable *> callbacks;
    int argc;
    char **argv;
    po::options_description generic;
    po::options_description hidden;
    std::string configfilename;
  };

  /**
   * @brief A Class which can be configured by options.
   *        The option handler takes care of this.
   */
  class Optionable
  {
  public:
    /**
     * @brief Adds the options of the object to the options description.
     * @param config The options description to which the options are added.
     */
    virtual void addOptionsCallback(po::options_description &config) = 0;

    Optionable() : m_option_handler(nullptr) {}

    /**
     * Constructs the object which is configurable by options.
     * The option_handler takes care of this construction.
     */
    Optionable(const std::shared_ptr<OptionHandler> &option_handler)
    {
      m_option_handler = option_handler;
      option_handler->addAddOptionsCallback(*this);
    }

    virtual ~Optionable()
    {
      if (m_option_handler)
      {
        m_option_handler->removeAddOptionsCallback(*this);
      }
    }

  protected:
    /** The corresponding optionHandler.*/
    std::shared_ptr<OptionHandler> m_option_handler;
  };

  inline void OptionHandler::
      parse(bool checkunused)
  {
    po::options_description config("Configuration");
    std::string simulationname;
    config.add_options()("Simulation,S", po::value<std::string>(&simulationname)->default_value("Simulation"), "Simulation");

    for (std::size_t i = 0; i < callbacks.size(); ++i)
    {
      if (callbacks[i])
      {
        callbacks[i]->addOptionsCallback(config);
      }
    }

    if (adddouble)
    {
      config.add_options()(doublename.c_str(), po::value<double>(&doublevalue)->default_value(doubledefaultvalue), doublename.c_str());
    }

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);
    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);
    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, argv).options(cmdline_options).allow_unregistered().run();
    po::store(parsed, vm);
    std::ifstream ifs(configfilename.c_str());
    if (checkunused)
    {
      po::store(po::parse_config_file(ifs, config_file_options), vm);
    }
    else
    {
      po::store(po::parse_config_file(ifs, config_file_options, true), vm);
    }
    po::notify(vm);

    if (checkunused)
    {
      if (vm.count("help"))
      {
        std::cout << generic << config << "\n";
        exit(0);
      }

      std::vector<std::string> to_pass_further = collect_unrecognized(parsed.options, po::include_positional);
      if (!to_pass_further.empty())
      {
        std::cout << "Unrecognized options:\n";
        for (std::size_t i = 0; i < to_pass_further.size(); i++)
        {
          std::cout << to_pass_further[i] << " ";
        }
        std::cout << std::endl;
      }
    }
  }

}
// end namespace hlangana

#endif /** __Hlangana_OptionHandler_HPP__ */