# ==============================================================================
#                Flujo: Copyright Valiant Scientific
#
# Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
# ==============================================================================

INCLUDE(TribitsLibraryMacros)

TRIBITS_CONFIGURE_FILE(Flujo_config.hpp)

TRIBITS_LIBRARY(Flujo)

TRIBITS_LIBRARY_TEST(Flujo_test)

