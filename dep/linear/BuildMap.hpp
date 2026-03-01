#ifndef BuildMap_HPP
#define BuildMap_HPP


// Wrapper for MPI_Comm
#include <Teuchos_DefaultMpiComm.hpp>

// Teuchos reference counted pointers
#include <Teuchos_RCP.hpp>

// Tpetra Map
#include <Tpetra_Map.hpp>

// Parallel mesh
#include <ParGrid2D.hpp>


/**
 * \brief A class for building the maps from a parallel mesh
 * 
 * Creates three maps that are needed in the simulation. The first
 * is a map that includes overlapping with other procs. The second
 * is a map that includes ghosted global vertices for each proc, &
 * the third is a one to one map for forming linear alg. objects.
 */
class BuildMap
{

public:

  //! Default Ctor
  BuildMap(const Teuchos::RCP<Teuchos::MpiComm<int>>& tcomm, const Teuchos::RCP<ParGrid2D>& par_grid);

  //! Dtor
  virtual ~BuildMap() {}

  //! Builds a map with shared interface vertices
  void buildStartingMap();

  //! Builds a map with larger overlap including ghosted vertices
  void buildGhostedMap();
  
  //! Builds a one to one map
  void buildOneToOneMap();

  //! Returns the pointer to the starting map
  const Teuchos::RCP<Tpetra::Map<int,int>>& getStartingMap() const;

  //! Returns the pointer to the ghosted map
  const Teuchos::RCP<Tpetra::Map<int,int>>& getGhostedMap() const;

  //! Returns the pointer to the final one to one map
  const Teuchos::RCP<Tpetra::Map<int,int>>& getOneToOneMap() const;

  //! Returns the shared pointer to a constant starting map
  const Teuchos::RCP<const Tpetra::Map<int,int>>& getConstStartingMap() const;

  //! Returns the shared pointer to a constant ghosted map
  const Teuchos::RCP<const Tpetra::Map<int,int>>& getConstGhostedMap() const;

  //! Returns the shared pointer to a constant one to one map
  const Teuchos::RCP<const Tpetra::Map<int,int>>& getConstOneToOneMap() const;
  
  //! Build all maps
  inline void buildAllMaps()
  {
    buildStartingMap();
    buildGhostedMap();
    buildOneToOneMap();
  }

  //! number of owned points
  int getNumberOfOwnedPoints() const;

private:

  //! The main parallel communicator
  const Teuchos::RCP<Teuchos::MpiComm<int>> m_comm;

  //! The mail parallel mesh
  const Teuchos::RCP<ParGrid2D> m_grid;

  //! The starting map
  Teuchos::RCP<Tpetra::Map<int,int>> m_start_map;

  //! The ghosted map
  Teuchos::RCP<Tpetra::Map<int,int>> m_ghosted_map;

  //! The one to one map
  Teuchos::RCP<Tpetra::Map<int,int>> m_one_to_one_map;

  //! The constant starting map
  Teuchos::RCP<const Tpetra::Map<int,int>> m_const_start_map;

  //! The constant ghosted map
  Teuchos::RCP<const Tpetra::Map<int,int>> m_const_ghosted_map;

  //! The constant one to one map
  Teuchos::RCP<const Tpetra::Map<int,int>> m_const_one_to_one_map;

};


#endif
