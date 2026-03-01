
#include "BuildMap.hpp"


BuildMap::BuildMap(const Teuchos::RCP< Teuchos::MpiComm< int > >& tcomm, 
                   const Teuchos::RCP< ParGrid2D >& par_grid) : m_comm(tcomm), m_grid(par_grid)
{
}

void BuildMap::buildStartingMap()
{
  int numOfGlobalPoints = m_grid->m_points.size();
  int numOfLocalPoints = m_grid->m_points_ind_loc.size();
  auto arrayLocalGlobIndexes = m_grid->m_points_ind_loc.data();
  m_start_map = Teuchos::rcp(new Tpetra::Map<int,int>(numOfGlobalPoints,arrayLocalGlobIndexes,numOfLocalPoints,0,m_comm));
  m_const_start_map = Teuchos::rcp(new const Tpetra::Map<int,int>(numOfGlobalPoints,arrayLocalGlobIndexes,numOfLocalPoints,0,m_comm));
  
// Constructor for Tpetra::Map :
//      Map (const global_size_t numGlobalElements,
//           const GlobalOrdinal indexList[],
//           const LocalOrdinal indexListSize,
//           const GlobalOrdinal indexBase,
//           const Teuchos::RCP<const Teuchos::Comm<int> >& comm);
}

void BuildMap::buildGhostedMap()
{
  int numOfGlobalPoints = m_grid->m_points.size();
  int numOfLocalPoints = m_grid->m_points_ind.size();
  auto arrayLocalGlobIndexes = m_grid->m_points_ind.data();
  m_ghosted_map = Teuchos::rcp(new Tpetra::Map<int,int>(numOfGlobalPoints,arrayLocalGlobIndexes,numOfLocalPoints,0,m_comm));
  m_const_ghosted_map = Teuchos::rcp(new const Tpetra::Map<int,int>(numOfGlobalPoints,arrayLocalGlobIndexes,numOfLocalPoints,0,m_comm));
}

void BuildMap::buildOneToOneMap()
{
  if (m_start_map.is_null())
  {
    std::cerr << "Error in BuildMap::buildOneToOneMap, starting map not built. Call BuildMap::buildStartingMap first." << std::endl;
    exit(1);
  }
  
  if (m_start_map->isOneToOne())
  {
    m_one_to_one_map = m_start_map;
    m_const_one_to_one_map = m_const_start_map;
  } else {
    const Tpetra::Map<int,int> rawStartMap = *(m_start_map.get());
    Teuchos::RCP<const Tpetra::Map<int,int>> tmpStartMap = Teuchos::rcp(&rawStartMap,false);
    auto tmpOneToOneMap = Tpetra::createOneToOne(tmpStartMap);
    Tpetra::Map<int,int> rawOneToOneMap = *(tmpOneToOneMap.get());
    m_one_to_one_map = Teuchos::rcp(&rawOneToOneMap,false);

    // form the const 1-1 map
    m_const_one_to_one_map = Tpetra::createOneToOne(m_const_start_map);

  }
}

const Teuchos::RCP< Tpetra::Map< int, int > >& BuildMap::getStartingMap() const
{
  return m_start_map;
}

const Teuchos::RCP< Tpetra::Map< int, int > >& BuildMap::getGhostedMap() const
{
  return m_ghosted_map;
}

const Teuchos::RCP< Tpetra::Map< int, int > >& BuildMap::getOneToOneMap() const
{
  return m_one_to_one_map;
}

const Teuchos::RCP<const Tpetra::Map<int,int>>& BuildMap::getConstStartingMap() const
{
  return m_const_start_map;
}

const Teuchos::RCP<const Tpetra::Map<int,int>>& BuildMap::getConstGhostedMap() const
{
  return m_const_ghosted_map;
}

const Teuchos::RCP<const Tpetra::Map<int,int>>& BuildMap::getConstOneToOneMap() const
{
  return m_const_one_to_one_map;
}

int BuildMap::getNumberOfOwnedPoints() const
{
  return m_grid->m_points_ind_loc.size();
}
