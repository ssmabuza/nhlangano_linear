
#ifndef ExportArray_HPP
#define ExportArray_HPP

#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

/**
 * \brief A class that exports data from an array object in one map to another
 */
template <typename ArrayType>
class ExportArray
{
public:

  //! Ctor
  ExportArray(const Teuchos::RCP<const Tpetra::Map<int,int>>& source, const Teuchos::RCP<const Tpetra::Map<int,int>>& target);

  //! Dtor
  ~ExportArray() {}

  //! main method for exporting data
  void exportData(const Teuchos::RCP<ArrayType>& sourceMat, Teuchos::RCP<ArrayType>& targetMat) const;
  
  //! get the main export object
  const Teuchos::RCP<Tpetra::Export<int,int>>& getExporter() const;

private:

  Teuchos::RCP<Tpetra::Export<int,int>> m_exporter;

};


#endif
