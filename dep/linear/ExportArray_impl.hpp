
#ifndef ExportArray_Impl_HPP
#define ExportArray_Impl_HPP

#include "ExportArray.hpp"

template <typename ArrayType>
ExportArray<ArrayType>::ExportArray(const Teuchos::RCP<const Tpetra::Map<int,int>>& source, 
                                    const Teuchos::RCP<const Tpetra::Map<int,int>>& target)
{
  m_exporter = Teuchos::rcp(new Tpetra::Export<int,int>(source,target));
}

template <typename ArrayType>
void ExportArray<ArrayType>::exportData(const Teuchos::RCP<ArrayType>& sourceMat, 
                                        Teuchos::RCP<ArrayType>& targetMat) const
{
  targetMat->doExport(*sourceMat, *m_exporter, Tpetra::CombineMode::INSERT);
}

template <typename ArrayType>
const Teuchos::RCP<Tpetra::Export<int,int>>& ExportArray<ArrayType>::getExporter() const
{
  return m_exporter;
}

#endif

