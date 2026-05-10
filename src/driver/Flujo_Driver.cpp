// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#include "Flujo_Driver.hpp"

#include <map>
#include <string>

#include <Shards_CellTopology.hpp>

#include <Teuchos_TestForException.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Panzer_GlobalData.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_WorksetContainer.hpp"

#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace flujo {

Driver::Driver(const Teuchos::RCP<const Teuchos::MpiComm<int>>& comm,
               const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
               const Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>>& cm_factory,
               const Teuchos::RCP<const panzer::BCStrategyFactory>& bc_factory)
    : comm_(comm),
      eqset_factory_(eqset_factory),
      cm_factory_(cm_factory),
      bc_factory_(bc_factory) {}

void Driver::setup(const Teuchos::RCP<Teuchos::ParameterList>& input_params) {
  TEUCHOS_TEST_FOR_EXCEPTION(input_params == Teuchos::null, std::logic_error,
                             "Flujo::Driver::setup: input_params is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(comm_ == Teuchos::null, std::logic_error,
                             "Flujo::Driver::setup: MPI communicator is null.");

  Teuchos::ParameterList& mesh_pl = input_params->sublist("Mesh");
  Teuchos::ParameterList& assembly_pl = input_params->sublist("Assembly");
  Teuchos::ParameterList& block_to_physics_pl = input_params->sublist("Block ID to Physics ID Mapping");

  Teuchos::RCP<Teuchos::ParameterList> physics_blocks_pl =
      Teuchos::rcp(new Teuchos::ParameterList(input_params->sublist("Physics Blocks")));

  global_data_ = panzer::createGlobalData();

  const std::size_t workset_size = static_cast<std::size_t>(assembly_pl.get<int>("Workset Size", 2000));
  const int default_integration_order = assembly_pl.get<int>("Default Integration Order", -1);
  const bool build_transient_support = assembly_pl.get<bool>("Build Transient Support", false);

  panzer_stk::ModelEvaluatorFactory<double> panzer_me_helpers;

  mesh_factory_ = panzer_me_helpers.buildSTKMeshFactory(mesh_pl);
  TEUCHOS_TEST_FOR_EXCEPTION(mesh_factory_ == Teuchos::null, std::runtime_error,
                             "Flujo::Driver::setup: mesh factory could not be constructed "
                             "(check Mesh/Source and mesh sublists).");

  mesh_ = mesh_factory_->buildUncommitedMesh(*(comm_->getRawMpiComm()));

  std::map<std::string, std::string> block_ids_to_physics_ids;
  panzer::buildBlockIdToPhysicsIdMap(block_ids_to_physics_ids, block_to_physics_pl);

  std::map<std::string, Teuchos::RCP<const shards::CellTopology>> block_ids_to_cell_topo;
  for (const auto& kv : block_ids_to_physics_ids) {
    block_ids_to_cell_topo[kv.first] = mesh_->getCellTopology(kv.first);
    TEUCHOS_TEST_FOR_EXCEPTION(block_ids_to_cell_topo[kv.first] == Teuchos::null, std::logic_error,
                             "Flujo::Driver::setup: missing cell topology for element block \""
                                 << kv.first << "\".");
  }

  std::vector<std::string> tangent_param_names;
  panzer::buildPhysicsBlocks(block_ids_to_physics_ids, block_ids_to_cell_topo, physics_blocks_pl,
                             default_integration_order, workset_size, eqset_factory_, global_data_,
                             build_transient_support, physics_blocks_, tangent_param_names);

  panzer_me_helpers.finalizeMeshConstruction(*mesh_factory_, physics_blocks_, *comm_, *mesh_);

  conn_manager_ = Teuchos::rcp(new panzer_stk::STKConnManager(mesh_));

  std::string field_order = assembly_pl.get<std::string>("Field Order", "");
  if (field_order.empty() || !panzer::BlockedDOFManagerFactory::requiresBlocking(field_order)) {
    panzer::DOFManagerFactory dof_factory;
    global_indexer_ =
        dof_factory.buildGlobalIndexer(comm_->getRawMpiComm(), physics_blocks_, conn_manager_, field_order);
  } else {
    panzer::BlockedDOFManagerFactory dof_factory;
    global_indexer_ =
        dof_factory.buildGlobalIndexer(comm_->getRawMpiComm(), physics_blocks_, conn_manager_, field_order);
  }

  Teuchos::RCP<panzer_stk::WorksetFactory> wkst_factory = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh_));
  workset_container_ = Teuchos::rcp(new panzer::WorksetContainer);
  workset_container_->setFactory(wkst_factory);
  for (std::size_t i = 0; i < physics_blocks_.size(); ++i) {
    workset_container_->setNeeds(physics_blocks_[i]->elementBlockID(), physics_blocks_[i]->getWorksetNeeds());
  }
  workset_container_->setWorksetSize(workset_size);
  workset_container_->setGlobalIndexer(global_indexer_);
}

void Driver::solve() const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR: Driver::solve not implemented!");
}

}  // namespace flujo
