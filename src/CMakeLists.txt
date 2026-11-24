# Core library
add_library(protein_folding_core
    core/ProteinStateSpace.cpp
    core/ProbabilisticMotionValidator.cpp
    core/StateExporter.cpp
    core/ConfigLoader.cpp
    planners/PlannerRegistry.cpp
    planners/QuadtreePlanner.cpp
    rules/StateRuleRegistry.cpp
    energy/EnergyEvaluatorRegistry.cpp
    protein/ProteinStructure.cpp
    protein/BackboneReconstruction.cpp
    models/ProteinModel.cpp
    projection/Projection2D.cpp
    projection/ProjectionBounds.cpp
    quadtree/AdaptiveQuadtree.cpp
    experiments/ExperimentRunner.cpp
)

target_link_libraries(protein_folding_core
    ${OMPL_LIBRARIES}
    ${Boost_LIBRARIES}
)

# Main executable for running experiments
add_executable(protein_planner
    main.cpp
)

target_link_libraries(protein_planner
    protein_folding_core
    ${OMPL_LIBRARIES}
    ${Boost_LIBRARIES}
)

# Experiment runner executable
add_executable(experiment_runner
    experiments/experiment_runner_main.cpp
)

target_link_libraries(experiment_runner
    protein_folding_core
    ${OMPL_LIBRARIES}
    ${Boost_LIBRARIES}
)
