# Core library
add_library(protein_folding_core
    core/ProteinStateSpace.cpp
    core/ProbabilisticMotionValidator.cpp
    core/StateExporter.cpp
    core/ConfigLoader.cpp
    planners/PlannerRegistry.cpp
    rules/StateRuleRegistry.cpp
    energy/EnergyEvaluatorRegistry.cpp
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
