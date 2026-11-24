# Test executable
add_executable(test_protein_folding
    test_main.cpp
    test_state_space.cpp
    test_motion_validator.cpp
    test_planner_integration.cpp
    test_storage.cpp
    test_export.cpp
)

target_link_libraries(test_protein_folding
    protein_folding_core
    ${OMPL_LIBRARIES}
    ${Boost_LIBRARIES}
)

# Register tests with CTest
add_test(NAME StateSpaceTests COMMAND test_protein_folding state_space)
add_test(NAME MotionValidatorTests COMMAND test_protein_folding motion_validator)
add_test(NAME PlannerIntegrationTests COMMAND test_protein_folding planner)
add_test(NAME StorageTests COMMAND test_protein_folding storage)
add_test(NAME ExportTests COMMAND test_protein_folding export)
add_test(NAME AllTests COMMAND test_protein_folding all)
