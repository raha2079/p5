# Test executable
add_executable(test_protein_folding
    test_main.cpp
    test_state_space.cpp
    test_motion_validator.cpp
    test_planner_integration.cpp
    test_storage.cpp
    test_export.cpp
    test_protein_reconstruction.cpp
    test_projection_bounds.cpp
    test_native_sanity.cpp
    test_reconstruction_full.cpp
    test_geometry_verification.cpp
)

target_link_libraries(test_protein_folding
    protein_folding_core
    ${OMPL_LIBRARIES}
    ${Boost_LIBRARIES}
)

# Target-steering focused tests
add_executable(test_target_steering
    test_target_steering.cpp
)

target_link_libraries(test_target_steering
    protein_folding_core
    ${OMPL_LIBRARIES}
    ${Boost_LIBRARIES}
)

# Simple target-steering trace test
add_executable(test_target_steering_simple
    test_target_steering_simple.cpp
)

target_link_libraries(test_target_steering_simple
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
add_test(NAME TargetSteeringTests COMMAND test_target_steering)
