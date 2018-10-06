add_executable(run_echelle_tests
        basic_test.cpp)

target_link_libraries(run_echelle_tests gtest gtest_main)
target_link_libraries(run_echelle_tests echellelib)