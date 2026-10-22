# CMake generated Testfile for 
# Source directory: /home/artur/inzyn/quantum_decaf/tests
# Build directory: /home/artur/inzyn/quantum_decaf/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[test_dist]=] "/home/artur/inzyn/quantum_decaf/tests/test_dist")
set_tests_properties([=[test_dist]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/artur/inzyn/quantum_decaf/cmake/util.cmake;8;add_test;/home/artur/inzyn/quantum_decaf/tests/CMakeLists.txt;8;add_logging_test;/home/artur/inzyn/quantum_decaf/tests/CMakeLists.txt;0;")
add_test([=[test_integral]=] "/home/artur/inzyn/quantum_decaf/tests/test_integral")
set_tests_properties([=[test_integral]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/artur/inzyn/quantum_decaf/cmake/util.cmake;8;add_test;/home/artur/inzyn/quantum_decaf/tests/CMakeLists.txt;8;add_logging_test;/home/artur/inzyn/quantum_decaf/tests/CMakeLists.txt;0;")
add_test([=[test_lap]=] "/home/artur/inzyn/quantum_decaf/tests/test_lap")
set_tests_properties([=[test_lap]=] PROPERTIES  _BACKTRACE_TRIPLES "/home/artur/inzyn/quantum_decaf/cmake/util.cmake;8;add_test;/home/artur/inzyn/quantum_decaf/tests/CMakeLists.txt;8;add_logging_test;/home/artur/inzyn/quantum_decaf/tests/CMakeLists.txt;0;")
subdirs("modules")
