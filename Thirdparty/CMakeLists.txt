cmake_minimum_required(VERSION 3.10)

message(CHECK_START "Configuring third party libraries")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Add third party libraries
add_subdirectory(RcdMathLib)

message(CHECK_PASS "done")
list(POP_BACK CMAKE_MESSAGE_INDENT)