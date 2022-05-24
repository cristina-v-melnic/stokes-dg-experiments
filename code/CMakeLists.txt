cmake_minimum_required(VERSION 3.9)

find_package(parmoon)
add_executable(nse2d "${CMAKE_CURRENT_SOURCE_DIR}/NSE2D_ParMooN.C")
target_link_libraries(nse2d PUBLIC parmoon::parmoon_2d_SEQUENTIAL)
