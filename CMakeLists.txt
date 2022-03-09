cmake_minimum_required(VERSION 3.19)
project(graphbolt)

set(CMAKE_CXX_STANDARD 14)

add_definitions(-DEDGELONG)

include_directories(core)
include_directories(core/common)
include_directories(core/graph)
include_directories(core/graphBolt)
include_directories(tools/common)
include_directories(apps)

add_executable(SNAPtoAdjConverter tools/converters/SNAPtoAdjConverter.C)
add_executable(PageRank apps/PageRank.C)
add_executable(SSSP apps/SSSP.C)
add_executable(LabelPropagation apps/LabelPropagation.C)
add_executable(BFS apps/BFS.C)
add_executable(CF apps/CF.C)
add_executable(COEM apps/COEM.C)
