add_library(truss_fem_lib STATIC
  Types.h
  Matrix.h
  Node.h
  Element.h Element.cpp
  BoundaryCondition.h
  NodalLoad.h
  TrussModel.h TrussModel.cpp
  Assembler.h Assembler.cpp
  BoundaryConditionHandler.h BoundaryConditionHandler.cpp
  Solver.h
  GaussianSolver.h GaussianSolver.cpp
  JsonParser.h JsonParser.cpp
  PostProcessor.h PostProcessor.cpp
)
target_include_directories(truss_fem_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(truss_fem_lib PUBLIC nlohmann_json::nlohmann_json)

add_executable(truss_solver main.cpp)
target_link_libraries(truss_solver PRIVATE truss_fem_lib)
