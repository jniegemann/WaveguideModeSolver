project(WaveguideModeSolver)

  add_executable(ModeSolver WaveguideModeSolverOverlap.cpp ClockTickCounter.cpp)
  target_link_libraries(ModeSolver ${HERMES_COMMON_LIBRARY} ${HERMES_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
