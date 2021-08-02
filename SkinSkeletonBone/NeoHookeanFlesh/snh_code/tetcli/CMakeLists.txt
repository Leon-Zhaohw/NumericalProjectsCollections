set(Sources
  tetcli.cpp
)

add_executable(tetcli ${Sources})

target_link_libraries(tetcli cubesim)
