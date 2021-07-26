set(Sources
  hexcli.cpp
)

add_executable(hexcli ${Sources})

target_link_libraries(hexcli cubesim)
