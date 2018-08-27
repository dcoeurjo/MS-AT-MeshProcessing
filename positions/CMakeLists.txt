
set(target_base "positions")

add_executable(at-${target_base}-mesh
    options.cpp
    main.cpp
    Viewer.cpp
    )

target_link_libraries(at-${target_base}-mesh
    dec
    ${Boost_LIBRARIES}
    )

set_target_properties(at-${target_base}-mesh
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
