add_executable(cell_ordering
    cell_ordering.cpp)
target_link_libraries(cell_ordering
    OpenFOAM
    finiteVolume
    Pstream
    )
