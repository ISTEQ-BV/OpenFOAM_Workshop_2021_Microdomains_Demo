add_executable(of_transport
    of_transport.cpp
    lsgrad.hpp lsgrad.cpp
    microdomain.hpp microdomain.cpp
    meshdata.hpp meshdata.cpp
    )
target_link_libraries(of_transport
    fvOptions
    OpenFOAM
    finiteVolume
    Pstream
    meshTools
    )

add_executable(test_lsgrad
    test_lsgrad.cpp
    lsgrad.hpp lsgrad.cpp
    meshdata.hpp meshdata.cpp
    )
target_link_libraries(test_lsgrad PUBLIC
    Pstream
    OpenFOAM
    finiteVolume
    meshTools
    )
