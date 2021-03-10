# shared library をつくる
add_library(scf
    SHARED
        test.cpp
        molecule.cpp
        iwl.cpp
        basisset.cpp
        hf.cpp
)

# libraryのincludeファイル用

target_include_directories(scf
    PUBLIC 
        ${PROJECT_SOURCE_DIR}
        "/usr/local/libint/2.6.0/include"
        "/usr/local/include/eigen3"
)

target_link_libraries(scf 
    PUBLIC "/usr/local/libint/2.6.0/lib/libint2.a"
)
