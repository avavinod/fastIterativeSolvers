file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS
    matrixUtils/*.cpp
    CSRMatrix/*.cpp
    MSRMatrix/*.cpp
    gmresUtils/*.cpp
)

add_library(CSRMatrix STATIC ${LIB_SOURCES})
add_library(MSRMatrix STATIC ${LIB_SOURCES})
add_library(matrixUtils STATIC ${LIB_SOURCES})
add_library(gmresUtils STATIC ${LIB_SOURCES})

target_include_directories(CSRMatrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(MSRMatrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(matrixUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(gmresUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
