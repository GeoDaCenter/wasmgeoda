set(project "jsgeoda")

project(${project})
cmake_minimum_required(VERSION 3.12)

set (DEP_LIBRARY_PATH ${CMAKE_SOURCE_DIR}/3rd_party)

message(STATUS "dep library path: ${DEP_LIBRARY_PATH}")
set (RAPIDJSON_INCLUDE_DIR ${DEP_LIBRARY_PATH}/rapidjson-1.1.0/include)

set(src
		src/jsgeoda.cpp
		src/geojson.cpp
		src/weights/GeodaWeight.cpp
		src/weights/GalWeight.cpp
		src/weights/GwtWeight.cpp
		src/weights/PolysToContigWeights.cpp
		src/weights/PointsToContigWeights.cpp
		src/gda_weights.cpp
)

set(exports 
	  _print_json
	  _get_num_obs
	  _get_map_type
	  _queen_weights
)


# process exported functions
set(exports_string "")
list(JOIN exports "," exports_string)

# set compiler and flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TARGET_NAME ${project})
if (USE_EMCC)
	message("MY_VARIABLE evaluates to True")
	set(EXE_NAME ${project}.html)
	SET(CMAKE_C_COMPILER emcc)
	SET(CMAKE_CPP_COMPILER em++)
	SET(CMAKE_CXX_COMPILER em++)
	set( CMAKE_CXX_FLAGS "-s BINARYEN_TRAP_MODE='clamp' -s ASSERTIONS=1 -s ALLOW_MEMORY_GROWTH=1 -s ERROR_ON_UNDEFINED_SYMBOLS=0 --bind -s EXPORTED_FUNCTIONS=\"[${exports_string}]\" -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'"  )
endif()

# specify the project
add_executable(${TARGET_NAME} ${src})

if (USE_EMCC)
	target_compile_definitions(${TARGET_NAME} PRIVATE __JSGEODA__=1)
endif()

target_include_directories(${TARGET_NAME} PUBLIC ${RAPIDJSON_INCLUDE_DIR} )

###########################################################
#
#  googletest testing framework
#
###########################################################
# git clone https://github.com/google/googletest.git
ADD_SUBDIRECTORY (googletest)
enable_testing()

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

################################
# Unit Tests
################################
# Add test cpp file
set (TEST_SOURCE_FILES ${src}
		test/test_weights.cpp
)
add_executable( runUnitTests ${TEST_SOURCE_FILES} )

target_include_directories(runUnitTests PUBLIC
		${gtest_SOURCE_DIR}/include
		${gmock_SOURCE_DIR}/include
		${RAPIDJSON_INCLUDE_DIR} )

# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main)