set(project "jsgeoda")

project(${project})
cmake_minimum_required(VERSION 3.12)

set (DEP_LIBRARY_PATH ${CMAKE_SOURCE_DIR}/3rd_party)

message(STATUS "dep library path: ${DEP_LIBRARY_PATH}")
set (RAPIDJSON_INCLUDE_DIR ${DEP_LIBRARY_PATH}/rapidjson-1.1.0/include)
set (BOOST_INCLUDE_DIR ${DEP_LIBRARY_PATH}/boost-1.70.0)

set(src
		src/jsgeoda.cpp
		src/geojson.cpp
		src/weights/GeodaWeight.cpp
		src/weights/GalWeight.cpp
		src/weights/GwtWeight.cpp
		src/weights/PolysToContigWeights.cpp
		src/weights/VoronoiUtils.cpp
		src/gda_weights.cpp
		src/gda_clustering.cpp
		src/gda_sa.cpp
		src/GenGeomAlgs.cpp
		src/SpatialIndAlgs.cpp
		src/GenUtils.cpp
		src/sa/LISA.cpp
		src/sa/UniG.cpp
		src/sa/UniGstar.cpp
		src/sa/UniGeary.cpp
		src/sa/UniJoinCount.cpp
		src/sa/UniLocalMoran.cpp
		src/clustering/cluster.cpp
		src/clustering/maxp.cpp
		src/clustering/maxp_wrapper.cpp
		src/clustering/redcap.cpp
		src/clustering/redcap_wrapper.cpp
)

set(exports 
	  _print_json
)


# process exported functions
set(exports_string "")
list(JOIN exports "," exports_string)

# set compiler and flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TARGET_NAME ${project})
if (USE_EMCC)
	set(TARGET_NAME ${project}.html)
	message("MY_VARIABLE evaluates to ${TARGET_NAME}")
	SET(CMAKE_C_COMPILER /Users/xunli/Downloads/emsdk/fastcomp/emscripten/emcc)
	SET(CMAKE_CPP_COMPILER /Users/xunli/Downloads/emsdk/fastcomp/emscripten/em++)
	SET(CMAKE_CXX_COMPILER /Users/xunli/Downloads/emsdk/fastcomp/emscripten/em++)
	set( CMAKE_CXX_FLAGS "-s DISABLE_EXCEPTION_CATCHING=0 -s TOTAL_MEMORY=181207040 -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=8 -s USE_BOOST_HEADERS=1 -s BINARYEN_TRAP_MODE='clamp' --bind -s EXPORTED_FUNCTIONS=\"[${exports_string}]\" -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'"  )
endif()

# specify the project
add_executable(${TARGET_NAME} ${src})

if (USE_EMCC)
	target_compile_definitions(${TARGET_NAME} PRIVATE __JSGEODA__=1)
	target_compile_definitions(${TARGET_NAME} PRIVATE EMCC_DEBUG=1)
else()
	target_include_directories(${TARGET_NAME} PUBLIC ${BOOST_INCLUDE_DIR} )
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
		test/test_clustering.cpp
		test/test_localSA.cpp
)
add_executable( runUnitTests ${TEST_SOURCE_FILES} )

target_include_directories(runUnitTests PUBLIC
		${gtest_SOURCE_DIR}/include
		${gmock_SOURCE_DIR}/include
		${RAPIDJSON_INCLUDE_DIR} )

if (USE_EMCC)
else()
	target_include_directories(runUnitTests PUBLIC ${BOOST_INCLUDE_DIR} )
endif()
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main)