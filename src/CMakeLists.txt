set(project "jsgeoda")

project(${project})
cmake_minimum_required(VERSION 3.12)

set (DEP_LIBRARY_PATH ${CMAKE_SOURCE_DIR}/../3rd_party)

message(STATUS "dep library path: ${DEP_LIBRARY_PATH}")
set (RAPIDJSON_INCLUDE_DIR ${DEP_LIBRARY_PATH}/rapidjson-1.1.0/include)

set(src 
	  jsgeoda.cpp
	  geojson.cpp
	  weights/GeodaWeight.cpp
	  weights/GalWeight.cpp
	  weights/GwtWeight.cpp
	  weights/PolysToContigWeights.cpp
	  gda_weights.cpp
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

SET(CMAKE_C_COMPILER emcc)
SET(CMAKE_CPP_COMPILER em++)
SET(CMAKE_CXX_COMPILER em++)
set( CMAKE_CXX_FLAGS "-s BINARYEN_TRAP_MODE='clamp' -s ASSERTIONS=1 -s ALLOW_MEMORY_GROWTH=1 -s ERROR_ON_UNDEFINED_SYMBOLS=0 --bind -s EXPORTED_FUNCTIONS=\"[${exports_string}]\" -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'"  )

# specify the project
add_executable(${project}.html ${src})

target_compile_definitions(${project}.html PRIVATE __JSGEODA__=1)

target_include_directories(${project}.html PUBLIC ${GEOJSON_INCLUDE_DIR} )
target_include_directories(${project}.html PUBLIC ${GEOMETRYHPP_INCLUDE_DIR} )
target_include_directories(${project}.html PUBLIC ${VARIANT_INCLUDE_DIR} )
target_include_directories(${project}.html PUBLIC ${RAPIDJSON_INCLUDE_DIR} )
