set(project "jsgeoda")

cmake_minimum_required(VERSION 3.12)
project(${project} VERSION "0.0.6")

# process exported functions
set(exports _print_json)
set(exports_string "")
list(JOIN exports "," exports_string)

# set compiler and flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TARGET_NAME ${project})
message("MY_VARIABLE evaluates to ${TARGET_NAME}")

# emcc version 1.39.0  (commit e047fe4c1ecfae6ba471ca43f2f630b79516706b)
SET(CMAKE_C_COMPILER emcc)
SET(CMAKE_CPP_COMPILER em++)
SET(CMAKE_CXX_COMPILER em++)

#set( CMAKE_CXX_FLAGS "-O3 -s TOTAL_MEMORY=1812070400 -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=16 -s USE_BOOST_HEADERS=1 -s BINARYEN_TRAP_MODE='clamp' --bind -s EXPORTED_FUNCTIONS=\"[${exports_string}]\" -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'"  )
#set( CMAKE_CXX_FLAGS "-s DISABLE_EXCEPTION_CATCHING=0 -s TOTAL_MEMORY=1812070400 -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=8 -s USE_BOOST_HEADERS=1 -s BINARYEN_TRAP_MODE='clamp' --bind -s EXPORTED_FUNCTIONS=\"[${exports_string}]\" -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'"  )
#set( CMAKE_CXX_FLAGS "-O2 -s ALLOW_MEMORY_GROWTH=1 -s TOTAL_MEMORY=200MB -s USE_PTHREADS=0 -s USE_BOOST_HEADERS=1 -s BINARYEN_TRAP_MODE='clamp' --bind -s EXPORTED_FUNCTIONS=\"[${exports_string}]\" -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\", \"cwrap\"]'"  )
set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG --bind -s ASSERTIONS=0 -s WASM=1 -s ALLOW_MEMORY_GROWTH=1")
#set(CMAKE_CXX_FLAGS "-g --bind -s ASSERTIONS=1 -s WASM=1 -s ALLOW_MEMORY_GROWTH=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s NO_FILESYSTEM=1 -s USE_BOOST_HEADERS=1")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s FILESYSTEM=1 -s FORCE_FILESYSTEM=1")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D MEMFS")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=6")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"ccall\"]'")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s BINARYEN_TRAP_MODE='clamp'")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s EXPORTED_FUNCTIONS=\"[${exports_string}]\"")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s MODULARIZE=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --pre-js ${CMAKE_SOURCE_DIR}/src/pre.js")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s EXPORT_NAME='jsgeoda'")

set(src
		src/jsgeoda.cpp
		src/jsgeoda_breaks.cpp
		src/jsgeoda_weights.cpp
		src/geojson.cpp
		libgeoda_src/knn/ANN.cpp
		libgeoda_src/knn/kd_tree.cpp
		libgeoda_src/knn/kd_split.cpp
		libgeoda_src/knn/kd_util.cpp
		libgeoda_src/knn/kd_search.cpp
		libgeoda_src/knn/kd_dump.cpp
		libgeoda_src/knn/kd_fix_rad_search.cpp
		libgeoda_src/knn/kd_pr_search.cpp
		libgeoda_src/weights/GeodaWeight.cpp
		libgeoda_src/weights/GalWeight.cpp
		libgeoda_src/weights/GwtWeight.cpp
		libgeoda_src/weights/PolysToContigWeights.cpp
		libgeoda_src/weights/VoronoiUtils.cpp
		libgeoda_src/gda_weights.cpp
		libgeoda_src/gda_clustering.cpp
		libgeoda_src/gda_data.cpp
		libgeoda_src/gda_sa.cpp
		libgeoda_src/GenGeomAlgs.cpp
		libgeoda_src/SpatialIndAlgs.cpp
		libgeoda_src/GenUtils.cpp
		libgeoda_src/sa/LISA.cpp
		libgeoda_src/sa/UniG.cpp
		libgeoda_src/sa/UniGstar.cpp
		libgeoda_src/sa/UniGeary.cpp
		libgeoda_src/sa/UniJoinCount.cpp
		libgeoda_src/sa/UniLocalMoran.cpp
		libgeoda_src/sa/MultiJoinCount.cpp
		libgeoda_src/sa/MultiGeary.cpp
		libgeoda_src/clustering/cluster.cpp
		libgeoda_src/clustering/maxp_wrapper.cpp
		libgeoda_src/clustering/redcap_wrapper.cpp
		libgeoda_src/clustering/azp_wrapper.cpp
		libgeoda_src/clustering/schc_wrapper.cpp
		libgeoda_src/clustering/redcap.cpp
		libgeoda_src/clustering/azp.cpp
		libgeoda_src/clustering/DorlingCartogram.cpp
		src/jsgeoda_lisa.cpp
		src/jsgeoda_clustering.cpp )



# specify the project
add_executable(${TARGET_NAME} ${src})

target_compile_definitions(${TARGET_NAME} PRIVATE __JSGEODA__=1)
target_compile_definitions(${TARGET_NAME} PRIVATE __NO_THREAD__=1)
target_compile_definitions(${TARGET_NAME} PRIVATE EMCC_DEBUG=0)

target_include_directories(${TARGET_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/3rd_party/rapidjson-1.1.0/include)
