cd build
emcmake cmake .. -DUSE_EMCC=True -DCMAKE_C_COMPILER=emcc -DCMAKE_CXX_COMPILER=em++ -DCMAKE_CPP_COMPILER=em++
emmake make jsgeoda.html

