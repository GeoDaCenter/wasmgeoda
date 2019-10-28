cd build
cmake .. -DUSE_EMCC=True -DCMAKE_C_COMPILER=emcc -DCMAKE_CXX_COMPILER=em++ -DCMAKE_CPP_COMPILER=em++ -DCMAKE_CXX_STANDARD=14
make jsgeoda.html

mv jsgeoda.js ../web/gen/
mv jsgeoda.wasm ../web/gen/