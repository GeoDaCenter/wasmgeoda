rm build/ -rf
mkdir build
cd build
cmake ../src -DUSE_EMCC=True
make

mv jsgeoda.js ../web/gen/
mv jsgeoda.wasm ../web/gen/