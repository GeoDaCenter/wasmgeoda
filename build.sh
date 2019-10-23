rm build/ -rf
mkdir build
cd build
cmake ../src
make

mv jsgeoda.js ../web/gen/
mv jsgeoda.wasm ../web/gen/