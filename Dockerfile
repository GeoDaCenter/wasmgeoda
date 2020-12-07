FROM trzeci/emscripten:1.39.17-upstream

WORKDIR /tmp/wasmgeoda

COPY . .

WORKDIR /tmp/wasmgeoda/build

RUN cmake .. -DUSE_EMCC=True -DCMAKE_C_COMPILER=emcc -DCMAKE_CXX_COMPILER=em++ -DCMAKE_CPP_COMPILER=em++
RUN make jsgeoda.html