emcc -O2 -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=8 -o test_pthreads.js test_pthreads.c
em++ -O2 -s USE_PTHREADS=1 -std=c++11 -s -o test_multithreads.js test_multithreads.cpp
