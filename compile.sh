#g++ $1/*.h $1/*.cpp -o $2 -std=c++11 -Ofast -O3 -m64 -Wconversion -Wall
#g++ $1/*.h $1/*.cpp -o $2 -O3 -m64 -fpermissive # -Wconversion -Wall -Wmisleading-indentation -Wunused-variable -Wint-in-bool-context -Wconversion -Wdelete-non-virtual-dtor -Wdeprecated-declarations -Wsign-compare
icpc $1/*.h $1/*.cpp -o $2 -std=c++14 -Ofast -no-ansi-alias -m64
