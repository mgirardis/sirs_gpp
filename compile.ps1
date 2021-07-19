param([String]$inputDir,[String]$outputFile,[String]$extraParam="")

#g++ $inputDir\*.h $inputDir\*.cpp -o $outputFile -std=c++11 -Ofast -O3 -m64 -Wconversion -Wall $extraParam
g++ $inputDir\*.h $inputDir\*.cpp -o $outputFile -std=c++11 -Ofast -O3 -m64 $extraParam
