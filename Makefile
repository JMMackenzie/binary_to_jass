all:
	g++ --std=c++14 -msse4 -DCI_FORCEINLINE -O3 -Wall bin_to_jass.cpp compress_qmx.c -o bin_to_jass
	g++ --std=c++14 -msse4 -DCI_FORCEINLINE -O3 -Wall bin_to_jass_disk.cpp compress_qmx.c -o bin_to_jass_lowmemory

