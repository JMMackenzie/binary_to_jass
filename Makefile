all:
	g++ --std=c++14 -mavx2 -DCI_FORCEINLINE -O3 -Wall bin_to_jass.cpp compress_qmx.cpp compress_integer_elias_gamma_simd.cpp -o bin_to_jass
	g++ --std=c++14 -mavx2 -DCI_FORCEINLINE -O3 -Wall bin_to_jass_disk.cpp compress_qmx.cpp compress_integer_elias_gamma_simd.cpp -o bin_to_jass_lowmemory

