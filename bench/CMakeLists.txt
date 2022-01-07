find_package(benchmark REQUIRED)

function(jaro_winkler_add_benchmark NAME SOURCE)
	add_executable(bench_${NAME} ${SOURCE})
	target_link_libraries(bench_${NAME} ${PROJECT_NAME})
	target_link_libraries(bench_${NAME} benchmark::benchmark)
endfunction()

# todo
#jaro_winkler_add_benchmark(fuzz bench-fuzz.cpp)
