find_package(Catch2 3 REQUIRED)

function(jaro_winkler_add_test NAME SOURCE)
	add_executable(test_${NAME} ${SOURCE})
	target_link_libraries(test_${NAME} ${PROJECT_NAME})
	target_link_libraries(test_${NAME} Catch2::Catch2WithMain )
	add_test(NAME ${NAME} COMMAND test_${NAME})
endfunction()

jaro_winkler_add_test(jaro-winkler tests-jaro-winkler.cpp)
