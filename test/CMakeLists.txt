include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY  https://github.com/catchorg/Catch2.git
        GIT_TAG         v3.0.0-preview4
)

FetchContent_MakeAvailable(Catch2)

function(jaro_winkler_add_test NAME SOURCE)
	add_executable(test_${NAME} ${SOURCE})
	target_link_libraries(test_${NAME} ${PROJECT_NAME})
	target_link_libraries(test_${NAME} Catch2::Catch2WithMain)
	add_test(NAME ${NAME} COMMAND test_${NAME})
endfunction()

jaro_winkler_add_test(jaro-winkler tests-jaro-winkler.cpp)
