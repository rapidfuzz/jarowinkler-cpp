find_package(Catch2 3 QUIET)
if (Catch2_FOUND)
    message("Using system supplied version of Catch2")
else()
    message("Using FetchContent to load Catch2")
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY  https://github.com/catchorg/Catch2.git
        GIT_TAG         v3.0.1
    )
    FetchContent_MakeAvailable(Catch2)
endif()

function(jaro_winkler_add_test NAME SOURCE)
	add_executable(test_${NAME} ${SOURCE})
	target_link_libraries(test_${NAME} ${PROJECT_NAME})
	target_link_libraries(test_${NAME} Catch2::Catch2WithMain)
	add_test(NAME ${NAME} COMMAND test_${NAME})
endfunction()

jaro_winkler_add_test(jaro-winkler tests-jaro-winkler.cpp)
