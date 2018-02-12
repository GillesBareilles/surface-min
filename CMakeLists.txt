cmake_minimum_required(VERSION 2.6)



#Configuration du projet

project(MyProject)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

set(CMAKE_BUILD_TYPE Debug)

add_definitions(-std=c++11 -Wall -g -O3 -pedantic) #pedantic...

#Configuration de l'exécutable

file(

        GLOB_RECURSE

        source_files

        src/*

)

add_executable(

        my_exe

        ${source_files}

)
