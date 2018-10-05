cmake_minimum_required(VERSION 3.5)
project(filter)

#execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOTINCS)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOTINCS} -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -m64 -I/home/jsvirzi/root/include -std=c++11")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -m64 -I/home/jsvirzi/root/include -std=c++11")

# manually throw the ROOT cflags in
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -m64")

include_directories($ENV{ROOTSYS}/include)
include_directories($ENV{HOME}/projects/ComputerVision/scv/include)
include_directories($ENV{HOME}/projects/ComputerVision/LaneDetection/include)

#execute_process(COMMAND pkg-config --libs opencv OUTPUT_VARIABLE OPENCVLIBS)
#execute_process(COMMAND pkg-config --cflags opencv OUTPUT_VARIABLE OPENCVINCS)

include_directories(/usr/local/include)

link_libraries(-L$ENV{ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic)

set(SOURCE_FILES src/filter_response.cpp)
add_executable(filter ${SOURCE_FILES})