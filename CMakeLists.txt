cmake_minimum_required(VERSION 2.6)
file(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}/CMake" p)
list(APPEND CMAKE_MODULE_PATH "${p}") #For old Imagine++
list(APPEND CMAKE_SYSTEM_FRAMEWORK_PATH /Library/Frameworks) #Mac, why not auto?
find_package(Imagine REQUIRED)

project(MOPSI)
add_executable(MOPSI main.cpp lucas_kanade.cpp lucas_kanade.h horn_schunck.cpp horn_schunck.h functions.cpp functions.h videos.cpp videos.h)
ImagineUseModules(MOPSI Graphics LinAlg Images)
