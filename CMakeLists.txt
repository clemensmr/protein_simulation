include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}/Include"
  "${CMAKE_CURRENT_SOURCE_DIR}")


# add executable targets

add_executable(prot
                            prot.cpp
                            protdata.cpp
                            protdata.hpp)
#target_link_libraries(readdx ahmed)