# Copyright (C) 2025 Andrew D Smith
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

find_package(GSL REQUIRED)
find_package(HTSLIB REQUIRED)

file(GLOB cpp_files "*.cpp")

# Gather all the object files that will be put in the static library
# and prepare to compile them.
set(LIBRARY_OBJECTS "")
foreach(cpp_file ${cpp_files})
  get_filename_component(BASE_NAME ${cpp_file} NAME_WE)
  add_library(${BASE_NAME} OBJECT ${cpp_file})
  target_link_libraries(${BASE_NAME} PUBLIC bamxx smithlab_cpp)
  target_link_libraries(${BASE_NAME} PRIVATE GSL::gsl HTSLIB::HTSLIB)
  target_include_directories(${BASE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
