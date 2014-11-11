project (BoWSLAM)

file(GLOB BoWSLAM_SRC
    "./BoWSLAM_Main/*.cpp"
)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ../cmake)
message(STATUS "CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")

find_package(BOARD 0.9)

if(BOARD_FOUND)
  include_directories(BOARD_INCLUDE_DIR)
  add_executable(BoWSLAM ${BoWSLAM_SRC})
  target_link_libraries(BoWSLAM board)
else()
  message(STATUS "NOT BUILDING BoWSLAM!! libboard must be installed to build BoWSLAM (sudo ln -s /home/XX/lib/board/include/Board.h /usr/include/Board.h; sudo ln -s /home/XX/lib/board/lib/libboard.so /usr/lib/libboard.so; sudo ln -s /home/XX/lib/board/include/board/ /usr/include/board)")
endif()