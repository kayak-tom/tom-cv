project (timing)

file(GLOB time_SRC
    "./time/*.*"
)

add_library(timing ${time_SRC})