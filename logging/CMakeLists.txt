project (logging)

file(GLOB logging_SRC
    "./logging/*.*"
)

add_library(logging ${logging_SRC})