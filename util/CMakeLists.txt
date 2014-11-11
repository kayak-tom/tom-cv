project (util)

file(GLOB util_SRC
    "./util/*.*"
)

add_library(util ${util_SRC})