project (params)

file(GLOB params_SRC
    "./params/*.*"
)

add_library(params ${params_SRC})