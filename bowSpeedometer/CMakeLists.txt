project (bowSpeedometer)

file(GLOB scale_SRC
    "./scale/*.*"
)

add_library(bowSpeedometer ${scale_SRC})