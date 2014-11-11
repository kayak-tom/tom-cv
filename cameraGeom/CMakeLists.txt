project (cameraGeom)

file(GLOB geom_SRC
    "./geom/*.*"
)

add_library(cameraGeom ${geom_SRC})