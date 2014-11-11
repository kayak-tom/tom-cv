project (imageSource)

file(GLOB imageSource_SRC
    "./imageSource/*.*"
)

add_library(imageSource ${imageSource_SRC})