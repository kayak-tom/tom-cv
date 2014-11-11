project (image)

file(GLOB image_SRC
    "./image/*.*"
)

add_library(image ${image_SRC})