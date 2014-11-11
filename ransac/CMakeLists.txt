project (ransac)

file(GLOB ransac_SRC
    "./ransac/*.*"
)

add_library(ransac ${ransac_SRC})