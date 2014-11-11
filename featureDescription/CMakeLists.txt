project (featureDescription)

file(GLOB description_SRC
    "./description/*.*"
)

add_library(featureDescription ${description_SRC})