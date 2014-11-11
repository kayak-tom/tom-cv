project (bagOfWords)

file(GLOB bow_SRC
    "./bow/*.*"
)

add_library(bagOfWords ${bow_SRC})