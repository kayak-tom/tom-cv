project (featureMatching)

file(GLOB featureMatching_SRC
    "./featureMatching/*.*"
)

add_library(featureMatching ${featureMatching_SRC})