project (featureExtract)

file(GLOB featureExtract_SRC
    "./featureExtract/*.*"
    "./featureExtract/*/*.*"
)

add_library(featureExtract ${featureExtract_SRC})