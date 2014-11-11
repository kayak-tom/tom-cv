project (mosaicing)

file(GLOB mosaicing_SRC
    "./MosaicingSoftware/mosaicing/*.cpp"
)

add_executable(mosaicing ${mosaicing_SRC})