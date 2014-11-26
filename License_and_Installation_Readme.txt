LICENSE: 

Essential matrix estimation code (in 'ransac' library) is based on Henrik Stewenius's MATLAB code (http://vis.uky.edu/~stewe/FIVEPOINT/) and is for academic use only. Please see the notice in the file calibrated_fivepoint_helper.c.

All other code is free to use/modify/redistribute for any purpose. Please send me an email if you find it useful, or if anything doesn't work, or when you find bugs! tom.botterill - at - grcnz.com. 


PRE-REQUISITES:

1) boost (I use 1.54.0, works with other versions, libraries need to be built)
2) openCV (latest version 2.4.9 recomended)
3) Eigen 3.1 or later
4) cmake
(5) for BoWSLAM, 'board' library is required)

BUILD IN LINUX:

Use cmake:

cd to the tom-cv directory.

mkdir build
cd build
cmake ..
make -j8

(This will build everything, and will put the executables under 'build', and libraries in RelWithDebInfo (default cmake build type).

BUILD IN WINDOWS: (sorry it is complicated! I don't know any better way to make cmake work in Windows)

1) Download+install "Github for Windows" and "Cmake for Windows"

2) Use Github for Windows to clone the tom-cv repo
https://github.com/kayak-tom/tom-cv

3) **Edit the top few lines of tom-cv/CMakeLists.txt to set the location of opencv, Eigen, Boost.**

4) Run cmake on the tom-cv directory (click "Configure" then "Generate")

5) Visual studio project files should be created. [You may need to add OpenCV (everything) and Boost libraries (system, filesystem, boost) to the project settings]



DOCUMENTATION:

Online at www.hilandtom.com/tombotterill/code


REFERENCING:

My publications (with links and bibtex) are listed at www.hilandtom.com/tombotterill

My Bag-of-Words scheme is described in: "Speeded-up Bag-of-Words algorithm for robot localisation through scene recognition" and used in "A Bag-of-Words Speedometer for Single Camera SLAM" and "An Integrated IMU, GNSS and Image Recognition Sensor for Pedestrian Navigation"

BaySAC is described in "New Conditional Sampling Strategies for Speeded-Up RANSAC"

Mosaicing is described in "Real-time aerial image mosaicing"
