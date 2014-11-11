LICENSE: 

Essential matrix estimation code (in 'ransac' library) is based on Henrik Stewenius's MATLAB code (http://vis.uky.edu/~stewe/FIVEPOINT/) and is for academic use only. Please see the notice in the file calibrated_fivepoint_helper.c.

All other code is free to use/modify/redistribute for any purpose. Please send me an email if you find it useful, or if anything doesn't work, or when you find bugs! tom.botterill - at - grcnz.com. 


PRE-REQUISITES:

Most components require:
1) boost (I use 1.39.0, works with other versions, libraries need to be built)
2) openCV (latest version 2.0 recomended, needed for homography estimation, N-point Essential matrix estimation (N>=8) and anything using images)
3) Eigen (a fast header-only matrix library, much faster if the convergence threshhold for eigenvector computation is raised and the complex eigenvector code is removed)


BUILD:

Project files for the following IDEs are included. You may have to fix some absolute paths, especially under Linux.

WINDOWS:
Visual Studio 2008 (open util.sln)
Visual Studio 2005 (MODIFY all the .sln and .vcproj files to change Version="9.00" to Version="7.00")

LINUX:
1) Make 3 symlinks in the workspace directory: ln -s PATH_TO_YOUR_OPENCV_DIR opencv; ln -s PATH_TO_YOUR_EIGEN_DIR opencv; ln -s PATH_TO_YOUR_BOOST_DIR boost

Eclipse CDT (Galileo)

The makefiles supplied are just the ones generated from Eclipse but should work. E.g. "cd util/Debug; make" will build libutil.a with debugging info, etc.


DOCUMENTATION:

Online at www.hilandtom.com/tombotterill/code


REFERENCING:

My publications (with links and bibtex) are listed at www.hilandtom.com/tombotterill

My Bag-of-Words scheme is described in: "Speeded-up Bag-of-Words algorithm for robot localisation through scene recognition" and used in "A Bag-of-Words Speedometer for Single Camera SLAM" and "An Integrated IMU, GNSS and Image Recognition Sensor for Pedestrian Navigation"

BaySAC is described in "New Conditional Sampling Strategies for Speeded-Up RANSAC"

Another paper on Visual SLAM (BoWSLAM) that uses all of this code is waiting for review; a video is online at www.hilandtom.com/tombotterill
