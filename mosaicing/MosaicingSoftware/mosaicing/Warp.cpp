/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

/*!  \file Warp.cpp \brief This file is derived from cvWarp.cpp in OpenCV. It contains a much faster warp implementation.
*/

#include "util/exception.h"
#include "Warp.h"
//#include "util/convert.h"

//#define WIN32
//#include "../src/_cv.h"
//#include "/home/tom/workspace/opencv/src/cv/_cv.h"
//#include "GRCException.h"
#include <iostream>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include "util/opencv.h"

#ifndef USE_OLD_OPENCV
#include "opencv2/opencv.hpp"
//#include "opencv2/core/internal.hpp" //deprecated now... TODO fix me
//#include "opencv2/core/private.hpp"
#else
#define __BEGIN__       {
#define __END__         goto exit; exit: ; }
#define __CLEANUP__
#define EXIT            goto exit
#endif

// Constants for double too:
const double icv8x64dTab_cv[] =
{
    -256.0, -255.0, -254.0, -253.0, -252.0, -251.0, -250.0, -249.0,
    -248.0, -247.0, -246.0, -245.0, -244.0, -243.0, -242.0, -241.0,
    -240.0, -239.0, -238.0, -237.0, -236.0, -235.0, -234.0, -233.0,
    -232.0, -231.0, -230.0, -229.0, -228.0, -227.0, -226.0, -225.0,
    -224.0, -223.0, -222.0, -221.0, -220.0, -219.0, -218.0, -217.0,
    -216.0, -215.0, -214.0, -213.0, -212.0, -211.0, -210.0, -209.0,
    -208.0, -207.0, -206.0, -205.0, -204.0, -203.0, -202.0, -201.0,
    -200.0, -199.0, -198.0, -197.0, -196.0, -195.0, -194.0, -193.0,
    -192.0, -191.0, -190.0, -189.0, -188.0, -187.0, -186.0, -185.0,
    -184.0, -183.0, -182.0, -181.0, -180.0, -179.0, -178.0, -177.0,
    -176.0, -175.0, -174.0, -173.0, -172.0, -171.0, -170.0, -169.0,
    -168.0, -167.0, -166.0, -165.0, -164.0, -163.0, -162.0, -161.0,
    -160.0, -159.0, -158.0, -157.0, -156.0, -155.0, -154.0, -153.0,
    -152.0, -151.0, -150.0, -149.0, -148.0, -147.0, -146.0, -145.0,
    -144.0, -143.0, -142.0, -141.0, -140.0, -139.0, -138.0, -137.0,
    -136.0, -135.0, -134.0, -133.0, -132.0, -131.0, -130.0, -129.0,
    -128.0, -127.0, -126.0, -125.0, -124.0, -123.0, -122.0, -121.0,
    -120.0, -119.0, -118.0, -117.0, -116.0, -115.0, -114.0, -113.0,
    -112.0, -111.0, -110.0, -109.0, -108.0, -107.0, -106.0, -105.0,
    -104.0, -103.0, -102.0, -101.0, -100.0,  -99.0,  -98.0,  -97.0,
     -96.0,  -95.0,  -94.0,  -93.0,  -92.0,  -91.0,  -90.0,  -89.0,
     -88.0,  -87.0,  -86.0,  -85.0,  -84.0,  -83.0,  -82.0,  -81.0,
     -80.0,  -79.0,  -78.0,  -77.0,  -76.0,  -75.0,  -74.0,  -73.0,
     -72.0,  -71.0,  -70.0,  -69.0,  -68.0,  -67.0,  -66.0,  -65.0,
     -64.0,  -63.0,  -62.0,  -61.0,  -60.0,  -59.0,  -58.0,  -57.0,
     -56.0,  -55.0,  -54.0,  -53.0,  -52.0,  -51.0,  -50.0,  -49.0,
     -48.0,  -47.0,  -46.0,  -45.0,  -44.0,  -43.0,  -42.0,  -41.0,
     -40.0,  -39.0,  -38.0,  -37.0,  -36.0,  -35.0,  -34.0,  -33.0,
     -32.0,  -31.0,  -30.0,  -29.0,  -28.0,  -27.0,  -26.0,  -25.0,
     -24.0,  -23.0,  -22.0,  -21.0,  -20.0,  -19.0,  -18.0,  -17.0,
     -16.0,  -15.0,  -14.0,  -13.0,  -12.0,  -11.0,  -10.0,   -9.0,
      -8.0,   -7.0,   -6.0,   -5.0,   -4.0,   -3.0,   -2.0,   -1.0,
       0.0,    1.0,    2.0,    3.0,    4.0,    5.0,    6.0,    7.0,
       8.0,    9.0,   10.0,   11.0,   12.0,   13.0,   14.0,   15.0,
      16.0,   17.0,   18.0,   19.0,   20.0,   21.0,   22.0,   23.0,
      24.0,   25.0,   26.0,   27.0,   28.0,   29.0,   30.0,   31.0,
      32.0,   33.0,   34.0,   35.0,   36.0,   37.0,   38.0,   39.0,
      40.0,   41.0,   42.0,   43.0,   44.0,   45.0,   46.0,   47.0,
      48.0,   49.0,   50.0,   51.0,   52.0,   53.0,   54.0,   55.0,
      56.0,   57.0,   58.0,   59.0,   60.0,   61.0,   62.0,   63.0,
      64.0,   65.0,   66.0,   67.0,   68.0,   69.0,   70.0,   71.0,
      72.0,   73.0,   74.0,   75.0,   76.0,   77.0,   78.0,   79.0,
      80.0,   81.0,   82.0,   83.0,   84.0,   85.0,   86.0,   87.0,
      88.0,   89.0,   90.0,   91.0,   92.0,   93.0,   94.0,   95.0,
      96.0,   97.0,   98.0,   99.0,  100.0,  101.0,  102.0,  103.0,
     104.0,  105.0,  106.0,  107.0,  108.0,  109.0,  110.0,  111.0,
     112.0,  113.0,  114.0,  115.0,  116.0,  117.0,  118.0,  119.0,
     120.0,  121.0,  122.0,  123.0,  124.0,  125.0,  126.0,  127.0,
     128.0,  129.0,  130.0,  131.0,  132.0,  133.0,  134.0,  135.0,
     136.0,  137.0,  138.0,  139.0,  140.0,  141.0,  142.0,  143.0,
     144.0,  145.0,  146.0,  147.0,  148.0,  149.0,  150.0,  151.0,
     152.0,  153.0,  154.0,  155.0,  156.0,  157.0,  158.0,  159.0,
     160.0,  161.0,  162.0,  163.0,  164.0,  165.0,  166.0,  167.0,
     168.0,  169.0,  170.0,  171.0,  172.0,  173.0,  174.0,  175.0,
     176.0,  177.0,  178.0,  179.0,  180.0,  181.0,  182.0,  183.0,
     184.0,  185.0,  186.0,  187.0,  188.0,  189.0,  190.0,  191.0,
     192.0,  193.0,  194.0,  195.0,  196.0,  197.0,  198.0,  199.0,
     200.0,  201.0,  202.0,  203.0,  204.0,  205.0,  206.0,  207.0,
     208.0,  209.0,  210.0,  211.0,  212.0,  213.0,  214.0,  215.0,
     216.0,  217.0,  218.0,  219.0,  220.0,  221.0,  222.0,  223.0,
     224.0,  225.0,  226.0,  227.0,  228.0,  229.0,  230.0,  231.0,
     232.0,  233.0,  234.0,  235.0,  236.0,  237.0,  238.0,  239.0,
     240.0,  241.0,  242.0,  243.0,  244.0,  245.0,  246.0,  247.0,
     248.0,  249.0,  250.0,  251.0,  252.0,  253.0,  254.0,  255.0,
     256.0,  257.0,  258.0,  259.0,  260.0,  261.0,  262.0,  263.0,
     264.0,  265.0,  266.0,  267.0,  268.0,  269.0,  270.0,  271.0,
     272.0,  273.0,  274.0,  275.0,  276.0,  277.0,  278.0,  279.0,
     280.0,  281.0,  282.0,  283.0,  284.0,  285.0,  286.0,  287.0,
     288.0,  289.0,  290.0,  291.0,  292.0,  293.0,  294.0,  295.0,
     296.0,  297.0,  298.0,  299.0,  300.0,  301.0,  302.0,  303.0,
     304.0,  305.0,  306.0,  307.0,  308.0,  309.0,  310.0,  311.0,
     312.0,  313.0,  314.0,  315.0,  316.0,  317.0,  318.0,  319.0,
     320.0,  321.0,  322.0,  323.0,  324.0,  325.0,  326.0,  327.0,
     328.0,  329.0,  330.0,  331.0,  332.0,  333.0,  334.0,  335.0,
     336.0,  337.0,  338.0,  339.0,  340.0,  341.0,  342.0,  343.0,
     344.0,  345.0,  346.0,  347.0,  348.0,  349.0,  350.0,  351.0,
     352.0,  353.0,  354.0,  355.0,  356.0,  357.0,  358.0,  359.0,
     360.0,  361.0,  362.0,  363.0,  364.0,  365.0,  366.0,  367.0,
     368.0,  369.0,  370.0,  371.0,  372.0,  373.0,  374.0,  375.0,
     376.0,  377.0,  378.0,  379.0,  380.0,  381.0,  382.0,  383.0,
     384.0,  385.0,  386.0,  387.0,  388.0,  389.0,  390.0,  391.0,
     392.0,  393.0,  394.0,  395.0,  396.0,  397.0,  398.0,  399.0,
     400.0,  401.0,  402.0,  403.0,  404.0,  405.0,  406.0,  407.0,
     408.0,  409.0,  410.0,  411.0,  412.0,  413.0,  414.0,  415.0,
     416.0,  417.0,  418.0,  419.0,  420.0,  421.0,  422.0,  423.0,
     424.0,  425.0,  426.0,  427.0,  428.0,  429.0,  430.0,  431.0,
     432.0,  433.0,  434.0,  435.0,  436.0,  437.0,  438.0,  439.0,
     440.0,  441.0,  442.0,  443.0,  444.0,  445.0,  446.0,  447.0,
     448.0,  449.0,  450.0,  451.0,  452.0,  453.0,  454.0,  455.0,
     456.0,  457.0,  458.0,  459.0,  460.0,  461.0,  462.0,  463.0,
     464.0,  465.0,  466.0,  467.0,  468.0,  469.0,  470.0,  471.0,
     472.0,  473.0,  474.0,  475.0,  476.0,  477.0,  478.0,  479.0,
     480.0,  481.0,  482.0,  483.0,  484.0,  485.0,  486.0,  487.0,
     488.0,  489.0,  490.0,  491.0,  492.0,  493.0,  494.0,  495.0,
     496.0,  497.0,  498.0,  499.0,  500.0,  501.0,  502.0,  503.0,
     504.0,  505.0,  506.0,  507.0,  508.0,  509.0,  510.0,  511.0,
};
const float icv8x64fTab_cv[] =
{
    -256.0f, -255.0f, -254.0f, -253.0f, -252.0f, -251.0f, -250.0f, -249.0f,
    -248.0f, -247.0f, -246.0f, -245.0f, -244.0f, -243.0f, -242.0f, -241.0f,
    -240.0f, -239.0f, -238.0f, -237.0f, -236.0f, -235.0f, -234.0f, -233.0f,
    -232.0f, -231.0f, -230.0f, -229.0f, -228.0f, -227.0f, -226.0f, -225.0f,
    -224.0f, -223.0f, -222.0f, -221.0f, -220.0f, -219.0f, -218.0f, -217.0f,
    -216.0f, -215.0f, -214.0f, -213.0f, -212.0f, -211.0f, -210.0f, -209.0f,
    -208.0f, -207.0f, -206.0f, -205.0f, -204.0f, -203.0f, -202.0f, -201.0f,
    -200.0f, -199.0f, -198.0f, -197.0f, -196.0f, -195.0f, -194.0f, -193.0f,
    -192.0f, -191.0f, -190.0f, -189.0f, -188.0f, -187.0f, -186.0f, -185.0f,
    -184.0f, -183.0f, -182.0f, -181.0f, -180.0f, -179.0f, -178.0f, -177.0f,
    -176.0f, -175.0f, -174.0f, -173.0f, -172.0f, -171.0f, -170.0f, -169.0f,
    -168.0f, -167.0f, -166.0f, -165.0f, -164.0f, -163.0f, -162.0f, -161.0f,
    -160.0f, -159.0f, -158.0f, -157.0f, -156.0f, -155.0f, -154.0f, -153.0f,
    -152.0f, -151.0f, -150.0f, -149.0f, -148.0f, -147.0f, -146.0f, -145.0f,
    -144.0f, -143.0f, -142.0f, -141.0f, -140.0f, -139.0f, -138.0f, -137.0f,
    -136.0f, -135.0f, -134.0f, -133.0f, -132.0f, -131.0f, -130.0f, -129.0f,
    -128.0f, -127.0f, -126.0f, -125.0f, -124.0f, -123.0f, -122.0f, -121.0f,
    -120.0f, -119.0f, -118.0f, -117.0f, -116.0f, -115.0f, -114.0f, -113.0f,
    -112.0f, -111.0f, -110.0f, -109.0f, -108.0f, -107.0f, -106.0f, -105.0f,
    -104.0f, -103.0f, -102.0f, -101.0f, -100.0f,  -99.0f,  -98.0f,  -97.0f,
     -96.0f,  -95.0f,  -94.0f,  -93.0f,  -92.0f,  -91.0f,  -90.0f,  -89.0f,
     -88.0f,  -87.0f,  -86.0f,  -85.0f,  -84.0f,  -83.0f,  -82.0f,  -81.0f,
     -80.0f,  -79.0f,  -78.0f,  -77.0f,  -76.0f,  -75.0f,  -74.0f,  -73.0f,
     -72.0f,  -71.0f,  -70.0f,  -69.0f,  -68.0f,  -67.0f,  -66.0f,  -65.0f,
     -64.0f,  -63.0f,  -62.0f,  -61.0f,  -60.0f,  -59.0f,  -58.0f,  -57.0f,
     -56.0f,  -55.0f,  -54.0f,  -53.0f,  -52.0f,  -51.0f,  -50.0f,  -49.0f,
     -48.0f,  -47.0f,  -46.0f,  -45.0f,  -44.0f,  -43.0f,  -42.0f,  -41.0f,
     -40.0f,  -39.0f,  -38.0f,  -37.0f,  -36.0f,  -35.0f,  -34.0f,  -33.0f,
     -32.0f,  -31.0f,  -30.0f,  -29.0f,  -28.0f,  -27.0f,  -26.0f,  -25.0f,
     -24.0f,  -23.0f,  -22.0f,  -21.0f,  -20.0f,  -19.0f,  -18.0f,  -17.0f,
     -16.0f,  -15.0f,  -14.0f,  -13.0f,  -12.0f,  -11.0f,  -10.0f,   -9.0f,
      -8.0f,   -7.0f,   -6.0f,   -5.0f,   -4.0f,   -3.0f,   -2.0f,   -1.0f,
       0.0f,    1.0f,    2.0f,    3.0f,    4.0f,    5.0f,    6.0f,    7.0f,
       8.0f,    9.0f,   10.0f,   11.0f,   12.0f,   13.0f,   14.0f,   15.0f,
      16.0f,   17.0f,   18.0f,   19.0f,   20.0f,   21.0f,   22.0f,   23.0f,
      24.0f,   25.0f,   26.0f,   27.0f,   28.0f,   29.0f,   30.0f,   31.0f,
      32.0f,   33.0f,   34.0f,   35.0f,   36.0f,   37.0f,   38.0f,   39.0f,
      40.0f,   41.0f,   42.0f,   43.0f,   44.0f,   45.0f,   46.0f,   47.0f,
      48.0f,   49.0f,   50.0f,   51.0f,   52.0f,   53.0f,   54.0f,   55.0f,
      56.0f,   57.0f,   58.0f,   59.0f,   60.0f,   61.0f,   62.0f,   63.0f,
      64.0f,   65.0f,   66.0f,   67.0f,   68.0f,   69.0f,   70.0f,   71.0f,
      72.0f,   73.0f,   74.0f,   75.0f,   76.0f,   77.0f,   78.0f,   79.0f,
      80.0f,   81.0f,   82.0f,   83.0f,   84.0f,   85.0f,   86.0f,   87.0f,
      88.0f,   89.0f,   90.0f,   91.0f,   92.0f,   93.0f,   94.0f,   95.0f,
      96.0f,   97.0f,   98.0f,   99.0f,  100.0f,  101.0f,  102.0f,  103.0f,
     104.0f,  105.0f,  106.0f,  107.0f,  108.0f,  109.0f,  110.0f,  111.0f,
     112.0f,  113.0f,  114.0f,  115.0f,  116.0f,  117.0f,  118.0f,  119.0f,
     120.0f,  121.0f,  122.0f,  123.0f,  124.0f,  125.0f,  126.0f,  127.0f,
     128.0f,  129.0f,  130.0f,  131.0f,  132.0f,  133.0f,  134.0f,  135.0f,
     136.0f,  137.0f,  138.0f,  139.0f,  140.0f,  141.0f,  142.0f,  143.0f,
     144.0f,  145.0f,  146.0f,  147.0f,  148.0f,  149.0f,  150.0f,  151.0f,
     152.0f,  153.0f,  154.0f,  155.0f,  156.0f,  157.0f,  158.0f,  159.0f,
     160.0f,  161.0f,  162.0f,  163.0f,  164.0f,  165.0f,  166.0f,  167.0f,
     168.0f,  169.0f,  170.0f,  171.0f,  172.0f,  173.0f,  174.0f,  175.0f,
     176.0f,  177.0f,  178.0f,  179.0f,  180.0f,  181.0f,  182.0f,  183.0f,
     184.0f,  185.0f,  186.0f,  187.0f,  188.0f,  189.0f,  190.0f,  191.0f,
     192.0f,  193.0f,  194.0f,  195.0f,  196.0f,  197.0f,  198.0f,  199.0f,
     200.0f,  201.0f,  202.0f,  203.0f,  204.0f,  205.0f,  206.0f,  207.0f,
     208.0f,  209.0f,  210.0f,  211.0f,  212.0f,  213.0f,  214.0f,  215.0f,
     216.0f,  217.0f,  218.0f,  219.0f,  220.0f,  221.0f,  222.0f,  223.0f,
     224.0f,  225.0f,  226.0f,  227.0f,  228.0f,  229.0f,  230.0f,  231.0f,
     232.0f,  233.0f,  234.0f,  235.0f,  236.0f,  237.0f,  238.0f,  239.0f,
     240.0f,  241.0f,  242.0f,  243.0f,  244.0f,  245.0f,  246.0f,  247.0f,
     248.0f,  249.0f,  250.0f,  251.0f,  252.0f,  253.0f,  254.0f,  255.0f,
     256.0f,  257.0f,  258.0f,  259.0f,  260.0f,  261.0f,  262.0f,  263.0f,
     264.0f,  265.0f,  266.0f,  267.0f,  268.0f,  269.0f,  270.0f,  271.0f,
     272.0f,  273.0f,  274.0f,  275.0f,  276.0f,  277.0f,  278.0f,  279.0f,
     280.0f,  281.0f,  282.0f,  283.0f,  284.0f,  285.0f,  286.0f,  287.0f,
     288.0f,  289.0f,  290.0f,  291.0f,  292.0f,  293.0f,  294.0f,  295.0f,
     296.0f,  297.0f,  298.0f,  299.0f,  300.0f,  301.0f,  302.0f,  303.0f,
     304.0f,  305.0f,  306.0f,  307.0f,  308.0f,  309.0f,  310.0f,  311.0f,
     312.0f,  313.0f,  314.0f,  315.0f,  316.0f,  317.0f,  318.0f,  319.0f,
     320.0f,  321.0f,  322.0f,  323.0f,  324.0f,  325.0f,  326.0f,  327.0f,
     328.0f,  329.0f,  330.0f,  331.0f,  332.0f,  333.0f,  334.0f,  335.0f,
     336.0f,  337.0f,  338.0f,  339.0f,  340.0f,  341.0f,  342.0f,  343.0f,
     344.0f,  345.0f,  346.0f,  347.0f,  348.0f,  349.0f,  350.0f,  351.0f,
     352.0f,  353.0f,  354.0f,  355.0f,  356.0f,  357.0f,  358.0f,  359.0f,
     360.0f,  361.0f,  362.0f,  363.0f,  364.0f,  365.0f,  366.0f,  367.0f,
     368.0f,  369.0f,  370.0f,  371.0f,  372.0f,  373.0f,  374.0f,  375.0f,
     376.0f,  377.0f,  378.0f,  379.0f,  380.0f,  381.0f,  382.0f,  383.0f,
     384.0f,  385.0f,  386.0f,  387.0f,  388.0f,  389.0f,  390.0f,  391.0f,
     392.0f,  393.0f,  394.0f,  395.0f,  396.0f,  397.0f,  398.0f,  399.0f,
     400.0f,  401.0f,  402.0f,  403.0f,  404.0f,  405.0f,  406.0f,  407.0f,
     408.0f,  409.0f,  410.0f,  411.0f,  412.0f,  413.0f,  414.0f,  415.0f,
     416.0f,  417.0f,  418.0f,  419.0f,  420.0f,  421.0f,  422.0f,  423.0f,
     424.0f,  425.0f,  426.0f,  427.0f,  428.0f,  429.0f,  430.0f,  431.0f,
     432.0f,  433.0f,  434.0f,  435.0f,  436.0f,  437.0f,  438.0f,  439.0f,
     440.0f,  441.0f,  442.0f,  443.0f,  444.0f,  445.0f,  446.0f,  447.0f,
     448.0f,  449.0f,  450.0f,  451.0f,  452.0f,  453.0f,  454.0f,  455.0f,
     456.0f,  457.0f,  458.0f,  459.0f,  460.0f,  461.0f,  462.0f,  463.0f,
     464.0f,  465.0f,  466.0f,  467.0f,  468.0f,  469.0f,  470.0f,  471.0f,
     472.0f,  473.0f,  474.0f,  475.0f,  476.0f,  477.0f,  478.0f,  479.0f,
     480.0f,  481.0f,  482.0f,  483.0f,  484.0f,  485.0f,  486.0f,  487.0f,
     488.0f,  489.0f,  490.0f,  491.0f,  492.0f,  493.0f,  494.0f,  495.0f,
     496.0f,  497.0f,  498.0f,  499.0f,  500.0f,  501.0f,  502.0f,  503.0f,
     504.0f,  505.0f,  506.0f,  507.0f,  508.0f,  509.0f,  510.0f,  511.0f,
};

/////////////////////////////// IPP warpaffine functions /////////////////////////////////
//Possibly very slightly faster than below...
template<typename FLOAT>
inline FLOAT toFloat(int i);

template<>
inline double toFloat<double>(int i)
{
    return icv8x64dTab_cv[i+256];
}

template<>
inline float toFloat<float>(int i)
{
    return icv8x64fTab_cv[i+256];
}

/*inline double toFloat(int i)
{
    return (double)i;
}*/


//Use a macro to define inverse point multiplications -- lots of 0s
#define INV_MAT_MUL(i, x, y) \
double ddest_x##i = A_inv11*x + A_inv12*y + A_inv13;\
double ddest_y##i = A_inv21*x + A_inv22*y + A_inv23;\
double inv_ddest_t##i = 1.0/(A_inv31*x + A_inv32*y + A_inv33);\
int dest_x##i = cvRound(inv_ddest_t##i*ddest_x##i);\
int dest_y##i = cvRound(inv_ddest_t##i*ddest_y##i);

const BB getBB(const double * inv_mat, const CvSize dsize, const CvSize ssize)
{
	double A_inv11 = (double)(inv_mat[0]), A_inv12 = (double)(inv_mat[1]), A_inv13 = (double)(inv_mat[2]);
	double A_inv21 = (double)(inv_mat[3]), A_inv22 = (double)(inv_mat[4]), A_inv23 = (double)(inv_mat[5]);
	double A_inv31 = (double)(inv_mat[6]), A_inv32 = (double)(inv_mat[7]), A_inv33 = (double)(inv_mat[8]);
	INV_MAT_MUL(1, 0.0, 0.0);
	INV_MAT_MUL(2, 0.0, ssize.height);
	INV_MAT_MUL(3, ssize.width, 0.0);
	INV_MAT_MUL(4, ssize.width, ssize.height);
	int minX = std::min<int>(std::min<int>(dest_x1, dest_x2), std::min<int>(dest_x3, dest_x4));
	int minY = std::min<int>(std::min<int>(dest_y1, dest_y2), std::min<int>(dest_y3, dest_y4));
	int maxX = std::max<int>(std::max<int>(dest_x1, dest_x2), std::max<int>(dest_x3, dest_x4));
	int maxY = std::max<int>(std::max<int>(dest_y1, dest_y2), std::max<int>(dest_y3, dest_y4));
	if(minX < 0)
		minX = 0;

	if(minY < 0)
		minY = 0;

	if(maxX > dsize.width)
		maxX = dsize.width;

	if(maxY > dsize.height)
		maxY = dsize.height;

	return BB(minX, minY, maxX, maxY);
}

//! Faster version of various internal OpenCV warping functions
/*! Nearest-neighbour interpolation supported. A bounding box limits the area we iterate over. Division (perspective transformation)
 *  is implemented with a binomial expansion. Floating-point maths done with doubles (faster than floats). Repeated arithmatic, 
 *  filling of outliers, and unnecessary logic removed. Also implements affine warping.
*/
template<int NUM_CHANNELS, int INTERP_METHOD, bool AFFINE>
class CFastWarp
{
public:

	static void warpInBB(const CvSize & ssize, uchar *dst, const int dststep, const double *matrix, const uchar *src, const int step, const BB & bb)
	{
        double A31 = 0, A32 = 0;
        if(!AFFINE)
            A31 = (double)(matrix[6]), A32 = (double)(matrix[7]);
        double dNearAffine = fabs(A31)+fabs(A32);

        //float is a little faster, particularly for division
        if(dNearAffine < 0.00005 ) //Often have an ID transform
        	warpInBB_int<true, true, float>(ssize, dst, dststep, matrix, src, step, bb);
        else if(dNearAffine < 0.0008 ) //about 1/1000=mosaic max size
        	warpInBB_int<true, AFFINE, float>(ssize, dst, dststep, matrix, src, step, bb);
        else
        	warpInBB_int<false, AFFINE, float>(ssize, dst, dststep, matrix, src, step, bb);

	}

	template<bool FAST_DIVIDE, bool IS_AFFINE, typename FLOAT>
    static void warpInBB_int(const CvSize & ssize, uchar *dst, const int dststep, const double *matrix, const uchar *src, const int step, const BB & bb)
    {
        unsigned uMaxSourceWidth = (unsigned )((ssize.width));
        unsigned uMaxSourceHeight = (unsigned )((ssize.height));
        //FLOAT x_cached = 0.;
        FLOAT x_inv_cached = 0.; //Remember what x and 1/x were last time
        //const FLOAT DELTA_MAX = 1e-2;
        // For linear interp.
        unsigned uMaxSourceX = (unsigned )((ssize.width - 1));
        unsigned uMaxSourceY = (unsigned )((ssize.height - 1));
        unsigned uMaxSourceWidthPlus1 = (unsigned )((ssize.width + 1));
        unsigned uMaxSourceHeightPlus1 = (unsigned )((ssize.height + 1));
        int maxSourceX = (ssize.width - 1);
        int maxSourceY = (ssize.height - 1);
        ////
        //const int MAX_TIME_BETWEEN_EXACT_DIVISIONS = 2000;
        //int maxTimeUntilNextDiv = MAX_TIME_BETWEEN_EXACT_DIVISIONS;

        dst += dststep * bb.minY;

        FLOAT A11 = (FLOAT)(matrix[0]), A12 = (FLOAT)(matrix[1]), A13 = (FLOAT)(matrix[2]);
        FLOAT A21 = (FLOAT)(matrix[3]), A22 = (FLOAT)(matrix[4]), A23 = (FLOAT)(matrix[5]);
        FLOAT A31 = 0;
        FLOAT A32 = 0, A33 = 1;
        if(!IS_AFFINE)
            A31 = (FLOAT)(matrix[6]), A32 = (FLOAT)(matrix[7]), A33 = (FLOAT)(matrix[8]);

        FLOAT xs0_init = A12 * bb.minY + A13;
        FLOAT ys0_init = A22 * bb.minY + A23;
        FLOAT ws_init = A32 * bb.minY + A33; //about 1 if A32 about 0

        DEBUGONLY(int nFailures = 0, nSuccesses = 0);

		FLOAT dMinX = (FLOAT)bb.minX;

		uchar * dst_rowStart = dst+bb.minX*NUM_CHANNELS;

        for ( int y = bb.minY; y < bb.maxY; y++ )
		{
        	//if(dst[bb.minX*NUM_CHANNELS] && dst[bb.maxX*NUM_CHANNELS])
        		//continue; //assume line filled. TODO binary chop to find minX and maxX

			FLOAT xs0 = xs0_init + dMinX*A11;
			FLOAT ys0 = ys0_init + dMinX*A21;
			FLOAT ws = ws_init + dMinX*A31;

			x_inv_cached = 1.0f / ws;
			//x_cached = ws;
			uchar* dspPx = dst_rowStart;

			for ( int x = bb.minX; x < bb.maxX; x++, dspPx += NUM_CHANNELS)
			{
				//uchar* dspPx = dst + x*NUM_CHANNELS;

				//Do division:
				//if we're close enough then binomial expand around previous soln
				if(!IS_AFFINE)
				{
					//FLOAT x_minus_x_old = ws-x_cached;

					//We're close enough (todo: if we do this lots then may need to stop errors accumulating??)
					//Gen. binomial expansion (1+x)^-1 ~= 1-x
					//
					//	given x^-1 want (x+d)^-1
					//
					//	(x+d)^-1 = ((x*x^-1)*(x+d))^-1 = ((x)*(1+x^-1*d))^-1 = x^-1 * (1+x^-1*d)^-1
					//	 = x^-1 * (1 - x^-1*d) when x^-1*d small
					//

					if(FAST_DIVIDE)
					{
						/*if(fabs(x_minus_x_old - A31) > 1e-8)
						{
							cout << "x_minus_x_old=" << x_minus_x_old  << endl;
							cout << "A31=" << A31 << endl;
						}*/
						x_inv_cached = x_inv_cached*((FLOAT)(1.) - x_inv_cached*A31);
						//x_cached = ws;

						DEBUGONLY(
						FLOAT true_inv = (1.0f / ws);
						if(fabs((FLOAT)(x_inv_cached - true_inv)*std::max<FLOAT>(xs0, ys0)) > 0.5) //more than 1/2 a px out
							nFailures++;
						else
							nSuccesses++);
					}
					else
						x_inv_cached = 1.0f / ws;
				}

				FLOAT xs = IS_AFFINE ? xs0 : xs0*x_inv_cached;
				FLOAT ys = IS_AFFINE ? ys0 : ys0*x_inv_cached;
				int ixs = cvRound(xs);//Coordinates of this dest px in source image
				int iys = cvRound(ys);

				if(INTERP_METHOD == CV_INTER_LINEAR)
				{
					FLOAT a = xs - ixs;
					FLOAT b = ys - iys;
					FLOAT p0, p1;
					// (unsigned) casts mean can do 0 <= x < val
					if ( (unsigned)ixs < uMaxSourceX && (unsigned)iys < uMaxSourceY ) {
						const uchar* ptr = src + step*iys + ixs*NUM_CHANNELS;
						for ( int k = 0; k < NUM_CHANNELS; k++ ) {
							p0 = toFloat<FLOAT>(ptr[k]) + a * (toFloat<FLOAT>(ptr[k+NUM_CHANNELS]) - toFloat<FLOAT>(ptr[k]));
							p1 = toFloat<FLOAT>(ptr[k+step]) + a * (toFloat<FLOAT>(ptr[k+NUM_CHANNELS+step]) - toFloat<FLOAT>(ptr[k+step]));
							FLOAT val = p0 + b*(p1 - p0);
							if(val<0.) val = 0.;
							if(val>255.) val = 255.;
							dspPx[k] = (uchar)cvRound(val);
						}
					} else if ( (unsigned)(ixs+1) < uMaxSourceWidthPlus1 && (unsigned)(iys+1) < uMaxSourceHeightPlus1) {
						int x0 = ((unsigned)(ixs) < uMaxSourceWidth ? (ixs) : (ixs) < 0 ? 0 : maxSourceX);
						int y0 = ((unsigned)(iys) < uMaxSourceHeight ? (iys) : (iys) < 0 ? 0 : maxSourceY);
						int x1 = ((unsigned)(ixs + 1) < uMaxSourceWidth ? (ixs + 1) : (ixs + 1) < 0 ? 0 : maxSourceX);
						int y1 = ((unsigned)(iys + 1) < uMaxSourceHeight ? (iys + 1) : (iys + 1) < 0 ? 0 : maxSourceY);

						//These are the pointers to the 4 corners for lin. interp.
						const uchar* ptr0, *ptr1, *ptr2, *ptr3;
						ptr0 = src + y0*step + x0*NUM_CHANNELS;//Repeated arith.
						ptr1 = src + y0*step + x1*NUM_CHANNELS;
						ptr2 = src + y1*step + x0*NUM_CHANNELS;
						ptr3 = src + y1*step + x1*NUM_CHANNELS;
						for ( int k = 0; k < NUM_CHANNELS; k++ )
						{
							//bilinear interpolation
							p0 = toFloat<FLOAT>(ptr0[k]) + a * (toFloat<FLOAT>(ptr1[k]) - toFloat<FLOAT>(ptr0[k]));
							p1 = toFloat<FLOAT>(ptr2[k]) + a * (toFloat<FLOAT>(ptr3[k]) - toFloat<FLOAT>(ptr2[k]));
							FLOAT val = p0 + b*(p1 - p0);
							if(val<0.) val = 0.;
							if(val>255.) val = 255.;
							dspPx[k] = (uchar)cvRound(val);
						}
					}
				}
				else
				{
					// (unsigned) casts mean can do 0 <= x < val
					if ( (unsigned)ixs < uMaxSourceWidth && (unsigned)iys < uMaxSourceHeight ) {
						const uchar* ptr = src + step*iys + ixs*NUM_CHANNELS;

						//for ( k = 0; k < NUM_CHANNELS; k++ ) {
						//    dst[x*NUM_CHANNELS+k] = ptr[k];
					   // }
						//expand loop
						//uchar* dspPx = dst + x*NUM_CHANNELS;
						if(NUM_CHANNELS == 4)
						{
							const int * pSrcInt = (int *)(void *)ptr;
							int * pDestInt = (int *)(void *)dspPx;
							*pDestInt = *pSrcInt;
						}
						else
						{
							dspPx[0] = ptr[0];
							if(NUM_CHANNELS == 3)
							{
								dspPx[1] = ptr[1];
								dspPx[2] = ptr[2];
							}
						}
					}
				}

				xs0 += A11;
				ys0 += A21;
				if(!FAST_DIVIDE || IS_DEBUG)
					ws += A31;
			}

			dst_rowStart += dststep;

			xs0_init += A12;
			ys0_init += A22;
			ws_init += A32;
		}
        DEBUGONLY( if(nFailures > 0)
        	cout << "Approx. division has failed " << (nFailures*100)/(nSuccesses+nFailures) << "% of the time\nA31=" << A31 << ", A32=" << A32 << " ws_init + dMinX*A31=" << ws_init + dMinX*A31 << endl);
    }
    /*static void warpInBB(const CvSize & ssize, uchar *dst, const int dststep, const double *matrix, const uchar *src, const int step, const int cn, const BB & bb)
    {
        unsigned uMaxSourceWidth = (unsigned )((ssize.width));
        unsigned uMaxSourceHeight = (unsigned )((ssize.height));
        double x_cached = 0., x_inv_cached = 0.; //Remember what x and 1/x were last time
        const double DELTA_MAX = 1e-2;
        // For linear interp.
        unsigned uMaxSourceX = (unsigned )((ssize.width - 1));
        unsigned uMaxSourceY = (unsigned )((ssize.height - 1));
        unsigned uMaxSourceWidthPlus1 = (unsigned )((ssize.width + 1));
        unsigned uMaxSourceHeightPlus1 = (unsigned )((ssize.height + 1));
        int maxSourceX = (ssize.width - 1);
        int maxSourceY = (ssize.height - 1);
        ////
        const int MAX_TIME_BETWEEN_EXACT_DIVISIONS = 2000;
        int maxTimeUntilNextDiv = MAX_TIME_BETWEEN_EXACT_DIVISIONS;

        dst += dststep * bb.minY;

        double A11 = (double)(matrix[0]), A12 = (double)(matrix[1]), A13 = (double)(matrix[2]);
        double A21 = (double)(matrix[3]), A22 = (double)(matrix[4]), A23 = (double)(matrix[5]);
        double A31 = 0, A32 = 0, A33 = 1;
        if(!AFFINE)
            A31 = (double)(matrix[6]), A32 = (double)(matrix[7]), A33 = (double)(matrix[8]);

        double xs0_init = A12 * bb.minY + A13;
        double ys0_init = A22 * bb.minY + A23;
        double ws_init = A32 * bb.minY + A33;

        for ( int y = bb.minY; y < bb.maxY; y++ )
		{
			double dMinX = (double)bb.minX;
			double xs0 = xs0_init + dMinX*A11;
			double ys0 = ys0_init + dMinX*A21;
			double ws = ws_init + dMinX*A31;

			for ( int x = bb.minX; x < bb.maxX; x++)
			{
				//Do division:
				//if we're close enough then binomial expand around previous soln
				if(!AFFINE)
				{
					double x_minus_x_old = ws-x_cached;
					if(maxTimeUntilNextDiv==0 || fabs(x_minus_x_old) > DELTA_MAX)
					{
						x_inv_cached = 1. / ws;
						maxTimeUntilNextDiv = MAX_TIME_BETWEEN_EXACT_DIVISIONS;
					}
					else
					{
						//We're close enough (todo: if we do this lots then may need to stop errors accumulating??)
						//Gen. binomial expansion (1+x)^-1 ~= 1-x
						//
						//	given x^-1 want (x+d)^-1
						//
						//	(x+d)^-1 = ((x*x^-1)*(x+d))^-1 = ((x)*(1+x^-1*d))^-1 = x^-1 * (1+x^-1*d)^-1
						//	 = x^-1 * (1 - x^-1*d) when x^-1*d small
						//
						x_inv_cached = x_inv_cached*(1. - x_inv_cached*x_minus_x_old);
						maxTimeUntilNextDiv--;
					}
					x_cached = ws;
					DEBUGONLY(
					double true_inv = (1. / ws);
					if(fabs((double)(x_inv_cached - true_inv)) > 0.0015)
					    THROW("Approx. division has failed"));
				}

				double xs = AFFINE ? xs0 : xs0*x_inv_cached;
				double ys = AFFINE ? ys0 : ys0*x_inv_cached;
				int ixs = cvRound(xs);//Coordinates of this dest px in source image
				int iys = cvRound(ys);

				if(INTERP_METHOD == CV_INTER_LINEAR)
				{
					double a = xs - ixs;
					double b = ys - iys;
					double p0, p1;
					// (unsigned) casts mean can do 0 <= x < val
					if ( (unsigned)ixs < uMaxSourceX && (unsigned)iys < uMaxSourceY ) {
						const uchar* ptr = src + step*iys + ixs*cn;
						for ( int k = 0; k < cn; k++ ) {
							p0 = toFloat(ptr[k]) + a * (toFloat(ptr[k+cn]) - toFloat(ptr[k]));
							p1 = toFloat(ptr[k+step]) + a * (toFloat(ptr[k+cn+step]) - toFloat(ptr[k+step]));
							double val = p0 + b*(p1 - p0);
							if(val<0.) val = 0.;
							if(val>255.) val = 255.;
							dst[x*cn+k] = (uchar)cvRound(val);
						}
					} else if ( (unsigned)(ixs+1) < uMaxSourceWidthPlus1 && (unsigned)(iys+1) < uMaxSourceHeightPlus1) {
						int x0 = ((unsigned)(ixs) < uMaxSourceWidth ? (ixs) : (ixs) < 0 ? 0 : maxSourceX);
						int y0 = ((unsigned)(iys) < uMaxSourceHeight ? (iys) : (iys) < 0 ? 0 : maxSourceY);
						int x1 = ((unsigned)(ixs + 1) < uMaxSourceWidth ? (ixs + 1) : (ixs + 1) < 0 ? 0 : maxSourceX);
						int y1 = ((unsigned)(iys + 1) < uMaxSourceHeight ? (iys + 1) : (iys + 1) < 0 ? 0 : maxSourceY);

						//These are the pointers to the 4 corners for lin. interp.
						const uchar* ptr0, *ptr1, *ptr2, *ptr3;
						ptr0 = src + y0*step + x0*cn;//Repeated arith.
						ptr1 = src + y0*step + x1*cn;
						ptr2 = src + y1*step + x0*cn;
						ptr3 = src + y1*step + x1*cn;
						for ( int k = 0; k < NUM_CHANNELS; k++ )
						{
							//bilinear interpolation
							p0 = toFloat(ptr0[k]) + a * (toFloat(ptr1[k]) - toFloat(ptr0[k]));
							p1 = toFloat(ptr2[k]) + a * (toFloat(ptr3[k]) - toFloat(ptr2[k]));
							double val = p0 + b*(p1 - p0);
							if(val<0.) val = 0.;
							if(val>255.) val = 255.;
							dst[x*cn+k] = (uchar)cvRound(val);
						}
					}
				}
				else
				{
					// (unsigned) casts mean can do 0 <= x < val
					if ( (unsigned)ixs < uMaxSourceWidth && (unsigned)iys < uMaxSourceHeight ) {
						const uchar* ptr = src + step*iys + ixs*cn;
						////Copy 3 chars at once by casting to int: Could cause access violation if the last 3 bytes are in pos 2 or 3
						// *((int *)(dst + x*cn)) = *((const int *)ptr);
						//for ( k = 0; k < cn; k++ ) {
						//    dst[x*cn+k] = ptr[k];
					   // }
						//expand loop
						uchar* dspPx = dst + x*cn;
						dspPx[0] = ptr[0];
						if(NUM_CHANNELS == 3)
						{
							dspPx[1] = ptr[1];
							dspPx[2] = ptr[2];
						}
					}
				}

				xs0 += A11;
				ys0 += A21;
				ws += A31;
			}

			dst += dststep;

			xs0_init += A12;
			ys0_init += A22;
			ws_init += A32;
		}
    }*/

    static void warp_BB_8u(const uchar *src, int step, CvSize ssize, uchar *dst, int dststep, CvSize dsize, const double *matrix, const double *inv_mat)
    {
        step /= sizeof (src[0]);
        dststep /= sizeof (dst[0]);
        //Calculate the positions in the dest image of the 4 corners of the source im:
        //inv_mat*p_src = p_dest
        double AinvData[9];
        if(AFFINE)
		{
			double Adata[9];
			CvMat A = cvMat(3, 3, CV_64FC1, Adata);
			CvMat Ainv = cvMat(3, 3, CV_64FC1, AinvData);

			for(int i=0;i<6;i++)
				Adata[i] = matrix[i];
			Adata[6] = 0; Adata[7]= 0; Adata[8] = 1;

			cvInv(&A, &Ainv);

			inv_mat = AinvData;
		}
        BB bbAll(getBB(inv_mat, dsize, ssize));

        //const int NUM_THREADS = 2; //todo only 2 atm
        BB bbFirstThread;
        BB bb2ndThread;
        bbAll.split(bbFirstThread, bb2ndThread);

        boost::thread oneThread(boost::bind(warpInBB, boost::ref(ssize), dst, dststep, matrix, src, step, boost::ref(bbFirstThread)));
        warpInBB(ssize, dst, dststep, matrix, src, step, bb2ndThread);
        oneThread.join();
	}
};

void
WarpAffine( const CvArr* srcarr, CvArr* dstarr, const CvMat* matrix,
              int flags )
{
    //static CvFuncTable bilin_tab;

    //CV_FUNCNAME( "cvWarpAffine" );

    CvMat srcstub, *src = (CvMat*)srcarr;
    CvMat dststub, *dst = (CvMat*)dstarr;
    int type, depth, cn;//, *ofs = 0;
    double src_matrix[6], dst_matrix[6];
    int method = flags & 3;
    CvMat srcAb = cvMat( 2, 3, CV_64F, src_matrix ),
          dstAb = cvMat( 2, 3, CV_64F, dst_matrix ),
          A, b, invA, invAb;

    //CvWarpAffineFunc func;
    CvSize ssize, dsize;

    void (CV_CDECL *func)(const uchar *,int,CvSize,uchar *,int,CvSize,const double *,const double *) = 0;

     src = cvGetMat( srcarr, &srcstub );
     dst = cvGetMat( dstarr, &dststub );
    type = CV_MAT_TYPE(src->type);
    depth = CV_MAT_DEPTH(type);
    cn = CV_MAT_CN(type);
    if( cn > 4 )
        THROW( "CV_BadNumChannels" );

    if( !CV_ARE_TYPES_EQ( src, dst ))
        THROW( "CV_StsUnmatchedFormats" );

    if( !CV_IS_MAT(matrix) || CV_MAT_CN(matrix->type) != 1 ||
        CV_MAT_DEPTH(matrix->type) < CV_32F || matrix->rows != 2 || matrix->cols != 3 )
        THROW(  "Transformation matrix should be 2x3 floating-point single-channel matrix" );

    if( flags & CV_WARP_INVERSE_MAP )
        cvConvertScale( matrix, &dstAb );
    else
    {
        // [R|t] -> [R^-1 | -(R^-1)*t]
        cvConvertScale( matrix, &srcAb );
        cvGetCols( &srcAb, &A, 0, 2 );
        cvGetCol( &srcAb, &b, 2 );
        cvGetCols( &dstAb, &invA, 0, 2 );
        cvGetCol( &dstAb, &invAb, 2 );
        cvInvert( &A, &invA, CV_SVD );
        cvGEMM( &invA, &b, -1, 0, 0, &invAb );
    }


    ssize = cvSize(src->width, src->height);// cvGetMatSize(src);
    dsize = cvSize(dst->width, dst->height);//cvGetMatSize(dst);

    /*if( icvWarpAffineBack_8u_C1R_p && MIN( ssize.width, dsize.width ) >= 4 &&
        MIN( ssize.height, dsize.height ) >= 4 )
    {
        CvWarpAffineBackIPPFunc ipp_func =
            type == CV_8UC1 ? icvWarpAffineBack_8u_C1R_p :
            type == CV_8UC3 ? icvWarpAffineBack_8u_C3R_p :
            type == CV_8UC4 ? icvWarpAffineBack_8u_C4R_p :
            type == CV_32FC1 ? icvWarpAffineBack_32f_C1R_p :
            type == CV_32FC3 ? icvWarpAffineBack_32f_C3R_p :
            type == CV_32FC4 ? icvWarpAffineBack_32f_C4R_p : 0;

        if( ipp_func && CV_INTER_NN <= method && method <= CV_INTER_AREA )
        {
            int srcstep = src->step ? src->step : CV_STUB_STEP;
            int dststep = dst->step ? dst->step : CV_STUB_STEP;
            CvRect srcroi = {0, 0, ssize.width, ssize.height};
            CvRect dstroi = {0, 0, dsize.width, dsize.height};

            // this is not the most efficient way to fill outliers
            if( flags & CV_WARP_FILL_OUTLIERS )
                cvSet( dst, fillval );

            if( ipp_func( src->data.ptr, ssize, srcstep, srcroi,
                          dst->data.ptr, dststep, dstroi,
                          dstAb.data.db, 1 << method ) >= 0 )
                EXIT;
        }
    }

    cvScalarToRawData( &fillval, fillbuf, CV_MAT_TYPE(src->type), 0 );
    ofs = (int*)cvStackAlloc( dst->cols*2*sizeof(ofs[0]) );
    for( k = 0; k < dst->cols; k++ )
    {
        ofs[2*k] = CV_FLT_TO_FIX( dst_matrix[0]*k, ICV_WARP_SHIFT );
        ofs[2*k+1] = CV_FLT_TO_FIX( dst_matrix[3]*k, ICV_WARP_SHIFT );
    }*/

    if(cn == 1)
    {
        if(method == CV_INTER_NN)
            func = CFastWarp<1, CV_INTER_NN, true>::warp_BB_8u;
        else
            func = CFastWarp<1, CV_INTER_LINEAR, true>::warp_BB_8u;
    }
    else if(cn == 3)
    {
        if(method == CV_INTER_NN)
            func = CFastWarp<3, CV_INTER_NN, true>::warp_BB_8u;
        else
            func = CFastWarp<3, CV_INTER_LINEAR, true>::warp_BB_8u;
    }
    else if(cn == 4)
    {
        if(method == CV_INTER_NN)
            func = CFastWarp<3, CV_INTER_NN, true>::warp_BB_8u;
        else
            func = CFastWarp<3, CV_INTER_LINEAR, true>::warp_BB_8u;
    }

    if( !func )
        THROW( "CV_StsUnsupportedFormat No Affine warp function available" );

    func( src->data.ptr, src->step, ssize, dst->data.ptr,
                     dst->step, dsize, dst_matrix, 0 );

}


/****************************************************************************************\
*                                    WarpPerspective                                     *
\****************************************************************************************/
/*
#define ICV_DEF_WARP_PERSPECTIVE_BILINEAR_FUNC( flavor, arrtype, load_macro, cast_macro )\
static CvStatus CV_STDCALL                                                  \
icvWarpPerspective_Bilinear_##flavor##_CnR(                                 \
    const arrtype* src, int step, CvSize ssize,                             \
    arrtype* dst, int dststep, CvSize dsize,                                \
    const double* matrix, int cn,                                           \
    const arrtype* fillval )                                                \
{                                                                           \
    int x, y, k;                                                            \
    float A11 = (float)matrix[0], A12 = (float)matrix[1], A13 = (float)matrix[2];\
    float A21 = (float)matrix[3], A22 = (float)matrix[4], A23 = (float)matrix[5];\
    float A31 = (float)matrix[6], A32 = (float)matrix[7], A33 = (float)matrix[8];\
                                                                            \
    step /= sizeof(src[0]);                                                 \
    dststep /= sizeof(dst[0]);                                              \
                                                                            \
    for( y = 0; y < dsize.height; y++, dst += dststep )                     \
    {                                                                       \
        float xs0 = A12*y + A13;                                            \
        float ys0 = A22*y + A23;                                            \
        float ws = A32*y + A33;                                             \
                                                                            \
        for( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 )\
        {                                                                   \
            float inv_ws = 1.f/ws;                                          \
            float xs = xs0*inv_ws;                                          \
            float ys = ys0*inv_ws;                                          \
            int ixs = cvFloor(xs);                                          \
            int iys = cvFloor(ys);                                          \
            float a = xs - ixs;                                             \
            float b = ys - iys;                                             \
            float p0, p1;                                                   \
                                                                            \
            if( (unsigned)ixs < (unsigned)(ssize.width - 1) &&              \
                (unsigned)iys < (unsigned)(ssize.height - 1) )              \
            {                                                               \
                const arrtype* ptr = src + step*iys + ixs*cn;               \
                                                                            \
                for( k = 0; k < cn; k++ )                                   \
                {                                                           \
                    p0 = load_macro(ptr[k]) +                               \
                        a * (load_macro(ptr[k+cn]) - load_macro(ptr[k]));   \
                    p1 = load_macro(ptr[k+step]) +                          \
                        a * (load_macro(ptr[k+cn+step]) -                   \
                             load_macro(ptr[k+step]));                      \
                    dst[x*cn+k] = (arrtype)cast_macro(p0 + b*(p1 - p0));    \
                }                                                           \
            }                                                               \
            else if( (unsigned)(ixs+1) < (unsigned)(ssize.width+1) &&       \
                     (unsigned)(iys+1) < (unsigned)(ssize.height+1))        \
            {                                                               \
                int x0 = ICV_WARP_CLIP_X( ixs );                            \
                int y0 = ICV_WARP_CLIP_Y( iys );                            \
                int x1 = ICV_WARP_CLIP_X( ixs + 1 );                        \
                int y1 = ICV_WARP_CLIP_Y( iys + 1 );                        \
                const arrtype* ptr0, *ptr1, *ptr2, *ptr3;                   \
                                                                            \
                ptr0 = src + y0*step + x0*cn;                               \
                ptr1 = src + y0*step + x1*cn;                               \
                ptr2 = src + y1*step + x0*cn;                               \
                ptr3 = src + y1*step + x1*cn;                               \
                                                                            \
                for( k = 0; k < cn; k++ )                                   \
                {                                                           \
                    p0 = load_macro(ptr0[k]) +                              \
                        a * (load_macro(ptr1[k]) - load_macro(ptr0[k]));    \
                    p1 = load_macro(ptr2[k]) +                              \
                        a * (load_macro(ptr3[k]) - load_macro(ptr2[k]));    \
                    dst[x*cn+k] = (arrtype)cast_macro(p0 + b*(p1 - p0));    \
                }                                                           \
            }                                                               \
            else if( fillval )                                              \
                for( k = 0; k < cn; k++ )                                   \
                    dst[x*cn+k] = fillval[k];                               \
        }                                                                   \
    }                                                                       \
                                                                            \
    return CV_OK;                                                           \
}


#define ICV_WARP_SCALE_ALPHA(x) ((x)*(1./(ICV_WARP_MASK+1)))*/

//ICV_DEF_WARP_PERSPECTIVE_BILINEAR_FUNC( 8u, uchar, CV_8TO32F, cvRound ) //TB: This is the macro function defn, expoands the macro above
//ICV_DEF_WARP_PERSPECTIVE_BILINEAR_FUNC( 16u, ushort, CV_NOP, cvRound )
//ICV_DEF_WARP_PERSPECTIVE_BILINEAR_FUNC( 32f, float, CV_NOP, CV_NOP )

//These expand to:
//static CvStatus  icvWarpPerspective_Bilinear_8u_CnR( const uchar* src, int step, CvSize ssize, uchar* dst, int dststep, CvSize dsize, const double* matrix, int cn, const uchar* fillval ) { int x, y, k; float A11 = (float)matrix[0], A12 = (float)matrix[1], A13 = (float)matrix[2]; float A21 = (float)matrix[3], A22 = (float)matrix[4], A23 = (float)matrix[5]; float A31 = (float)matrix[6], A32 = (float)matrix[7], A33 = (float)matrix[8]; step /= sizeof(src[0]); dststep /= sizeof(dst[0]); for( y = 0; y < dsize.height; y++, dst += dststep ) { float xs0 = A12*y + A13; float ys0 = A22*y + A23; float ws = A32*y + A33; for( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) { float inv_ws = 1.f/ws; float xs = xs0*inv_ws; float ys = ys0*inv_ws; int ixs = cvFloor(xs); int iys = cvFloor(ys); float a = xs - ixs; float b = ys - iys; float p0, p1; if( (unsigned)ixs < (unsigned)(ssize.width - 1) && (unsigned)iys < (unsigned)(ssize.height - 1) ) { const uchar* ptr = src + step*iys + ixs*cn; for( k = 0; k < cn; k++ ) { p0 = icv8x32fTab_cv[(ptr[k])+256] + a * (icv8x32fTab_cv[(ptr[k+cn])+256] - icv8x32fTab_cv[(ptr[k])+256]); p1 = icv8x32fTab_cv[(ptr[k+step])+256] + a * (icv8x32fTab_cv[(ptr[k+cn+step])+256] - icv8x32fTab_cv[(ptr[k+step])+256]); dst[x*cn+k] = (uchar)cvRound(p0 + b*(p1 - p0)); } } else if( (unsigned)(ixs+1) < (unsigned)(ssize.width+1) && (unsigned)(iys+1) < (unsigned)(ssize.height+1)) { int x0 = ((unsigned)(ixs) < (unsigned)ssize.width ? (ixs) : (ixs) < 0 ? 0 : ssize.width - 1); int y0 = ((unsigned)(iys) < (unsigned)ssize.height ? (iys) : (iys) < 0 ? 0 : ssize.height - 1); int x1 = ((unsigned)(ixs + 1) < (unsigned)ssize.width ? (ixs + 1) : (ixs + 1) < 0 ? 0 : ssize.width - 1); int y1 = ((unsigned)(iys + 1) < (unsigned)ssize.height ? (iys + 1) : (iys + 1) < 0 ? 0 : ssize.height - 1); const uchar* ptr0, *ptr1, *ptr2, *ptr3; ptr0 = src + y0*step + x0*cn; ptr1 = src + y0*step + x1*cn; ptr2 = src + y1*step + x0*cn; ptr3 = src + y1*step + x1*cn; for( k = 0; k < cn; k++ ) { p0 = icv8x32fTab_cv[(ptr0[k])+256] + a * (icv8x32fTab_cv[(ptr1[k])+256] - icv8x32fTab_cv[(ptr0[k])+256]); p1 = icv8x32fTab_cv[(ptr2[k])+256] + a * (icv8x32fTab_cv[(ptr3[k])+256] - icv8x32fTab_cv[(ptr2[k])+256]); dst[x*cn+k] = (uchar)cvRound(p0 + b*(p1 - p0)); } } else if( fillval ) for( k = 0; k < cn; k++ ) dst[x*cn+k] = fillval[k]; } } return CV_OK; }
//static CvStatus  icvWarpPerspective_Bilinear_16u_CnR( const ushort* src, int step, CvSize ssize, ushort* dst, int dststep, CvSize dsize, const double* matrix, int cn, const ushort* fillval ) { int x, y, k; float A11 = (float)matrix[0], A12 = (float)matrix[1], A13 = (float)matrix[2]; float A21 = (float)matrix[3], A22 = (float)matrix[4], A23 = (float)matrix[5]; float A31 = (float)matrix[6], A32 = (float)matrix[7], A33 = (float)matrix[8]; step /= sizeof(src[0]); dststep /= sizeof(dst[0]); for( y = 0; y < dsize.height; y++, dst += dststep ) { float xs0 = A12*y + A13; float ys0 = A22*y + A23; float ws = A32*y + A33; for( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) { float inv_ws = 1.f/ws; float xs = xs0*inv_ws; float ys = ys0*inv_ws; int ixs = cvFloor(xs); int iys = cvFloor(ys); float a = xs - ixs; float b = ys - iys; float p0, p1; if( (unsigned)ixs < (unsigned)(ssize.width - 1) && (unsigned)iys < (unsigned)(ssize.height - 1) ) { const ushort* ptr = src + step*iys + ixs*cn; for( k = 0; k < cn; k++ ) { p0 = (ptr[k]) + a * ((ptr[k+cn]) - (ptr[k])); p1 = (ptr[k+step]) + a * ((ptr[k+cn+step]) - (ptr[k+step])); dst[x*cn+k] = (ushort)cvRound(p0 + b*(p1 - p0)); } } else if( (unsigned)(ixs+1) < (unsigned)(ssize.width+1) && (unsigned)(iys+1) < (unsigned)(ssize.height+1)) { int x0 = ((unsigned)(ixs) < (unsigned)ssize.width ? (ixs) : (ixs) < 0 ? 0 : ssize.width - 1); int y0 = ((unsigned)(iys) < (unsigned)ssize.height ? (iys) : (iys) < 0 ? 0 : ssize.height - 1); int x1 = ((unsigned)(ixs + 1) < (unsigned)ssize.width ? (ixs + 1) : (ixs + 1) < 0 ? 0 : ssize.width - 1); int y1 = ((unsigned)(iys + 1) < (unsigned)ssize.height ? (iys + 1) : (iys + 1) < 0 ? 0 : ssize.height - 1); const ushort* ptr0, *ptr1, *ptr2, *ptr3; ptr0 = src + y0*step + x0*cn; ptr1 = src + y0*step + x1*cn; ptr2 = src + y1*step + x0*cn; ptr3 = src + y1*step + x1*cn; for( k = 0; k < cn; k++ ) { p0 = (ptr0[k]) + a * ((ptr1[k]) - (ptr0[k])); p1 = (ptr2[k]) + a * ((ptr3[k]) - (ptr2[k])); dst[x*cn+k] = (ushort)cvRound(p0 + b*(p1 - p0)); } } else if( fillval ) for( k = 0; k < cn; k++ ) dst[x*cn+k] = fillval[k]; } } return CV_OK; }
//static CvStatus  icvWarpPerspective_Bilinear_32f_CnR( const float* src, int step, CvSize ssize, float* dst, int dststep, CvSize dsize, const double* matrix, int cn, const float* fillval ) { int x, y, k; float A11 = (float)matrix[0], A12 = (float)matrix[1], A13 = (float)matrix[2]; float A21 = (float)matrix[3], A22 = (float)matrix[4], A23 = (float)matrix[5]; float A31 = (float)matrix[6], A32 = (float)matrix[7], A33 = (float)matrix[8]; step /= sizeof(src[0]); dststep /= sizeof(dst[0]); for( y = 0; y < dsize.height; y++, dst += dststep ) { float xs0 = A12*y + A13; float ys0 = A22*y + A23; float ws = A32*y + A33; for( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) { float inv_ws = 1.f/ws; float xs = xs0*inv_ws; float ys = ys0*inv_ws; int ixs = cvFloor(xs); int iys = cvFloor(ys); float a = xs - ixs; float b = ys - iys; float p0, p1; if( (unsigned)ixs < (unsigned)(ssize.width - 1) && (unsigned)iys < (unsigned)(ssize.height - 1) ) { const float* ptr = src + step*iys + ixs*cn; for( k = 0; k < cn; k++ ) { p0 = (ptr[k]) + a * ((ptr[k+cn]) - (ptr[k])); p1 = (ptr[k+step]) + a * ((ptr[k+cn+step]) - (ptr[k+step])); dst[x*cn+k] = (float)(p0 + b*(p1 - p0)); } } else if( (unsigned)(ixs+1) < (unsigned)(ssize.width+1) && (unsigned)(iys+1) < (unsigned)(ssize.height+1)) { int x0 = ((unsigned)(ixs) < (unsigned)ssize.width ? (ixs) : (ixs) < 0 ? 0 : ssize.width - 1); int y0 = ((unsigned)(iys) < (unsigned)ssize.height ? (iys) : (iys) < 0 ? 0 : ssize.height - 1); int x1 = ((unsigned)(ixs + 1) < (unsigned)ssize.width ? (ixs + 1) : (ixs + 1) < 0 ? 0 : ssize.width - 1); int y1 = ((unsigned)(iys + 1) < (unsigned)ssize.height ? (iys + 1) : (iys + 1) < 0 ? 0 : ssize.height - 1); const float* ptr0, *ptr1, *ptr2, *ptr3; ptr0 = src + y0*step + x0*cn; ptr1 = src + y0*step + x1*cn; ptr2 = src + y1*step + x0*cn; ptr3 = src + y1*step + x1*cn; for( k = 0; k < cn; k++ ) { p0 = (ptr0[k]) + a * ((ptr1[k]) - (ptr0[k])); p1 = (ptr2[k]) + a * ((ptr3[k]) - (ptr2[k])); dst[x*cn+k] = (float)(p0 + b*(p1 - p0)); } } else if( fillval ) for( k = 0; k < cn; k++ ) dst[x*cn+k] = fillval[k]; } } return CV_OK; }

inline double invertApprox(double x)
{
    const double DELTA_MAX = 1e-1;
    double x_minus_1 = x-1.;
    if(x_minus_1 < DELTA_MAX && x_minus_1 > -DELTA_MAX)
    {
        //Taylor expansion (1+x)^-1 ~= 1-x
        return 1.-x_minus_1;
    }
    return 1./x;
}

//Watch out for statics, not really threadsafe
inline double invertApproxCache(double x)
{
    const double DELTA_MAX = 1e-2;

    static double x_cached=0., x_inv_cached=0.; //Remember what x and 1/x were last time

    //and if we're close enough then taylor expand around previous soln
    double x_minus_x_old = x-x_cached;
    if(x_minus_x_old < DELTA_MAX && x_minus_x_old > -DELTA_MAX)
    {
        //We're close enough (todo: if we do this lots then may need to stop errors accumulating??)
        //Taylor expansion (1+x)^-1 ~= 1-x
        /*
            given x^-1 want (x+d)^-1

            (x+d)^-1 = ((x*x^-1)*(x+d))^-1 = ((x)*(1+x^-1*d))^-1 = x^-1 * (1+x^-1*d)^-1
             = x^-1 * (1 - x^-1*d) when x^-1*d small
        */
        x_inv_cached = x_inv_cached*(1. - x_inv_cached*x_minus_x_old);
    }
    else
    {
        x_inv_cached = 1. / x;
    }

    x_cached = x;

    return x_inv_cached;
}

/*static CvStatus icvWarpPerspective_NN_BB_8u_C1( const uchar* src, int step, CvSize ssize, uchar* dst, int dststep, CvSize dsize, const double* matrix, int cn, const double* inv_mat ) {

    double A11 = (double)matrix[0], A12 = (double)matrix[1], A13 = (double)matrix[2];
    double A21 = (double)matrix[3], A22 = (double)matrix[4], A23 = (double)matrix[5];
    double A31 = (double)matrix[6], A32 = (double)matrix[7], A33 = (double)matrix[8];
    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);

    unsigned uMaxSourceWidth = (unsigned)(ssize.width);
    unsigned uMaxSourceHeight = (unsigned)(ssize.height);

    //Calculate the positions in the dest image of the 4 corners of the source im:
    //inv_mat*p_src = p_dest
    double A_inv11 = (double)inv_mat[0], A_inv12 = (double)inv_mat[1], A_inv13 = (double)inv_mat[2];
    double A_inv21 = (double)inv_mat[3], A_inv22 = (double)inv_mat[4], A_inv23 = (double)inv_mat[5];
    double A_inv31 = (double)inv_mat[6], A_inv32 = (double)inv_mat[7], A_inv33 = (double)inv_mat[8];

    INV_MAT_MUL(1, 0.0, 0.0);
    INV_MAT_MUL(2, 0.0, ssize.height);
    INV_MAT_MUL(3, ssize.width, 0.0);
    INV_MAT_MUL(4, ssize.width, ssize.height);

    int minX = std::min<int>(std::min<int>(dest_x1, dest_x2), std::min<int>(dest_x3, dest_x4));
    int minY = std::min<int>(std::min<int>(dest_y1, dest_y2), std::min<int>(dest_y3, dest_y4));
    int maxX = std::max<int>(std::max<int>(dest_x1, dest_x2), std::max<int>(dest_x3, dest_x4));
    int maxY = std::max<int>(std::max<int>(dest_y1, dest_y2), std::max<int>(dest_y3, dest_y4));

    if(minX < 0) minX = 0;
    if(minY < 0) minY = 0;
    if(maxX > dsize.width) maxX = dsize.width;
    if(maxY > dsize.height) maxY = dsize.height;

    double x_cached=0., x_inv_cached=0.; //Remember what x and 1/x were last time
    const double DELTA_MAX = 1e-2;

    const int MAX_TIME_BETWEEN_EXACT_DIVISIONS = 5000;
    int maxTimeUntilNextDiv=MAX_TIME_BETWEEN_EXACT_DIVISIONS;

    dst += dststep*minY;

    double xs0_init = A12*minY + A13;
    double ys0_init = A22*minY + A23;
    double ws_init = A32*minY + A33;

    for ( int y = minY; y < maxY; y++ )
    {
        double dMinX = (double)minX;
        double xs0 = xs0_init + dMinX*A11;
        double ys0 = ys0_init + dMinX*A21;
        double ws = ws_init + dMinX*A31;

        for ( int x = minX; x < maxX; x++)
        {
            //Do division:
            //if we're close enough then binomial expand around previous soln
            double x_minus_x_old = ws-x_cached;
            if(maxTimeUntilNextDiv > 0 && x_minus_x_old < DELTA_MAX && x_minus_x_old > -DELTA_MAX)
            {
                //We're close enough (todo: if we do this lots then may need to stop errors accumulating??)
                //Gen. binomial expansion (1+x)^-1 ~= 1-x
                / *
                    given x^-1 want (x+d)^-1

                    (x+d)^-1 = ((x*x^-1)*(x+d))^-1 = ((x)*(1+x^-1*d))^-1 = x^-1 * (1+x^-1*d)^-1
                     = x^-1 * (1 - x^-1*d) when x^-1*d small
                * /
                x_inv_cached = x_inv_cached*(1. - x_inv_cached*x_minus_x_old);
                maxTimeUntilNextDiv--;
            }
            else
            {
                x_inv_cached = 1. / ws;
                maxTimeUntilNextDiv = MAX_TIME_BETWEEN_EXACT_DIVISIONS;
            }
            x_cached = ws;

            //double true_inv = (1. / ws);
            //if(fabs((double)(x_inv_cached - true_inv)) > 0.001)
            //    throw new grc::GRCException("Approx. division has failed");

            double xs = xs0*x_inv_cached;
            double ys = ys0*x_inv_cached;
            int ixs = cvRound(xs);//Coordinates of this dest px in source image
            int iys = cvRound(ys);

            // (unsigned) casts mean can do 0 <= x < val
            if ( (unsigned)ixs < uMaxSourceWidth && (unsigned)iys < uMaxSourceHeight ) {
                const uchar* ptr = src + step*iys + ixs*cn;
                ////Copy 3 chars at once by casting to int: Could cause access violation if the last 3 bytes are in pos 2 or 3
                //*((int *)(dst + x*cn)) = *((const int *)ptr);
                //for ( k = 0; k < cn; k++ ) {
                //    dst[x*cn+k] = ptr[k];
               // }
                //expand loop
                uchar* dspPx = dst + x*cn;
                dspPx[0] = ptr[0];
            }

            xs0 += A11;
            ys0 += A21;
            ws += A31;
        }

        dst += dststep;

        xs0_init += A12;
        ys0_init += A22;
        ws_init += A32;

    } return CV_OK;
}

static CvStatus icvWarpAffine_NN_BB_8u_C3( const uchar* src, int step, CvSize ssize, uchar* dst, int dststep, CvSize dsize, const double* matrix, int cn, const double* inv_mat ) {

    double A11 = (double)matrix[0], A12 = (double)matrix[1], A13 = (double)matrix[2];
    double A21 = (double)matrix[3], A22 = (double)matrix[4], A23 = (double)matrix[5];

    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);

    unsigned uMaxSourceWidth = (unsigned)(ssize.width);
    unsigned uMaxSourceHeight = (unsigned)(ssize.height);

    //Calculate the positions in the dest image of the 4 corners of the source im:
    //first need a full invertable 3x3 mat:

    double AinvData[9], Adata[9];
    CvMat A = cvMat(3, 3, CV_64FC1, Adata);
    CvMat Ainv = cvMat(3, 3, CV_64FC1, AinvData);

    for(int i=0;i<6;i++)
        Adata[i] = matrix[i];
    Adata[6] = 0; Adata[7]= 0; Adata[8] = 1;

    cvInv(&A, &Ainv);

    //inv_mat*p_src = p_dest
    double A_inv11 = (double)AinvData[0], A_inv12 = (double)AinvData[1], A_inv13 = (double)AinvData[2];
    double A_inv21 = (double)AinvData[3], A_inv22 = (double)AinvData[4], A_inv23 = (double)AinvData[5];
    double A_inv31 = (double)AinvData[6], A_inv32 = (double)AinvData[7], A_inv33 = (double)AinvData[8];

    INV_MAT_MUL(1, 0.0, 0.0);
    INV_MAT_MUL(2, 0.0, ssize.height);
    INV_MAT_MUL(3, ssize.width, 0.0);
    INV_MAT_MUL(4, ssize.width, ssize.height);

    int minX = std::min<int>(std::min<int>(dest_x1, dest_x2), std::min<int>(dest_x3, dest_x4));
    int minY = std::min<int>(std::min<int>(dest_y1, dest_y2), std::min<int>(dest_y3, dest_y4));
    int maxX = std::max<int>(std::max<int>(dest_x1, dest_x2), std::max<int>(dest_x3, dest_x4));
    int maxY = std::max<int>(std::max<int>(dest_y1, dest_y2), std::max<int>(dest_y3, dest_y4));

    if(minX < 0) minX = 0;
    if(minY < 0) minY = 0;
    if(maxX > dsize.width) maxX = dsize.width;
    if(maxY > dsize.height) maxY = dsize.height;

    dst += dststep*minY;

    double xs0_init = A12*minY + A13;
    double ys0_init = A22*minY + A23;

    for ( int y = minY; y < maxY; y++ )
    {
        double dMinX = (double)minX;
        double xs0 = xs0_init + dMinX*A11;
        double ys0 = ys0_init + dMinX*A21;

        for ( int x = minX; x < maxX; x++)
        {
            double xs = xs0;
            double ys = ys0;
            int ixs = cvRound(xs);//Coordinates of this dest px in source image
            int iys = cvRound(ys);

            // (unsigned) casts mean can do 0 <= x < val
            if ( (unsigned)ixs < uMaxSourceWidth && (unsigned)iys < uMaxSourceHeight ) {
                const uchar* ptr = src + step*iys + ixs*cn;
                //expand loop
                uchar* dspPx = dst + x*cn;
                dspPx[0] = ptr[0];
                dspPx[1] = ptr[1];
                dspPx[2] = ptr[2];
            }

            xs0 += A11;
            ys0 += A21;
        }

        dst += dststep;

        xs0_init += A12;
        ys0_init += A22;

    } return CV_OK;
}

static CvStatus icvWarpAffine_NN_BB_8u_C1( const uchar* src, int step, CvSize ssize, uchar* dst, int dststep, CvSize dsize, const double* matrix, int cn, const double* inv_mat ) {

    double A11 = (double)matrix[0], A12 = (double)matrix[1], A13 = (double)matrix[2];
    double A21 = (double)matrix[3], A22 = (double)matrix[4], A23 = (double)matrix[5];

    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);

    unsigned uMaxSourceWidth = (unsigned)(ssize.width);
    unsigned uMaxSourceHeight = (unsigned)(ssize.height);

    //Calculate the positions in the dest image of the 4 corners of the source im:
    //first need a full invertable 3x3 mat:

    double AinvData[9], Adata[9];
    CvMat A = cvMat(3, 3, CV_64FC1, Adata);
    CvMat Ainv = cvMat(3, 3, CV_64FC1, AinvData);

    for(int i=0;i<6;i++)
        Adata[i] = matrix[i];
    Adata[6] = 0; Adata[7]= 0; Adata[8] = 1;

    cvInv(&A, &Ainv);

    //inv_mat*p_src = p_dest
    double A_inv11 = (double)AinvData[0], A_inv12 = (double)AinvData[1], A_inv13 = (double)AinvData[2];
    double A_inv21 = (double)AinvData[3], A_inv22 = (double)AinvData[4], A_inv23 = (double)AinvData[5];
    double A_inv31 = (double)AinvData[6], A_inv32 = (double)AinvData[7], A_inv33 = (double)AinvData[8];

    INV_MAT_MUL(1, 0.0, 0.0);
    INV_MAT_MUL(2, 0.0, ssize.height);
    INV_MAT_MUL(3, ssize.width, 0.0);
    INV_MAT_MUL(4, ssize.width, ssize.height);

    int minX = std::min<int>(std::min<int>(dest_x1, dest_x2), std::min<int>(dest_x3, dest_x4));
    int minY = std::min<int>(std::min<int>(dest_y1, dest_y2), std::min<int>(dest_y3, dest_y4));
    int maxX = std::max<int>(std::max<int>(dest_x1, dest_x2), std::max<int>(dest_x3, dest_x4));
    int maxY = std::max<int>(std::max<int>(dest_y1, dest_y2), std::max<int>(dest_y3, dest_y4));

    if(minX < 0) minX = 0;
    if(minY < 0) minY = 0;
    if(maxX > dsize.width) maxX = dsize.width;
    if(maxY > dsize.height) maxY = dsize.height;

    dst += dststep*minY;

    double xs0_init = A12*minY + A13;
    double ys0_init = A22*minY + A23;

    for ( int y = minY; y < maxY; y++ )
    {
        double dMinX = (double)minX;
        double xs0 = xs0_init + dMinX*A11;
        double ys0 = ys0_init + dMinX*A21;

        for ( int x = minX; x < maxX; x++)
        {
            double xs = xs0;
            double ys = ys0;
            int ixs = cvRound(xs);//Coordinates of this dest px in source image
            int iys = cvRound(ys);

            // (unsigned) casts mean can do 0 <= x < val
            if ( (unsigned)ixs < uMaxSourceWidth && (unsigned)iys < uMaxSourceHeight ) {
                const uchar* ptr = src + step*iys + ixs*cn;
                uchar* dspPx = dst + x*cn;
                dspPx[0] = ptr[0];
            }

            xs0 += A11;
            ys0 += A21;
        }

        dst += dststep;

        xs0_init += A12;
        ys0_init += A22;

    } return CV_OK;
}

//expand inner loop instead of unsafe cast
//go back to incrementing xs0 etc.
//test division

/ * Approx NN version with inlined taylor series division (no statics)--Yes this is faster that the inlining * /
CvStatus icvWarpPerspective_NN_8u_CnR( const uchar* src, int step, CvSize ssize, uchar* dst, int dststep, CvSize dsize, const double* matrix, int cn, const uchar* fillval ) {

    double A11 = (double)matrix[0], A12 = (double)matrix[1], A13 = (double)matrix[2];
    double A21 = (double)matrix[3], A22 = (double)matrix[4], A23 = (double)matrix[5];
    double A31 = (double)matrix[6], A32 = (double)matrix[7], A33 = (double)matrix[8];
    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);

    unsigned uMaxSourceWidth = (unsigned)(ssize.width);
    unsigned uMaxSourceHeight = (unsigned)(ssize.height);

    double x_cached=0., x_inv_cached=0.; //Remember what x and 1/x were last time

    const double DELTA_MAX = 1e-2; // These vals work well with an error usually much less than 0.001
    const int MAX_TIME_BETWEEN_EXACT_DIVISIONS = 5000;

    int maxTimeUntilNextDiv=MAX_TIME_BETWEEN_EXACT_DIVISIONS;
    for ( int y = 0; y < dsize.height; y++, dst += dststep ) {
        double xs0 = A12*y + A13;
        double ys0 = A22*y + A23;
        double ws = A32*y + A33;
        for ( int x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) {

            //Do division:
            //if we're close enough then binomial expand around previous soln
            double x_minus_x_old = ws-x_cached;
            if(maxTimeUntilNextDiv > 0 && x_minus_x_old < DELTA_MAX && x_minus_x_old > -DELTA_MAX)
            {
                //We're close enough (todo: if we do this lots then may need to stop errors accumulating??)
                //Gen. binomial expansion (1+x)^-1 ~= 1-x
                / *
                    given x^-1 want (x+d)^-1

                    (x+d)^-1 = ((x*x^-1)*(x+d))^-1 = ((x)*(1+x^-1*d))^-1 = x^-1 * (1+x^-1*d)^-1
                     = x^-1 * (1 - x^-1*d) when x^-1*d small
                * /
                x_inv_cached = x_inv_cached*(1. - x_inv_cached*x_minus_x_old);
                maxTimeUntilNextDiv--;
            }
            else
            {
                x_inv_cached = 1. / ws;
                maxTimeUntilNextDiv = MAX_TIME_BETWEEN_EXACT_DIVISIONS;
            }
            x_cached = ws;

            //double true_inv = (1. / ws);
            //if(fabs((double)(x_inv_cached - true_inv)) > 0.001)
            //    throw new grc::GRCException("TEST APPROX DIVISION: Approx. division has failed");

            double xs = xs0*x_inv_cached;
            double ys = ys0*x_inv_cached;
            int ixs = cvRound(xs);//Coordinates of this dest px in source image
            int iys = cvRound(ys);

            // (unsigned) casts mean can do 0 <= x < val
            if ( (unsigned)ixs < uMaxSourceWidth && (unsigned)iys < uMaxSourceHeight ) {
                const uchar* ptr = src + step*iys + ixs*cn;
                ////Copy 3 chars at once by casting to int: this works, could cause access violation if the last 3 bytes are in pos 2 or 3
                //*((int *)(dst + x*cn)) = *((const int *)ptr);
                //for ( k = 0; k < cn; k++ ) {
                //    dst[x*cn+k] = ptr[k];
               // }
                //expand loop
                dst[x*cn] = ptr[0];
                dst[x*cn+1] = ptr[1];
                dst[x*cn+2] = ptr[2];
            }
        }
    } return CV_OK;
}*/

/* Bilinear version with faster division */
// isn't fully optimised (no bounding box, repeated arithmatic, division needs checked)--see above
/***static CvStatus icvWarpPerspective_Bilinear_8u_CnR( const uchar* src, int step, CvSize ssize, uchar* dst, int dststep, CvSize dsize, const double* matrix, int cn, const uchar* fillval ) {
    int x, y, k;
    double A11 = (double)matrix[0], A12 = (double)matrix[1], A13 = (double)matrix[2];
    double A21 = (double)matrix[3], A22 = (double)matrix[4], A23 = (double)matrix[5];
    double A31 = (double)matrix[6], A32 = (double)matrix[7], A33 = (double)matrix[8];
    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);

    unsigned uMaxSourceX = (unsigned)(ssize.width - 1);
    unsigned uMaxSourceY = (unsigned)(ssize.height - 1);
    unsigned uMaxSourceWidth = (unsigned)(ssize.width);
    unsigned uMaxSourceHeight = (unsigned)(ssize.height);
    unsigned uMaxSourceWidthPlus1 = (unsigned)(ssize.width+1);
    unsigned uMaxSourceHeightPlus1 = (unsigned)(ssize.height+1);
    int maxSourceX = (ssize.width - 1);
    int maxSourceY = (ssize.height - 1);

    double x_cached=0., x_inv_cached=0.; //Remember what x and 1/x were last time

    const double DELTA_MAX = 1e-2; // These vals work well with an error usually much less than 0.001
    const int MAX_TIME_BETWEEN_EXACT_DIVISIONS = 5000;

    int maxTimeUntilNextDiv=MAX_TIME_BETWEEN_EXACT_DIVISIONS;

    for ( y = 0; y < dsize.height; y++, dst += dststep ) {
        double xs0 = A12*y + A13;
        double ys0 = A22*y + A23;
        double ws = A32*y + A33;
        for ( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) {
            //Do division:
            //if we're close enough then binomial expand around previous soln
            double x_minus_x_old = ws-x_cached;
            if(maxTimeUntilNextDiv > 0 && x_minus_x_old < DELTA_MAX && x_minus_x_old > -DELTA_MAX)
            {
                x_inv_cached = x_inv_cached*(1. - x_inv_cached*x_minus_x_old);
                maxTimeUntilNextDiv--;
            }
            else
            {
                x_inv_cached = 1. / ws;
                maxTimeUntilNextDiv = MAX_TIME_BETWEEN_EXACT_DIVISIONS;
            }
            x_cached = ws;

            double xs = xs0*x_inv_cached;
            double ys = ys0*x_inv_cached;

            int ixs = cvRound(xs);
            int iys = cvRound(ys);
            double a = xs - ixs;
            double b = ys - iys;
            double p0, p1;
            // (unsigned) casts mean can do 0 <= x < val
            if ( (unsigned)ixs < uMaxSourceX && (unsigned)iys < uMaxSourceY ) {
                const uchar* ptr = src + step*iys + ixs*cn;
                for ( k = 0; k < cn; k++ ) {
                    p0 = toFloat(ptr[k]) + a * (toFloat(ptr[k+cn]) - toFloat(ptr[k]));
                    p1 = toFloat(ptr[k+step]) + a * (toFloat(ptr[k+cn+step]) - toFloat(ptr[k+step]));
                    dst[x*cn+k] = (uchar)cvRound(p0 + b*(p1 - p0));
                }
            } else if ( (unsigned)(ixs+1) < uMaxSourceWidthPlus1 && (unsigned)(iys+1) < uMaxSourceHeightPlus1) {
                int x0 = ((unsigned)(ixs) < uMaxSourceWidth ? (ixs) : (ixs) < 0 ? 0 : maxSourceX);
                int y0 = ((unsigned)(iys) < uMaxSourceHeight ? (iys) : (iys) < 0 ? 0 : maxSourceY);
                int x1 = ((unsigned)(ixs + 1) < uMaxSourceWidth ? (ixs + 1) : (ixs + 1) < 0 ? 0 : maxSourceX);
                int y1 = ((unsigned)(iys + 1) < uMaxSourceHeight ? (iys + 1) : (iys + 1) < 0 ? 0 : maxSourceY);

                //These are the pointers to the 4 corners for lin. interp.
                const uchar* ptr0, *ptr1, *ptr2, *ptr3;
                ptr0 = src + y0*step + x0*cn;//Repeated arith.
                ptr1 = src + y0*step + x1*cn;
                ptr2 = src + y1*step + x0*cn;
                ptr3 = src + y1*step + x1*cn;
                for ( k = 0; k < cn; k++ )
                {
                    //bilinear interpolation
                    p0 = toFloat(ptr0[k]) + a * (toFloat(ptr1[k]) - toFloat(ptr0[k]));
                    p1 = toFloat(ptr2[k]) + a * (toFloat(ptr3[k]) - toFloat(ptr2[k]));
                    dst[x*cn+k] = (uchar)cvRound(p0 + b*(p1 - p0));
                }
            }// else if ( fillval ) for ( k = 0; k < cn; k++ ) dst[x*cn+k] = fillval[k];
        }
    } return CV_OK;
}***/

/* Bilinear version with faster division */
/*
static CvStatus  icvWarpPerspective_Bilinear_8u_CnR( const uchar* src, int step, CvSize ssize, uchar* dst, int dststep, CvSize dsize, const double* matrix, int cn, const uchar* fillval ) {
    int x, y, k;
    double A11 = (double)matrix[0], A12 = (double)matrix[1], A13 = (double)matrix[2];
    double A21 = (double)matrix[3], A22 = (double)matrix[4], A23 = (double)matrix[5];
    double A31 = (double)matrix[6], A32 = (double)matrix[7], A33 = (double)matrix[8];
    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);
    for ( y = 0; y < dsize.height; y++, dst += dststep ) {
        double xs0 = A12*y + A13;
        double ys0 = A22*y + A23;
        double ws = A32*y + A33;
        for ( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) {
            double inv_ws = invertApproxCache(ws); //1./ws; //ws usually 1.0000xxxx
            double xs = xs0*inv_ws;
            double ys = ys0*inv_ws;
            int ixs = cvFloor(xs);
            int iys = cvFloor(ys);
            double a = xs - ixs;
            double b = ys - iys;
            double p0, p1;
            if ( (unsigned)ixs < (unsigned)(ssize.width - 1) && (unsigned)iys < (unsigned)(ssize.height - 1) ) { //Repeated ssize.height - 1 arithmatic
                const uchar* ptr = src + step*iys + ixs*cn;
                for ( k = 0; k < cn; k++ ) {
                    p0 = toFloat(ptr[k]) + a * (toFloat(ptr[k+cn]) - toFloat(ptr[k]));
                    p1 = toFloat(ptr[k+step]) + a * (toFloat(ptr[k+cn+step]) - toFloat(ptr[k+step]));
                    dst[x*cn+k] = (uchar)cvRound(p0 + b*(p1 - p0));
                }
            } else if ( (unsigned)(ixs+1) < (unsigned)(ssize.width+1) && (unsigned)(iys+1) < (unsigned)(ssize.height+1)) {
                int x0 = ((unsigned)(ixs) < (unsigned)ssize.width ? (ixs) : (ixs) < 0 ? 0 : ssize.width - 1);
                int y0 = ((unsigned)(iys) < (unsigned)ssize.height ? (iys) : (iys) < 0 ? 0 : ssize.height - 1);
                int x1 = ((unsigned)(ixs + 1) < (unsigned)ssize.width ? (ixs + 1) : (ixs + 1) < 0 ? 0 : ssize.width - 1);
                int y1 = ((unsigned)(iys + 1) < (unsigned)ssize.height ? (iys + 1) : (iys + 1) < 0 ? 0 : ssize.height - 1);

                //These are the pointers to the 4 corners for lin. interp.
                const uchar* ptr0, *ptr1, *ptr2, *ptr3;
                ptr0 = src + y0*step + x0*cn;//Repeated arith.
                ptr1 = src + y0*step + x1*cn;
                ptr2 = src + y1*step + x0*cn;
                ptr3 = src + y1*step + x1*cn;
                for ( k = 0; k < cn; k++ ) {
                    //bilinear interpolation
                    p0 = toFloat(ptr0[k]) + a * (toFloat(ptr1[k]) - toFloat(ptr0[k]));
                    p1 = toFloat(ptr2[k]) + a * (toFloat(ptr3[k]) - toFloat(ptr2[k]));
                    dst[x*cn+k] = (uchar)cvRound(p0 + b*(p1 - p0));
                }
            } else if ( fillval ) for ( k = 0; k < cn; k++ ) dst[x*cn+k] = fillval[k];
        }
    } return CV_OK;
}
/ *
static CvStatus  icvWarpPerspective_Bilinear_8u_CnR( const uchar* src, int step, CvSize ssize, uchar* dst, int dststep, CvSize dsize, const double* matrix, int cn, const uchar* fillval ) {
    int x, y, k;
    double A11 = (double)matrix[0], A12 = (double)matrix[1], A13 = (double)matrix[2];
    double A21 = (double)matrix[3], A22 = (double)matrix[4], A23 = (double)matrix[5];
    double A31 = (double)matrix[6], A32 = (double)matrix[7], A33 = (double)matrix[8];
    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);
    for ( y = 0; y < dsize.height; y++, dst += dststep ) {
        double xs0 = A12*y + A13;
        double ys0 = A22*y + A23;
        double ws = A32*y + A33;
        for ( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) {
            double inv_ws = 1./ws; //ws usually 1.0000x
            double xs = xs0*inv_ws;
            double ys = ys0*inv_ws;
            int ixs = cvFloor(xs);
            int iys = cvFloor(ys);
            double a = xs - ixs;
            double b = ys - iys;
            double p0, p1;
            if ( (unsigned)ixs < (unsigned)(ssize.width - 1) && (unsigned)iys < (unsigned)(ssize.height - 1) ) { //Repeated ssize.height - 1 arithmatic
                const uchar* ptr = src + step*iys + ixs*cn;
                for ( k = 0; k < cn; k++ ) {
                    p0 = toFloat(ptr[k]) + a * (toFloat(ptr[k+cn]) - toFloat(ptr[k]));
                    p1 = toFloat(ptr[k+step]) + a * (toFloat(ptr[k+cn+step]) - toFloat(ptr[k+step]));
                    dst[x*cn+k] = (uchar)cvRound(p0 + b*(p1 - p0));
                }
            } else if ( (unsigned)(ixs+1) < (unsigned)(ssize.width+1) && (unsigned)(iys+1) < (unsigned)(ssize.height+1)) {
                int x0 = ((unsigned)(ixs) < (unsigned)ssize.width ? (ixs) : (ixs) < 0 ? 0 : ssize.width - 1);
                int y0 = ((unsigned)(iys) < (unsigned)ssize.height ? (iys) : (iys) < 0 ? 0 : ssize.height - 1);
                int x1 = ((unsigned)(ixs + 1) < (unsigned)ssize.width ? (ixs + 1) : (ixs + 1) < 0 ? 0 : ssize.width - 1);
                int y1 = ((unsigned)(iys + 1) < (unsigned)ssize.height ? (iys + 1) : (iys + 1) < 0 ? 0 : ssize.height - 1);

                //These are the pointers to the 4 corners for lin. interp.
                const uchar* ptr0, *ptr1, *ptr2, *ptr3;
                ptr0 = src + y0*step + x0*cn;//Repeated arith.
                ptr1 = src + y0*step + x1*cn;
                ptr2 = src + y1*step + x0*cn;
                ptr3 = src + y1*step + x1*cn;
                for ( k = 0; k < cn; k++ ) {
                    //bilinear interpolation
                    p0 = toFloat(ptr0[k]) + a * (toFloat(ptr1[k]) - toFloat(ptr0[k]));
                    p1 = toFloat(ptr2[k]) + a * (toFloat(ptr3[k]) - toFloat(ptr2[k]));
                    dst[x*cn+k] = (uchar)cvRound(p0 + b*(p1 - p0));
                }
            } else if ( fillval ) for ( k = 0; k < cn; k++ ) dst[x*cn+k] = fillval[k];
        }
    } return CV_OK;
}
static CvStatus  icvWarpPerspective_Bilinear_16u_CnR( const ushort* src, int step, CvSize ssize, ushort* dst, int dststep, CvSize dsize, const double* matrix, int cn, const ushort* fillval ) {
    int x, y, k;
    float A11 = (float)matrix[0], A12 = (float)matrix[1], A13 = (float)matrix[2];
    float A21 = (float)matrix[3], A22 = (float)matrix[4], A23 = (float)matrix[5];
    float A31 = (float)matrix[6], A32 = (float)matrix[7], A33 = (float)matrix[8];
    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);
    for ( y = 0; y < dsize.height; y++, dst += dststep ) {
        float xs0 = A12*y + A13;
        float ys0 = A22*y + A23;
        float ws = A32*y + A33;
        for ( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) {
            float inv_ws = 1.f/ws;
            float xs = xs0*inv_ws;
            float ys = ys0*inv_ws;
            int ixs = cvFloor(xs);
            int iys = cvFloor(ys);
            float a = xs - ixs;
            float b = ys - iys;
            float p0, p1;
            if ( (unsigned)ixs < (unsigned)(ssize.width - 1) && (unsigned)iys < (unsigned)(ssize.height - 1) ) {
                const ushort* ptr = src + step*iys + ixs*cn;
                for ( k = 0; k < cn; k++ ) {
                    p0 = (ptr[k]) + a * ((ptr[k+cn]) - (ptr[k]));
                    p1 = (ptr[k+step]) + a * ((ptr[k+cn+step]) - (ptr[k+step]));
                    dst[x*cn+k] = (ushort)cvRound(p0 + b*(p1 - p0));
                }
            } else if ( (unsigned)(ixs+1) < (unsigned)(ssize.width+1) && (unsigned)(iys+1) < (unsigned)(ssize.height+1)) {
                int x0 = ((unsigned)(ixs) < (unsigned)ssize.width ? (ixs) : (ixs) < 0 ? 0 : ssize.width - 1);
                int y0 = ((unsigned)(iys) < (unsigned)ssize.height ? (iys) : (iys) < 0 ? 0 : ssize.height - 1);
                int x1 = ((unsigned)(ixs + 1) < (unsigned)ssize.width ? (ixs + 1) : (ixs + 1) < 0 ? 0 : ssize.width - 1);
                int y1 = ((unsigned)(iys + 1) < (unsigned)ssize.height ? (iys + 1) : (iys + 1) < 0 ? 0 : ssize.height - 1);
                const ushort* ptr0, *ptr1, *ptr2, *ptr3;
                ptr0 = src + y0*step + x0*cn;
                ptr1 = src + y0*step + x1*cn;
                ptr2 = src + y1*step + x0*cn;
                ptr3 = src + y1*step + x1*cn;
                for ( k = 0; k < cn; k++ ) {
                    p0 = (ptr0[k]) + a * ((ptr1[k]) - (ptr0[k]));
                    p1 = (ptr2[k]) + a * ((ptr3[k]) - (ptr2[k]));
                    dst[x*cn+k] = (ushort)cvRound(p0 + b*(p1 - p0));
                }
            } else if ( fillval ) for ( k = 0; k < cn; k++ ) dst[x*cn+k] = fillval[k];
        }
    } return CV_OK;
}
static CvStatus  icvWarpPerspective_Bilinear_32f_CnR( const float* src, int step, CvSize ssize, float* dst, int dststep, CvSize dsize, const double* matrix, int cn, const float* fillval ) {
    int x, y, k;
    float A11 = (float)matrix[0], A12 = (float)matrix[1], A13 = (float)matrix[2];
    float A21 = (float)matrix[3], A22 = (float)matrix[4], A23 = (float)matrix[5];
    float A31 = (float)matrix[6], A32 = (float)matrix[7], A33 = (float)matrix[8];
    step /= sizeof(src[0]);
    dststep /= sizeof(dst[0]);
    for ( y = 0; y < dsize.height; y++, dst += dststep ) {
        float xs0 = A12*y + A13;
        float ys0 = A22*y + A23;
        float ws = A32*y + A33;
        for ( x = 0; x < dsize.width; x++, xs0 += A11, ys0 += A21, ws += A31 ) {
            float inv_ws = 1.f/ws;
            float xs = xs0*inv_ws;
            float ys = ys0*inv_ws;
            int ixs = cvFloor(xs);
            int iys = cvFloor(ys);
            float a = xs - ixs;
            float b = ys - iys;
            float p0, p1;
            if ( (unsigned)ixs < (unsigned)(ssize.width - 1) && (unsigned)iys < (unsigned)(ssize.height - 1) ) {
                const float* ptr = src + step*iys + ixs*cn;
                for ( k = 0; k < cn; k++ ) {
                    p0 = (ptr[k]) + a * ((ptr[k+cn]) - (ptr[k]));
                    p1 = (ptr[k+step]) + a * ((ptr[k+cn+step]) - (ptr[k+step]));
                    dst[x*cn+k] = (float)(p0 + b*(p1 - p0));
                }
            } else if ( (unsigned)(ixs+1) < (unsigned)(ssize.width+1) && (unsigned)(iys+1) < (unsigned)(ssize.height+1)) {
                int x0 = ((unsigned)(ixs) < (unsigned)ssize.width ? (ixs) : (ixs) < 0 ? 0 : ssize.width - 1);
                int y0 = ((unsigned)(iys) < (unsigned)ssize.height ? (iys) : (iys) < 0 ? 0 : ssize.height - 1);
                int x1 = ((unsigned)(ixs + 1) < (unsigned)ssize.width ? (ixs + 1) : (ixs + 1) < 0 ? 0 : ssize.width - 1);
                int y1 = ((unsigned)(iys + 1) < (unsigned)ssize.height ? (iys + 1) : (iys + 1) < 0 ? 0 : ssize.height - 1);
                const float* ptr0, *ptr1, *ptr2, *ptr3;
                ptr0 = src + y0*step + x0*cn;
                ptr1 = src + y0*step + x1*cn;
                ptr2 = src + y1*step + x0*cn;
                ptr3 = src + y1*step + x1*cn;
                for ( k = 0; k < cn; k++ ) {
                    p0 = (ptr0[k]) + a * ((ptr1[k]) - (ptr0[k]));
                    p1 = (ptr2[k]) + a * ((ptr3[k]) - (ptr2[k]));
                    dst[x*cn+k] = (float)(p0 + b*(p1 - p0));
                }
            } else if ( fillval ) for ( k = 0; k < cn; k++ ) dst[x*cn+k] = fillval[k];
        }
    } return CV_OK;
}


typedef CvStatus (CV_STDCALL * CvWarpPerspectiveFunc)(
    const void* src, int srcstep, CvSize ssize,
    void* dst, int dststep, CvSize dsize,
    const double* matrix, int cn, const void* fillval );

/ *static void icvInitWarpPerspectiveTab( CvFuncTable* bilin_tab, int method, int channels  )
{
    if(method == CV_INTER_LINEAR)
    {
        if(channels == 1)
        {
            //std::cout << "WARNING: linear interpolation with 1 channel not implemented properly\n";
            bilin_tab->fn_2d[CV_8U] = icvWarpPerspective_Bilinear_8u_CnR;
        }
        else if(channels == 3)
            bilin_tab->fn_2d[CV_8U] = (void*)icvWarpPerspective_Bilinear_8u_CnR;
        else
            bilin_tab->fn_2d[CV_8U] = 0;
    }
    else
    {
        if(channels == 3)
            bilin_tab->fn_2d[CV_8U] = (void*)icvWarpPerspective_NN_BB_8u_C3;
        else if(channels == 1)
            bilin_tab->fn_2d[CV_8U] = (void*)icvWarpPerspective_NN_BB_8u_C1;
        else bilin_tab->fn_2d[CV_8U] = 0;
    }

    bilin_tab->fn_2d[CV_16U] = 0;//(void*)icvWarpPerspective_Bilinear_16u_CnR;
    bilin_tab->fn_2d[CV_32F] = 0;//(void*)icvWarpPerspective_Bilinear_32f_CnR;
}*/


/////////////////////////// IPP warpperspective functions ////////////////////////////////

/*icvWarpPerspectiveBack_8u_C1R_t icvWarpPerspectiveBack_8u_C1R_p = 0;
icvWarpPerspectiveBack_8u_C3R_t icvWarpPerspectiveBack_8u_C3R_p = 0;
icvWarpPerspectiveBack_8u_C4R_t icvWarpPerspectiveBack_8u_C4R_p = 0;
icvWarpPerspectiveBack_32f_C1R_t icvWarpPerspectiveBack_32f_C1R_p = 0;
icvWarpPerspectiveBack_32f_C3R_t icvWarpPerspectiveBack_32f_C3R_p = 0;
icvWarpPerspectiveBack_32f_C4R_t icvWarpPerspectiveBack_32f_C4R_p = 0;

icvWarpPerspective_8u_C1R_t icvWarpPerspective_8u_C1R_p = 0;
icvWarpPerspective_8u_C3R_t icvWarpPerspective_8u_C3R_p = 0;
icvWarpPerspective_8u_C4R_t icvWarpPerspective_8u_C4R_p = 0;
icvWarpPerspective_32f_C1R_t icvWarpPerspective_32f_C1R_p = 0;
icvWarpPerspective_32f_C3R_t icvWarpPerspective_32f_C3R_p = 0;
icvWarpPerspective_32f_C4R_t icvWarpPerspective_32f_C4R_p = 0;

typedef CvStatus (CV_STDCALL * CvWarpPerspectiveBackIPPFunc)
( const void* src, CvSize srcsize, int srcstep, CvRect srcroi,
  void* dst, int dststep, CvRect dstroi,
  const double* coeffs, int interpolation );*/

//////////////////////////////////////////////////////////////////////////////////////////

void
WarpPerspective( const CvArr* srcarr, CvArr* dstarr,
                   const CvMat* matrix, int flags )
{
    //CV_FUNCNAME( "WarpPerspective" );

    CvMat srcstub, *src = (CvMat*)srcarr;
    CvMat dststub, *dst = (CvMat*)dstarr;
    int type, depth, channels;
    int method = flags & 3;
    double src_matrix[9], dst_matrix[9];
    //double fillbuf[4];
    CvMat A = cvMat( 3, 3, CV_64F, src_matrix ),
          invA = cvMat( 3, 3, CV_64F, dst_matrix );
    CvSize ssize, dsize;

    void (CV_CDECL *func)(const uchar *,int,CvSize,uchar *,int,CvSize,const double *,const double *) = 0;

    if( method == CV_INTER_AREA )
        THROW("Area interpolation not supported" );

    //Must always init fn tab
    src = cvGetMat( srcarr, &srcstub );
    dst = cvGetMat( dstarr, &dststub );

    type = CV_MAT_TYPE(src->type);
    depth = CV_MAT_DEPTH(type);
    channels = CV_MAT_CN(type);
    //int channels = CV_MAT_CN(CV_MAT_TYPE(src->type));
    //icvInitWarpPerspectiveTab( &bilin_tab, method, cn );

    if( !CV_ARE_TYPES_EQ( src, dst ))
        THROW( "CV_StsUnmatchedFormats" );

    if( !CV_IS_MAT(matrix) || CV_MAT_CN(matrix->type) != 1 ||
        CV_MAT_DEPTH(matrix->type) < CV_32F || matrix->rows != 3 || matrix->cols != 3 )
        THROW( "Transformation matrix should be 3x3 floating-point single-channel matrix" );

    if( flags & CV_WARP_INVERSE_MAP )
        cvConvertScale( matrix, &invA );
    else
    {
        cvConvertScale( matrix, &A );
        cvInvert( &A, &invA, CV_SVD );
    }

    if( channels > 4 )
        THROW( "CV_BadNumChannels" );

    if( flags & CV_WARP_FILL_OUTLIERS )
        THROW( "CV_StsBadArg Outlier fiulling no longer implemented" );

    ssize = cvSize(src->width, src->height);// cvGetMatSize(src);
    dsize = cvSize(dst->width, dst->height);//cvGetMatSize(dst);

    /*if( icvWarpPerspectiveBack_8u_C1R_p )
    {
        CvWarpPerspectiveBackIPPFunc ipp_func =
            type == CV_8UC1 ? icvWarpPerspectiveBack_8u_C1R_p :
            type == CV_8UC3 ? icvWarpPerspectiveBack_8u_C3R_p :
            type == CV_8UC4 ? icvWarpPerspectiveBack_8u_C4R_p :
            type == CV_32FC1 ? icvWarpPerspectiveBack_32f_C1R_p :
            type == CV_32FC3 ? icvWarpPerspectiveBack_32f_C3R_p :
            type == CV_32FC4 ? icvWarpPerspectiveBack_32f_C4R_p : 0;

        if( ipp_func && CV_INTER_NN <= method && method <= CV_INTER_AREA &&
            MIN(ssize.width,ssize.height) >= 4 && MIN(dsize.width,dsize.height) >= 4 )
        {
            int srcstep = src->step ? src->step : CV_STUB_STEP;
            int dststep = dst->step ? dst->step : CV_STUB_STEP;
            CvStatus status;
            CvRect srcroi = {0, 0, ssize.width, ssize.height};
            CvRect dstroi = {0, 0, dsize.width, dsize.height};

            // this is not the most efficient way to fill outliers
            if( flags & CV_WARP_FILL_OUTLIERS )
                cvSet( dst, fillval );

            status = ipp_func( src->data.ptr, ssize, srcstep, srcroi,
                               dst->data.ptr, dststep, dstroi,
                               invA.data.db, 1 << method );
            if( status >= 0 )
                EXIT;

            ipp_func = type == CV_8UC1 ? icvWarpPerspective_8u_C1R_p :
                type == CV_8UC3 ? icvWarpPerspective_8u_C3R_p :
                type == CV_8UC4 ? icvWarpPerspective_8u_C4R_p :
                type == CV_32FC1 ? icvWarpPerspective_32f_C1R_p :
                type == CV_32FC3 ? icvWarpPerspective_32f_C3R_p :
                type == CV_32FC4 ? icvWarpPerspective_32f_C4R_p : 0;

            if( ipp_func )
            {
                if( flags & CV_WARP_INVERSE_MAP )
                    cvInvert( &invA, &A, CV_SVD );

                status = ipp_func( src->data.ptr, ssize, srcstep, srcroi,
                               dst->data.ptr, dststep, dstroi,
                               A.data.db, 1 << method );
                if( status >= 0 )
                    EXIT;
            }
        }
    }

    cvScalarToRawData( &fillval, fillbuf, CV_MAT_TYPE(src->type), 0 );*/
    if(depth != CV_8U)
        THROW( "CV_StsUnsupportedFormat Only 8U images supported" );

    //func = 0;

    if(channels == 1)
    {
        if(method == CV_INTER_NN)
            func = CFastWarp<1, CV_INTER_NN, false>::warp_BB_8u;
        else
            func = CFastWarp<1, CV_INTER_LINEAR, false>::warp_BB_8u;
    }
    else if(channels == 3)
    {
        if(method == CV_INTER_NN)
            func = CFastWarp<3, CV_INTER_NN, false>::warp_BB_8u;
        else
            func = CFastWarp<3, CV_INTER_LINEAR, false>::warp_BB_8u;
    }
    else if(channels == 4)
    {
        if(method == CV_INTER_NN)
            func = CFastWarp<4, CV_INTER_NN, false>::warp_BB_8u;
        else
            func = CFastWarp<4, CV_INTER_LINEAR, false>::warp_BB_8u;
    }

    if( !func )
        THROW( "CV_StsUnsupportedFormat No perspective warp fn" );

     func( src->data.ptr, src->step, ssize, dst->data.ptr,
                     dst->step, dsize, dst_matrix,
                     src_matrix ); //use the now-unused void pointer param to pass in inverse matrix (for bounding box)

}


/* Calculates coefficients of perspective transformation
 * which maps (xi,yi) to (ui,vi), (i=1,2,3,4):
 *
 *      c00*xi + c01*yi + c02
 * ui = ---------------------
 *      c20*xi + c21*yi + c22
 *
 *      c10*xi + c11*yi + c12
 * vi = ---------------------
 *      c20*xi + c21*yi + c22
 *
 * Coefficients are calculated by solving linear system:
 * / x0 y0  1  0  0  0 -x0*u0 -y0*u0 \ /c00\ /u0\
 * | x1 y1  1  0  0  0 -x1*u1 -y1*u1 | |c01| |u1|
 * | x2 y2  1  0  0  0 -x2*u2 -y2*u2 | |c02| |u2|
 * | x3 y3  1  0  0  0 -x3*u3 -y3*u3 |.|c10|=|u3|,
 * |  0  0  0 x0 y0  1 -x0*v0 -y0*v0 | |c11| |v0|
 * |  0  0  0 x1 y1  1 -x1*v1 -y1*v1 | |c12| |v1|
 * |  0  0  0 x2 y2  1 -x2*v2 -y2*v2 | |c20| |v2|
 * \  0  0  0 x3 y3  1 -x3*v3 -y3*v3 / \c21/ \v3/
 *
 * where:
 *   cij - matrix coefficients, c22 = 1
 CV_IMPL*/
 CvMat*
cvGetPerspectiveTransform( const CvPoint2D32f* src,
                          const CvPoint2D32f* dst,
                          CvMat* matrix )
{
    //CV_FUNCNAME( "cvGetPerspectiveTransform" );

    double a[8][8];
    double b[8], x[9];

    CvMat A = cvMat( 8, 8, CV_64FC1, a );
    CvMat B = cvMat( 8, 1, CV_64FC1, b );
    CvMat X = cvMat( 8, 1, CV_64FC1, x );

    int i;

    if( !src || !dst || !matrix )
        THROW("CV_StsNullPtr" );

    for( i = 0; i < 4; ++i )
    {
        a[i][0] = a[i+4][3] = src[i].x;
        a[i][1] = a[i+4][4] = src[i].y;
        a[i][2] = a[i+4][5] = 1;
        a[i][3] = a[i][4] = a[i][5] =
        a[i+4][0] = a[i+4][1] = a[i+4][2] = 0;
        a[i][6] = -src[i].x*dst[i].x;
        a[i][7] = -src[i].y*dst[i].x;
        a[i+4][6] = -src[i].x*dst[i].y;
        a[i+4][7] = -src[i].y*dst[i].y;
        b[i] = dst[i].x;
        b[i+4] = dst[i].y;
    }

    cvSolve( &A, &B, &X, CV_SVD );
    x[8] = 1;

    X = cvMat( 3, 3, CV_64FC1, x );
    cvConvert( &X, matrix );

    return matrix;
}

/* Calculates coefficients of affine transformation
 * which maps (xi,yi) to (ui,vi), (i=1,2,3):
 *
 * ui = c00*xi + c01*yi + c02
 *
 * vi = c10*xi + c11*yi + c12
 *
 * Coefficients are calculated by solving linear system:
 * / x0 y0  1  0  0  0 \ /c00\ /u0\
 * | x1 y1  1  0  0  0 | |c01| |u1|
 * | x2 y2  1  0  0  0 | |c02| |u2|
 * |  0  0  0 x0 y0  1 | |c10| |v0|
 * |  0  0  0 x1 y1  1 | |c11| |v1|
 * \  0  0  0 x2 y2  1 / |c12| |v2|
 *
 * where:
 *   cij - matrix coefficients
 CV_IMPL*/
 CvMat*
cvGetAffineTransform( const CvPoint2D32f * src, const CvPoint2D32f * dst, CvMat * map_matrix )
{
    //CV_FUNCNAME( "cvGetAffineTransform" );

    CvMat mA, mX, mB;
    double A[6*6];
    double B[6];
	double x[6];
    int i;

    cvInitMatHeader(&mA, 6, 6, CV_64F, A);
    cvInitMatHeader(&mB, 6, 1, CV_64F, B);
	cvInitMatHeader(&mX, 6, 1, CV_64F, x);

	if( !src || !dst || !map_matrix )
        THROW("CV_StsNullPtr" );

    for( i = 0; i < 3; i++ )
    {
        int j = i*12;
        int k = i*12+6;
        A[j] = A[k+3] = src[i].x;
        A[j+1] = A[k+4] = src[i].y;
        A[j+2] = A[k+5] = 1;
        A[j+3] = A[j+4] = A[j+5] = 0;
        A[k] = A[k+1] = A[k+2] = 0;
        B[i*2] = dst[i].x;
        B[i*2+1] = dst[i].y;
    }
    cvSolve(&mA, &mB, &mX);

    mX = cvMat( 2, 3, CV_64FC1, x );
	cvConvert( &mX, map_matrix );

    return map_matrix;
}

/* End of file. */
