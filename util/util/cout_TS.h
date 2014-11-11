#pragma once
#ifndef COUT_TS_H_
#define COUT_TS_H_

class CTSOut : public std::ostringstream
{
public:
    ~CTSOut()
    {
        std::cout << std::ostringstream::str() << endl;
    }
};
#endif //COUT_TS_H_;
