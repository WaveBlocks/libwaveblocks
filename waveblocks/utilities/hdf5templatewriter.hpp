#pragma once

//need H5 Cpp header
#include "H5Cpp.h"
//using strings
#include <string>
#include <list>
#include <Eigen/Core>

#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"

namespace waveblocks
{
    namespace utilities
    {

    //using namespaces for convenience
    using namespace H5;

    //define complex struct used to write back
    struct ctype{ //our complex datatype
        double real=0.;
        double imag=0.;
    }instanceof;

    //template Dimension for writing
    template<int D>
    class hdf5writertemplate
    {
        //declare friend for easy use
        friend struct ctype;
    public:
        hdf5writertemplate(std::string name):filename_(name),mytype_(sizeof(instanceof)),file_(filename_,H5F_ACC_TRUNC)
        {
            //fixed writing type for easy overlap with python interface
            mytype_.insertMember( "r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
            mytype_.insertMember( "i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);

            //register datafields which should be written

        }

    private:
        H5std_string filename_;
        CompType mytype_;
        H5File file_;
        //std::list<std::string> wlist;
    };


    }
}
