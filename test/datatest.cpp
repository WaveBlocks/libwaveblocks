#include <string>
#include <iostream>
#include <vector>

#include "gtest/gtest.h"

#include "H5Cpp.h"

#define abstol 1e-6

using namespace H5;

int global_argc;
char** global_argv;

namespace test {

    class TestHDF : public ::testing::Test
    {
        protected:
        TestHDF()
        {
            if(global_argc<=2)
            {
            std::cout<<"Not enough arguments given!\n";
            std::abort();
            }
            else
            {
                std::string filename_1 = global_argv[1];
                std::string filename_2 = global_argv[2];
                if(filename_1.find("py")!=std::string::npos)
                {
                    pyname = filename_1;
                    cppname = filename_2;
                }
                else if(filename_1.find("cpp")!=std::string::npos)
                {
                    cppname = filename_1;
                    pyname = filename_2;
                }
                else
                {
                    std::cout<<"Failed to open. In argv *cpp.hdf5 and *py.hdf5 file needed";
                    std::abort();
                }
                cppfile=H5File(cppname,H5F_ACC_RDONLY);
                pyfile=H5File(pyname,H5F_ACC_RDONLY);
            }
        }
        virtual ~TestHDF()
        {
            cppfile.close();
            pyfile.close();
        }
        void SetUp()
        {
            Group block_0 = cppfile.openGroup("/datablock_0");
            Attribute apacket = block_0.openAttribute("packet");
            Attribute aenergy = block_0.openAttribute("energy");
            Attribute anorm = block_0.openAttribute("norm");
            Attribute adt_cpp = block_0.openAttribute("dt");
            apacket.read(PredType::NATIVE_INT,&bool_packet);
            aenergy.read(PredType::NATIVE_INT,&bool_energies);
            anorm.read(PredType::NATIVE_INT,&bool_norm);
            adt_cpp.read(PredType::NATIVE_DOUBLE,&dt_cpp);
            apacket.close();
            aenergy.close();
            anorm.close();
            adt_cpp.close();
            block_0.close();

            block_0 = pyfile.openGroup("/datablock_0");
            Attribute adt_py = block_0.openAttribute("ext:dt");
            StrType mystrtype=adt_py.getStrType();
            H5std_string buff_dt_py="";
            adt_py.read(mystrtype,buff_dt_py);
            dt_py=std::strtod(buff_dt_py.c_str(),NULL);
            mystrtype.close();
            adt_py.close();
            block_0.close();

        }
        void TearDown()
        {

        }
        H5File cppfile;
        H5File pyfile;
        H5std_string cppname;
        H5std_string pyname;
        int bool_packet;
        int bool_energies;
        int bool_norm;
        double dt_cpp;
        double dt_py;
    };

    TEST_F(TestHDF,Testpacket)
    {
        ASSERT_TRUE(bool_packet) << "Packet data missing cannot compare. Abort!";
        EXPECT_EQ(dt_cpp,dt_py);

        DataSet its1 = cppfile.openDataSet("/datablock_0/wavepacket/timegrid");
        DataSet its2 = pyfile.openDataSet("/datablock_0/wavepacket/timegrid");
        DataSpace itspace1 = its1.getSpace();
        DataSpace itspace2 = its2.getSpace();

        int itrank1 = itspace1.getSimpleExtentNdims();
        int itrank2 = itspace1.getSimpleExtentNdims();

        hsize_t* itdim1 = new hsize_t[itrank1];
        hsize_t* itdim2 = new hsize_t[itrank2];

        itspace1.getSimpleExtentDims(itdim1);
        itspace2.getSimpleExtentDims(itdim2);

        DataType type1 = its1.getDataType();
        DataType type2 = its1.getDataType();

//        IntType it1 = its1.getIntType();
//        IntType it2 = its2.getIntType();

//        int* index_array_cpp = new int[itdim1[0]];
//        int* index_array_py = new int[itdim2[0]];

//        hsize_t* itelem=new hsize_t[itrank1];
//        itelem[0]=itdim1[0];

//        DataSpace itelemspace(itrank1,itelem);
//        hsize_t* itstart1=new hsize_t[itrank1];
//        itstart1[0]=0;
//        hsize_t* itcount1= new hsize_t[itrank1];
//        itcount1[0]=itdim1[0];
//        hsize_t* itstride1=new hsize_t[itrank1];
//        itstride1[0]=1;
//        hsize_t* itblock1=new hsize_t[itrank1];
//        itblock1[0]=1;

//        itelemspace.selectHyperslab(H5S_SELECT_SET,itcount1,itstart1,itstride1,itblock1);
//        itspace1.selectHyperslab(H5S_SELECT_SET,itcount1,itstart1,itstride1,itblock1);
//        itspace2.selectHyperslab(H5S_SELECT_SET,itcount1,itstart1,itstride1,itblock1);

//        its1.read(index_array_cpp,it1,itelemspace,itspace1);
//        its2.read(index_array_py,it2,itelemspace,itspace2);

//        std::vector<int*> res;
//        int acc[2]={-1,-1};
//        for(unsigned int k=0;k<itdim1[0];++k)
//        {
//            int index_i_cpp=index_array_cpp[k];
//            for(unsigned int p=0;p<itdim2[0];++p)
//            {
//                int index_j_py=index_array_py[p];
//                int temp = index_j_py*(dt_py/dt_cpp);
//                if(index_i_cpp==temp)
//                {
//                    acc[0]=index_i_cpp;
//                    acc[1]=index_j_py;
//                    res.push_back(acc);
//                }
//            }
//        }
//        if(res.empty())
//        {
//            res.push_back(acc);
//        }
//        delete [] itdim1;
//        delete [] itdim2;
//        delete [] index_array_cpp;
//        delete [] index_array_py;
//        delete [] itstart1;
//        delete [] itcount1;
//        delete [] itstride1;
//        delete [] itblock1;

//        itelemspace.close();
//        itspace1.close();
//        itspace2.close();
//        it1.close();
//        it2.close();
        its1.close();
        its2.close();

    }

    TEST_F(TestHDF,Testenergies)
    {
        ASSERT_TRUE(bool_energies) << "Energy data missing cannot compare. Abort!";
    }

    TEST_F(TestHDF,Testnorm)
    {
        ASSERT_TRUE(bool_norm) << "Norm data missing cannot compare. Abort!";
    }

}
int main(int argc,char* argv[])
{
    global_argc=argc;
    global_argv=argv;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

