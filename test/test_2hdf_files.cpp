#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>
#include <cassert>

//include GTest
#include "gtest/gtest.h"

#include "H5Cpp.h"

#define abstol 1e-6

//global variables for argv
int global_argc;
char** global_argv;
//namespace for h5
using namespace H5;

/**
 * @brief Testclass needed for GTest interface
 *
 * This class need to be derived from ::testing::Test in gtest.h and has to define two functions.
 * SetUp and BreakDown. SetUp for initialization of variables/parameters for different testfissures.
 * BreakDown for cleanup of variables/parameters if needed.
 */
class Test2files: public ::testing::Test
{
    public:

    /**
     * @brief SetUp for TEST_F
     *
     * Initialize variables/parameters
     */
    void SetUp()
    {
    //read filenames
    assert(global_argc==3);
    //char* this_file=global_argv[0];
    std::string filename_1 = global_argv[1];
    std::string filename_2 = global_argv[2];
    if(filename_1.find("py")!=std::string::npos)
    {
        pythonname = filename_1;
        cppname = filename_2;
    }
    else if(filename_1.find("cpp")!=std::string::npos)
    {
        cppname = filename_1;
        pythonname = filename_2;
    }
    else
    {
        FAIL()<< "Failed to compare. In argv *cpp.hdf5 and *py.hdf5 file needed";
    }
    cppfile=H5File(cppname,H5F_ACC_RDONLY);
    pyfile=H5File(pythonname,H5F_ACC_RDONLY);

    datablock_cpp=std::make_shared<Group>(cppfile.openGroup(datablock_group));
    datablock_py=std::make_shared<Group>(pyfile.openGroup(datablock_group));
    apacket_cpp=datablock_cpp->openAttribute("packet");
    aenergy_cpp=datablock_cpp->openAttribute("energy");
    anorm_cpp=datablock_cpp->openAttribute("norm");
    adt_cpp=datablock_cpp->openAttribute("dt");
    adt_py=datablock_py->openAttribute("ext:dt");
    apacket_cpp.read(PredType::NATIVE_INT,&bool_packet);
    apacket_cpp.close();
    aenergy_cpp.read(PredType::NATIVE_INT,&bool_energy);
    aenergy_cpp.close();
    anorm_cpp.read(PredType::NATIVE_INT,&bool_norm);
    anorm_cpp.close();
    adt_cpp.read(PredType::NATIVE_DOUBLE,&dt_cpp);
    adt_cpp.close();

    StrType mystrtype=adt_py.getStrType();
    adt_py.read(mystrtype,buff_dt_py);
    dt_py=std::strtod(buff_dt_py.c_str(),NULL);
    mystrtype.close();
    adt_py.close();

    mytype=CompType(sizeof(instanceof));
    mytype.insertMember( "r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
    mytype.insertMember( "i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);

    datablock_cpp->close();
    datablock_py->close();

    datasetQpath=datablock_group+wavepacket_group+Pi_group+"/Q";
    datasetPpath=datablock_group+wavepacket_group+Pi_group+"/P";
    datasetqpath=datablock_group+wavepacket_group+Pi_group+"/q";
    datasetppath=datablock_group+wavepacket_group+Pi_group+"/p";
    datasetcpath=datablock_group+wavepacket_group+coeffs_group+"/c_0";
    timepathpacket=datablock_group+wavepacket_group+"/timegrid";
    datasetekinpath=datablock_group+energies_group+"/ekin";
    timepathekin=datablock_group+energies_group+"/timegrid_ekin";
    datasetepotpath=datablock_group+energies_group+"/epot";
    timepathepot=datablock_group+energies_group+"/timegrid_epot";
    datasetnormpath=datablock_group+norm_group+"/norm";
    timepathnorm=datablock_group+norm_group+"/timegrid";
    }
    /**
     * @brief BreakDown for TEST_F
     *
     * Cleanup for variables/parameters. In our case empty.
     */
    void BreakDown()
    {
        cppfile.close();
        pyfile.close();
    }
    struct ctype{ //our complex datatype
        double real=0.;
        double imag=0.;
    } instanceof;

    H5std_string cppname;//!<filename of cpp file
    H5File cppfile;//!<instance of cpp file
    H5std_string pythonname;//!<filename of py file
    H5File pyfile;//!<instance of py file
    CompType mytype;//!<predefined type used to write and read data
    H5std_string datablock_group="/datablock_0";//!<name for datablock group. DEFAULT:/datablock_0
    H5std_string energies_group="/energies";//!<name for energies group. DEFAULT:/energies
    H5std_string coeffs_group="/coefficients";//!<name for coefficients group. DEFAULT:/coefficients
    H5std_string norm_group="/norm";//!<name for norm group. DEFAULT:/norms
    H5std_string Pi_group="/Pi";//!<name for Pi group. DEFAULT:/Pi
    H5std_string wavepacket_group="/wavepacket";//!<name for wavepacket group. DEFAULT:/wavepacket
    std::shared_ptr<Group> datablock_cpp;//!<datablockgroup for cppfile identifier
    std::shared_ptr<Group> datablock_py;//!<datablockgroup for pyfile identifier
    H5std_string datasetQpath;//!<string for path to dataset Q
    H5std_string datasetPpath;//!<string for path to dataset P
    H5std_string datasetqpath;//!<string for path to dataset p
    H5std_string datasetppath;//!<string for path to dataset q
    H5std_string datasetcpath;//!<string for path to dataset c_0
    H5std_string datasetekinpath;//!<string for path to dataset ekin
    H5std_string datasetepotpath;//!<string for path to dataset epot
    H5std_string datasetnormpath;//!<string for path to dataset norm
    Attribute apacket_cpp;//!<attribute identifier for packet for cpp file
    Attribute aenergy_cpp;//!<attribute identifier for energy for cpp file
    Attribute anorm_cpp;//!<attribute identifier for norm for cpp file
    Attribute adt_cpp;//!<attribute identifier for dt for cpp file
    Attribute adt_py;//!<attribute identifier for dt for py file
    int bool_packet;//!<bool if packet is written in cpp file
    int bool_energy;//!<bool if energy is written in cpp file
    int bool_norm;//!<bool if norm is written in cpp file
    double dt_cpp;//!<timestep in cpp file
    double dt_py=0.01;//!<timestep in py file
    H5std_string buff_dt_py="";//!<stringbuffer needed for reading Attribute adt_py
    H5std_string timepathpacket;
    H5std_string timepathnorm;
    H5std_string timepathekin;
    H5std_string timepathepot;
};

TEST_F(Test2files,TestdatasetQ)
{
    std::cout<<"**************************0";

    if(bool_packet)
    {
        //expect RANK=3

        //***********************injection
        std::vector<int*> res;

        //expect RANK=1
        DataSet its1 = cppfile.openDataSet(timepathpacket);

        //segfault dataspace close()
        DataSet its2 = pyfile.openDataSet("/datablock_0/wavepacket/timegrid");


        IntType it1 = its1.getIntType();
        IntType it2 = its2.getIntType();

        DataSpace itspace1 = its1.getSpace();
        DataSpace itspace2 = its2.getSpace();

        int itrank1 = itspace1.getSimpleExtentNdims();
        int itrank2 = itspace1.getSimpleExtentNdims();

        hsize_t* itdim1 = new hsize_t[itrank1];
        hsize_t* itdim2 = new hsize_t[itrank2];

        itspace1.getSimpleExtentDims(itdim1);
        itspace2.getSimpleExtentDims(itdim2);

        int* index_array_cpp = new int[itdim1[0]];
        int* index_array_py = new int[itdim2[0]];

        hsize_t* itelem=new hsize_t[itrank1];
        itelem[0]=itdim1[0];

        DataSpace itelemspace(itrank1,itelem);
        hsize_t* itstart1=new hsize_t[itrank1];
        itstart1[0]=0;
        hsize_t* itcount1= new hsize_t[itrank1];
        itcount1[0]=1;
        hsize_t* itstride1=new hsize_t[itrank1];
        itstride1[0]=1;
        hsize_t* itblock1=new hsize_t[itrank1];
        itblock1[0]=1;

        itelemspace.selectHyperslab(H5S_SELECT_SET,itcount1,itstart1,itstride1,itblock1);

        for(unsigned int k=0;k<itdim1[0];++k)
        {
            hsize_t* it_dyn_start=new hsize_t[itrank1];
            it_dyn_start[0]=k;
            itspace1.selectHyperslab(H5S_SELECT_SET,itcount1,it_dyn_start,itstride1,itblock1);
            itspace2.selectHyperslab(H5S_SELECT_SET,itcount1,it_dyn_start,itstride1,itblock1);

            its1.read(&index_array_cpp[k],it1,itelemspace,itspace1);
            its2.read(&index_array_py[k],it2,itelemspace,itspace2);

        }

        int acc[2]={-1,-1};
        for(unsigned int k=0;k<itdim1[0];++k)
        {
            int index_i_cpp=index_array_cpp[k];
            for(unsigned int p=0;p<itdim2[0];++p)
            {
                int index_j_py=index_array_py[p];
                int temp = index_j_py*(dt_py/dt_cpp);
                if(index_i_cpp==temp)
                {
                    acc[0]=index_i_cpp;
                    acc[1]=index_j_py;
                    res.push_back(acc);
                }
            }
        }
        if(res.empty())
        {
            res.push_back(acc);
        }
//        its1.close();
//        its2.close();
        //*******************************************

        DataSet ds1 = cppfile.openDataSet(datasetQpath);
        DataSet ds2 = pyfile.openDataSet(datasetQpath);

        CompType comp1=ds1.getCompType(); //should be equal to mytype
        CompType comp2=ds2.getCompType(); //should be equal to mytype
        ASSERT_EQ(comp1,comp2);
        ASSERT_EQ(comp1,mytype);

        DataSpace dspace1 = ds1.getSpace();
        DataSpace dspace2 = ds2.getSpace();

        int rank1 = dspace1.getSimpleExtentNdims();
        int rank2 = dspace1.getSimpleExtentNdims();
        ASSERT_EQ(rank1,rank2);

        hsize_t* dim1 = new hsize_t[rank1];
        hsize_t* dim2 = new hsize_t[rank2];

        dspace1.getSimpleExtentDims(dim1);
        dspace2.getSimpleExtentDims(dim2);

        ASSERT_EQ(dim1[0],dim2[0]);
        ASSERT_EQ(dim1[1],dim2[1]);
        ASSERT_EQ(dim1[2],dim2[2]);
        int RANK3=3;
        ASSERT_EQ(rank1,RANK3);

        hsize_t* elem=new hsize_t[rank1];
        elem[0]=1;
        elem[1]=dim1[1];
        elem[2]=dim1[2];
        ctype* outdat1 = new ctype[dim1[1]*dim1[2]];
        ctype* outdat2 = new ctype[dim2[1]*dim2[2]];

        DataSpace elemspace(rank1,elem);
        hsize_t* start1=new hsize_t[rank1];
        start1[0]=0;
        start1[1]=0;
        start1[2]=0;
        hsize_t* count1= new hsize_t[rank1];
        count1[0]=1;
        count1[1]=dim1[1];
        count1[2]=dim1[2];
        hsize_t* stride1=new hsize_t[rank1];
        stride1[0]=1;
        stride1[1]=1;
        stride1[2]=1;
        hsize_t* block1=new hsize_t[rank1];
        block1[0]=1;
        block1[1]=1;
        block1[2]=1;

        elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

        if(res.empty())
        {
            FAIL()<<"No matching timepoint.Abort!";
        }
        for(auto index_element:res)
        {
            if(index_element[0]==-1)
            {
                FAIL()<<"No matching timepoint.Abort!";
            }

            hsize_t* start2 = new hsize_t[rank1];
            start2[0]=index_element[0];
            start2[1]=0;
            start2[2]=0;
            hsize_t* count2 = new hsize_t[rank1];
            count2[0]=1;
            count2[1]=dim1[1];
            count2[2]=dim1[2];
            hsize_t* stride2 = new hsize_t[rank1];
            stride2[0]=1;
            stride2[1]=1;
            stride2[2]=1;
            hsize_t* block2 = new hsize_t[rank1];
            block2[0]=1;
            block2[1]=1;
            block2[2]=1;
            hsize_t* start3 = new hsize_t[rank1];
            start3[0]=index_element[1];
            start3[1]=0;
            start3[2]=0;

            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start3, stride2, block2);
            ds1.read(outdat1,comp1,elemspace,dspace1);
            ds2.read(outdat2,comp2,elemspace,dspace2);

            EXPECT_NEAR(outdat1[0].real,outdat2[0].real,abstol);
            EXPECT_NEAR(outdat1[1].real,outdat2[1].real,abstol);
            EXPECT_NEAR(outdat1[2].real,outdat2[2].real,abstol);
            EXPECT_NEAR(outdat1[3].real,outdat2[3].real,abstol);

            EXPECT_NEAR(outdat1[0].imag,outdat2[0].imag,abstol);
            EXPECT_NEAR(outdat1[1].imag,outdat2[1].imag,abstol);
            EXPECT_NEAR(outdat1[2].imag,outdat2[2].imag,abstol);
            EXPECT_NEAR(outdat1[3].imag,outdat2[3].imag,abstol);
        }

//        for(unsigned int k=0;k<dim1[0];++k)
//        {
//            hsize_t* start2 = new hsize_t[rank1];
//            start2[0]=k;
//            start2[1]=0;
//            start2[2]=0;
//            hsize_t* count2 = new hsize_t[rank1];
//            count2[0]=1;
//            count2[1]=dim1[1];
//            count2[2]=dim1[2];
//            hsize_t* stride2 = new hsize_t[rank1];
//            stride2[0]=1;
//            stride2[1]=1;
//            stride2[2]=1;
//            hsize_t* block2 = new hsize_t[rank1];
//            block2[0]=1;
//            block2[1]=1;
//            block2[2]=1;

//            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
//            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
//            ds1.read(outdat1,comp1,elemspace,dspace1);
//            ds2.read(outdat2,comp2,elemspace,dspace2);

//            EXPECT_NEAR(outdat1[0].real,outdat2[0].real,abstol);
//            EXPECT_NEAR(outdat1[1].real,outdat2[1].real,abstol);
//            EXPECT_NEAR(outdat1[2].real,outdat2[2].real,abstol);
//            EXPECT_NEAR(outdat1[3].real,outdat2[3].real,abstol);

//            EXPECT_NEAR(outdat1[0].imag,outdat2[0].imag,abstol);
//            EXPECT_NEAR(outdat1[1].imag,outdat2[1].imag,abstol);
//            EXPECT_NEAR(outdat1[2].imag,outdat2[2].imag,abstol);
//            EXPECT_NEAR(outdat1[3].imag,outdat2[3].imag,abstol);
//        }
    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2files,TestdatasetP)
{
    if(bool_packet)
    {
        DataSet ds1 = cppfile.openDataSet(datasetPpath);
        DataSet ds2 = pyfile.openDataSet(datasetPpath);

        CompType comp1=ds1.getCompType(); //should be equal to mytype
        CompType comp2=ds2.getCompType(); //should be equal to mytype
        ASSERT_EQ(comp1,comp2);
        ASSERT_EQ(comp1,mytype);

        DataSpace dspace1 = ds1.getSpace();
        DataSpace dspace2 = ds2.getSpace();

        int rank1 = dspace1.getSimpleExtentNdims();
        int rank2 = dspace1.getSimpleExtentNdims();
        ASSERT_EQ(rank1,rank2);

        hsize_t* dim1 = new hsize_t[rank1];
        hsize_t* dim2 = new hsize_t[rank2];

        dspace1.getSimpleExtentDims(dim1);
        dspace2.getSimpleExtentDims(dim2);

        ASSERT_EQ(dim1[0],dim2[0]);
        ASSERT_EQ(dim1[1],dim2[1]);
        ASSERT_EQ(dim1[2],dim2[2]);
        int RANK3=3;
        ASSERT_EQ(rank1,RANK3);

        hsize_t* elem=new hsize_t[rank1];
        elem[0]=1;
        elem[1]=dim1[1];
        elem[2]=dim1[2];
        ctype* outdat1 = new ctype[dim1[1]*dim1[2]];
        ctype* outdat2 = new ctype[dim2[1]*dim2[2]];

        DataSpace elemspace(rank1,elem);
        hsize_t* start1=new hsize_t[rank1];
        start1[0]=0;
        start1[1]=0;
        start1[2]=0;
        hsize_t* count1= new hsize_t[rank1];
        count1[0]=1;
        count1[1]=dim1[1];
        count1[2]=dim1[2];
        hsize_t* stride1=new hsize_t[rank1];
        stride1[0]=1;
        stride1[1]=1;
        stride1[2]=1;
        hsize_t* block1=new hsize_t[rank1];
        block1[0]=1;
        block1[1]=1;
        block1[2]=1;

        elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

        for(unsigned int k=0;k<dim1[0];++k)
        {
            hsize_t* start2 = new hsize_t[rank1];
            start2[0]=k;
            start2[1]=0;
            start2[2]=0;
            hsize_t* count2 = new hsize_t[rank1];
            count2[0]=1;
            count2[1]=dim1[1];
            count2[2]=dim1[2];
            hsize_t* stride2 = new hsize_t[rank1];
            stride2[0]=1;
            stride2[1]=1;
            stride2[2]=1;
            hsize_t* block2 = new hsize_t[rank1];
            block2[0]=1;
            block2[1]=1;
            block2[2]=1;

            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            ds1.read(outdat1,comp1,elemspace,dspace1);
            ds2.read(outdat2,comp2,elemspace,dspace2);

            EXPECT_NEAR(outdat1[0].real,outdat2[0].real,abstol);
            EXPECT_NEAR(outdat1[1].real,outdat2[1].real,abstol);
            EXPECT_NEAR(outdat1[2].real,outdat2[2].real,abstol);
            EXPECT_NEAR(outdat1[3].real,outdat2[3].real,abstol);

            EXPECT_NEAR(outdat1[0].imag,outdat2[0].imag,abstol);
            EXPECT_NEAR(outdat1[1].imag,outdat2[1].imag,abstol);
            EXPECT_NEAR(outdat1[2].imag,outdat2[2].imag,abstol);
            EXPECT_NEAR(outdat1[3].imag,outdat2[3].imag,abstol);
        }
    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2files,Testdatasetq)
{
    if(bool_packet)
    {
        DataSet ds1 = cppfile.openDataSet(datasetqpath);
        DataSet ds2 = pyfile.openDataSet(datasetqpath);

        CompType comp1=ds1.getCompType(); //should be equal to mytype
        CompType comp2=ds2.getCompType(); //should be equal to mytype
        ASSERT_EQ(comp1,comp2);
        ASSERT_EQ(comp1,mytype);

        DataSpace dspace1 = ds1.getSpace();
        DataSpace dspace2 = ds2.getSpace();

        int rank1 = dspace1.getSimpleExtentNdims();
        int rank2 = dspace1.getSimpleExtentNdims();
        ASSERT_EQ(rank1,rank2);

        hsize_t* dim1 = new hsize_t[rank1];
        hsize_t* dim2 = new hsize_t[rank2];

        dspace1.getSimpleExtentDims(dim1);
        dspace2.getSimpleExtentDims(dim2);

        ASSERT_EQ(dim1[0],dim2[0]);
        ASSERT_EQ(dim1[1],dim2[1]);
        ASSERT_EQ(dim1[2],dim2[2]);
        int RANK3=3;
        ASSERT_EQ(rank1,RANK3);

        hsize_t* elem=new hsize_t[rank1];
        elem[0]=1;
        elem[1]=dim1[1];
        elem[2]=dim1[2];
        ctype* outdat1 = new ctype[dim1[1]*dim1[2]];
        ctype* outdat2 = new ctype[dim2[1]*dim2[2]];

        DataSpace elemspace(rank1,elem);
        hsize_t* start1=new hsize_t[rank1];
        start1[0]=0;
        start1[1]=0;
        start1[2]=0;
        hsize_t* count1= new hsize_t[rank1];
        count1[0]=1;
        count1[1]=dim1[1];
        count1[2]=dim1[2];
        hsize_t* stride1=new hsize_t[rank1];
        stride1[0]=1;
        stride1[1]=1;
        stride1[2]=1;
        hsize_t* block1=new hsize_t[rank1];
        block1[0]=1;
        block1[1]=1;
        block1[2]=1;

        elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

        for(unsigned int k=0;k<dim1[0];++k)
        {
            hsize_t* start2 = new hsize_t[rank1];
            start2[0]=k;
            start2[1]=0;
            start2[2]=0;
            hsize_t* count2 = new hsize_t[rank1];
            count2[0]=1;
            count2[1]=dim1[1];
            count2[2]=dim1[2];
            hsize_t* stride2 = new hsize_t[rank1];
            stride2[0]=1;
            stride2[1]=1;
            stride2[2]=1;
            hsize_t* block2 = new hsize_t[rank1];
            block2[0]=1;
            block2[1]=1;
            block2[2]=1;

            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            ds1.read(outdat1,comp1,elemspace,dspace1);
            ds2.read(outdat2,comp2,elemspace,dspace2);

            EXPECT_NEAR(outdat1[0].real,outdat2[0].real,abstol);
            EXPECT_NEAR(outdat1[1].real,outdat2[1].real,abstol);

            EXPECT_NEAR(outdat1[0].imag,outdat2[0].imag,abstol);
            EXPECT_NEAR(outdat1[1].imag,outdat2[1].imag,abstol);
        }
    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2files,Testdatasetp)
{
    if(bool_packet)
    {
        DataSet ds1 = cppfile.openDataSet(datasetppath);
        DataSet ds2 = pyfile.openDataSet(datasetppath);

        CompType comp1=ds1.getCompType(); //should be equal to mytype
        CompType comp2=ds2.getCompType(); //should be equal to mytype
        ASSERT_EQ(comp1,comp2);
        ASSERT_EQ(comp1,mytype);

        DataSpace dspace1 = ds1.getSpace();
        DataSpace dspace2 = ds2.getSpace();

        int rank1 = dspace1.getSimpleExtentNdims();
        int rank2 = dspace1.getSimpleExtentNdims();
        ASSERT_EQ(rank1,rank2);

        hsize_t* dim1 = new hsize_t[rank1];
        hsize_t* dim2 = new hsize_t[rank2];

        dspace1.getSimpleExtentDims(dim1);
        dspace2.getSimpleExtentDims(dim2);

        ASSERT_EQ(dim1[0],dim2[0]);
        ASSERT_EQ(dim1[1],dim2[1]);
        ASSERT_EQ(dim1[2],dim2[2]);
        int RANK3=3;
        ASSERT_EQ(rank1,RANK3);

        hsize_t* elem=new hsize_t[rank1];
        elem[0]=1;
        elem[1]=dim1[1];
        elem[2]=dim1[2];
        ctype* outdat1 = new ctype[dim1[1]*dim1[2]];
        ctype* outdat2 = new ctype[dim2[1]*dim2[2]];

        DataSpace elemspace(rank1,elem);
        hsize_t* start1=new hsize_t[rank1];
        start1[0]=0;
        start1[1]=0;
        start1[2]=0;
        hsize_t* count1= new hsize_t[rank1];
        count1[0]=1;
        count1[1]=dim1[1];
        count1[2]=dim1[2];
        hsize_t* stride1=new hsize_t[rank1];
        stride1[0]=1;
        stride1[1]=1;
        stride1[2]=1;
        hsize_t* block1=new hsize_t[rank1];
        block1[0]=1;
        block1[1]=1;
        block1[2]=1;

        elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

        for(unsigned int k=0;k<dim1[0];++k)
        {
            hsize_t* start2 = new hsize_t[rank1];
            start2[0]=k;
            start2[1]=0;
            start2[2]=0;
            hsize_t* count2 = new hsize_t[rank1];
            count2[0]=1;
            count2[1]=dim1[1];
            count2[2]=dim1[2];
            hsize_t* stride2 = new hsize_t[rank1];
            stride2[0]=1;
            stride2[1]=1;
            stride2[2]=1;
            hsize_t* block2 = new hsize_t[rank1];
            block2[0]=1;
            block2[1]=1;
            block2[2]=1;

            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            ds1.read(outdat1,comp1,elemspace,dspace1);
            ds2.read(outdat2,comp2,elemspace,dspace2);

            EXPECT_NEAR(outdat1[0].real,outdat2[0].real,abstol);
            EXPECT_NEAR(outdat1[1].real,outdat2[1].real,abstol);

            EXPECT_NEAR(outdat1[0].imag,outdat2[0].imag,abstol);
            EXPECT_NEAR(outdat1[1].imag,outdat2[1].imag,abstol);
        }
    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2files,Testdatasetcoefficients)
{
    if(bool_packet)
    {
        //expect rank of coeffs always to be 2
        DataSet ds1 = cppfile.openDataSet(datasetcpath);
        DataSet ds2 = pyfile.openDataSet(datasetcpath);

        CompType comp1=ds1.getCompType(); //should be equal to mytype
        CompType comp2=ds2.getCompType(); //should be equal to mytype
        ASSERT_EQ(comp1,comp2);
        ASSERT_EQ(comp1,mytype);

        DataSpace dspace1 = ds1.getSpace();
        DataSpace dspace2 = ds2.getSpace();

        int rank1 = dspace1.getSimpleExtentNdims();
        int rank2 = dspace1.getSimpleExtentNdims();
        ASSERT_EQ(rank1,rank2);
        int RANK2=2;
        ASSERT_EQ(rank1,RANK2);

        hsize_t* dim1 = new hsize_t[rank1];
        hsize_t* dim2 = new hsize_t[rank2];

        dspace1.getSimpleExtentDims(dim1);
        dspace2.getSimpleExtentDims(dim2);

        ASSERT_EQ(dim1[0],dim2[0]);
        ASSERT_EQ(dim1[1],dim2[1]);

        hsize_t* elem= new hsize_t[rank1];
        elem[0]=1;
        elem[1]=dim1[1];
        DataSpace elemspace(rank1,elem);

        ctype* outdat1 = new ctype[dim1[1]];
        ctype* outdat2 = new ctype[dim2[1]];

        hsize_t* start1 = new hsize_t[rank1];
        start1[0]=0;
        start1[1]=0;
        hsize_t* count1 = new hsize_t[rank1];
        count1[0]=1;
        count1[1]=dim1[1];
        hsize_t* stride1=new hsize_t[rank1];
        stride1[0]=1;
        stride1[1]=1;
        hsize_t* block1=new hsize_t[rank1];
        block1[0]=1;
        block1[1]=1;

        elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

        for(unsigned int k=0;k<dim1[0];++k)
        {
            hsize_t* start2 = new hsize_t[rank1];
            start2[0]=k;
            start2[1]=0;
            hsize_t* count2 = new hsize_t[rank1];
            count2[0]=1;
            count2[1]=dim1[1];
            hsize_t* stride2 = new hsize_t[rank1];
            stride2[0]=1;
            stride2[1]=1;
            hsize_t* block2 = new hsize_t[rank1];
            block2[0]=1;
            block2[1]=1;

            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            ds1.read(outdat1,comp1,elemspace,dspace1);
            ds2.read(outdat2,comp2,elemspace,dspace2);

            for(unsigned int i=0;i<dim1[1];++i)
            {
            EXPECT_NEAR(outdat1[i].real,outdat2[i].real,abstol);
            EXPECT_NEAR(outdat1[i].imag,outdat2[i].imag,abstol);
            }

        }
    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2files,DISABLED_Testdatasetekin)
{
    if(bool_energy)
    {
        //expect rank of ekin always be 2
        DataSet ds1 = cppfile.openDataSet(datasetekinpath);
        DataSet ds2 = pyfile.openDataSet(datasetekinpath);

        DataType tp1 = ds1.getDataType();
        DataType tp2 = ds1.getDataType();

        DataSpace dspace1 = ds1.getSpace();
        DataSpace dspace2 = ds2.getSpace();

        int rank1 = dspace1.getSimpleExtentNdims();
        int rank2 = dspace1.getSimpleExtentNdims();
        ASSERT_EQ(rank1,rank2);

        hsize_t* dim1 = new hsize_t[rank1];
        hsize_t* dim2 = new hsize_t[rank2];

        dspace1.getSimpleExtentDims(dim1);
        dspace2.getSimpleExtentDims(dim2);

        ASSERT_EQ(dim1[0],dim2[0]);
        ASSERT_EQ(dim1[1],dim2[1]);
        int RANK2=2;
        ASSERT_EQ(rank1,RANK2);

        hsize_t* elem=new hsize_t[rank1];
        elem[0]=1;
        elem[1]=dim1[1];
        double* outdat1 = new double[dim1[1]];
        double* outdat2 = new double[dim2[1]];

        DataSpace elemspace(rank1,elem);
        hsize_t* start1=new hsize_t[rank1];
        start1[0]=0;
        start1[1]=0;
        hsize_t* count1= new hsize_t[rank1];
        count1[0]=1;
        count1[1]=dim1[1];
        hsize_t* stride1=new hsize_t[rank1];
        stride1[0]=1;
        stride1[1]=1;
        hsize_t* block1=new hsize_t[rank1];
        block1[0]=1;
        block1[1]=1;

        elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

        for(unsigned int k=0;k<dim1[0];++k)
        {
            hsize_t* start2 = new hsize_t[rank1];
            start2[0]=k;
            start2[1]=0;
            hsize_t* count2 = new hsize_t[rank1];
            count2[0]=1;
            count2[1]=dim1[1];
            hsize_t* stride2 = new hsize_t[rank1];
            stride2[0]=1;
            stride2[1]=1;
            hsize_t* block2 = new hsize_t[rank1];
            block2[0]=1;
            block2[1]=1;

            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            ds1.read(outdat1,tp1,elemspace,dspace1);
            ds2.read(outdat2,tp2,elemspace,dspace2);

            EXPECT_NEAR(outdat1[0],outdat2[0],abstol);
        }

    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2files,DISABLED_Testdatasetepot)
{
    if(bool_energy)
    {
        //expect rank of epot always be 2
        DataSet ds1 = cppfile.openDataSet(datasetepotpath);
        DataSet ds2 = pyfile.openDataSet(datasetepotpath);

        DataType tp1 = ds1.getDataType();
        DataType tp2 = ds1.getDataType();

        DataSpace dspace1 = ds1.getSpace();
        DataSpace dspace2 = ds2.getSpace();

        int rank1 = dspace1.getSimpleExtentNdims();
        int rank2 = dspace1.getSimpleExtentNdims();
        ASSERT_EQ(rank1,rank2);

        hsize_t* dim1 = new hsize_t[rank1];
        hsize_t* dim2 = new hsize_t[rank2];

        dspace1.getSimpleExtentDims(dim1);
        dspace2.getSimpleExtentDims(dim2);

        ASSERT_EQ(dim1[0],dim2[0]);
        ASSERT_EQ(dim1[1],dim2[1]);
        int RANK2=2;
        ASSERT_EQ(rank1,RANK2);

        hsize_t* elem=new hsize_t[rank1];
        elem[0]=1;
        elem[1]=dim1[1];
        double* outdat1 = new double[dim1[1]];
        double* outdat2 = new double[dim2[1]];

        DataSpace elemspace(rank1,elem);
        hsize_t* start1=new hsize_t[rank1];
        start1[0]=0;
        start1[1]=0;
        hsize_t* count1= new hsize_t[rank1];
        count1[0]=1;
        count1[1]=dim1[1];
        hsize_t* stride1=new hsize_t[rank1];
        stride1[0]=1;
        stride1[1]=1;
        hsize_t* block1=new hsize_t[rank1];
        block1[0]=1;
        block1[1]=1;

        elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

        for(unsigned int k=0;k<dim1[0];++k)
        {
            hsize_t* start2 = new hsize_t[rank1];
            start2[0]=k;
            start2[1]=0;
            hsize_t* count2 = new hsize_t[rank1];
            count2[0]=1;
            count2[1]=dim1[1];
            hsize_t* stride2 = new hsize_t[rank1];
            stride2[0]=1;
            stride2[1]=1;
            hsize_t* block2 = new hsize_t[rank1];
            block2[0]=1;
            block2[1]=1;

            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            ds1.read(outdat1,tp1,elemspace,dspace1);
            ds2.read(outdat2,tp2,elemspace,dspace2);

            EXPECT_NEAR(outdat1[0],outdat2[0],abstol);
        }

    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2files,DISABLED_Testdatasetnorm)
{
    if(bool_norm)
    {
        //expect rank of norm always be 2
        DataSet ds1 = cppfile.openDataSet(datasetnormpath);
        DataSet ds2 = pyfile.openDataSet(datasetnormpath);

        DataType tp1 = ds1.getDataType();
        DataType tp2 = ds1.getDataType();

        DataSpace dspace1 = ds1.getSpace();
        DataSpace dspace2 = ds2.getSpace();

        int rank1 = dspace1.getSimpleExtentNdims();
        int rank2 = dspace1.getSimpleExtentNdims();
        ASSERT_EQ(rank1,rank2);

        hsize_t* dim1 = new hsize_t[rank1];
        hsize_t* dim2 = new hsize_t[rank2];

        dspace1.getSimpleExtentDims(dim1);
        dspace2.getSimpleExtentDims(dim2);

        ASSERT_EQ(dim1[0],dim2[0]);
        ASSERT_EQ(dim1[1],dim2[1]);
        int RANK2=2;
        ASSERT_EQ(rank1,RANK2);

        hsize_t* elem=new hsize_t[rank1];
        elem[0]=1;
        elem[1]=dim1[1];
        double* outdat1 = new double[dim1[1]];
        double* outdat2 = new double[dim2[1]];

        DataSpace elemspace(rank1,elem);
        hsize_t* start1=new hsize_t[rank1];
        start1[0]=0;
        start1[1]=0;
        hsize_t* count1= new hsize_t[rank1];
        count1[0]=1;
        count1[1]=dim1[1];
        hsize_t* stride1=new hsize_t[rank1];
        stride1[0]=1;
        stride1[1]=1;
        hsize_t* block1=new hsize_t[rank1];
        block1[0]=1;
        block1[1]=1;

        elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

        for(unsigned int k=0;k<dim1[0];++k)
        {
            hsize_t* start2 = new hsize_t[rank1];
            start2[0]=k;
            start2[1]=0;
            hsize_t* count2 = new hsize_t[rank1];
            count2[0]=1;
            count2[1]=dim1[1];
            hsize_t* stride2 = new hsize_t[rank1];
            stride2[0]=1;
            stride2[1]=1;
            hsize_t* block2 = new hsize_t[rank1];
            block2[0]=1;
            block2[1]=1;

            dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            dspace2.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
            ds1.read(outdat1,tp1,elemspace,dspace1);
            ds2.read(outdat2,tp2,elemspace,dspace2);

            EXPECT_NEAR(outdat1[0],outdat2[0],abstol);
        }

    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}
/**
 * @brief main
 * @param argc
 * @param argv
 * @return run all google defined tests
 *
 * Uses global variables to use argv in google SetUp function.
 */
int main(int argc,char* argv[])
{
    global_argc=argc;
    global_argv=argv;
    ::testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}
//void compute_time_matching(H5std_string timegridpath_cpp,H5std_string timegridpath_py)
//{
//    std::vector<int*> res;

//    //expect RANK=1
//    DataSet its1 = cppfile.openDataSet(timegridpath_cpp);
//    DataSet its2 = pyfile.openDataSet(timegridpath_py);

//    IntType it1 = its1.getIntType();
//    IntType it2 = its2.getIntType();

//    DataSpace dspace1 = its1.getSpace();
//    DataSpace dspace2 = its2.getSpace();

//    int rank1 = dspace1.getSimpleExtentNdims();
//    int rank2 = dspace1.getSimpleExtentNdims();
//    assert(rank1==1);
//    assert(rank2==1);
//    hsize_t* dim1 = new hsize_t[rank1];
//    hsize_t* dim2 = new hsize_t[rank2];

//    dspace1.getSimpleExtentDims(dim1);
//    dspace2.getSimpleExtentDims(dim2);
//    assert(dim1[0]==dim2[0]);

//    int* index_array_cpp = new int[dim1[0]];
//    int* index_array_py = new int[dim2[0]];

//    hsize_t* elem=new hsize_t[rank1];
//    elem[0]=dim1[0];

//    DataSpace elemspace(rank1,elem);
//    hsize_t* start1=new hsize_t[rank1];
//    start1[0]=0;
//    hsize_t* count1= new hsize_t[rank1];
//    count1[0]=dim1[0];
//    hsize_t* stride1=new hsize_t[rank1];
//    stride1[0]=1;
//    hsize_t* block1=new hsize_t[rank1];
//    block1[0]=1;

//    elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);
//    dspace1.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);
//    dspace2.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

//    its1.read(index_array_cpp,it1,elemspace,dspace1);
//    its2.read(index_array_py,it2,elemspace,dspace2);

//    int acc[2]={-1,-1};
//    for(unsigned int k=0;k<dim1[0];++k)
//    {
//        int index_i_cpp=index_array_cpp[k];
//        for(unsigned int p=0;p<dim2[0];++p)
//        {
//            int index_j_py=index_array_py[p];
//            int temp = index_j_py*(dt_py/dt_cpp);
//            if(index_i_cpp==temp)
//            {
//                acc[0]=index_i_cpp;
//                acc[1]=index_j_py;
//                res.push_back(acc);
//            }
//        }
//    }
//    if(res.empty())
//    {
//        res.push_back(acc);
//    }

//}
