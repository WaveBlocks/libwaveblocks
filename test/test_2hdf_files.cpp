#include <iostream>
#include <string>
#include <cmath>
#include <memory>

//include GTest
#include "gtest/gtest.h"

#include "H5Cpp.h"

#define abstol 1e-6

//global variables for argv
int global_argc;
char** global_argv;
//namespace for h5
using namespace H5;

class Test2files: public ::testing::Test
{
    public:
    void SetUp()
    {
    //read filenames
    assert(global_argc==3);
    //char* this_file=global_argv[0];
    char* filename_1 = global_argv[1];
    char* filename_2 = global_argv[2];
    filename1 = filename_1;
    filename2 = filename_2;
    mytype=CompType(sizeof(instanceof));
    mytype.insertMember( "r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
    mytype.insertMember( "i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);

    file1=H5File(filename1,H5F_ACC_RDONLY);
    file2=H5File(filename2,H5F_ACC_RDONLY);

    datablock_file1=std::make_shared<Group>(file1.openGroup(datablock_group));
    datablock_file2=std::make_shared<Group>(file2.openGroup(datablock_group));
    apacket1=datablock_file1->openAttribute("packet");
    aenergy1=datablock_file1->openAttribute("energy");
    anorm1=datablock_file1->openAttribute("norm");



    datasetQpath=datablock_group+wavepacket_group+Pi_group+"/Q";
    datasetPpath=datablock_group+wavepacket_group+Pi_group+"/P";
    datasetqpath=datablock_group+wavepacket_group+Pi_group+"/q";
    datasetppath=datablock_group+wavepacket_group+Pi_group+"/p";
    datasetcpath=datablock_group+wavepacket_group+coeffs_group+"/c_0";
    datasetekinpath=datablock_group+energies_group+"/ekin";
    datasetepotpath=datablock_group+energies_group+"/epot";
    datasetnormspath=datablock_group+norms_group;//TODO
    }
    void BreakDown()
    {
    
    }
    struct ctype{ //our complex datatype
        double real=0.;
        double imag=0.;
    } instanceof;

    H5std_string filename1;
    H5File file1;
    H5std_string filename2;
    H5File file2;
    CompType mytype;
    H5std_string datablock_group="/datablock_0";
    H5std_string energies_group="/energies";
    H5std_string coeffs_group="/coefficients";
    H5std_string norms_group="/norms";
    H5std_string Pi_group="/Pi";
    H5std_string wavepacket_group="/wavepacket";
    std::shared_ptr<Group> datablock_file1;
    std::shared_ptr<Group> datablock_file2;

    H5std_string datasetQpath;
    H5std_string datasetPpath;
    H5std_string datasetqpath;
    H5std_string datasetppath;
    H5std_string datasetcpath;
    H5std_string datasetekinpath;
    H5std_string datasetepotpath;
    H5std_string datasetnormspath;
    Attribute apacket1;
    Attribute aenergy1;
    Attribute anorm1;
    Attribute apacket2;
    Attribute aenergy2;
    Attribute anorm2;
};


TEST_F(Test2files,TestdatasetQ)
{
    //expect RANK=3
    DataSet ds1 = file1.openDataSet(datasetQpath);
    DataSet ds2 = file2.openDataSet(datasetQpath);

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    DataSpace dspace1 = ds1.getSpace();
    DataSpace dspace2 = ds2.getSpace();

    int rank1 = dspace1.getSimpleExtentNdims();
    int rank2 = dspace1.getSimpleExtentNdims();
    EXPECT_EQ(rank1,rank2);

    hsize_t* dim1 = new hsize_t[rank1];
    hsize_t* dim2 = new hsize_t[rank2];

    dspace1.getSimpleExtentDims(dim1);
    dspace2.getSimpleExtentDims(dim2);

    ASSERT_EQ(dim1[0],dim2[0]);
    ASSERT_EQ(dim1[1],dim2[1]);
    ASSERT_EQ(dim1[2],dim2[2]);
    int RANK3=3;
    EXPECT_EQ(rank1,RANK3);

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

TEST_F(Test2files,TestdatasetP)
{

    DataSet ds1 = file1.openDataSet(datasetPpath);
    DataSet ds2 = file2.openDataSet(datasetPpath);

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    DataSpace dspace1 = ds1.getSpace();
    DataSpace dspace2 = ds2.getSpace();

    int rank1 = dspace1.getSimpleExtentNdims();
    int rank2 = dspace1.getSimpleExtentNdims();
    EXPECT_EQ(rank1,rank2);


    hsize_t* dim1 = new hsize_t[rank1];
    hsize_t* dim2 = new hsize_t[rank2];

    dspace1.getSimpleExtentDims(dim1);
    dspace2.getSimpleExtentDims(dim2);

    ASSERT_EQ(dim1[0],dim2[0]);
    ASSERT_EQ(dim1[1],dim2[1]);
    ASSERT_EQ(dim1[2],dim2[2]);
    int RANK3=3;
    EXPECT_EQ(rank1,RANK3);

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

TEST_F(Test2files,Testdatasetq)
{

    DataSet ds1 = file1.openDataSet(datasetqpath);
    DataSet ds2 = file2.openDataSet(datasetqpath);

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    DataSpace dspace1 = ds1.getSpace();
    DataSpace dspace2 = ds2.getSpace();

    int rank1 = dspace1.getSimpleExtentNdims();
    int rank2 = dspace1.getSimpleExtentNdims();
    EXPECT_EQ(rank1,rank2);

    hsize_t* dim1 = new hsize_t[rank1];
    hsize_t* dim2 = new hsize_t[rank2];

    dspace1.getSimpleExtentDims(dim1);
    dspace2.getSimpleExtentDims(dim2);

    ASSERT_EQ(dim1[0],dim2[0]);
    ASSERT_EQ(dim1[1],dim2[1]);
    ASSERT_EQ(dim1[2],dim2[2]);
    int RANK3=3;
    EXPECT_EQ(rank1,RANK3);

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

TEST_F(Test2files,Testdatasetp)
{

    DataSet ds1 = file1.openDataSet(datasetppath);
    DataSet ds2 = file2.openDataSet(datasetppath);

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    DataSpace dspace1 = ds1.getSpace();
    DataSpace dspace2 = ds2.getSpace();

    int rank1 = dspace1.getSimpleExtentNdims();
    int rank2 = dspace1.getSimpleExtentNdims();

    hsize_t* dim1 = new hsize_t[rank1];
    hsize_t* dim2 = new hsize_t[rank2];

    dspace1.getSimpleExtentDims(dim1);
    dspace2.getSimpleExtentDims(dim2);

    ASSERT_EQ(dim1[0],dim2[0]);
    ASSERT_EQ(dim1[1],dim2[1]);
    ASSERT_EQ(dim1[2],dim2[2]);
    int RANK3=3;
    EXPECT_EQ(rank1,RANK3);

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

TEST_F(Test2files,Testdatasetcoefficients)
{
    //expect rank of coeffs always to be 2
    DataSet ds1 = file1.openDataSet(datasetcpath);
    DataSet ds2 = file2.openDataSet(datasetcpath);

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    DataSpace dspace1 = ds1.getSpace();
    DataSpace dspace2 = ds2.getSpace();

    int rank1 = dspace1.getSimpleExtentNdims();
    int rank2 = dspace1.getSimpleExtentNdims();
    EXPECT_EQ(rank1,rank2);
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

int main(int argc,char* argv[])
{
    global_argc=argc;
    global_argv=argv;
    ::testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}
