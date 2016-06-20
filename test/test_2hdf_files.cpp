
#include <iostream>
#include <string>
#include <cmath>

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
    char* this_file=global_argv[0];
    char* filename_1 = global_argv[1];
    char* filename_2 = global_argv[2];
    filename1 = filename_1;
    filename2 = filename_2;
    mytype=CompType(sizeof(instanceof));
    mytype.insertMember( "r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
    mytype.insertMember( "i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);


    file1=H5File(filename1,H5F_ACC_RDONLY);
    file2=H5File(filename2,H5F_ACC_RDONLY);

//    std::cout<<"Objects:"<<file1.getNumObjs()<<'\n';
//    std::cout<<"Objects:"<<file2.getNumObjs()<<'\n';

//    std::cout<<"Attributes:"<<file1.getNumAttrs()<<'\n';
//    std::cout<<"Attributes:"<<file2.getNumAttrs()<<'\n';

    datasetQpath="/datablock_0/wavepacket/Pi/Q";
    datasetPpath="/datablock_0/wavepacket/Pi/P";
    datasetqpath="/datablock_0/wavepacket/Pi/q";
    datasetppath="/datablock_0/wavepacket/Pi/p";
    datasetcpath="/datablock_0/wavepacket/coefficients/c_0";
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
    H5std_string datasetQpath;
    H5std_string datasetPpath;
    H5std_string datasetqpath;
    H5std_string datasetppath;
    H5std_string datasetcpath;
};


TEST_F(Test2files,TestdatasetQ)
{

    DataSet ds1 = file1.openDataSet(datasetQpath);
    DataSet ds2 = file2.openDataSet(datasetQpath);

    H5T_class_t c1 = ds1.getTypeClass();//gives CompType
    H5T_class_t c2 = ds2.getTypeClass();//gives CompType

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    EXPECT_EQ(c1,c2);

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
    hsize_t elem[3]={1,dim1[1],dim1[2]};
    ctype* outdat1 = new ctype[dim1[1]*dim1[2]];
    ctype* outdat2 = new ctype[dim2[1]*dim2[2]];

    DataSpace elemspace(RANK3,elem);
    hsize_t start1[3]={0,0,0};
    hsize_t count1[3]={1,2,2};
    hsize_t stride1[3]={1,1,1};
    hsize_t block1[3]={1,1,1};

    elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

    for(unsigned int k=0;k<dim1[0];++k)
    {
        hsize_t start2[3];
        hsize_t count2[3]={1,2,2};
        hsize_t stride2[3]={1,1,1};
        hsize_t block2[3]={1,1,1};
        start2[0]=k;
        start2[1]=0;
        start2[2]=0;

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

    H5T_class_t c1 = ds1.getTypeClass();//gives CompType
    H5T_class_t c2 = ds2.getTypeClass();//gives CompType

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    EXPECT_EQ(c1,c2);

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
    hsize_t elem[3]={1,dim1[1],dim1[2]};
    ctype* outdat1 = new ctype[dim1[1]*dim1[2]];
    ctype* outdat2 = new ctype[dim2[1]*dim2[2]];

    DataSpace elemspace(RANK3,elem);
    hsize_t start1[3]={0,0,0};
    hsize_t count1[3]={1,2,2};
    hsize_t stride1[3]={1,1,1};
    hsize_t block1[3]={1,1,1};

    elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

    for(unsigned int k=0;k<dim1[0];++k)
    {
        hsize_t start2[3];
        hsize_t count2[3]={1,2,2};
        hsize_t stride2[3]={1,1,1};
        hsize_t block2[3]={1,1,1};
        start2[0]=k;
        start2[1]=0;
        start2[2]=0;

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

    H5T_class_t c1 = ds1.getTypeClass();//gives CompType
    H5T_class_t c2 = ds2.getTypeClass();//gives CompType

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    EXPECT_EQ(c1,c2);

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
    hsize_t elem[3]={1,dim1[1],dim1[2]};
    ctype* outdat1 = new ctype[dim1[1]*dim1[2]];
    ctype* outdat2 = new ctype[dim2[1]*dim2[2]];

    DataSpace elemspace(RANK3,elem);
    hsize_t start1[3]={0,0,0};
    hsize_t count1[3]={1,2,1};
    hsize_t stride1[3]={1,1,1};
    hsize_t block1[3]={1,1,1};

    elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

    for(unsigned int k=0;k<dim1[0];++k)
    {
        hsize_t start2[3];
        hsize_t count2[3]={1,2,1};
        hsize_t stride2[3]={1,1,1};
        hsize_t block2[3]={1,1,1};
        start2[0]=k;
        start2[1]=0;
        start2[2]=0;

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

    H5T_class_t c1 = ds1.getTypeClass();//gives CompType
    H5T_class_t c2 = ds2.getTypeClass();//gives CompType

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    EXPECT_EQ(c1,c2);

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
    hsize_t elem[3]={1,dim1[1],dim1[2]};
    ctype* outdat1 = new ctype[dim1[1]*dim1[2]];
    ctype* outdat2 = new ctype[dim2[1]*dim2[2]];

    DataSpace elemspace(RANK3,elem);
    hsize_t start1[3]={0,0,0};
    hsize_t count1[3]={1,2,1};
    hsize_t stride1[3]={1,1,1};
    hsize_t block1[3]={1,1,1};

    elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

    for(unsigned int k=0;k<dim1[0];++k)
    {
        hsize_t start2[3];
        hsize_t count2[3]={1,2,1};
        hsize_t stride2[3]={1,1,1};
        hsize_t block2[3]={1,1,1};
        start2[0]=k;
        start2[1]=0;
        start2[2]=0;

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
    DataSet ds1 = file1.openDataSet(datasetcpath);
    DataSet ds2 = file2.openDataSet(datasetcpath);

    H5T_class_t c1 = ds1.getTypeClass();//gives CompType
    H5T_class_t c2 = ds2.getTypeClass();//gives CompType

    CompType comp1=ds1.getCompType(); //should be equal to mytype
    CompType comp2=ds2.getCompType(); //should be equal to mytype

    EXPECT_EQ(c1,c2);

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
    int RANK2=2;
    hsize_t elem[2]={1,dim1[1]};
    ctype* outdat1 = new ctype[dim1[1]];
    ctype* outdat2 = new ctype[dim2[1]];

    DataSpace elemspace(RANK2,elem);
    hsize_t start1[2]={0,0};
    hsize_t count1[2]={1,16};
    hsize_t stride1[2]={1,1};
    hsize_t block1[2]={1,1};

    elemspace.selectHyperslab(H5S_SELECT_SET,count1,start1,stride1,block1);

    for(unsigned int k=0;k<dim1[0];++k)
    {
        hsize_t start2[2];
        hsize_t count2[2]={1,16};
        hsize_t stride2[2]={1,1};
        hsize_t block2[2]={1,1};
        start2[0]=k;
        start2[1]=0;

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
