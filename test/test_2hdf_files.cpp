#include <iostream>
#include <string>
#include <vector>

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
 * This class need to be derived from ::testing::Test in gtest.h. This class is used afterwards in a testfissure TEST_F
 * to make test with the given initalized variables/parameters.
 */
class Test2HDFfiles: public ::testing::Test
{
    public:
    /**
     * @brief Test2HDFfiles
     *
     * Constructor used to open the two give files over argv.
     */
    Test2HDFfiles()
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
    /**
     * @brief ~Test2HDFfiles
     *
     * Destructor to close used files.
     */
    virtual ~Test2HDFfiles()
    {
        cppfile.close();
        pyfile.close();
    }
    /**
     * @brief SetUp for TEST_F
     *
     * Initialize strings, ints and doubles.
     */
    void SetUp()
    {
        Group block_0 = cppfile.openGroup("/datablock_0");
        Attribute apacket = block_0.openAttribute("packet");
        Attribute aenergy = block_0.openAttribute("energy");
        Attribute anorm = block_0.openAttribute("norm");
        Attribute adt_cpp = block_0.openAttribute("dt");
        apacket.read(PredType::NATIVE_INT,&bool_packet);
        aenergy.read(PredType::NATIVE_INT,&bool_energy);
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

        mytype=CompType(sizeof(instanceof));
        mytype.insertMember( "r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
        mytype.insertMember( "i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);

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
     * Cleanup for used Type
     */
    void BreakDown()
    {
        mytype.close();
    }
    /**
     * @brief The ctype struct
     *
     * Our used struct for complex datatypes compatible with python interface.
     */
    struct ctype{
        double real=0.;
        double imag=0.;
    } instanceof;
    /**
     * @brief compute time matching for two timegrid paths
     * @param res
     * @param timegridpath_cpp
     * @param timegridpath_py
     *
     * The output is a matching vector with minimal length of one.
     */
    void compute_time_matching(std::vector<int*>& res,H5std_string timegridpath_cpp,H5std_string timegridpath_py)
    {
        //expect RANK=1
        DataSet its1 = cppfile.openDataSet(timegridpath_cpp);
        DataSet its2 = pyfile.openDataSet(timegridpath_py);

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

        hsize_t* it_dyn_start=new hsize_t[itrank1];

        for(unsigned int k=0;k<itdim1[0];++k)
        {
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
    }
    /**
     * @brief Test dataset Q
     * @param res time matching
     *
     * Test dataset Q over all elements in the time matching.
     */
    void Test_Q(std::vector<int*>& res)
    {
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
            ADD_FAILURE()<<"No matching timepoint for Q.Abort!";
        }
        else if(res.data()[0][0]==-1)
        {
            ADD_FAILURE()<<"No matching timepoint for P.Abort!";
        }
        else
        {
            std::cout<<res.size()<<" matching timepoints found in Q\n";
            hsize_t* start2 = new hsize_t[rank1];
            hsize_t* count2 = new hsize_t[rank1];
            hsize_t* stride2 = new hsize_t[rank1];
            hsize_t* block2 = new hsize_t[rank1];
            hsize_t* start3 = new hsize_t[rank1];
            for(auto index_element:res)
            {
                start2[0]=index_element[0];
                start2[1]=0;
                start2[2]=0;
                count2[0]=1;
                count2[1]=dim1[1];
                count2[2]=dim1[2];
                stride2[0]=1;
                stride2[1]=1;
                stride2[2]=1;
                block2[0]=1;
                block2[1]=1;
                block2[2]=1;
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
        }
    }
    /**
     * @brief Test dataset P
     * @param res time matching
     *
     * Test dataset P over all elements in the time matching.
     */
    void Test_P(std::vector<int*>& res)
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

        if(res.empty())
        {
            ADD_FAILURE()<<"No matching timepoint for P.Abort!";
        }
        else if(res.data()[0][0]==-1)
        {
            ADD_FAILURE()<<"No matching timepoint for P.Abort!";
        }
        else
        {
            std::cout<<res.size()<<" matching timepoints found in P\n";
            hsize_t* start2 = new hsize_t[rank1];
            hsize_t* count2 = new hsize_t[rank1];
            hsize_t* stride2 = new hsize_t[rank1];
            hsize_t* block2 = new hsize_t[rank1];
            hsize_t* start3 = new hsize_t[rank1];
            for(auto index_element:res)
            {
                start2[0]=index_element[0];
                start2[1]=0;
                start2[2]=0;
                count2[0]=1;
                count2[1]=dim1[1];
                count2[2]=dim1[2];
                stride2[0]=1;
                stride2[1]=1;
                stride2[2]=1;
                block2[0]=1;
                block2[1]=1;
                block2[2]=1;
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
        }
    }
    /**
     * @brief Test dataset q
     * @param res time matching
     *
     * Test dataset q over all elements in the time matching.
     */
    void Test_q(std::vector<int*>& res)
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

        if(res.empty())
        {
            ADD_FAILURE()<<"No matching timepoint for q.Abort!";
        }
        else if(res.data()[0][0]==-1)
        {
            ADD_FAILURE()<<"No matching timepoint for q.Abort!";
        }
        else
        {
            std::cout<<res.size()<<" matching timepoints found in q\n";
            for(auto index_element:res)
            {
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

                EXPECT_NEAR(outdat1[0].imag,outdat2[0].imag,abstol);
                EXPECT_NEAR(outdat1[1].imag,outdat2[1].imag,abstol);
            }
        }
    }
    /**
     * @brief Test dataset p
     * @param res time matching
     *
     * Test dataset p over all elements in the time matching.
     */
    void Test_p(std::vector<int*>& res)
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

        if(res.empty())
        {
            ADD_FAILURE()<<"No matching timepoint for p.Abort!";
        }
        else if(res.data()[0][0]==-1)
        {
            ADD_FAILURE()<<"No matching timepoint for p.Abort!";
        }
        else
        {
            std::cout<<res.size()<<" matching timepoints found in p\n";
            for(auto index_element:res)
            {
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

                EXPECT_NEAR(outdat1[0].imag,outdat2[0].imag,abstol);
                EXPECT_NEAR(outdat1[1].imag,outdat2[1].imag,abstol);
            }
        }
    }
    /**
     * @brief Test dataset coefficients
     * @param res time matching
     *
     * Test dataset coefficients over all elements in the time matching.
     */
    void Test_coefficients(std::vector<int*>& res)
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

        if(res.empty())
        {
            ADD_FAILURE()<<"No matching timepoint for coefficients.Abort!";
        }
        else if(res.data()[0][0]==-1)
        {
            ADD_FAILURE()<<"No matching timepoint for coefficients.Abort!";
        }
        else
        {
            std::cout<<res.size()<<" matching timepoints found in coefficients\n";
            for(auto index_element:res)
            {
                hsize_t* start2 = new hsize_t[rank1];
                start2[0]=index_element[0];
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
                hsize_t* start3 = new hsize_t[rank1];
                start3[0]=index_element[1];
                start3[1]=0;

                dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
                dspace2.selectHyperslab(H5S_SELECT_SET, count2, start3, stride2, block2);
                ds1.read(outdat1,comp1,elemspace,dspace1);
                ds2.read(outdat2,comp2,elemspace,dspace2);

                for(unsigned int i=0;i<dim1[1];++i)
                {
                EXPECT_NEAR(outdat1[i].real,outdat2[i].real,abstol);
                EXPECT_NEAR(outdat1[i].imag,outdat2[i].imag,abstol);
                }
            }
        }
    }
    /**
     * @brief Test dataset ekin
     * @param res time matching
     *
     * Test dataset ekin over all elements in the time matching.
     */
    void Test_ekin(std::vector<int*>& res)
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

        if(res.empty())
        {
            ADD_FAILURE()<<"No matching timepoint for ekin.Abort!";
        }
        else if(res.data()[0][0]==-1)
        {
            ADD_FAILURE()<<"No matching timepoint for ekin.Abort!";
        }
        else
        {
            std::cout<<res.size()<<" matching timepoints found in ekin\n";
            for(auto index_element:res)
            {
                hsize_t* start2 = new hsize_t[rank1];
                start2[0]=index_element[0];
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
                hsize_t* start3 = new hsize_t[rank1];
                start3[0]=index_element[1];
                start3[1]=0;
                dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
                dspace2.selectHyperslab(H5S_SELECT_SET, count2, start3, stride2, block2);
                ds1.read(outdat1,tp1,elemspace,dspace1);
                ds2.read(outdat2,tp2,elemspace,dspace2);

                EXPECT_NEAR(outdat1[0],outdat2[0],abstol);
            }
        }
    }
    /**
     * @brief Test dataset epot
     * @param res time matching
     *
     * Test dataset epot over all elements in the time matching.
     */
    void Test_epot(std::vector<int*>& res)
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

        if(res.empty())
        {
            ADD_FAILURE()<<"No matching timepoint for epot.Abort!";
        }
        else if(res.data()[0][0]==-1)
        {
            ADD_FAILURE()<<"No matching timepoint for epot.Abort!";
        }
        else
        {
            hsize_t* start2 = new hsize_t[rank1];
            hsize_t* count2 = new hsize_t[rank1];
            hsize_t* stride2 = new hsize_t[rank1];
            hsize_t* block2 = new hsize_t[rank1];
            hsize_t* start3 = new hsize_t[rank1];

            std::cout<<res.size()<<" matching timepoints found in epot\n";
            for(auto index_element:res)
            {
                start2[0]=index_element[0];
                start2[1]=0;
                count2[0]=1;
                count2[1]=dim1[1];
                stride2[0]=1;
                stride2[1]=1;
                block2[0]=1;
                block2[1]=1;
                start3[0]=index_element[1];
                start3[1]=0;

                dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
                dspace2.selectHyperslab(H5S_SELECT_SET, count2, start3, stride2, block2);
                ds1.read(outdat1,tp1,elemspace,dspace1);
                ds2.read(outdat2,tp2,elemspace,dspace2);

                EXPECT_NEAR(outdat1[0],outdat2[0],abstol);
            }
//            delete [] start2;
//            delete [] count2;
//            delete [] stride2;
//            delete [] block2;
//            delete [] start3;
        }
//        delete [] dim1;
//        delete [] dim2;
//        delete [] elem;
//        delete [] start1;
//        delete [] stride1;
//        delete [] block1;
//        delete [] count1;
//        ds1.close();
//        ds2.close();
//        tp1.close();
//        tp2.close();
//        dspace1.close();
//        dspace2.close();
    }

    H5std_string cppname;//!<filename of cpp file
    H5File cppfile;//!<instance of cpp file
    H5std_string pyname;//!<filename of py file
    H5File pyfile;//!<instance of py file
    CompType mytype;//!<predefined type used to write and read data
    H5std_string datablock_group="/datablock_0";//!<name for datablock group. DEFAULT:/datablock_0
    H5std_string energies_group="/energies";//!<name for energies group. DEFAULT:/energies
    H5std_string coeffs_group="/coefficients";//!<name for coefficients group. DEFAULT:/coefficients
    H5std_string norm_group="/norm";//!<name for norm group. DEFAULT:/norms
    H5std_string Pi_group="/Pi";//!<name for Pi group. DEFAULT:/Pi
    H5std_string wavepacket_group="/wavepacket";//!<name for wavepacket group. DEFAULT:/wavepacket
    H5std_string datasetQpath;//!<string for path to dataset Q
    H5std_string datasetPpath;//!<string for path to dataset P
    H5std_string datasetqpath;//!<string for path to dataset p
    H5std_string datasetppath;//!<string for path to dataset q
    H5std_string datasetcpath;//!<string for path to dataset c_0
    H5std_string datasetekinpath;//!<string for path to dataset ekin
    H5std_string datasetepotpath;//!<string for path to dataset epot
    H5std_string datasetnormpath;//!<string for path to dataset norm
    int bool_packet;//!<bool if packet is written in cpp file
    int bool_energy;//!<bool if energy is written in cpp file
    int bool_norm;//!<bool if norm is written in cpp file
    double dt_cpp;//!<timestep in cpp file
    double dt_py;//!<timestep in py file
    H5std_string buff_dt_py="";//!<stringbuffer needed for reading Attribute adt_py
    H5std_string timepathpacket;//!<string for path timegrid packet
    H5std_string timepathnorm;//!<string for path timegrid norm
    H5std_string timepathekin;//!<string for path timegrid ekin
    H5std_string timepathepot;//!<string for path timegrid epot
};

TEST_F(Test2HDFfiles,Testpacket)
{
    if(bool_packet)
    {
        std::vector<int*> res;

        compute_time_matching(res,timepathpacket,timepathpacket);

        Test_Q(res);
        Test_q(res);
        Test_P(res);
        Test_p(res);
        Test_coefficients(res);
    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2HDFfiles,DISABLED_Testenergies)
{
    if(bool_energy)
    {
        std::vector<int*> res;

        compute_time_matching(res,timepathekin,timepathekin);

        Test_ekin(res);

        res.clear();

        compute_time_matching(res,timepathepot,timepathepot);

        Test_epot(res);
    }
    else
    {
        ADD_FAILURE()<<"Cannot compare. Data missing";
    }
}

TEST_F(Test2HDFfiles,DISABLED_Testnorm)
{
    if(bool_norm)
    {
        std::vector<int*> res;

        compute_time_matching(res,timepathnorm,timepathnorm);

        //Test_norm(res);

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

//        if(res.empty())
//        {
//            ADD_FAILURE()<<"No matching timepoint for norm.Abort!";
//        }
//        else if(res.data()[0][0]==-1)
//        {
//            ADD_FAILURE()<<"No matching timepoint for norm.Abort!";
//        }
//        else
//        {
//            std::cout<<res.size()<<" matching timepoints found in norm\n";
//            for(auto index_element:res)
//            {
            for(unsigned int k=0;k<dim1[0];++k)
            {
                hsize_t* start2 = new hsize_t[rank1];
                start2[0]=k;//index_element[0];
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
                hsize_t* start3 = new hsize_t[rank1];
                start3[0]=k;//index_element[1];
                start3[1]=0;

                dspace1.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
                dspace2.selectHyperslab(H5S_SELECT_SET, count2, start3, stride2, block2);
                ds1.read(outdat1,tp1,elemspace,dspace1);
                ds2.read(outdat2,tp2,elemspace,dspace2);

                EXPECT_NEAR(outdat1[0],outdat2[0],abstol);
            }
//            }
//        }
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
