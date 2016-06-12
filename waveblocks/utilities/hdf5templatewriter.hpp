#pragma once

//need H5 Cpp header
#include "H5Cpp.h"
//using strings
#include <string>

#include <Eigen/Core>

#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"

namespace waveblocks
{
    namespace utilities
    {

    //using namespaces for convenience
    using namespace H5;

    /**
     * \brief The ctype struct for writing complex numbers
     *
     * The ctype struct is used for defining an H5:CompType which is used to write complex numbers
     * in a simplified manner which is compatible with the python HDF interface. Also the default
     * value for this type is defined as 0 + 0*i which is used for the allocation of memory in the
     * HDF format
     */
    struct ctype{
        double real=0.;
        double imag=0.;
    }instanceof;

    /**
     * \brief Our HDF5 writer class
     *
     * This class is templated on the Dimension D which is also used to also define the dimensionality
     * of the HaWpBasisVector
     */
    template<int D>
    class hdf5writertemplate
    {
        //declare friend for easy use
        friend struct ctype;
    public:
        /**
         *\brief Construction of hdf5 writer class with a name
         *\param name for construction
         *
         * Construction of hdf5 writer from std::string \param name where also implicitly the definition
         * of our H5:CompType is hidden. Also the H5::File gets constructed by the \param name where
         * it is open with H5F_ACC_TRUNC value
         * */
        hdf5writertemplate(std::string name):filename_(name),mytype_(sizeof(instanceof)),file_(filename_,H5F_ACC_TRUNC)
        {
            //fixed writing type for easy overlap with python interface
            mytype_.insertMember( "r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
            mytype_.insertMember( "i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);

            //set up chunk dimension
            set_chunk_dim();

            //set up elem space
            set_elem_space();

            //set group structure
            set_group_structure();

            //
            select_elem_hyperslabs();
        }
        /**
         * @brief sets bool value for writing the packet
         * @param flag
         */
        void set_write_packet(bool flag) const
        {
            wlist:packet=flag;
        }
        /**
         * @brief sets bool value for writing the energies
         * @param flag
         */
        void set_write_energy(bool flag) const
        {
            wlist:energy=flag;
        }
        /**
         * @brief sets bool value for writing coefficients
         * @param flag
         */
        void set_write_coefficients(bool flag) const
        {
            wlist:coefficients=flag;
        }
        /**
         * @brief sets bool value for writing timegrid
         * @param flag
         */
        void set_write_timegrid(bool flag) const
        {
            wlist:timegrid=flag;
        }
        /**
         * @brief set string for root group
         * @param name
         */
        void set_datablockstring(std::string name)
        {
            datablock_string=name;
        }
        /**
         * @brief set string for wavepacket group
         * @param name
         */
        void set_wavepacketstring(std::string name)
        {
            wavepacket_group_string=name;
        }
        /**
         * @brief set string for packet group
         * @param name
         */
        void set_packetgroupstring(std::string name)
        {
            packet_group_string=name;
        }
        /**
         * @brief set string for coefficient group
         * @param name
         */
        void set_coefficientgroupstring(std::string name)
        {
        coefficient_group_string=name;
        }
        /**
         * @brief set string for energy group
         * @param name
         */
        void set_energygroupstring(std::string name)
        {
        energy_group_string=name;
        }
        /**
         * \brief Set up chunk dimension for the written variables
         *
         * The HDF Interface needs that the PropList exits with chunked dimension otherwise
         * it is not possible to extend this particular DataSet
         */
        void set_chunk_dim(void) const
        {
            constexpr int dim = D; //alias for D
            if(packet)
            {
                hsize_t chunk_dims1[3]={1,dim,1};
                plist_qp.setChunk(RANK3,chunk_dims1);
                plist_qp.setFillValue(mytype_,&instanceof);

                hsize_t chunk_dims2[3]={1,dim,dim};
                plist_QP.setChunk(RANK3,chunk_dims2);
                plist_QP.setFillValue(mytype_,&instanceof);

                hsize_t chunk_dims3[2]={1,dim-1};
                plist_S.setChunk(RANK2,chunk_dims3);
                plist_S.setFillValue(mytype_,&instanceof);
            }
            if(energy)
            {
                hsize_t chunk_dims4[2]={1,3};
                plist_energy.setChunk(RANK2,chunk_dims4);
                plist_energy.setFillValue(PredType::NATIVE_DOUBLE,&dref);
            }
            if(coefficients)
            {
                hsize_t chunk_dims5[2]={1,16};
                plist_c.setChunk(RANK2,chunk_dims5);
                plist_c.setFillValue(mytype_,&instanceof);

            }
            if(timegrid)
            {
                hsize_t chunk_dims6[1]={1};
                plist_time.setChunk(RANK1,chunk_dims6);
                plist_time.setFillValue(PredType::NATIVE_DOUBLE,&dref);
            }
        }
        /**
         *\brief Set up elem space for HDF interface
         *
         *The HDF Interface needs a DataSpace for elements written from data to the file.
         *A dataspace is constructet from a rank identifier and corresponding hsize_t dimension array
         * which has to be of size rank
         */
        void set_elem_space(void)
        {
            constexpr int dim = D;
            if(packet)
            {
                qpelem[0]=1;
                qpelem[1]=dim;
                qpelem[2]=1;
                DataSpace t1(RANK3,qpelem);
                qpelemspace=t1;
                QPelem[0]=1;
                QPelem[1]=dim;
                QPelem[2]=dim;
                DataSpace t2(RANK3,QPelem);
                QPelemspace=t2;
                Selem[0]=1;
                Selem[1]=1;
                DataSpace t3(RANK2,Selem);
                Selemspace=t3;
            }
            if(energy)
            {
                energyelem[0]=1;
                energyelem[1]=3;
                DataSpace t4(RANK2,energyelem);
                energyelemspace=t4;
            }
            if(coefficients)
            {
                celem[0]=1;
                celem[1]=16;
                DataSpace t5(RANK2,celem);
                celemspace=t5;
            }
            if(timegrid)
            {
                timeelem[0]=1;
                DataSpace t6(RANK1,timeelem);
                timelemspace=t6;
            }
        }
        /**
         * @brief set up group structure in file
         */
        void set_group_structure(void)
        {
            H5std_string a1=datablock_string;
            gblock = std::make_shared<Group>(file_.createGroup(a1));
            if(packet&&coefficients)
            {
                H5std_string a2=(datablock_string+wavepacket_group_string);
                gpacket = std::make_shared<Group>(file_.createGroup(a2));
                H5std_string a3=(datablock_string+wavepacket_group_string+packet_group_string);
                gPi = std::make_shared<Group>(file_.createGroup(a3));
                H5std_string a4=(datablock_string+wavepacket_group_string+coefficient_group_string);
                gcoefficient = std::make_shared<Group>(file_.createGroup(a4));
            }
            else if(packet)
            {
                H5std_string a2=(datablock_string+wavepacket_group_string);
                gpacket = std::make_shared<Group>(file_.createGroup(a2));
                H5std_string a3=(datablock_string+wavepacket_group_string+packet_group_string);
                gPi = std::make_shared<Group>(file_.createGroup(a3));
            }
            else if(coefficients)
            {
                H5std_string a2=(datablock_string+wavepacket_group_string);
                gpacket = std::make_shared<Group>(file_.createGroup(a2));
                H5std_string a4=(datablock_string+wavepacket_group_string+coefficient_group_string);
                gcoefficient = std::make_shared<Group>(file_.createGroup(a4));
            }

            if(energy)
            {
                H5std_string a5=(datablock_string+energy_group_string);
                genergy = std::make_shared<Group>(file_.createGroup(a5));
            }
            if(timegrid)
            {
                H5std_string a6=(datablock_string+timegrid_group_string);
                gtimegrid = std::make_shared<Group>(file_.createGroup(a6));
            }
        }
        /**
         * \brief From HDF interface select element hyperslabs
         *
         * These elements are selected such that in every timestep it can be written to file
         */
        void select_elem_hyperslabs(void)
        {
            constexpr int dim = D;

            if(packet)
            {
                //qp
                hsize_t count1[3]={1,dim,1};
                hsize_t start1[3]={0,0,0};
                hsize_t stride1[3]={1,1,1};
                hsize_t block1[3]={1,1,1};
                //QP
                hsize_t count2[3]={1,dim,dim};
                hsize_t start2[3]={0,0,0};
                hsize_t stride2[3]={1,1,1};
                hsize_t block2[3]={1,1,1};
                //S
                hsize_t count3[2]={1,1};
                hsize_t start3[2]={0,0};
                hsize_t stride3[2]={1,1};
                hsize_t block3[2]={1,1};

                //qp
                qpelemspace.selectHyperslab(H5S_SELECT_SET, count1, start1, stride1, block1);
                //QP
                QPelemspace.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
                //Selem
                Selemspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);

            }
            if(coefficients)
            {
                //coefficients
                hsize_t count4[2]={1,16};
                hsize_t start4[2]={0,0};
                hsize_t stride4[2]={1,1};
                hsize_t block4[2]={1,1};
                celemspace.selectHyperslab(H5S_SELECT_SET, count4, start4, stride4, block4);
            }
            if(energy)
            {
                //energies
                hsize_t count5[2]={1,3};
                hsize_t start5[2]={0,0};
                hsize_t stride5[2]={1,1};
                hsize_t block5[2]={1,1};
                energyelemspace.selectHyperslab(H5S_SELECT_SET, count5, start5, stride5, block5);
            }
            if(timegrid)
            {
                //timegrid
                hsize_t count6[1]={1};
                hsize_t start6[1]={0};
                hsize_t stride6[1]={1};
                hsize_t block6[1]={1};
                timelemspace.selectHyperslab(H5S_SELECT_SET, count6, start6, stride6, block6);
            }
        }
        /**
         * @brief increase index
         */
        void increase_index(void)
        {
            current_index+=1;
        }
        /**
         * @brief decrease index
         */
        void decrease_index(void)
        {
            current_index-=1;
        }
        /**
         * @brief allocate space for datasets
         */
        void allocate_datasets(void)
        {
            if(packet)
            {
                qs=std::make_shared<DataSet>(file_.createDataSet(q,mytype_,qspace,plist_qp));
                ps=std::make_shared<DataSet>(file_.createDataSet(p,mytype_,qspace,plist_qp));
                Qs=std::make_shared<DataSet>(file_.createDataSet(Q,mytype_,qspace,plist_QP));
                Ps=std::make_shared<DataSet>(file_.createDataSet(P,mytype_,qspace,plist_QP));
                Ss=std::make_shared<DataSet>(file_.createDataSet(S,mytype_,Sspace,plist_S));
            }
            if(coefficients)
            {
                coeffs=std::make_shared<DataSet>(file_.createDataSet(c,mytype_,cspace,plist_c));
            }
            if(energy)
            {
                energys=std::make_shared<DataSet>(file_.createDataSet(energies,PredType::NATIVE_DOUBLE,energyspace,plist_energy));
            }
            if(timegrid)
            {
                times=std::make_shared<DataSet>(file_.createDataSet(time,PredType::NATIVE_DOUBLE,timespace,plist_time));
            }
        }

    private:
        H5std_string filename_;/**<identifier for filename*/
        CompType mytype_;/**<Declaration of H5:CompType member */
        H5File file_; /**<H5File member */
        enum wlist:bool{packet=1,energy=0,coefficients=1,timegrid=0};/**<enum for storing writing bool values*/

        double dref=0.;/**<fillvalue for energys for allocation */
        H5std_string packet_group_string="/Pi";/**<String for H5Group to save packet to. Default:Pi*/
        H5std_string datablock_string="/datablock_0";/**<String for H5Group for datablock.default. Datablock_0*/
        H5std_string coefficient_group_string="/coefficients";/**<String for H5Group of coefficients. Default:coefficients*/
        H5std_string wavepacket_group_string="/wavepacket";/**<String for H5Group for packet and coefficients. Default:wavepacket*/
        H5std_string energy_group_string="/energies";/**<String for H5Group for energies. Default:energies*/
        H5std_string timegrid_group_string="/timegrid";/**<String for H5Group for timegrid. Default:timegrid*/

        DSetCreatPropList plist_qp;/**<PropList for packet.q() packet.p()*/
        DSetCreatPropList plist_QP;/**<PropList for packet.Q() packet.P()*/
        DSetCreatPropList plist_S;/**<PropList for packet.S()*/
        DSetCreatPropList plist_energy;/**<PropList for energies*/
        DSetCreatPropList plist_c;/**<PropList for coefficients*/
        DSetCreatPropList plist_time;/**<PropList for timegrid*/

        const int RANK1=1;/**<rank 1 identifier */
        const int RANK2=2;/**<rank 2 identifier */
        const int RANK3=3;/**<rank 3 identifier */

        const hsize_t maxdims1[1]={H5S_UNLIMITED};/**<max dim identifier for rank1 for extension */
        const hsize_t maxdims2[2]={H5S_UNLIMITED,H5S_UNLIMITED};/**<max dim identifier for rank2 for extension */
        const hsize_t maxdims3[3]={H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};/**<max dim identifier for rank3 for extension */

        hsize_t exqp[3];/**<extenstion for packet.q() packet.p() */
        hsize_t exQP[3];/**<extenstion for packet.Q() packet.P() */
        hsize_t exS[2];/**<extenstion for packet.S()*/
        hsize_t exenergy[2];/**<extenstion for energies*/
        hsize_t exc[2];/**<extenstion for coefficients*/
        hsize_t extime[1];/**<extenstion for timegrid*/

        hsize_t qpelem[3]; /**<size of q,p element written from program to file needed by HDF interface */
        DataSpace qpelemspace;/**<space of q,p element written from program to file needed by HDF interface */
        hsize_t QPelem[3];/**<size of Q,P element written from program to file needed by HDF interface */
        DataSpace QPelemspace;/**<space of Q,P element written from program to file needed by HDF interface */
        hsize_t Selem[2];/**<size of S element written from program to file needed by HDF interface */
        DataSpace Selemspace;/**<space of S element written from program to file needed by HDF interface */
        hsize_t energyelem[2];/**<size of energy element written from program to file needed by HDF interface */
        DataSpace energyelemspace;/**<space of energy element written from program to file needed by HDF interface */
        hsize_t timeelem[1];/**<size of timegrid element written from program to file needed by HDF interface */
        DataSpace timelemspace;/**<space of timegrid element written from program to file needed by HDF interface */
        hsize_t celem[2];/**<size of coefficient element written from program to file needed by HDF interface */
        DataSpace celemspace;/**<space of coefficient element written from program to file needed by HDF interface */

        int current_index=1;/**<current index used to determine position in file in time dimension*/

        std::shared_ptr<Group> gblock;/**<group for datablock */
        std::shared_ptr<Group> gpacket;/**<group for packet */
        std::shared_ptr<Group> gPi;/**<group for matrices in packet */
        std::shared_ptr<Group> gcoefficient;/**<group for coefficients in packet */
        std::shared_ptr<Group> genergy;/**<group for energies */
        std::shared_ptr<Group> gtimegrid;/**<group for timegrid */

        H5std_string q="q";/**<name for packet.q()*/
        DataSpace qspace;/**<space for packet.q() in file*/
        std::shared_ptr<DataSet> qs;/**<dataset for packet.q() in file*/
        H5std_string p="p";/**<name for packet.p()*/
        DataSpace pspace;/**<space for packet.p() in file*/
        std::shared_ptr<DataSet> ps;/**<dataset for packet.p() in file*/
        H5std_string Q="Q";/**<name for packet.Q()*/
        DataSpace Qspace;/**<space for packet.Q() in file*/
        std::shared_ptr<DataSet> Qs;/**<dataset for packet.Q() in file*/
        H5std_string P="P";/**<name for packet.P()*/
        DataSpace Pspace;/**<space for packet.P() in file*/
        std::shared_ptr<DataSet> Ps;/**<dataset for packet.P() in file*/
        H5std_string S="S";/**<name for packet.S()*/
        DataSpace Sspace;/**<space for packet.S() in file*/
        std::shared_ptr<DataSet> Ss;/**<dataset for packet.S() in file*/
        H5std_string c="c_0";/**<name for coefficients*/
        DataSpace cspace;/**<space for coefficients in file*/
        std::shared_ptr<DataSet> coeffs;/**<dataset for coefficients in file*/
        H5std_string time="timegrid";/**<name for timegrid*/
        DataSpace timespace;/**<space for timegrid in file*/
        std::shared_ptr<DataSet> times;/**<dataset for timegrid in file*/
        H5std_string energies="energies";/**<name for energies*/
        DataSpace energyspace;/**<space for energies in file*/
        std::shared_ptr<DataSet> energys;/**<dataset for energies in file*/
    };


    }
}
