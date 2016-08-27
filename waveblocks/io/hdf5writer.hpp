#pragma once

//need H5 Cpp header for HDF interface
#include "H5Cpp.h"

#include <string>
#include <map>
#include <cstdio>
#include <vector>

#include <Eigen/Core>

#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"
#include "waveblocks/types.hpp"
#include "waveblocks/wavepackets/shapes/tiny_multi_index.hpp"

namespace waveblocks
{
    namespace io
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
     * of the HaWpBasisVector. Also after ctor is called there has to follow a prestructuring<MultiIndex>()
     * for correct construction
     */
    template<int D>
    class hdf5writer
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
        hdf5writer(std::string name):filename_(name),mytype_(sizeof(instanceof)),file_(filename_,H5F_ACC_TRUNC)
        {
            //fixed writing type for easy overlap with python interface
            mytype_.insertMember( "r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
            mytype_.insertMember( "i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);
        }
        /**
         * \brief runtime function to evaluate number of coefficients
         * \param packet to deduce number of coefficients
         * Used in prestructuring<MultiIndex>
         */
        template<class MultiIndex>
        void set_coeff_dim(waveblocks::wavepackets::ScalarHaWp<D,MultiIndex> packet)
        {
            Eigen::Matrix<complex_t, Eigen::Dynamic, 1> coeffs=utilities::PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packet);
            int a=coeffs.rows();
            int b=coeffs.cols();
            csize=a*b;
        }
        /**
         * \brief prestructure after knowing bool values and packet
         * \param packet used for evaluation number of coefficients
         * \param dt to save timestepsize as attribute
         */
        template<class MultiIndex>
        void prestructuring(waveblocks::wavepackets::ScalarHaWp<D,MultiIndex> packet,double dt)
        {
            //get number of coefficients
            set_coeff_dim(packet);

            //set group structure
            set_group_structure(dt);

            //set up chunk dimension
            set_chunk_dim();

            //set up space for file
            set_file_dataspace();

            //allocate datasets needs to be after select writespace
            allocate_datasets();

            //set up elementary spaces
            set_elem_space();

            //select space within elemtary space
            select_elem_hyperslabs();
        }
        /**
         * @brief sets bool value for writing the packet(matrices and coefficients)
         * @param flag
         */
        void set_write_packet(bool flag)
        {
            wrlist["packet"]=flag;
        }
        /**
         * @brief sets bool value for writing the energies
         * @param flag
         */
        void set_write_energies(bool flag)
        {
            wrlist["energy"]=flag;
        }
        /**
         * @brief sets bool value for writing norms
         * @param flag
         */
        void set_write_norm(bool flag)
        {
            wrlist["norm"]=flag;
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
        void set_energiesgroupstring(std::string name)
        {
            energies_group=name;
        }
        /**
         * @brief set string for norm group
         * @param name
         */
        void set_normsgroupstring(std::string name)
        {
            norms_group=name;
        }
        /**
         * @brief set string for observables group
         * @param name
         */
        void set_observablesgroupstring(std::string name)
        {
            observables_group=name;
        }
        /**
         * @brief set timestep for writing packet
         * @param timestep
         */
        void set_timestep_packet(int timestep)
        {
            timestepsize_packet=timestep;
        }
        /**
         * @brief set timestep for writing norms
         * @param timestep
         */
        void set_timestep_norm(int timestep)
        {
            timestepsize_norms=timestep;
        }
        /**
         * @brief set timestep for writing energies
         * @param timestep
         */
        void set_timestep_energies(int timestep)
        {
            timestepsize_epot=timestep;
            timestepsize_ekin=timestep;
        }
        /**
         * @brief set timestep for writing epot
         * @param timestep
         */
        void set_timestep_epot(int timestep)
        {
            timestepsize_epot=timestep;
        }
        /**
         * @brief set timestep for writing ekin
         * @param timestep
         */
        void set_timestep_ekin(int timestep)
        {
            timestepsize_ekin=timestep;
        }
        /**
         * \brief Set up chunk dimension for the written variables
         *
         * The HDF Interface needs that the PropList exits with chunked dimension otherwise
         * it is not possible to extend this particular DataSet
         */
        void set_chunk_dim(void)
        {
            constexpr int dim = D; //alias for D
            //time chunk dimension the same for all
            hsize_t chunk_dims6[]={1};
            plist_time.setChunk(RANK1,chunk_dims6);
            plist_time.setFillValue(PredType::NATIVE_INT,&iref);

            if(wrlist["packet"])
            {
                hsize_t chunk_dims1[]={1,dim,1};
                plist_qp.setChunk(RANK3,chunk_dims1);
                plist_qp.setFillValue(mytype_,&instanceof);

                hsize_t chunk_dims2[]={1,dim,dim};
                plist_QP.setChunk(RANK3,chunk_dims2);
                plist_QP.setFillValue(mytype_,&instanceof);

                hsize_t chunk_dims3[]={1,1,1};
                plist_S.setChunk(RANK3,chunk_dims3);
                plist_S.setFillValue(mytype_,&instanceof);

                hsize_t chunk_dims5[2];
                chunk_dims5[0]=1;
                chunk_dims5[1]=csize;
                plist_c.setChunk(RANK2,chunk_dims5);
                plist_c.setFillValue(mytype_,&instanceof);
            }
            if(wrlist["energy"])
            {
                hsize_t chunk_dims4[]={1,1};
                plist_energy.setChunk(RANK2,chunk_dims4);
                plist_energy.setFillValue(PredType::NATIVE_DOUBLE,&dref);

            }
            if(wrlist["norm"])
            {
                hsize_t chunk_dims6[]={1,1};
                plist_norms.setChunk(RANK2,chunk_dims6);
                plist_norms.setFillValue(PredType::NATIVE_DOUBLE,&dref);
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
            //timeelem space const for all timesteps
            timeelem[0]=1;
            DataSpace t6(RANK1,timeelem);
            timelemspace=t6;

            constexpr int dim = D;
            if(wrlist["packet"])
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
                Selem[2]=1;
                DataSpace t3(RANK3,Selem);
                Selemspace=t3;

                celem[0]=1;
                celem[1]=csize;
                DataSpace t5(RANK2,celem);
                celemspace=t5;
            }
            if(wrlist["energy"])
            {
                energyelem[0]=1;
                energyelem[1]=1;
                DataSpace t4(RANK2,energyelem);
                energyelemspace=t4;
            }
            if(wrlist["norm"])
            {
                normelem[0]=1;
                normelem[1]=1;
                DataSpace t6(RANK2,normelem);
                normelemspace=t6;
            }
        }
        /**
         * \brief set up group structure in file
         * \param dt for writing global timestep
         *
         * This is needed for the correct subtructure of all datasets. It also sets Attributes in datablock
         * to save which datasets were written and its correspondig global timestep.
         */
        void set_group_structure(double dt)
        {
            H5std_string a1=datablock_string;
            gblock = std::make_shared<Group>(file_.createGroup(a1));
            hsize_t q1[]={1};
            DataSpace s1(RANK1,q1);
            apacket=gblock->createAttribute("packet",PredType::NATIVE_INT,s1);
            aenergy=gblock->createAttribute("energy",PredType::NATIVE_INT,s1);
            anorm=gblock->createAttribute("norm",PredType::NATIVE_INT,s1);
            adt=gblock->createAttribute("dt",PredType::NATIVE_DOUBLE,s1);
            adt.write(PredType::NATIVE_DOUBLE,&dt);

            if(wrlist["packet"])
            {
                int ap=1;
                apacket.write(PredType::NATIVE_INT,&ap);
                H5std_string a2=(datablock_string+wavepacket_group_string);
                gpacket = std::make_shared<Group>(file_.createGroup(a2));
                H5std_string a3=(datablock_string+wavepacket_group_string+packet_group_string);
                gPi = std::make_shared<Group>(file_.createGroup(a3));
                H5std_string a4=(datablock_string+wavepacket_group_string+coefficient_group_string);
                gcoefficient = std::make_shared<Group>(file_.createGroup(a4));
            }
            else
            {
                int ap=0;
                apacket.write(PredType::NATIVE_INT,&ap);
            }
            H5std_string a7=datablock_string+observables_group;
            gobservables=std::make_shared<Group>(file_.createGroup(a7));
            if(wrlist["energy"])
            {
                int ae=1;
                aenergy.write(PredType::NATIVE_INT,&ae);
                H5std_string a5 = (datablock_string+observables_group+energies_group);
                genergy = std::make_shared<Group>(file_.createGroup(a5));
            }
            else
            {
                int ae=0;
                aenergy.write(PredType::NATIVE_INT,&ae);
            }
            if(wrlist["norm"])
            {
                int an=1;
                anorm.write(PredType::NATIVE_INT,&an);
                H5std_string a6 = (datablock_string+observables_group+norms_group);
                gnorms = std::make_shared<Group>(file_.createGroup(a6));
            }
            else
            {
                int an=0;
                anorm.write(PredType::NATIVE_INT,&an);
            }
        }
        /**
         * \brief From HDF interface select element hyperslabs
         *
         * These elements are selected such that in every timestep it can be written to file.
         * This is onlye called once because the elemantary spaces should never change.
         */
        void select_elem_hyperslabs(void)
        {
            constexpr int dim = D;

            hsize_t count6[]={1};
            hsize_t start6[]={0};
            hsize_t stride6[]={1};
            hsize_t block6[]={1};
            timelemspace.selectHyperslab(H5S_SELECT_SET, count6, start6, stride6, block6);

            if(wrlist["packet"])
            {
                //qp
                hsize_t count1[]={1,dim,1};
                hsize_t start1[]={0,0,0};
                hsize_t stride1[]={1,1,1};
                hsize_t block1[]={1,1,1};
                //QP
                hsize_t count2[]={1,dim,dim};
                hsize_t start2[]={0,0,0};
                hsize_t stride2[]={1,1,1};
                hsize_t block2[]={1,1,1};
                //S
                hsize_t count3[]={1,1,1};
                hsize_t start3[]={0,0,0};
                hsize_t stride3[]={1,1,1};
                hsize_t block3[]={1,1,1};

                //qp
                qpelemspace.selectHyperslab(H5S_SELECT_SET, count1, start1, stride1, block1);
                //QP
                QPelemspace.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
                //Selem also for adQ
                Selemspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);

                //coefficients
                hsize_t count4[2];
                count4[0]=1;
                count4[1]=csize;
                hsize_t start4[]={0,0};
                hsize_t stride4[]={1,1};
                hsize_t block4[]={1,1};
                celemspace.selectHyperslab(H5S_SELECT_SET, count4, start4, stride4, block4);
            }
            if(wrlist["energy"])
            {
                hsize_t count5[]={1,1};
                hsize_t start5[]={0,0};
                hsize_t stride5[]={1,1};
                hsize_t block5[]={1,1};
                energyelemspace.selectHyperslab(H5S_SELECT_SET, count5, start5, stride5, block5);
            }
            if(wrlist["norm"])
            {
                hsize_t count6[]={1,1};
                hsize_t start6[]={0,0};
                hsize_t stride6[]={1,1};
                hsize_t block6[]={1,1};
                normelemspace.selectHyperslab(H5S_SELECT_SET, count6, start6, stride6, block6);
            }
        }

        /**
         * @brief set dataspace used in file
         *
         * Initial construction of the spaces for all sets.
         * Also called only once because later with extensions it
         * has to be resett in every extension call with new sizes.
         */
        void set_file_dataspace(void)
        {
            constexpr int dim = D;
            if(wrlist["packet"])
            {
                //set up q and p
                hsize_t dim1[]={1,dim,1};
                DataSpace d1(RANK3,dim1,maxdims3);
                qspace=d1;
                pspace=d1;
                //set up Q and P
                hsize_t dim2[]={1,dim,dim};
                DataSpace d2(RANK3,dim2,maxdims3);
                Qspace=d2;
                Pspace=d2;
                //set up S and adQ
                hsize_t dim3[]={1,1,1};
                DataSpace d3(RANK3,dim3,maxdims3);
                Sspace=d3;
                adQspace=d3;

                //set up coefficients
                hsize_t dim5[2];
                dim5[0]=1;
                dim5[1]=csize;
                DataSpace d5(RANK2,dim5,maxdims2);
                cspace=d5;

                hsize_t dimt1[]={1};
                DataSpace dt1(RANK1,dimt1,maxdims1);
                timespace_packet=dt1;
            }
            if(wrlist["energy"])
            {
                hsize_t dim4[]={1,1};
                DataSpace d4(RANK2,dim4,maxdims2);
                energyspace_ekin=d4;
                energyspace_epot=d4;

                hsize_t dimt2[]={1};
                DataSpace dt2(RANK1,dimt2,maxdims1);
                timespace_ekin=dt2;
                timespace_epot=dt2;
            }
            if(wrlist["norm"])
            {
                hsize_t dim5[]={1,1};
                DataSpace d5(RANK2,dim5,maxdims2);
                normspace=d5;

                hsize_t dimt3[]={1};
                DataSpace dt3(RANK1,dimt3,maxdims1);
                timespace_norms=dt3;
            }
        }
        /**
         * @brief allocate space for datasets
         *
         * Allocate all datasets in file depending on the set bool values. This is called once
         * and uses shared_ptr for easy use no garbage collection needed.
         */
        void allocate_datasets(void)
        {
            if(wrlist["packet"])
            {
                H5std_string pack=datablock_string+wavepacket_group_string+packet_group_string;
                qs=std::make_shared<DataSet>(file_.createDataSet(pack+"/q",mytype_,qspace,plist_qp));
                ps=std::make_shared<DataSet>(file_.createDataSet(pack+"/p",mytype_,pspace,plist_qp));
                Qs=std::make_shared<DataSet>(file_.createDataSet(pack+"/Q",mytype_,Qspace,plist_QP));
                Ps=std::make_shared<DataSet>(file_.createDataSet(pack+"/P",mytype_,Pspace,plist_QP));
                Ss=std::make_shared<DataSet>(file_.createDataSet(pack+"/S",mytype_,Sspace,plist_S));
                adQs=std::make_shared<DataSet>(file_.createDataSet(pack+"/adQ",mytype_,adQspace,plist_S));
                H5std_string tpack=datablock_string+wavepacket_group_string;
                times_packet=std::make_shared<DataSet>(file_.createDataSet(tpack+"/timegrid",PredType::NATIVE_INT,timespace_packet,plist_time));

                H5std_string cff=datablock_string+wavepacket_group_string+coefficient_group_string;
                coeffs=std::make_shared<DataSet>(file_.createDataSet(cff+"/c_0",mytype_,cspace,plist_c));
            }
            if(wrlist["energy"])
            {
                H5std_string ename=datablock_string+observables_group+energies_group;
                H5std_string enameepot=ename+"/epot";
                energys_epot=std::make_shared<DataSet>(file_.createDataSet(enameepot,PredType::NATIVE_DOUBLE,energyspace_epot,plist_energy));
                H5std_string enameekin=ename+"/ekin";
                energys_ekin=std::make_shared<DataSet>(file_.createDataSet(enameekin,PredType::NATIVE_DOUBLE,energyspace_ekin,plist_energy));
                H5std_string engtime1=ename+"/timegrid_epot";
                times_epot=std::make_shared<DataSet>(file_.createDataSet(engtime1,PredType::NATIVE_INT,timespace_epot,plist_time));
                H5std_string engtime2=ename+"/timegrid_ekin";
                times_ekin=std::make_shared<DataSet>(file_.createDataSet(engtime2,PredType::NATIVE_INT,timespace_ekin,plist_time));
            }
            if(wrlist["norm"])
            {
                H5std_string tmm=datablock_string+observables_group+norms_group+"/norm";
                normss=std::make_shared<DataSet>(file_.createDataSet(tmm,PredType::NATIVE_DOUBLE,normspace,plist_norms));
                H5std_string tmt=datablock_string+observables_group+norms_group+"/timegrid";
                times_norms=std::make_shared<DataSet>(file_.createDataSet(tmt,PredType::NATIVE_INT,timespace_norms,plist_time));
            }
        }
        /**
         * @brief setup extension for norms
         *
         * Called after a subset was written.
         */
        void setup_extension_norm(void)
        {
            if(wrlist["norm"])
            {
                exnorms[0]=index_norm;
                exnorms[1]=1;
                ex_timegrid_norms[0]=index_norm;
            }
            else
            {
            std::cout<<"ERROR: setup_extension_norms called with bool=false\n";
            }
        }
        /**
         * @brief setup extension for epot
         *
         * Called after a subset was written.
         */
        void setup_extension_epot(void)
        {
            if(wrlist["energy"])
            {
                exepot[0]=index_epot;
                exepot[1]=1;
                ex_timegrid_epot[0]=index_epot;
            }
            else
            {
                std::cout<<"ERROR: setup_extension_epot called with bool=false\n";
            }
        }
        /**
         * @brief setup extension for ekin
         *
         * Called after a subset was written.
         */
        void setup_extension_ekin(void)
        {
            if(wrlist["energy"])
            {
                exekin[0]=index_ekin;
                exekin[1]=1;
                ex_timegrid_ekin[0]=index_ekin;
            }
            else
            {
                std::cout<<"ERROR: setup_extension_ekin called with bool=false\n";
            }
        }
        /**
         * @brief setup extension for packet
         *
         * Called after a subset was written.
         */
        void setup_extension_packet(void)
        {
            constexpr int dim=D;
            if(wrlist["packet"])
            {
                exqp[0]=index_packet;
                exqp[1]=dim;
                exqp[2]=1;
                exQP[0]=index_packet;
                exQP[1]=dim;
                exQP[2]=dim;
                exS[0]=index_packet;
                exS[1]=1;
                exS[2]=1;

                exc[0]=index_packet;
                exc[1]=csize;
                ex_timegrid_packet[0]=index_packet;
            }
            else
            {
                std::cout<<"ERROR: setup_extension_packet called with bool=false\n";
            }
        }
        /**
         * @brief extend dataset for next timestep for packet
         *
         * Called after a subset was written and extension is set.
         */
        void extend_dataset_packet(void)
        {
            if(wrlist["packet"])
            {
                //qp
                qs->extend(exqp);
                ps->extend(exqp);
                //QP
                Qs->extend(exQP);
                Ps->extend(exQP);
                //Senergies
                Ss->extend(exS);
                adQs->extend(exS);

                coeffs->extend(exc);
                times_packet->extend(ex_timegrid_packet);
            }
            else
            {
                std::cout<<"ERROR: extend_dataset_packet called with bool=false\n";
            }
        }
        /**
         * @brief extend dataset for next timestep for norms
         *
         * Called after a subset was written and extension is set.
         */
        void extend_dataset_norm(void)
        {
            if(wrlist["norm"])
            {
                normss->extend(exnorms);
                times_norms->extend(ex_timegrid_norms);
            }
            else
            {
                std::cout<<"ERROR: extend_dataset_norms called with bool=false\n";
            }
        }
        /**
         * @brief extend dataset for next timestep for epot
         *
         * Called after a subset was written and extension is set.
         */
        void extend_dataset_epot(void)
        {
            if(wrlist["energy"])
            {
                energys_epot->extend(exepot);
                times_epot->extend(ex_timegrid_epot);
            }
            else
            {
                std::cout<<"ERROR: extend_dataset_epot called with bool=false\n";
            }
        }
        /**
         * @brief extend dataset for next timestep for ekin
         *
         * Called after a subset was written and extension is set.
         */
        void extend_dataset_ekin(void)
        {
            if(wrlist["energy"])
            {
                energys_ekin->extend(exekin);
                times_ekin->extend(ex_timegrid_ekin);
            }
            else
            {
                std::cout<<"ERROR: extend_dataset_ekin called with bool=false\n";
            }
        }
        /**
         * @brief update ekin filespace
         *
         * Needs to be done after extending this dataset
         */
        void update_filespace_ekin(void)
        {
            if(wrlist["energy"])
            {
                energyspace_ekin=energys_ekin->getSpace();
                timespace_ekin=times_ekin->getSpace();
            }
            else
            {
                std::cout<<"ERROR: update_filespace_ekin called with bool=false\n";
            }
        }
        /**
         * @brief update epot filespace
         *
         * Needs to be done after extending this dataset
         */
        void update_filespace_epot(void)
        {
            if(wrlist["energy"])
            {
                energyspace_epot=energys_epot->getSpace();
                timespace_epot=times_epot->getSpace();
            }
            else
            {
                std::cout<<"ERROR: update_filespace_epot called with bool=false\n";
            }
        }
        /**
         * @brief update norms filespace
         *
         * Needs to be done after extending this dataset
         */
        void update_filespace_norm(void)
        {
            if(wrlist["norm"])
            {
                normspace=normss->getSpace();
                timespace_norms=times_norms->getSpace();
            }
            else
            {
                std::cout<<"ERROR: update_filespace_norms called with bool=false\n";
            }
        }
        /**
         * @brief update packet filespace
         *
         * Needs to be done after extending this dataset
         */
        void update_filespace_packet(void)
        {
            if(wrlist["packet"])
            {
                //qp
                qspace=qs->getSpace();
                pspace=ps->getSpace();
                //QP
                Qspace=Qs->getSpace();
                Pspace=Ps->getSpace();
                //S
                Sspace=Ss->getSpace();
                adQspace=adQs->getSpace();
                cspace=coeffs->getSpace();
                timespace_packet=times_packet->getSpace();
            }
            else
            {
                std::cout<<"ERROR: update_filespace_packet called with bool=false\n";
            }
        }
        /**
         * \brief transform Eigen::Matrix<complex_t,D,D> into std::vector<ctype>
         * \param newdata
         * \param mat
         *
         * Transforms an Eigen::Matrix<complex_t,D,D> into a writeable ctype*
         */
        void transform(std::vector<ctype>& newdata,Eigen::Matrix<complex_t,D,D> mat)
        {
            constexpr int dim=D;
            newdata.resize(dim*dim);
            std::complex<double>* tmp(mat.data());
            for(int p=0;p<dim*dim;++p)
            {
                newdata[p].real=tmp[p].real();
                newdata[p].imag=tmp[p].imag();
            }
        }
        /**
         * \brief transform Eigen::Matrix<real_t,D,1> into std::vector<ctype>
         * \param newdata
         * \param mat
         *
         * Transforms an Eigen::Matrix<real_t,D,1> into a writable std::vector<ctype>
         */
        void transform(std::vector<ctype>& newdata,Eigen::Matrix<real_t,D,1> mat)
        {
            constexpr int dim=D;
            newdata.resize(dim);
            double* tmp(mat.data());
            for(int p=0;p<dim;++p)
            {
                newdata[p].real=tmp[p];
                newdata[p].imag=0.;
            }
        }
        /**
         * \brief transform a std::complex variable into std::vector<ctype>
         * \param newdata
         * \param arg
         *
         * Transform a complex_t type into writable std::vector<ctype>
         */
        void transform(std::vector<ctype>& newdata,complex_t arg)
        {
            ctype newarg;
            newarg.real=arg.real();
            newarg.imag=arg.imag();
            newdata.push_back(newarg);
        }
        /**
         * \brief transform an Eigen::Matrix<complex_t,Eigen::Dynamic,1> to std::vector<ctype>
         * \param newdata
         * \param cmat
         *
         * Transform a Eigen::Matrix<complex_t,Eigen::Dynamic,1> into writable std::vector<ctype>
         */
        void transform(std::vector<ctype>& newdata,Eigen::Matrix<complex_t, Eigen::Dynamic, 1> cmat)
        {
            newdata.resize(csize);
            for(int q=0;q<csize;++q)
            {
                newdata[q].real=cmat[q].real();
                newdata[q].imag=cmat[q].imag();
            }
        }
        /**
         * \brief transform real_t to std::vector<ctype>
         * \param newdata
         * \param arg
         *
         * Transforms a real_t into writeable std::vector<ctype>
         */
        void transform(std::vector<ctype>& newdata,real_t arg)
        {
            ctype newarg;
            newarg.real=arg;
            newarg.imag=0.;
            newdata.push_back(newarg);
        }
        /**
         * \brief store ScalarHaWp<D,Multiindex> packet
         *
         * Uses the the transformer PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>to(packet)
         * and the packet.parameters() function to call their corresponding transform function which then
         * can be written to dataset.
         */
        template<class MultiIndex>
        void store_packet(const waveblocks::wavepackets::ScalarHaWp<D,MultiIndex>& packetto)
        {
            const auto& params = packetto.parameters();
            if(wrlist["packet"])
            {
                if(tindex_packet%timestepsize_packet==0)
                {
                    select_file_writespace_packet();

                    std::vector<ctype> myc;
                    transform(myc,utilities::PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packetto));
                    coeffs->write(myc.data(),mytype_,celemspace,cspace);

                    std::vector<ctype> myQ;
                    transform(myQ,params.Q());
                    Qs->write(myQ.data(),mytype_,QPelemspace,Qspace);

                    std::vector<ctype> myP;
                    transform(myP,params.P());
                    Ps->write(myP.data(),mytype_,QPelemspace,Pspace);

                    std::vector<ctype> myq;
                    transform(myq,params.q());
                    qs->write(myq.data(),mytype_,qpelemspace,qspace);

                    std::vector<ctype> myp;
                    transform(myp,params.p());
                    ps->write(myp.data(),mytype_,qpelemspace,pspace);

                    std::vector<ctype> myS;
                    transform(myS,params.S());
                    Ss->write(myS.data(),mytype_,Selemspace,Sspace);

                    std::vector<ctype> myadQ;
                    transform(myadQ,params.sdQ());
                    adQs->write(myadQ.data(),mytype_,Selemspace,adQspace);

                    times_packet->write(&tindex_packet,PredType::NATIVE_INT,timelemspace,timespace_packet);

                    advance_packet();
                }
                tindex_packet+=1;
            }
            else
            {
                std::cout<<"ERROR: store_packet called with bool=false\n";
            }
        }
        /**
         * @brief store energies in chosen timestep
         * @param epot_
         * @param ekin_
         */
        void store_energies(double epot_,double ekin_)
        {
            if(wrlist["energy"])
            {
                if(tindex_ekin%timestepsize_ekin==0)
                {
                    select_file_writespace_ekin();
                    energys_ekin->write(&ekin_,PredType::NATIVE_DOUBLE,energyelemspace,energyspace_ekin);
                    times_ekin->write(&tindex_ekin,PredType::NATIVE_INT,timelemspace,timespace_ekin);
                    advance_ekin();
                }
                tindex_ekin+=1;
                if(tindex_epot%timestepsize_epot==0)
                {
                    select_file_writespace_epot();
                    energys_epot->write(&epot_,PredType::NATIVE_DOUBLE,energyelemspace,energyspace_epot);
                    times_epot->write(&tindex_epot,PredType::NATIVE_INT,timelemspace,timespace_epot);
                    advance_epot();
                }
                tindex_epot+=1;
            }
            else
            {
                std::cout<<"ERROR: store_energies called with bool=false\n";
            }
        }
        /**
         * \brief store norms in chosen timestep
         *
         * Extendable for norms user is interested in.
         */
        template<class MultiIndex>
        void store_norm(const waveblocks::wavepackets::ScalarHaWp<D,MultiIndex>& packet)
        {
            if(wrlist["norm"])
            {
                if(tindex_norm%timestepsize_norms==0)
                {
                    double mynorm = waveblocks::observables::norm<D,MultiIndex>(packet);
                    select_file_writespace_norm();
                    normss->write(&mynorm,PredType::NATIVE_DOUBLE,normelemspace,normspace);
                    times_norms->write(&tindex_norm,PredType::NATIVE_INT,timelemspace,timespace_norms);
                    advance_norm();
                }
                tindex_norm+=1;
            }
            else
            {
                std::cout<<"ERROR: store_norms called with bool=false\n";
            }
        }
        /**
         * @brief select file writespace for packet
         *
         * Called always before packet is written to evaluate space to write.
         */
        void select_file_writespace_packet(void)
        {
            constexpr int dim=D;
            if(wrlist["packet"])
            {
                //qp
                hsize_t count1[]={1,dim,1};
                hsize_t start1[3];
                start1[0]=index_packet-1;
                start1[1]=0;
                start1[2]=0;
                hsize_t stride1[]={1,1,1};
                hsize_t block1[]={1,1,1};
                //QP
                hsize_t count2[]={1,dim,dim};
                hsize_t start2[3];
                start2[0]=index_packet-1;
                start2[1]=0;
                start2[2]=0;
                hsize_t stride2[]={1,1,1};
                hsize_t block2[]={1,1,1};
                //S
                hsize_t count3[]={1,1,1};
                hsize_t start3[3];
                start3[0]=index_packet-1;
                start3[1]=0;
                start3[2]=0;
                hsize_t stride3[]={1,1,1};
                hsize_t block3[]={1,1,1};
                //qp
                qspace.selectHyperslab(H5S_SELECT_SET, count1, start1, stride1, block1);
                pspace.selectHyperslab(H5S_SELECT_SET, count1, start1, stride1, block1);
                //QP
                Qspace.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
                Pspace.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
                //S
                Sspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);
                adQspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);

                //coefficients
                hsize_t count4[2];
                count4[0]=1;
                count4[1]=csize;
                hsize_t start4[2];
                start4[0]=index_packet-1;
                start4[1]=0;
                hsize_t stride4[]={1,1};
                hsize_t block4[]={1,1};
                cspace.selectHyperslab(H5S_SELECT_SET, count4, start4, stride4, block4);

                hsize_t countt1[]={1};
                hsize_t startt1[1];
                startt1[0]=index_packet-1;
                hsize_t stridet1[]={1};
                hsize_t blockt1[]={1};
                timespace_packet.selectHyperslab(H5S_SELECT_SET, countt1, startt1, stridet1, blockt1);
            }
            else
            {
                std::cout<<"ERROR: select_file_writespace_packet called with bool=false\n";
            }
        }
        /**
         * @brief advance writing position for packet and extends set
         *
         * Wrapper for index_packet++, set_up_extension_packet(), extend_dataset_packet()
         * and update_filespace_packet()
         */
        void advance_packet(void)
        {
            index_packet+=1;
            setup_extension_packet();
            extend_dataset_packet();
            update_filespace_packet();
        }
        /**
         * @brief select file writespace for ekin
         *
         * Called always before ekin is written to evaluate space to write.
         */
        void select_file_writespace_ekin(void)
        {
            if(wrlist["energy"])
            {
                hsize_t count5[]={1,1};
                hsize_t start5[2];
                start5[0]=index_ekin-1;
                start5[1]=0;
                hsize_t stride5[]={1,1};
                hsize_t block5[]={1,1};
                energyspace_ekin.selectHyperslab(H5S_SELECT_SET, count5, start5, stride5, block5);

                hsize_t countt2[]={1};
                hsize_t startt2[1];
                startt2[0]=index_ekin-1;
                hsize_t stridet2[]={1};
                hsize_t blockt2[]={1};
                timespace_ekin.selectHyperslab(H5S_SELECT_SET, countt2, startt2, stridet2, blockt2);
            }
            else
            {
                std::cout<<"ERROR: select_file_writespace_ekin called with bool=false\n";
            }
        }
        /**
         * @brief advance writing position for ekin and extends set
         *
         * Wrapper for index_ekin++, set_up_extension_ekin(), extend_dataset_ekin()
         * and update_filespace_ekin()
         */
        void advance_ekin(void)
        {
            index_ekin+=1;
            setup_extension_ekin();
            extend_dataset_ekin();
            update_filespace_ekin();
        }
        /**
         * @brief select file writespace for epot
         *
         * Called always before epot is written to evaluate space to write.
         */
        void select_file_writespace_epot(void)
        {
            if(wrlist["energy"])
            {
                hsize_t count5[]={1,1};
                hsize_t start5[2];
                start5[0]=index_epot-1;
                start5[1]=0;
                hsize_t stride5[]={1,1};
                hsize_t block5[]={1,1};
                energyspace_epot.selectHyperslab(H5S_SELECT_SET, count5, start5, stride5, block5);

                hsize_t countt3[]={1};
                hsize_t startt3[1];
                startt3[0]=index_epot-1;
                hsize_t stridet3[]={1};
                hsize_t blockt3[]={1};
                timespace_epot.selectHyperslab(H5S_SELECT_SET, countt3, startt3, stridet3, blockt3);
            }
            else
            {
                std::cout<<"ERROR: select_file_writespace_epot called with bool=false\n";
            }
        }
        /**
         * @brief advance writing position for epot and extends set
         *
         * Wrapper for index_epot++, set_up_extension_epot(), extend_dataset_epot()
         * and update_filespace_epot()
         */
        void advance_epot(void)
        {
            index_epot+=1;
            setup_extension_epot();
            extend_dataset_epot();
            update_filespace_epot();
        }
        /**
         * @brief select file writespace for norms
         *
         * Called always before norms are written to evaluate space to write.
         */
        void select_file_writespace_norm(void)
        {
            if(wrlist["norm"])
            {
                hsize_t count6[]={1,1};
                hsize_t start6[2];
                start6[0]=index_norm-1;
                start6[1]=0;
                hsize_t stride6[]={1,1};
                hsize_t block6[]={1,1};
                normspace.selectHyperslab(H5S_SELECT_SET, count6, start6, stride6, block6);

                hsize_t countt4[]={1};
                hsize_t startt4[1];
                startt4[0]=index_norm-1;
                hsize_t stridet4[]={1};
                hsize_t blockt4[]={1};
                timespace_norms.selectHyperslab(H5S_SELECT_SET, countt4, startt4, stridet4, blockt4);
            }
            else
            {
                std::cout<<"ERROR: select_file_writespace_norms called with bool=false\n";
            }
        }
        /**
         * @brief advance writing position for norms and extends set
         *
         * Wrapper for index_norm++, set_up_extension_norms(), extend_dataset_norms()
         * and update_filespace_norms()
         */
        void advance_norm(void)
        {
            index_norm+=1;
            setup_extension_norm();
            extend_dataset_norm();
            update_filespace_norm();
        }
        /**
         * @brief reverse last dataset extension for packet
         *
         * Last writing process extended set although it is the last one.
         * This just reverses the last extension so the dataset is exactly as big
         * as the number of timesteps written.
         */
        void cleanup_packet(void)
        {
            index_packet-=1;
            setup_extension_packet();
            extend_dataset_packet();
        }
        /**
         * @brief reverse last dataset extension for ekin
         *
         * Last writing process extended set although it is the last one.
         * This just reverses the last extension so the dataset is exactly as big
         * as the number of timesteps written.
         */
        void cleanup_ekin(void)
        {
            index_ekin-=1;
            setup_extension_ekin();
            extend_dataset_ekin();
        }
        /**
         * @brief reverse last dataset extension for epot
         *
         * Last writing process extended set although it is the last one.
         * This just reverses the last extension so the dataset is exactly as big
         * as the number of timesteps written.
         */
        void cleanup_epot(void)
        {
            index_epot-=1;
            setup_extension_epot();
            extend_dataset_epot();
        }
        /**
         * @brief reverse last dataset extension for norms
         *
         * Last writing process extended set although it is the last one.
         * This just reverses the last extension so the dataset is exactly as big
         * as the number of timesteps written.
         */
        void cleanup_norm(void)
        {
            index_norm-=1;
            setup_extension_norm();
            extend_dataset_norm();
        }
        /**
         * @brief last resize for exact length
         *
         * Wrapper for all cleanup calls
         */
        void poststructuring(void)
        {
            if(wrlist["packet"])
            {
                cleanup_packet();
            }
            if(wrlist["norm"])
            {
                cleanup_norm();
            }
            if(wrlist["energy"])
            {
                cleanup_ekin();
                cleanup_epot();
            }
        }

    private:
        H5std_string filename_;//!<placeholder for filename
        CompType mytype_;//!<declaration of H5:CompType member used for HDF interface to write ctype*
        H5File file_; //!<H5File placeholder
        std::map<std::string,bool> wrlist={{"packet",1},{"energy",0},{"norm",0}};//!<maps string to bool for constructing und writing defined variables
        Attribute apacket;//!<attribute to save bool packet in datablock_0
        Attribute aenergy;//!<attribute to save bool energy in datablock_0
        Attribute anorm;//!<attribute to save bool norm in datablock_0
        Attribute adt;//!<attribute to save dt in datablock_0
        double dref=0.;//!<fillvalue for energys for allocation
        int iref=0;//!<fillvalue for timegrids for allocation
        H5std_string packet_group_string="/Pi";//!<String for H5Group to save packet to. Default:Pi
        H5std_string datablock_string="/datablock_0";//!<String for H5Group for datablock.default. datablock_0
        H5std_string coefficient_group_string="/coefficients";//!<String for H5Group of coefficients. Default:coefficients
        H5std_string wavepacket_group_string="/wavepacket";//!<String for H5Group for packet and coefficients. Default:wavepacket
        H5std_string energies_group="/energies";//!<name for group energies Default:energies
        H5std_string norms_group="/norm";//!<name for group norms Default:norm
        H5std_string observables_group="/observables";//!<name for group observables Default:observables
        DSetCreatPropList plist_qp;//!<PropList for packet.q() packet.p()
        DSetCreatPropList plist_QP;//!<PropList for packet.Q() packet.P()
        DSetCreatPropList plist_S;//!<PropList for packet.S()
        DSetCreatPropList plist_energy;//!<PropList for energies
        DSetCreatPropList plist_c;//!<PropList for coefficients
        DSetCreatPropList plist_time;//!<PropList for timegrids
        DSetCreatPropList plist_norms;//!<PropList for norms
        const int RANK1=1;//!<rank 1 identifier
        const int RANK2=2;//!<rank 2 identifier
        const int RANK3=3;//!<rank 3 identifier
        const hsize_t maxdims1[1]={H5S_UNLIMITED};//!<max dim identifier for rank1 for extension
        const hsize_t maxdims2[2]={H5S_UNLIMITED,H5S_UNLIMITED};//!<max dim identifier for rank2 for extension
        const hsize_t maxdims3[3]={H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};//!<max dim identifier for rank3 for extension
        hsize_t exqp[3];//!<extension for packet.q() packet.p()
        hsize_t exQP[3];//!<extension for packet.Q() packet.P()
        hsize_t exS[3];//!<extension for packet.S() and packet.sdQ()
        hsize_t exc[2];//!<extension for coefficients
        hsize_t exekin[2];//!<extension for ekin
        hsize_t exepot[2];//!<extension for epot
        hsize_t exnorms[2];//!<extension for norms
        hsize_t ex_timegrid_norms[1];//!<extension for timegrid for norms
        hsize_t ex_timegrid_epot[1];//!<extension for timegrid for epot
        hsize_t ex_timegrid_ekin[1];//!<extension for timegrid for ekin
        hsize_t ex_timegrid_packet[1];//!<extension for timegrid for packet
        hsize_t qpelem[3]; //!<size of q,p element written from program to file needed by HDF interface
        DataSpace qpelemspace;//!<space of q,p element written from program to file needed by HDF interface
        hsize_t QPelem[3];//!<size of Q,P element written from program to file needed by HDF interface
        DataSpace QPelemspace;//!<space of Q,P element written from program to file needed by HDF interface
        hsize_t Selem[3];//!<size of S element written from program to file needed by HDF interface
        DataSpace Selemspace;//!<space of S element written from program to file needed by HDF interface
        hsize_t energyelem[2];//!<size of energy element written from program to file needed by HDF interface
        DataSpace energyelemspace;//!<space of energy element written from program to file needed by HDF interface
        hsize_t timeelem[1];//!<size of timegrid element written from program to file needed by HDF interface
        DataSpace timelemspace;//!<space of timegrid element written from program to file needed by HDF interface
        hsize_t celem[2];//!<size of coefficient element written from program to file needed by HDF interface
        DataSpace celemspace;//!<space of coefficient element written from program to file needed by HDF interface
        hsize_t normelem[2];//!<size of coefficient element written from program to file needed by HDF interface
        DataSpace normelemspace;//!<space of norm element written from program to file needed by HDF interface
        int index_packet=1;//!<index used for storing packet
        int index_norm=1;//!<index used for storing norm
        int index_ekin=1;//!<index used for storing ekin
        int index_epot=1;//!<index used for storing epot
        int timestepsize_norms=1;//!<timestepsize for norm timegrid
        int timestepsize_ekin=1;//!<timestepsize for ekin timegrid
        int timestepsize_epot=1;//!<timestepsize for epot timegrid
        int timestepsize_packet=1;//!<timestepsize for packet timegrid
        int tindex_ekin=0;//!<timeindex for modulo writing ekin
        int tindex_epot=0;//!<timeindex for modulo writing epot
        int tindex_norm=0;//!<timeindex for modulo writing norms
        int tindex_packet=0;//!<timeindex for modulo writing packet
        int csize;//!<runtime size coefficients
        std::shared_ptr<Group> gblock;//!<group for datablock
        std::shared_ptr<Group> gpacket;//!<group for packet
        std::shared_ptr<Group> gPi;//!<group for matrices in packet
        std::shared_ptr<Group> gcoefficient;//!<group for coefficients in packet
        std::shared_ptr<Group> genergy;//!<group for energies
        std::shared_ptr<Group> gnorms;//!<group for norms
        std::shared_ptr<Group> gobservables;//!<group for observables
        DataSpace qspace;//!<space for packet.q() in file
        std::shared_ptr<DataSet> qs;//!<dataset for packet.q() in file
        DataSpace pspace;//!<space for packet.p() in file
        std::shared_ptr<DataSet> ps;//!<dataset for packet.p() in file
        DataSpace Qspace;//!<space for packet.Q() in file
        std::shared_ptr<DataSet> Qs;//!<dataset for packet.Q() in file
        DataSpace Pspace;//!<space for packet.P() in file
        std::shared_ptr<DataSet> Ps;//!<dataset for packet.P() in file
        DataSpace Sspace;//!<space for packet.S() in file
        std::shared_ptr<DataSet> Ss;//!<dataset for packet.S() in file
        DataSpace adQspace;//!<space for packet.sdQ() in file
        std::shared_ptr<DataSet> adQs;//!<dataset for packet.sdQ() in file
        DataSpace cspace;//!<space for coefficients in file
        std::shared_ptr<DataSet> coeffs;//!<dataset for coefficients in file
        DataSpace timespace_packet;//!<space for timegrid for packet in file
        std::shared_ptr<DataSet> times_packet;//!<dataset for timegrid for packet in file
        DataSpace timespace_ekin;//!<space for timegrid ekin
        std::shared_ptr<DataSet> times_ekin;//!<dataset for timegrid for ekin in file
        DataSpace timespace_epot;//!<space for timegrid epot
        std::shared_ptr<DataSet> times_epot;//!<dataset for timegrid for epot in file
        DataSpace timespace_norms;//!<space for timegrid in file
        std::shared_ptr<DataSet> times_norms;//!<dataset for timegrid in file
        DataSpace energyspace_ekin;//!<space for ekin in file
        std::shared_ptr<DataSet> energys_ekin;//!<dataset for ekin in file
        DataSpace energyspace_epot;//!<space for epot in file
        std::shared_ptr<DataSet> energys_epot;//!<dataset for epot in file
        DataSpace normspace;//!<space for norms
        std::shared_ptr<DataSet> normss;//!<dataset for norms in file
    };
    }
}
