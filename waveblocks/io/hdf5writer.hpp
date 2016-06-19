#pragma once

//need H5 Cpp header
#include "H5Cpp.h"
//using strings
#include <string>
#include <map>
#include <cstdio>

#include <Eigen/Core>

#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"

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
     * @brief get size of coefficients in compile time
     * @return
     */
    constexpr int get_size_coefficients(void)
    {
        //TODO
        return 16;
    }

    /**
     * \brief Our HDF5 writer class
     *
     * This class is templated on the Dimension D which is also used to also define the dimensionality
     * of the HaWpBasisVector
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
         * \brief prestructure after knowing bool values
         */
        void prestructuring(void)
        {
            //set group structure
            set_group_structure();

            //set up chunk dimension
            set_chunk_dim();

            //set up space for file
            set_file_dataspace();

            //allocate datasets needs to be after select writespace
            allocate_datasets();

            //set up elem space
            set_elem_space();

            //select elemtary space
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
        void set_write_energy(bool flag)
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
         * @brief set timestep for writing packet
         * @param a
         */
        void set_timestep_packet(int a)
        {
            timestepsize_packet=a;
        }
        /**
         * @brief set timestep for writing norms
         * @param a
         */
        void set_timestep_norms(int a)
        {
            timestepsize_norms=a;
        }
        /**
         * @brief set timestep for writing epot
         * @param a
         */
        void set_timestep_epot(int a)
        {
            timestepsize_epot=a;
        }
        /**
         * @brief set timestep for writing ekin
         * @param a
         */
        void set_timestep_ekin(int a)
        {
            timestepsize_ekin=a;
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
            plist_time.setFillValue(PredType::NATIVE_DOUBLE,&dref);

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

                constexpr int cdim=get_size_coefficients();
                hsize_t chunk_dims5[]={1,cdim};
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
                //TODO
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
                constexpr int cdim=get_size_coefficients();
                celem[1]=cdim;
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
                //TODO
                normelem[0]=1;
                normelem[1]=1;
                DataSpace t6(RANK2,normelem);
                normelemspace=t6;
            }
        }
        /**
         * @brief set up group structure in file
         */
        void set_group_structure(void)
        {
            H5std_string a1=datablock_string;
            gblock = std::make_shared<Group>(file_.createGroup(a1));
            if(wrlist["packet"])
            {
                H5std_string a2=(datablock_string+wavepacket_group_string);
                gpacket = std::make_shared<Group>(file_.createGroup(a2));
                H5std_string a3=(datablock_string+wavepacket_group_string+packet_group_string);
                gPi = std::make_shared<Group>(file_.createGroup(a3));
                H5std_string a4=(datablock_string+wavepacket_group_string+coefficient_group_string);
                gcoefficient = std::make_shared<Group>(file_.createGroup(a4));
            }
            if(wrlist["energy"])
            {
                H5std_string a5 = (datablock_string+energies_group);
                genergy = std::make_shared<Group>(file_.createGroup(a5));
            }
            if(wrlist["norm"])
            {
                H5std_string a6 = (datablock_string+norms_group);
                gnorms = std::make_shared<Group>(file_.createGroup(a6));
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

                constexpr int cdim = get_size_coefficients();
                //coefficients
                hsize_t count4[]={1,cdim};
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
                //TODO
            }
        }

        /**
         * @brief set dataspace used in file
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

                constexpr int cdim = get_size_coefficients();
                //set up coefficients
                hsize_t dim5[]={1,cdim};
                DataSpace d5(RANK2,dim5,maxdims2);
                cspace=d5;
            }
            if(wrlist["energy"])
            {
                hsize_t dim4[]={1,1};
                DataSpace d4(RANK2,dim4,maxdims2);
                energyspace_ekin=d4;
                energyspace_epot=d4;

            }
            if(wrlist["norm"])
            {
                //TODO
            }
        }
        /**
         * @brief allocate space for datasets
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
                times_packet=std::make_shared<DataSet>(file_.createDataSet(pack+"/timegrid",PredType::NATIVE_DOUBLE,timespace_packet,plist_time));

                H5std_string cff=datablock_string+wavepacket_group_string+coefficient_group_string;
                coeffs=std::make_shared<DataSet>(file_.createDataSet(cff+"/c_0",mytype_,cspace,plist_c));
            }
            if(wrlist["energy"])
            {
                H5std_string ename=datablock_string+energies_group;
                H5std_string enameepot=ename+"/epot";
                energys_epot=std::make_shared<DataSet>(file_.createDataSet(enameepot,PredType::NATIVE_DOUBLE,energyspace_epot,plist_energy));
                H5std_string enameekin=ename+"/ekin";
                energys_ekin=std::make_shared<DataSet>(file_.createDataSet(enameekin,PredType::NATIVE_DOUBLE,energyspace_ekin,plist_energy));
                H5std_string engtime1=ename+"/timegrid_epot";
                times_epot=std::make_shared<DataSet>(file_.createDataSet(engtime1,PredType::NATIVE_DOUBLE,timespace_epot,plist_time));
                H5std_string engtime2=ename+"/timegrid_ekin";
                times_ekin=std::make_shared<DataSet>(file_.createDataSet(engtime2,PredType::NATIVE_DOUBLE,timespace_ekin,plist_time));
            }
            if(wrlist["norm"])
            {
                H5std_string tmm=datablock_string+norms_group+"/norms";
                normss=std::make_shared<DataSet>(file_.createDataSet(tmm,PredType::NATIVE_DOUBLE,normspace,plist_norms));
                H5std_string tmt=datablock_string+norms_group+"/timegrid";
                times_norms=std::make_shared<DataSet>(file_.createDataSet(tmt,PredType::NATIVE_DOUBLE,timespace_norms,plist_time));
            }
        }
        /**
         * @brief setup extension for norms
         */
        void setup_extension_norms(void)
        {
            //TODO
            exnorms[0]=index_norm;
            exnorms[1]=1;
        }
        /**
         * @brief setup extension for epot
         */
        void setup_extension_epot(void)
        {
            exepot[0]=index_epot;
        }
        /**
         * @brief setup extension for ekin
         */
        void setup_extension_ekin(void)
        {
            exekin[0]=index_ekin;
        }
        /**
         * @brief setup extension for packet
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

                constexpr int cdim =get_size_coefficients();
                exc[0]=index_packet;
                exc[1]=cdim;
            }
            else
            {
                std::cout<<"ERROR: setup_extension_packet called with bool=false\n";
            }
        }
        /**
         * @brief extend dataset for next timestep for packet
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
            }
            else
            {
                std::cout<<"ERROR: extend_dataset_packet called with bool=false\n";
            }
        }
        /**
         * @brief extend dataset for next timestep for norms
         */
        void extend_dataset_norms(void)
        {
            if(wrlist["norm"])
            {
                normss->extend(exnorms);
            }
            else
            {
                std::cout<<"ERROR: extend_dataset_norms called with bool=false\n";
            }
        }
        /**
         * @brief extend dataset for next timestep for epot
         */
        void extend_dataset_epot(void)
        {
            if(wrlist["energy"])
            {
                energys_epot->extend(exepot);
            }
            else
            {
                std::cout<<"ERROR: extend_dataset_epot called with bool=false\n";
            }
        }
        /**
         * @brief extend dataset for next timestep for ekin
         */
        void extend_dataset_ekin(void)
        {
            if(wrlist["energy"])
            {
                energys_ekin->extend(exekin);
            }
            else
            {
                std::cout<<"ERROR: extend_dataset_ekin called with bool=false\n";
            }
        }
        /**
         * @brief update ekin filespace
         *
         * Needs to be done after extending a dataset
         */
        void update_filespace_ekin(void)
        {
            if(wrlist["energy"])
            {
                energyspace_ekin=energys_ekin->getSpace();
            }
            else
            {
                std::cout<<"ERROR: update_filespace_ekin called with bool=false\n";
            }
        }
        /**
         * @brief update epot filespace
         *
         * Needs to be done after extending a dataset
         */
        void update_filespace_epot(void)
        {
            if(wrlist["energy"])
            {
                energyspace_epot=energys_epot->getSpace();
            }
            else
            {
                std::cout<<"ERROR: update_filespace_epot called with bool=false\n";
            }
        }
        /**
         * @brief update norms filespace
         *
         * Needs to be done after extending a dataset
         */
        void update_filespace_norms(void)
        {
            if(wrlist["norm"])
            {
                //TODO
            }
            else
            {
                std::cout<<"ERROR: update_filespace_norms called with bool=false\n";
            }
        }
        /**
         * @brief update packet filespace
         *
         * Needs to be done after extending a dataset
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
            }
            else
            {
                std::cout<<"ERROR: update_filespace_packet called with bool=false\n";
            }
        }
        /**
         * \brief transform Eigen::Matrix<complex_t,D,1> into ctype*
         * \param mat
         * \return ctype*
         *
         * TODO in transform fix memory leaks
         */
        ctype* transform(Eigen::Matrix<complex_t,D,D> mat)
        {
            constexpr int dim=D;
            ctype* newdat= new ctype[dim*dim];
            std::complex<double>* tmp(mat.data());
            for(int p=0;p<dim*dim;++p)
            {
                newdat[p].real=tmp[p].real();
                newdat[p].imag=tmp[p].imag();
            }
            return newdat;
        }
        /**
         * \brief transform Eigen::Matrix<real_t,D,1> into ctype*
         * \param mat
         * \return ctype*
         */
        ctype* transform(Eigen::Matrix<real_t,D,1> mat)
        {
            constexpr int dim=D;
            ctype* newdat = new ctype[dim];
            double* tmp(mat.data());
            for(int p=0;p<dim;++p)
            {
                newdat[p].real=tmp[p];
                newdat[p].imag=0.;
            }
            return newdat;
        }
        /**
         * \brief transform a std::complex variable into ctype*
         * \param arg
         * \return ctype*
         */
        ctype* transform(complex_t arg)
        {
            ctype* newarg = new ctype;
            newarg->real=arg.real();
            newarg->imag=arg.imag();
            return newarg;
        }
        /**
         * \brief transform an Eigen::Matrix<complex_t,Eigen::Dynamic,1> to ctype*
         * \param cmat
         * \return ctype*
         *
         * Be careful the dimension is hardcoded to a 16x1 matrix from coefficients
         * and will be written back as 1x16 matrix
         */
        ctype* transform(Eigen::Matrix<complex_t, Eigen::Dynamic, 1> cmat)
        {
            //int rowdim=cmat.rows();
            //int coldim=cmat.cols();
            //hardcoded because Eigen::dynamic is not const
            //need rowdim = 1 coldim = 16
            ctype* newdat=new ctype[16];

            for(int q=0;q<16;++q)
            {
                newdat[q].real=cmat[q].real();
                newdat[q].imag=cmat[q].imag();
            }
            return newdat;
        }
        /**
         * \brief transform real_t to ctype
         * \param arg
         * \return
         */
        ctype * transform(real_t arg)
        {
            ctype* newarg = new ctype;
            newarg->real=arg;
            newarg->imag=0.;
            return newarg;
        }
        /**
         * \brief store ScalarHaWp<D,Multiindex> packet
         *
         * Uses the the transformer PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>to(packet)
         * and the packet.parameters() function
         */
        template<class MultiIndex>
        void store_packet(const double& time_, const waveblocks::wavepackets::ScalarHaWp<D,MultiIndex>& packetto)
        {
            const auto& params = packetto.parameters();
            if(wrlist["packet"])
            {
                if((index_packet-1)%timestepsize_packet==0)
                {
                    select_file_writespace_packet();

                    ctype* myc=transform(utilities::PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packetto));
                    coeffs->write(myc,mytype_,celemspace,cspace);

                    ctype* myQw=transform(params.Q());
                    Qs->write(myQw,mytype_,QPelemspace,Qspace);
                    ctype* myP=transform(params.P());
                    Ps->write(myP,mytype_,QPelemspace,Pspace);

                    ctype* myq=transform(params.q());
                    qs->write(myq,mytype_,qpelemspace,qspace);
                    ctype* myp=transform(params.p());
                    ps->write(myp,mytype_,qpelemspace,pspace);

                    ctype* myS=transform(params.S());
                    Ss->write(myS,mytype_,Selemspace,Sspace);

                    ctype* myadQ=transform(params.sdQ());
                    adQs->write(myadQ,mytype_,Selemspace,adQspace);

                    advance_packet();
                }
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
                if((index_ekin-1)%timestepsize_ekin==0)
                {
                    select_file_writespace_ekin();
                    energys_ekin->write(&ekin_,PredType::NATIVE_DOUBLE,energyelemspace,energyspace_ekin);
                    advance_ekin();
                }
                if((index_epot-1)%timestepsize_epot==0)
                {
                    select_file_writespace_epot();
                    energys_epot->write(&epot_,PredType::NATIVE_DOUBLE,energyelemspace,energyspace_epot);
                    advance_epot();
                }
            }
            else
            {
                std::cout<<"ERROR: store_energies called with bool=false\n";
            }
        }
        /**
         * \brief store norms in chosen timestep
         */
        template<class MultiIndex>
        void store_norms(const waveblocks::wavepackets::ScalarHaWp<D,MultiIndex>& packet)
        {
            if(wrlist["norm"])
            {
                if((index_norm-1)%timestepsize_norms==0)
                {
                    double norms=0.5;
                    select_file_writespace_norms();
                    normss->write(&norms,PredType::NATIVE_DOUBLE,normelemspace,normspace);
                    advance_norms();
                }
            }
            else
            {
                std::cout<<"ERROR: store_norms called with bool=false\n";
            }
        }
        /**
         * @brief select file writespace for packet
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

                constexpr int cdim = get_size_coefficients();
                //coefficients
                hsize_t count4[]={1,cdim};
                hsize_t start4[2];
                start4[0]=index_packet-1;
                start4[1]=0;
                hsize_t stride4[]={1,1};
                hsize_t block4[]={1,1};
                cspace.selectHyperslab(H5S_SELECT_SET, count4, start4, stride4, block4);
            }
            else
            {
                std::cout<<"ERROR: select_file_writespace_packet called with bool=false\n";
            }
        }
        /**
         * @brief advance writing position for packet and extends set
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
            }
            else
            {
                std::cout<<"ERROR: select_file_writespace_ekin called with bool=false\n";
            }
        }
        /**
         * @brief advance writing position for ekin and extends set
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
            }
            else
            {
                std::cout<<"ERROR: select_file_writespace_epot called with bool=false\n";
            }
        }
        /**
         * @brief advance writing position for epot and extends set
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
         */
        void select_file_writespace_norms(void)
        {
            if(wrlist["norm"])
            {
                //TODO
            }
            else
            {
                std::cout<<"ERROR: select_file_writespace_norms called with bool=false\n";
            }
        }
        /**
         * @brief advance writing position for norms and extends set
         */
        void advance_norms(void)
        {
            index_norm+=1;
            setup_extension_norms();
            extend_dataset_norms();
            update_filespace_norms();
        }
        /**
         * @brief reverse last dataset extension for packet
         */
        void cleanup_packet(void)
        {
            index_packet-=1;
            setup_extension_packet();
            extend_dataset_packet();
        }

    private:
        H5std_string filename_;///<identifier for filename
        CompType mytype_;///Declaration of H5:CompType member for HDF interface
        H5File file_; ///H5File member
        std::map<std::string,bool> wrlist={{"packet",1},{"energy",0},{"norm",0}};///map string->bool for constructing und writing defined variables

        double dref=0.;///fillvalue for energys for allocation
        H5std_string packet_group_string="/Pi";///String for H5Group to save packet to. Default:Pi
        H5std_string datablock_string="/datablock_0";///String for H5Group for datablock.default. datablock_0
        H5std_string coefficient_group_string="/coefficients";///String for H5Group of coefficients. Default:coefficients
        H5std_string wavepacket_group_string="/wavepacket";///String for H5Group for packet and coefficients. Default:wavepacket
        H5std_string energies_group="/energies";///name for group energies Default:energies
        H5std_string norms_group="/norms";///name for group norms Default:norms

        DSetCreatPropList plist_qp;///PropList for packet.q() packet.p()
        DSetCreatPropList plist_QP;///PropList for packet.Q() packet.P()
        DSetCreatPropList plist_S;///PropList for packet.S()
        DSetCreatPropList plist_energy;///PropList for energies
        DSetCreatPropList plist_c;///PropList for coefficients
        DSetCreatPropList plist_time;///PropList for timegrids
        DSetCreatPropList plist_norms;///PropList for norms


        const int RANK1=1;///rank 1 identifier
        const int RANK2=2;///rank 2 identifier
        const int RANK3=3;///rank 3 identifier

        const hsize_t maxdims1[1]={H5S_UNLIMITED};///max dim identifier for rank1 for extension
        const hsize_t maxdims2[2]={H5S_UNLIMITED,H5S_UNLIMITED};///max dim identifier for rank2 for extension
        const hsize_t maxdims3[3]={H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};///max dim identifier for rank3 for extension

        hsize_t exqp[3];///extension for packet.q() packet.p()
        hsize_t exQP[3];///extension for packet.Q() packet.P()
        hsize_t exS[3];///extension for packet.S() and packet.sdQ()
        hsize_t exc[2];///extension for coefficients
        hsize_t exekin[2];///extension for ekin
        hsize_t exepot[2];///extension for epot
        hsize_t exnorms[2];///extension for norms
        hsize_t ex_timegrid_norms[1];///extension for timegrid for norms
        hsize_t ex_timegrid_epot[1];///extension for timegrid for epot
        hsize_t ex_timegrid_ekin[1];///extension for timegrid for ekin
        hsize_t ex_timegrid__packet[1];///extension for timegrid for packet

        hsize_t qpelem[3]; ///size of q,p element written from program to file needed by HDF interface
        DataSpace qpelemspace;///space of q,p element written from program to file needed by HDF interface
        hsize_t QPelem[3];///size of Q,P element written from program to file needed by HDF interface
        DataSpace QPelemspace;///space of Q,P element written from program to file needed by HDF interface
        hsize_t Selem[3];///size of S element written from program to file needed by HDF interface
        DataSpace Selemspace;///space of S element written from program to file needed by HDF interface
        hsize_t energyelem[2];///size of energy element written from program to file needed by HDF interface
        DataSpace energyelemspace;///space of energy element written from program to file needed by HDF interface
        hsize_t timeelem[1];///size of timegrid element written from program to file needed by HDF interface
        DataSpace timelemspace;///space of timegrid element written from program to file needed by HDF interface
        hsize_t celem[2];///size of coefficient element written from program to file needed by HDF interface
        DataSpace celemspace;///space of coefficient element written from program to file needed by HDF interface
        hsize_t normelem[2];///size of coefficient element written from program to file needed by HDF interface
        DataSpace normelemspace;///space of norm element written from program to file needed by HDF interface

        int index_packet=1;///index for storing packet
        int index_norm=1;///index for storing norm
        int index_ekin=1;///index for storing ekin
        int index_epot=1;///index for storing epot

        int timestepsize_norms=1;///timestepsize for norm timegrid
        int timestepsize_ekin=1;///timestepsize for ekin timegrid
        int timestepsize_epot=1;///timestepsize for epot timegrid
        int timestepsize_packet=1;///timestepsize for packet timegrid

        std::shared_ptr<Group> gblock;///group for datablock
        std::shared_ptr<Group> gpacket;///group for packet
        std::shared_ptr<Group> gPi;///group for matrices in packet
        std::shared_ptr<Group> gcoefficient;///group for coefficients in packet
        std::shared_ptr<Group> genergy;///group for energies
        std::shared_ptr<Group> gnorms;//group for norms

        DataSpace qspace;///space for packet.q() in file
        std::shared_ptr<DataSet> qs;///dataset for packet.q() in file
        DataSpace pspace;///space for packet.p() in file
        std::shared_ptr<DataSet> ps;///dataset for packet.p() in file
        DataSpace Qspace;///space for packet.Q() in file
        std::shared_ptr<DataSet> Qs;///dataset for packet.Q() in file
        DataSpace Pspace;///space for packet.P() in file
        std::shared_ptr<DataSet> Ps;///dataset for packet.P() in file
        DataSpace Sspace;///space for packet.S() in file
        std::shared_ptr<DataSet> Ss;///dataset for packet.S() in file
        DataSpace adQspace;///space for packet.sdQ() in file
        std::shared_ptr<DataSet> adQs;///dataset for packet.sdQ() in file
        DataSpace cspace;///space for coefficients in file
        std::shared_ptr<DataSet> coeffs;///dataset for coefficients in file
        DataSpace timespace_packet;///space for timegrid for packet in file
        std::shared_ptr<DataSet> times_packet;///dataset for timegrid for packet in file
        DataSpace timespace_ekin;///space for timegrid ekin
        std::shared_ptr<DataSet> times_ekin;///dataset for timegrid for ekin in file
        DataSpace timespace_epot;///space for timegrid epot
        std::shared_ptr<DataSet> times_epot;///dataset for timegrid for epot in file
        DataSpace timespace_norms;///space for timegrid in file
        std::shared_ptr<DataSet> times_norms;///dataset for timegrid in file
        DataSpace energyspace_ekin;///space for ekin in file
        std::shared_ptr<DataSet> energys_ekin;///dataset for ekin in file
        DataSpace energyspace_epot;///space for epot in file
        std::shared_ptr<DataSet> energys_epot;///dataset for epot in file
        DataSpace normspace;///space for norms
        std::shared_ptr<DataSet> normss;///dataset for norms in file
    };
    }
}
