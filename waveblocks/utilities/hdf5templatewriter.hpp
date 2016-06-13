#pragma once

//need H5 Cpp header
#include "H5Cpp.h"
//using strings
#include <string>
#include <map>

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
         * @brief sets bool value for writing the packet
         * @param flag
         */
        void set_write_packet(bool flag)
        {
            writeflags[0]=flag;
            wrlist["packet"]=flag;
        }
        /**
         * @brief sets bool value for writing the energies
         * @param flag
         */
        void set_write_energy(bool flag)
        {
            writeflags[2]=flag;
            wrlist["energy"]=flag;
        }
        /**
         * @brief sets bool value for writing coefficients
         * @param flag
         */
        void set_write_coefficients(bool flag)
        {
            writeflags[1]=flag;
            wrlist["coefficients"]=flag;
        }
        /**
         * @brief sets bool value for writing timegrid
         * @param flag
         */
        void set_write_timegrid(bool flag)
        {
            writeflags[3]=flag;
            wrlist["timegrid"]=flag;
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
         * \brief Set up chunk dimension for the written variables
         *
         * The HDF Interface needs that the PropList exits with chunked dimension otherwise
         * it is not possible to extend this particular DataSet
         */
        void set_chunk_dim(void)
        {
            constexpr int dim = D; //alias for D
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
            }
            if(wrlist["energy"])
            {
                hsize_t chunk_dims4[]={1,3};
                plist_energy.setChunk(RANK2,chunk_dims4);
                plist_energy.setFillValue(PredType::NATIVE_DOUBLE,&dref);

            }
            if(wrlist["coefficients"])
            {
                hsize_t chunk_dims5[]={1,16};
                plist_c.setChunk(RANK2,chunk_dims5);
                plist_c.setFillValue(mytype_,&instanceof);

            }
            if(wrlist["timegrid"])
            {
                hsize_t chunk_dims6[]={1};
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
            }
            if(wrlist["energy"])
            {
                energyelem[0]=1;
                energyelem[1]=3;
                DataSpace t4(RANK2,energyelem);
                energyelemspace=t4;
            }
            if(wrlist["coefficients"])
            {
                celem[0]=1;
                celem[1]=16;
                DataSpace t5(RANK2,celem);
                celemspace=t5;
            }
            if(wrlist["timegrid"])
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
            if(wrlist["packet"] && wrlist["coefficients"])
            {
                H5std_string a2=(datablock_string+wavepacket_group_string);
                gpacket = std::make_shared<Group>(file_.createGroup(a2));
                H5std_string a3=(datablock_string+wavepacket_group_string+packet_group_string);
                gPi = std::make_shared<Group>(file_.createGroup(a3));
                H5std_string a4=(datablock_string+wavepacket_group_string+coefficient_group_string);
                gcoefficient = std::make_shared<Group>(file_.createGroup(a4));
            }
            else if(wrlist["packet"])
            {
                H5std_string a2=(datablock_string+wavepacket_group_string);
                gpacket = std::make_shared<Group>(file_.createGroup(a2));
                H5std_string a3=(datablock_string+wavepacket_group_string+packet_group_string);
                gPi = std::make_shared<Group>(file_.createGroup(a3));
            }
            else if(wrlist["coefficients"])
            {
                H5std_string a2=(datablock_string+wavepacket_group_string);
                gpacket = std::make_shared<Group>(file_.createGroup(a2));
                H5std_string a4=(datablock_string+wavepacket_group_string+coefficient_group_string);
                gcoefficient = std::make_shared<Group>(file_.createGroup(a4));
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
                //Selem
                Selemspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);
            }
            if(wrlist["coefficients"])
            {
                //coefficients
                hsize_t count4[]={1,16};
                hsize_t start4[]={0,0};
                hsize_t stride4[]={1,1};
                hsize_t block4[]={1,1};
                celemspace.selectHyperslab(H5S_SELECT_SET, count4, start4, stride4, block4);
            }
            if(wrlist["energy"])
            {
                hsize_t count5[]={1,3};
                hsize_t start5[]={0,0};
                hsize_t stride5[]={1,1};
                hsize_t block5[]={1,1};
                energyelemspace.selectHyperslab(H5S_SELECT_SET, count5, start5, stride5, block5);
            }
            if(wrlist["timegrid"])
            {
                hsize_t count6[]={1};
                hsize_t start6[]={0};
                hsize_t stride6[]={1};
                hsize_t block6[]={1};
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
                //set up S
                hsize_t dim3[]={1,1,1};
                DataSpace d3(RANK3,dim3,maxdims3);
                Sspace=d3;
            }
            if(wrlist["energy"])
            {
                hsize_t dim4[]={10,3};
                DataSpace d4(RANK2,dim4,maxdims2);
                energyspace=d4;

            }
            if(wrlist["coefficients"])
            {
                //set up coefficients
                hsize_t dim5[]={1,16};
                DataSpace d5(RANK2,dim5,maxdims2);
                cspace=d5;
            }
            if(wrlist["timegrid"])
            {
                hsize_t dim6[]={10};
                DataSpace d6(RANK1,dim6,maxdims1);
                timespace=d6;
            }
        }
        /**
         * @brief set the writespace(block) which is written to in file
         *
         * Is dependent on template int D and current_index
         */
        void select_file_writespace(void)
        {
            int tr=current_index-1;
            constexpr int dim=D;
            if(wrlist["packet"])
            {
                //qp
                hsize_t count1[]={1,dim,1};
                hsize_t start1[3];
                start1[0]=tr;
                start1[1]=0;
                start1[2]=0;
                hsize_t stride1[]={1,1,1};
                hsize_t block1[]={1,1,1};
                //QP
                hsize_t count2[]={1,dim,dim};
                hsize_t start2[3];
                start2[0]=tr;
                start2[1]=0;
                start2[2]=0;
                hsize_t stride2[]={1,1,1};
                hsize_t block2[]={1,1,1};
                //S
                hsize_t count3[]={1,1,1};
                hsize_t start3[3];
                start3[0]=tr;
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
            }
            if(wrlist["coefficients"])
            {
                //coefficients
                hsize_t count4[]={1,16};
                hsize_t start4[2];
                start4[0]=tr;
                start4[1]=0;
                hsize_t stride4[]={1,1};
                hsize_t block4[]={1,1};
                //coefficients
                cspace.selectHyperslab(H5S_SELECT_SET, count4, start4, stride4, block4);
            }
            if(wrlist["energy"])
            {
                hsize_t count5[]={1,3};
                hsize_t start5[2];
                start5[0]=tr;
                start5[1]=0;
                hsize_t stride5[]={1,1};
                hsize_t block5[]={1,1};
                energyspace.selectHyperslab(H5S_SELECT_SET, count5, start5, stride5, block5);

            }
            if(wrlist["timegrid"])
            {
                hsize_t count6[]={1};
                hsize_t start6[1];
                start6[0]=tr;
                hsize_t stride6[]={1};
                hsize_t block6[]={1};
                timespace.selectHyperslab(H5S_SELECT_SET, count6, start6, stride6, block6);
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
                qs=std::make_shared<DataSet>(file_.createDataSet(pack+q,mytype_,qspace,plist_qp));
                ps=std::make_shared<DataSet>(file_.createDataSet(pack+p,mytype_,pspace,plist_qp));
                Qs=std::make_shared<DataSet>(file_.createDataSet(pack+Q,mytype_,Qspace,plist_QP));
                Ps=std::make_shared<DataSet>(file_.createDataSet(pack+P,mytype_,Pspace,plist_QP));
                Ss=std::make_shared<DataSet>(file_.createDataSet(pack+S,mytype_,Sspace,plist_S));
            }
            if(wrlist["coefficients"])
            {
                H5std_string cff=datablock_string+wavepacket_group_string+coefficient_group_string;
                coeffs=std::make_shared<DataSet>(file_.createDataSet(cff+c,mytype_,cspace,plist_c));
            }
            if(wrlist["energy"])
            {
                H5std_string eng=datablock_string+energies;
                energys=std::make_shared<DataSet>(file_.createDataSet(eng,PredType::NATIVE_DOUBLE,energyspace,plist_energy));
            }
            if(wrlist["timegrid"])
            {
                H5std_string tmm=datablock_string+wavepacket_group_string+time;
                times=std::make_shared<DataSet>(file_.createDataSet(tmm,PredType::NATIVE_DOUBLE,timespace,plist_time));
            }
        }
        /**
         * @brief setup extensions for file for a timestep
         */
        void setup_extensions(void)
        {
            constexpr int dim=D;
            if(wrlist["packet"])
            {
                exqp[0]=current_index;
                exqp[1]=dim;
                exqp[2]=1;
                exQP[0]=current_index;
                exQP[1]=dim;
                exQP[2]=dim;
                exS[0]=current_index;
                exS[1]=1;
                exS[2]=1;
            }
            if(wrlist["coefficients"])
            {
                exc[0]=current_index;
                exc[1]=16;
            }
            if(wrlist["energy"])
            {
                exenergy[0]=current_index;
                exenergy[1]=3;
            }
            if(wrlist["timegrid"])
            {
                extime[0]=current_index;
            }
        }
        /**
         * @brief extend datasets for next timestep
         */
        void extend_datasets(void)
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
            }
            if(wrlist["coefficients"])
            {
                //coefficients
                coeffs->extend(exc);
            }
            if(wrlist["energy"])
            {
                energys->extend(exenergy);
            }
            if(wrlist["timegrid"])
            {
                times->extend(extime);
            }
        }
        /**
         * @brief update corresponding filespace
         *
         * Needs to be done after extending a dataset
         */
        void update_filespace(void)
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
            }
            if(wrlist["coefficients"])
            {
                cspace=coeffs->getSpace();
            }
            if(wrlist["energy"])
            {
                energyspace=energys->getSpace();
            }
            if(wrlist["timegrid"])
            {
                timespace=times->getSpace();
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
         * \brief store ScalarHaWp<D,Multiindex> packet
         *
         * Uses the the transformer PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>to(packet)
         * and the packet.parameters() function
         */
        template<class MultiIndex>
        void store_packet(const double& time_, const waveblocks::wavepackets::ScalarHaWp<D,MultiIndex>& packetto)
        {
            select_file_writespace();
            //PacketToCoefficients<Packet>::to(packet) get coefficient matrix;
            //auto arg=PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packet);
            //Eigen::Matrix<complex_t, Eigen::Dynamic, 1> cmat=PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packet);
            if(wrlist["coefficients"])
            {
                ctype* myc=transform(PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packetto));
                coeffs->write(myc,mytype_,celemspace,cspace);
            }
            //std::cout<<cmat.rows()<<'\t'<<cmat.cols()<<'\n';
            const auto& params = packetto.parameters();

        //  Eigen::Matrix<complex_t,D,D> myQ = params.Q();//complex_t = std::complex<double>
        //  Eigen::Matrix<complex_t,D,D> myP = params.P();

        //  std::complex<double>* tmp(myQ.data());
            if(wrlist["packet"])
            {
                ctype* myQ=transform(params.Q());
                Qs->write(myQ,mytype_,QPelemspace,Qspace);
                ctype* myP=transform(params.P());
                Ps->write(myP,mytype_,QPelemspace,Pspace);

                ctype* myq=transform(params.q());
                qs->write(myq,mytype_,qpelemspace,qspace);
                ctype* myp=transform(params.p());
                ps->write(myp,mytype_,qpelemspace,pspace);

                ctype* myS=transform(params.S());
                Ss->write(myS,mytype_,Selemspace,Sspace);
            }
            if(wrlist["timegrid"])
            {
                double* mtime = new double;
                //*mtime=time_;
                *mtime=1.*current_index-1;
                times->write(mtime,PredType::NATIVE_DOUBLE,timelemspace,timespace);
            }
        }
        /**
         * @brief store energies in a timestep
         * @param epot_
         * @param ekin_
         */
        void store_energies(double epot_,double ekin_)
        {
            if(wrlist["energy"])
            {
                double* en = new double[3];
                en[0]=epot_;
                en[1]=ekin_;
                en[2]=epot_+ekin_;
                energys->write(en,PredType::NATIVE_DOUBLE,energyelemspace,energyspace);
            }
        }
        /**
         * @brief advance witer after timestep
         */
        void advance(void)
        {
            increase_index();
            setup_extensions();
            extend_datasets();
            update_filespace();

        }
        /**
         * @brief reverse last dataset extension
         */
        void cleanup(void)
        {
            decrease_index();
            setup_extensions();
            extend_datasets();
        }

    private:
        H5std_string filename_;///<identifier for filename
        CompType mytype_;///Declaration of H5:CompType member
        H5File file_; ///H5File member
        bool writeflags[4]={true,true,false,false}; ///bool array for flags writer [0]=packet [1]=coefficients [2]=energy [3]=timegrid
        std::map<std::string,bool> wrlist={{"packet",1},{"coefficients",1},{"energy",0},{"timegrid",0}};

        double dref=0.;///fillvalue for energys for allocation
        H5std_string packet_group_string="/Pi";///String for H5Group to save packet to. Default:Pi
        H5std_string datablock_string="/datablock_0";///String for H5Group for datablock.default. datablock_0
        H5std_string coefficient_group_string="/coefficients";///String for H5Group of coefficients. Default:coefficients
        H5std_string wavepacket_group_string="/wavepacket";///String for H5Group for packet and coefficients. Default:wavepacket

        DSetCreatPropList plist_qp;///PropList for packet.q() packet.p()
        DSetCreatPropList plist_QP;///PropList for packet.Q() packet.P()
        DSetCreatPropList plist_S;///PropList for packet.S()
        DSetCreatPropList plist_energy;///PropList for energies
        DSetCreatPropList plist_c;///PropList for coefficients
        DSetCreatPropList plist_time;///PropList for timegrid

        const int RANK1=1;///rank 1 identifier
        const int RANK2=2;///rank 2 identifier
        const int RANK3=3;///rank 3 identifier

        const hsize_t maxdims1[1]={H5S_UNLIMITED};///max dim identifier for rank1 for extension
        const hsize_t maxdims2[2]={H5S_UNLIMITED,H5S_UNLIMITED};///max dim identifier for rank2 for extension
        const hsize_t maxdims3[3]={H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};///max dim identifier for rank3 for extension

        hsize_t exqp[3];///extenstion for packet.q() packet.p()
        hsize_t exQP[3];///extenstion for packet.Q() packet.P()
        hsize_t exS[3];///extenstion for packet.S()
        hsize_t exenergy[2];///extenstion for energies
        hsize_t exc[2];///extenstion for coefficients
        hsize_t extime[1];///extenstion for timegrid

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

        int current_index=1;///current index used to determine position in file in time dimension

        std::shared_ptr<Group> gblock;///group for datablock
        std::shared_ptr<Group> gpacket;///group for packet
        std::shared_ptr<Group> gPi;///group for matrices in packet
        std::shared_ptr<Group> gcoefficient;///group for coefficients in packet
        std::shared_ptr<Group> genergy;///group for energies
        std::shared_ptr<Group> gtimegrid;///group for timegrid

        H5std_string q="/q";///name for packet.q()
        DataSpace qspace;///space for packet.q() in file
        std::shared_ptr<DataSet> qs;///dataset for packet.q() in file
        H5std_string p="/p";///name for packet.p()
        DataSpace pspace;///space for packet.p() in file
        std::shared_ptr<DataSet> ps;///dataset for packet.p() in file
        H5std_string Q="/Q";///name for packet.Q()
        DataSpace Qspace;///space for packet.Q() in file
        std::shared_ptr<DataSet> Qs;///dataset for packet.Q() in file
        H5std_string P="/P";///name for packet.P()
        DataSpace Pspace;///space for packet.P() in file
        std::shared_ptr<DataSet> Ps;///dataset for packet.P() in file
        H5std_string S="/S";///name for packet.S()
        DataSpace Sspace;///space for packet.S() in file
        std::shared_ptr<DataSet> Ss;///dataset for packet.S() in file
        H5std_string c="/c_0";///name for coefficients
        DataSpace cspace;///space for coefficients in file
        std::shared_ptr<DataSet> coeffs;///dataset for coefficients in file
        H5std_string time="/timegrid";///name for timegrid
        DataSpace timespace;///space for timegrid in file
        std::shared_ptr<DataSet> times;///dataset for timegrid in file
        H5std_string energies="/energies";///name for energies
        DataSpace energyspace;///space for energies in file
        std::shared_ptr<DataSet> energys;///dataset for energies in file
    };


    }
}
