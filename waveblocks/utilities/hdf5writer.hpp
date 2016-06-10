//
// Created by flofrei on 01.06.16.
//
#pragma once

//need H5 Cpp header
#include "H5Cpp.h"
//using strings
#include <string>
#include <Eigen/Core>

#include "waveblocks/wavepackets/hawp_commons.hpp"
#include "waveblocks/wavepackets/hawp_paramset.hpp"


//utilities::PacketWriter<wavepackets::ScalarHaWp<D,MultiIndex >> wirter;


namespace waveblocks {
    namespace utilities {
    //using namespaces for convenience
    using namespace H5;

//template<int D>
class hdf5writer{
public:
    hdf5writer(std::string filename):filename_(filename),mytype(sizeof(instanceof)),file(filename_,H5F_ACC_TRUNC)
    {
        mytype.insertMember( "r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
        mytype.insertMember( "i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);

        //hardcode the chunk dimension to for qp
        hsize_t chunk_dims1[3]={1,2,1};
        plist_qp.setChunk(RANK3,chunk_dims1);
        plist_qp.setFillValue(mytype,&instanceof);

        //hardcode the chunk dimension to 1x2x2 for QP
        hsize_t chunk_dims2[3]={1,2,2};
        plist_QP.setChunk(RANK3,chunk_dims2);
        plist_QP.setFillValue(mytype,&instanceof);

        //hardcode the chunk dimension to for Sepotkintot
        hsize_t chunk_dims3[2]={1,1};
        plist_Sepotkintot.setChunk(RANK2,chunk_dims3);
        plist_Sepotkintot.setFillValue(mytype,&instanceof);

        //hardcode the chunk dimension for the coefficients
        hsize_t chunk_dims4[2]={1,16};
        plist_c.setChunk(RANK2,chunk_dims4);
        plist_c.setFillValue(mytype,&instanceof);

        groupdatablock="/datablock_0";
        groupwavepacket="/datablock_0/wavepacket";
        groupPi="/datablock_0/wavepacket/Pi";
        groupcoefficients="/datablock_0/wavepacket/coefficients";
        groupenergies="/datablock_0/energies";

        //names
        q="/datablock_0/wavepacket/Pi/q";
        p="/datablock_0/wavepacket/Pi/p";
        Q="/datablock_0/wavepacket/Pi/Q";
        P="/datablock_0/wavepacket/Pi/P";
        S="/datablock_0/wavepacket/Pi/S";
        c="/datablock_0/wavepacket/coefficients/c_0";
        epot="/datablock_0/energies/epot";
        ekin="/datablock_0/energies/ekin";
        etot="/datablock_0/energies/etot";
    }

    void setup_dataspaces(void)
    {
        //set up q and p
        hsize_t dim1[3]={1,2,1};
        DataSpace d1(RANK3,dim1,maxdims3);
        qspace=d1;
        pspace=d1;
        //set up Q and P
        hsize_t dim2[3]={1,2,2};
        DataSpace d2(RANK3,dim2,maxdims3);
        Qspace=d2;
        Pspace=d2;
        //set up S etot epot ekin
        hsize_t dim3[2]={1,1};
        DataSpace d3(RANK2,dim3,maxdims2);
        Sspace=d3;
        epotspace=d3;
        ekinspace=d3;
        etotspace=d3;
        //set up coefficients
        hsize_t dim4[2]={1,16};
        DataSpace d4(RANK2,dim4,maxdims2);
        cspace=d4;

        //qp
        hsize_t count1[3]={1,2,1};
        hsize_t start1[3]={0,0,0};
        hsize_t stride1[3]={1,1,1};
        hsize_t block1[3]={1,1,1};
        //QP
        hsize_t count2[3]={1,2,2};
        hsize_t start2[3]={0,0,0};
        hsize_t stride2[3]={1,1,1};
        hsize_t block2[3]={1,1,1};
        //Setotpotkin
        hsize_t count3[2]={1,1};
        hsize_t start3[2]={0,0};
        hsize_t stride3[2]={1,1};
        hsize_t block3[2]={1,1};
        //coefficients
        hsize_t count4[2]={1,16};
        hsize_t start4[2]={0,0};
        hsize_t stride4[2]={1,1};
        hsize_t block4[2]={1,1};

        //select hyperslabs for elements
        //qp
        DataSpace f1(RANK3,qpelem);
        qpelemspace=f1;
        qpelemspace.selectHyperslab(H5S_SELECT_SET, count1, start1, stride1, block1);
        //QP
        DataSpace f2(RANK3,QPelem);
        QPelemspace=f2;
        QPelemspace.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
        //Senergies
        DataSpace f3(RANK2,Senergyelem);
        Senergyelemspace=f3;
        Senergyelemspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);
        //coefficients
        DataSpace f4(RANK2,celem);
        celemspace=f4;
        celemspace.selectHyperslab(H5S_SELECT_SET, count4, start4, stride4, block4);
    }

    void setup_groups(void)
    {
    gblock = new Group(file.createGroup(groupdatablock));
    gpacket = new Group(file.createGroup(groupwavepacket));
    gPi = new Group(file.createGroup(groupPi));
    gc = new Group(file.createGroup(groupcoefficients));
    ge = new Group(file.createGroup(groupenergies));
    }

    void setup_datasets(void)
    {
        //qs = std::make_shared<DataSet>(file.createDataSet(q,mytype,qspace,plist_qp));
        qs = new DataSet(file.createDataSet(q,mytype,qspace,plist_qp));
        ps = new DataSet(file.createDataSet(p,mytype,pspace,plist_qp));
        Qs = new DataSet(file.createDataSet(Q,mytype,Qspace,plist_QP));
        Ps = new DataSet(file.createDataSet(P,mytype,Pspace,plist_QP));
        Ss = new DataSet(file.createDataSet(S,mytype,Sspace,plist_Sepotkintot));
        cs = new DataSet(file.createDataSet(c,mytype,cspace,plist_c));
        epots = new DataSet(file.createDataSet(epot,mytype,epotspace,plist_Sepotkintot));
        ekins = new DataSet(file.createDataSet(ekin,mytype,ekinspace,plist_Sepotkintot));
        etots = new DataSet(file.createDataSet(etot,mytype,etotspace,plist_Sepotkintot));
    }

    void select_writespace(void)
    {
        int tr=current_index-1;
        //qp
        hsize_t count1[3]={1,2,1};
        hsize_t start1[3];
        start1[0]=tr;
        start1[1]=0;
        start1[2]=0;
        hsize_t stride1[3]={1,1,1};
        hsize_t block1[3]={1,1,1};
        //QP
        hsize_t count2[3]={1,2,2};
        hsize_t start2[3];
        start2[0]=tr;
        start2[1]=0;
        start2[2]=0;
        hsize_t stride2[3]={1,1,1};
        hsize_t block2[3]={1,1,1};
        //Setotpotkin
        hsize_t count3[2]={1,1};
        hsize_t start3[2];
        start3[0]=tr;
        start3[1]=0;
        hsize_t stride3[2]={1,1};
        hsize_t block3[2]={1,1};
        //coefficients
        hsize_t count4[2]={1,16};
        hsize_t start4[2];
        start4[0]=tr;
        start4[1]=0;
        hsize_t stride4[2]={1,1};
        hsize_t block4[2]={1,1};

        //qp
        qspace.selectHyperslab(H5S_SELECT_SET, count1, start1, stride1, block1);
        pspace.selectHyperslab(H5S_SELECT_SET, count1, start1, stride1, block1);
        //QP
        Qspace.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
        Pspace.selectHyperslab(H5S_SELECT_SET, count2, start2, stride2, block2);
        //Senergies
        Sspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);
        epotspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);
        ekinspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);
        etotspace.selectHyperslab(H5S_SELECT_SET, count3, start3, stride3, block3);
        //coefficients
        cspace.selectHyperslab(H5S_SELECT_SET, count4, start4, stride4, block4);
    }

    void store_energies(const double& time,real_t epot_,real_t ekin_)
    {
        ctype* myepot = new ctype;
        myepot->real=epot_;
        myepot->imag=0.;
        ctype* myekin = new ctype;
        myekin->real=ekin_;
        myekin->imag=0.;
        real_t etot_=epot_+ekin_;
        ctype* myetot = new ctype;
        myetot->real=etot_;
        myetot->imag=0.;
        ekins->write(myekin,mytype,Senergyelemspace,ekinspace);
        epots->write(myepot,mytype,Senergyelemspace,epotspace);
        etots->write(myetot,mytype,Senergyelemspace,etotspace);
    }

    template<int D, class MultiIndex>
    void store_packet(const double& time, const waveblocks::wavepackets::ScalarHaWp<D,MultiIndex>& packet)
    {
    constexpr int dim=D;
    //PacketToCoefficients<Packet>::to(packet) get coefficient matrix;
    //auto arg=PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packet);
    //Eigen::Matrix<complex_t, Eigen::Dynamic, 1> cmat=PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packet);
    ctype* myc=transform(PacketToCoefficients<wavepackets::ScalarHaWp<D,MultiIndex>>::to(packet));
    cs->write(myc,mytype,celemspace,cspace);
    //std::cout<<cmat.rows()<<'\t'<<cmat.cols()<<'\n';
    const auto& params = packet.parameters();

//  Eigen::Matrix<complex_t,D,D> myQ = params.Q();//complex_t = std::complex<double>
//  Eigen::Matrix<complex_t,D,D> myP = params.P();

//  std::complex<double>* tmp(myQ.data());

    ctype* myQ=transform<D>(params.Q());
    Qs->write(myQ,mytype,QPelemspace,Qspace);
    ctype* myP=transform<D>(params.P());
    Ps->write(myP,mytype,QPelemspace,Pspace);

    ctype* myq=transform<D>(params.q());
    qs->write(myq,mytype,qpelemspace,qspace);
    ctype* myp=transform<D>(params.p());
    ps->write(myp,mytype,qpelemspace,pspace);

    ctype* myS=transform(params.S());
    Ss->write(myS,mytype,Senergyelemspace,Sspace);
    }

    void decrease_index(void)
    {
        current_index-=1;
    }

    void setup_extension(void)
    {
        exqp[0]=current_index;
        exqp[1]=2;
        exqp[2]=1;
        exQP[0]=current_index;
        exQP[1]=2;
        exQP[2]=2;
        exSepotkintot[0]=current_index;
        exSepotkintot[1]=1;
        exc[0]=current_index;
        exc[1]=16;
    }

    void extend_datasets(void)
    {
        //qp
        qs->extend(exqp);
        ps->extend(exqp);
        //QP
        Qs->extend(exQP);
        Ps->extend(exQP);
        //Senergies
        Ss->extend(exSepotkintot);
        epots->extend(exSepotkintot);
        ekins->extend(exSepotkintot);
        etots->extend(exSepotkintot);
        //coefficients
        cs->extend(exc);
    }

    void reduce_datasets(void)
    {
        decrease_index();
        setup_extension();
        extend_datasets();
    }

    void overwrite_dataspaces(void)
    {
        //qp
        qspace=qs->getSpace();
        pspace=ps->getSpace();
        //QP
        Qspace=Qs->getSpace();
        Pspace=Ps->getSpace();
        //Senergies
        Sspace=Ss->getSpace();
        epotspace=epots->getSpace();
        ekinspace=ekins->getSpace();
        etotspace=etots->getSpace();
        //coefficients
        cspace=cs->getSpace();
    }

    void increase_index(void)
    {
        current_index+=1;
    }

    void cleanup(void)
    {
        reduce_datasets();
        delete qs;
        delete ps;
        delete Qs;
        delete Ps;
        delete Ss;
        delete cs;
        delete epots;
        delete ekins;
        delete etots;

        delete gblock;
        delete gpacket;
        delete gc;
        delete ge;
        delete gPi;
    }

    //define the complex type which gets written to memory
    struct ctype{ //our complex datatype
        double real=0.;
        double imag=0.;
    };
    ctype instanceof;

    template<int D>
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
    //delete tmp;
    return newdat;
    }

    template<int D>
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
    //transforms a std::complex variable in a ctype variable
    ctype* transform(complex_t arg)
    {
    ctype* newarg = new ctype;
    newarg->real=arg.real();
    newarg->imag=arg.imag();
    return newarg;
    }

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

private:
    //filename
    H5std_string filename_;
    //hdf5 writabl complex type
    CompType mytype;
    //create file
    H5File file;
    //property list for chunking the file and extensions
    DSetCreatPropList plist_qp;
    DSetCreatPropList plist_QP;
    DSetCreatPropList plist_Sepotkintot;
    DSetCreatPropList plist_c;
    //Ranks
    const int RANK3=3;
    const int RANK2=2;
    //const int RANK1=1;
    //maxdims
    //hsize_t maxdims1[1]={H5S_UNLIMITED};
    hsize_t maxdims2[2]={H5S_UNLIMITED,H5S_UNLIMITED};
    hsize_t maxdims3[3]={H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};
    //extensions
    hsize_t exqp[3];
    hsize_t exQP[3];
    hsize_t exSepotkintot[2];
    hsize_t exc[2];

    //one element space for source
    hsize_t qpelem[3] = {1,2,1};
    DataSpace qpelemspace;
    hsize_t QPelem[3] = {1,2,2};
    DataSpace QPelemspace;
    hsize_t Senergyelem[2] = {1,1};
    DataSpace Senergyelemspace;
    hsize_t celem[2] ={1,16};
    DataSpace celemspace;

    //index
    int current_index=1;

    //packet
    //place holder for q native double dim 2x1 or ctype 2x1 imag==0
    H5std_string q;
    DataSpace qspace;
    //std::shared_ptr<DataSet> qs;
    DataSet* qs;
    //place holder for p native double dim 2x1 or ctype 2x1 imag==0
    H5std_string p;
    DataSpace pspace;
    DataSet* ps;
    //place holder for Q matrix ctype 2x2
    H5std_string Q;
    DataSpace Qspace;
    DataSet* Qs;
    //place holder for P matrix ctype 2x2
    H5std_string P;
    DataSpace Pspace;
    DataSet* Ps;
    //place holder for S scalar ctype dim 1x1
    H5std_string S;
    DataSpace Sspace;
    DataSet* Ss;

    //place holder for coefficients matrix PacketToCoefficients<Packet>::to(packet)
    H5std_string c; //ctype dim 16x1 alternative timescale x 16 matrix
    DataSpace cspace;
    DataSet* cs;

    //place holder for energies;
    H5std_string epot; //native double 1dim or ctype 1dim imag==0
    DataSpace epotspace;
    DataSet* epots;
    H5std_string ekin; //native double 1dim or ctype 1dim imag==0
    DataSpace ekinspace;
    DataSet* ekins;
    H5std_string etot; //native double 1dim or ctype 1dim imag==0
    DataSpace etotspace;
    DataSet* etots;

    H5std_string groupdatablock;
    Group* gblock;
    H5std_string groupwavepacket;
    Group* gpacket;
    H5std_string groupPi;
    Group* gPi;
    H5std_string groupcoefficients;
    Group* gc;
    H5std_string groupenergies;
    Group* ge;
};
    }
}
