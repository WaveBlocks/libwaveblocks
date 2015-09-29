#include <iostream>
#include <sstream>
#include <Eigen/Core>
#include "../hdf5/eigen3-hdf5.hpp"
#include "../types.hpp"

namespace waveblocks {
  namespace utilities {
    namespace packetWriter {
      template<class T>
      std::string toString(T in) {
        std::ostringstream strs;
        strs << in;
        return strs.str();
      }
      
     template<class Subtype, class Packet, class Phase>
      class Abstract {
      std::string file_name;
      H5::H5File file;
      public:
        Abstract(const std::string& file_name) : file_name(file_name), file(file_name, H5F_ACC_TRUNC) {}
        ~Abstract() {file.close();}

        void store_packet(const double& t, const Packet& packet, const Phase& S) {
          static_cast<Subtype*>(this)->store_packet_implementation(t,packet,S);
        }
        
        template<class M>
        void save_matrix(const M& mat, const std::string& matrix_name)
        {
            EigenHDF5::save(file, matrix_name, mat);
        }


        template<class M>
        void save_matrix(M& mat, const std::string& matrix_name)
        {
            EigenHDF5::load(file, matrix_name, mat);
        }

        template<class T>
        void save_scalar(const T& x, const std::string& scalar_name) {
            EigenHDF5::save_scalar_attribute(file,scalar_name,x);
        }

    };

    template<class Packet>
    class Standard;

    template<int D, class MultiIndex>
    struct Standard<ScalarHaWp<D,MultiIndex>> : public Abstract<Standard<ScalarHaWp<D,MultiIndex>>,ScalarHaWp<D,MultiIndex>, complex_t> {
      using Phase = complex_t;
      using Packet = ScalarHaWp<D,MultiIndex>;
      using Super = Abstract<Standard<ScalarHaWp<D,MultiIndex>>,ScalarHaWp<D,MultiIndex>, Phase>;
      Standard(const std::string& file_name) : Super(file_name) {}
    
      void store_packet_implementation(const double& time, const Packet& packet, const Phase& S) {
        
        std::string t = toString<double>(time);
        // coefficients
        Super::save_matrix(PacketToCoefficients<Packet>::to(packet), "coefficients@" + t);
        
        //
        const auto& params = packet.parameters();
        Super::save_matrix(params.q, "q@" + t);
        Super::save_matrix(params.p, "p@" + t);
        Super::save_matrix(params.Q, "Q@" + t);
        Super::save_matrix(params.P, "P@" + t);
        Super::save_scalar(S, "S@" + t);

        // eps
        Super::save_scalar(packet.eps(), "eps@" + t);
      }
    };
    }
    template<class Packet>
    using PacketWriter = packetWriter::Standard<Packet>;
  }
}
