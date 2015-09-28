#pragma once
#include "../../macros.hpp"
#include "../../types.hpp"
#include "../../utilities/evaluations.hpp"

namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace modules
    {
      namespace leadingLevelOwner {
       /**
       * \brief Abstract class for local quadratic evaluation
       *
       * A matrix potential inheriting an implementation of this module
       * can evaluate the local quadratic approximation of its' potential
       * elementwise
       *
       * This makes use of the CRTPattern
       *
       * \tparam Subtype The type extending this interface (used for static polymorphism)
       * \tparam Basis
       * Which basis (bases::Eigen or bases::Canonical) the potential is given in
       * \tparam N
       * Number of levels (dimension of square matrix when evaluated)
       * \tparam D
       * Dimension of argument space
       */
      template <class Owned>
      struct Standard {
      private:
        Owned owned;
      public:
        Owned& get_leading_level() {
          return owned;
        }
        const Owned& get_leading_level() const {
          return owned;
        }

        template<class... T>
        Standard(T... args) :owned(args...){}
      };


      }
      template<class Owned>
      using LeadingLevelOwner = leadingLevelOwner::Standard<Owned>;

    }
  }
}
