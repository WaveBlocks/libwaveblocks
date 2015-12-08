#pragma once

#include "../../types.hpp"
#include "../../utilities/evaluations.hpp"

namespace waveblocks
{
  namespace matrixPotentials
  {
    namespace modules
    {
      namespace leadingLevelOwner {

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
    }
  }
 template<class Owned>
      using LeadingLevelOwner = matrixPotentials::modules::leadingLevelOwner::Standard<Owned>;
}
