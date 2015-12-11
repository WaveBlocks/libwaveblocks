#pragma once

#include <chrono>


namespace waveblocks {
    namespace utilities {
        class Timer
        {
        public:
            double millis() const
            {
                return std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1,1000> > >(stop_ - start_).count();
            }

            double seconds() const
            {
                return std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1,1> > >(stop_ - start_).count();
            }

            void start()
            {
                start_ = std::chrono::high_resolution_clock::now();
            }

            void stop()
            {
                stop_ = std::chrono::high_resolution_clock::now();
            }

        private:
            std::chrono::high_resolution_clock::time_point start_;
            std::chrono::high_resolution_clock::time_point stop_;
        };
    }
}
