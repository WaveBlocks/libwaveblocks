#pragma once

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "../basic_types.hpp"
#include "../continuous_sqrt.hpp"


namespace waveblocks {
    namespace wavepackets {
        /**
         * \brief This class represents the Hagedorn parameter set \f$ \Pi = \{q, p, Q, P, S\} \f$.
         *
         * The first two parameters \f$ q \f$ and \f$ p \f$ are \f$ D \f$ dimensional
         * real-valued vectors. The second two \f$ Q \f$ and \f$ P \f$ are complex
         * \f$ D \times D \f$ matrices. The last parameter \f$ S \f$ is the global
         * complex phase.
         */
        template<dim_t D>
        struct HaWpParamSet
        {
        private:
            RMatrix<D,1> q_, p_;
            CMatrix<D,D> Q_, P_;
            complex_t S_;
            ContinuousSqrt<real_t> sqrt_detQ_;

        public:
            /** Construct a Hagedorn parameter set with default values.
             */
            HaWpParamSet()
                : q_(RMatrix<D,1>::Zero())
                , p_(RMatrix<D,1>::Zero())
                , Q_(CMatrix<D,D>::Identity())
                , P_(CMatrix<D,D>::Identity()*complex_t(0,1))
                , S_(complex_t(0,0))
                , sqrt_detQ_(1)
            { }

            /** Construct a Hagedorn parameter set by copying from another one.
             */
            HaWpParamSet(const HaWpParamSet &that)
                : q_(that.q_)
                , p_(that.p_)
                , Q_(that.Q_)
                , P_(that.P_)
                , S_(that.S_)
                , sqrt_detQ_(that.sqrt_detQ_)
            { }

            /** Construct a Hagedorn parameter set with explicit values.
             */
            HaWpParamSet(const RMatrix<D,1> &q,
                         const RMatrix<D,1> &p,
                         const CMatrix<D,D> &Q,
                         const CMatrix<D,D> &P,
                         const complex_t &S)
                : q_(q)
                , p_(p)
                , Q_(Q)
                , P_(P)
                , S_(S)
                , sqrt_detQ_(std::sqrt(Q_.determinant()))
            { }

            // Undocumented
            HaWpParamSet(const RMatrix<D,1> &q,
                         const RMatrix<D,1> &p,
                         const CMatrix<D,D> &Q,
                         const CMatrix<D,D> &P,
                         const complex_t &S,
                         ContinuousSqrt<real_t> sqrt_detQ)
                : q_(q)
                , p_(p)
                , Q_(Q)
                , P_(P)
                , S_(S)
                , sqrt_detQ_(sqrt_detQ)
            { }

            /** Construct a Hagedorn parameter set by assigning from another one.
             */
            HaWpParamSet &operator=(const HaWpParamSet &that)
            {
                q_ = that.q_;
                p_ = that.p_;
                Q_ = that.Q_;
                P_ = that.P_;
                S_ = that.S_;
                sqrt_detQ_ = that.sqrt_detQ_;
                return *this;
            }

            /** \brief Get the parameter \f$ q \f$ */
            inline RMatrix<D,1> const& q() const {return q_;}

            /** \brief Get the parameter \f$ p \f$ */
            inline RMatrix<D,1> const& p() const {return p_;}

            /** \brief Get the parameter \f$ Q \f$ */
            inline CMatrix<D,D> const& Q() const {return Q_;}

            /** \brief Get the parameter \f$ P \f$ */
            inline CMatrix<D,D> const& P() const {return P_;}

            /** \brief Get the parameter \f$ S \f$ */
            inline complex_t const& S() const {return S_;}

            // Undocumented
            inline complex_t const sdQ() const {return sqrt_detQ_();}

            /** \brief Set the parameter \f$ q \f$ */
            inline void q(const RMatrix<D,1>& qnew) {q_ = qnew;}

            /** \brief Set the parameter \f$ p \f$ */
            inline void p(const RMatrix<D,1>& pnew) {p_ = pnew;}

            /** \brief Set the parameter \f$ Q \f$ */
            inline void Q(const CMatrix<D,D>& Qnew) {Q_ = Qnew; resync();}

            /** \brief Set the parameter \f$ P \f$ */
            inline void P(const CMatrix<D,D>& Pnew) {P_ = Pnew;}

            /** \brief Set the parameter \f$ S \f$ */
            inline void S(const complex_t& Snew) {S_ = Snew;}

            /** \brief Update the parameter \f$ q \f$ */
            inline void updateq(const RMatrix<D,1>& qnew) {q_ += qnew;}

            /** \brief Update the parameter \f$ p \f$ */
            inline void updatep(const RMatrix<D,1>& pnew) {p_ += pnew;}

            /** \brief Update the parameter \f$ Q \f$ */
            inline void updateQ(const CMatrix<D,D>& Qnew) {Q_ += Qnew; resync();}

            /** \brief Update the parameter \f$ P \f$ */
            inline void updateP(const CMatrix<D,D>& Pnew) {P_ += Pnew;}

            /** \brief Update the parameter \f$ S \f$ */
            inline void updateS(const complex_t& Snew) {S_ += Snew;}

            /* Undocumented
             * Compute the continuous square root of \f$ \det Q \f$ after an update
             * of the \f$ Q \f$ parameter.
             */
            inline void resync() {sqrt_detQ_(Q_.determinant());}

            /**
             * Check the compatibility relations
             * \f[\begin{align}
             *      Q^{\texttt{H}} P - P^{\texttt{H}} Q & = 2i I \\
             *      Q^{\texttt{T}} P - P^{\texttt{T}} Q & = 0
             *    \end{align}\f]
             */
            bool compatible() const
            {
                const double tol = 1e-10;

                CMatrix<D,D> C1 = Q_.adjoint()*P_ - P_.adjoint()*Q_ - CMatrix<D,D>::Identity()*complex_t(0,2);
                CMatrix<D,D> C2 = Q_.transpose()*P_ - P_.transpose()*Q_;

                return C1.template lpNorm<1>() < tol && C2.template lpNorm<1>() < tol;
            }

            /**
             * Mix the two parameter sets \f$ \Pi_i \f$ and \f$ \Pi_j \f$
             * from the bra and the ket wavepackets \f$ \Phi\left[\Pi_i\right] \f$
             * and \f$ \Phi^\prime\left[\Pi_j\right] \f$.
             *
             * \param ket The parameter set \f$ \Pi_j \f$ from the ket part wavepacket.
             * \return The mixed parameters \f$ q_0 \f$ and \f$ Q_0 \f$.
             */
            std::pair< RMatrix<D,1>, RMatrix<D,D> >
            mix(const HaWpParamSet<D>& ket) const
            {
                // Mix the parameters
                CMatrix<D,D> Gbra = P_ * Q_.inverse();
                CMatrix<D,D> Gket = ket.P_ * ket.Q_.inverse();
                RMatrix<D,D> G = (Gket - Gbra.adjoint()).imag();
                RMatrix<D,1> g = (Gket*ket.q_ - Gbra.adjoint()*q_).imag();
                RMatrix<D,1> q0 = G.inverse() * g;
                RMatrix<D,D> Q0 = 0.5 * G;
                // We can not avoid the matrix root
                RMatrix<D,D> Qs = Q0.sqrt().inverse();
                return {q0,Qs};
            }
        };

        template<dim_t D>
        std::ostream &operator<<(std::ostream &out, const HaWpParamSet<D> &parameters)
        {
            Eigen::IOFormat CleanFmt(4, 0, ", ", "\n     ", "[", "]");

            out << "HaWpParamSet {\n";
            out << "  q: " << parameters.q().format(CleanFmt) << '\n';
            out << "  p: " << parameters.p().format(CleanFmt) << '\n';
            out << "  Q: " << parameters.Q().format(CleanFmt) << '\n';
            out << "  P: " << parameters.P().format(CleanFmt) << '\n';
            out << "  S: " << parameters.S() << '\n';
            out << "  sqrt(detQ): " << parameters.sdQ() << '\n';
            out << "  compatible(): " << (parameters.compatible() ? "yes" : "no") << '\n';
            out << "}" << std::endl;
            return out;
        }
    }
}
