#define that (*static_cast<Subtype *>(this))
#define GET_TYPES(NAME, S, B, N, D)                                            \
  using Basis = B<N, D>;                                                       \
  using Subtype = S<B, N, D>;                                                  \
  using potential_type = typename Basis::potential_type;                       \
  using jacobian_type = typename Basis::jacobian_type;                         \
  using hessian_type = typename Basis::hessian_type;                           \
  using local_quadratic_type = typename Basis::local_quadratic_type;           \
  using local_remainder_type = typename Basis::local_remainder_type;           \
  using potential_evaluation_type = typename Basis::potential_evaluation_type; \
  using jacobian_evaluation_type = typename Basis::jacobian_evaluation_type;   \
  using hessian_evaluation_type = typename Basis::hessian_evaluation_type;     \
  using potential_return_type = typename Basis::potential_return_type;         \
  using jacobian_return_type = typename Basis::jacobian_return_type;           \
  using hessian_return_type = typename Basis::hessian_return_type;             \
  using self_type = NAME<S, B, N, D>
