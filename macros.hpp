#define that (*static_cast<Subtype*>(this))
#define IMPORT_TYPES_FROM(B,N,D)                           \
  using potential_type = typename B<N,D>::potential_type;                       \
  using jacobian_type = typename B<N,D>::jacobian_type;                         \
  using hessian_type = typename B<N,D>::hessian_type;                           \
  using local_quadratic_type = typename B<N,D>::local_quadratic_type;           \
  using local_remainder_type = typename B<N,D>::local_remainder_type;           \
  using potential_evaluation_type = typename B<N,D>::potential_evaluation_type; \
  using jacobian_evaluation_type = typename B<N,D>::jacobian_evaluation_type;   \
  using hessian_evaluation_type = typename B<N,D>::hessian_evaluation_type;     \
  using potential_return_type = typename B<N,D>::potential_return_type;         \
  using jacobian_return_type = typename B<N,D>::jacobian_return_type;           \
  using hessian_return_type = typename B<N,D>::hessian_return_type             \
  
  
