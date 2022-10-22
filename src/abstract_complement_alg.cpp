
#include "abstract_complement_alg.hpp"

using namespace kofola;

/// equality operator
bool kofola::operator==(
  const abstract_complement_alg::mstate& lhs,
  const abstract_complement_alg::mstate& rhs)
{
  return lhs.eq(rhs);
}

/// disequality operator
bool kofola::operator!=(
  const abstract_complement_alg::mstate& lhs,
  const abstract_complement_alg::mstate& rhs)
{
  return !(lhs == rhs);
}

/// ordering relation
bool kofola::operator<(
  const abstract_complement_alg::mstate& lhs,
  const abstract_complement_alg::mstate& rhs)
{
  return lhs.lt(rhs);
}

std::ostream& kofola::operator<<(std::ostream& os, const abstract_complement_alg::mstate& ms)
{
  os << ms.to_string();
  return os;
}
