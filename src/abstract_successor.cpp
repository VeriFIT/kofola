#include "abstract_successor.hpp"

namespace kofola {
    bool operator==(const abstract_successor::mstate& lhs,
                    const abstract_successor::mstate& rhs) {
        return lhs.eq(rhs);
    }

    bool operator!=(const abstract_successor::mstate& lhs,
                    const abstract_successor::mstate& rhs) {
        return !(lhs == rhs);
    }

    bool operator<(const abstract_successor::mstate& lhs,
                   const abstract_successor::mstate& rhs) {
        return lhs.lt(rhs);
    }
}
