/**
 * @file abstract_successor.cpp
 * @author Ondrej Alexaj (xalexa09@stud.fit.vutbr.cz)
 * @brief neccessary methods for abstract_successor
 * @version 0.1
 * @date 2024-05-03
 * 
 * 
 */

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
