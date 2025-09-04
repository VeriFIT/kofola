#pragma once

namespace kofola
{
class waiting
{
private: // DATA MEMBERS

  std::set<std::set<unsigned>> states_;
  std::map<std::set<unsigned>, std::set<std::set<unsigned>>> trans_;
  std::map<std::set<unsigned>, std::set<std::set<unsigned>>> predecessors_;

public: // METHODS

  waiting(
    const std::set<std::set<unsigned>>&                                states,
    const std::map<std::set<unsigned>, std::set<std::set<unsigned>>>&  trans) :
      states_(states), trans_(trans)
  {
    std::map<std::set<unsigned>, std::set<std::set<unsigned>>> predecessors;

    for (auto state : states) {
      predecessors.insert({state, std::set<std::set<unsigned>>()});
    }

    for (auto src : states) {
      for (auto dst : trans.at(src)) {
        predecessors[dst].insert(src);
      }
    }

    predecessors_ = predecessors;
  }

  const std::set<std::set<unsigned>>& get_states() const { return states_; }

  const std::map<std::set<unsigned>, std::set<std::set<unsigned>>>&
  get_predecessors() const
  {
    return predecessors_;
  }

  friend std::ostream& operator<<(std::ostream& os, const waiting& wait)
  {
    os << "(states_: " + std::to_string(wait.states_);
    os << ", trans_: " + std::to_string(wait.trans_);
    os << ", predecessors_: " + std::to_string(wait.predecessors_);
    os << ")";

    return os;
  }
};
}
