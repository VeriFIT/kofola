// partial determinization algorithm for initial almost deterministic
// accepting components (IADACs)

#include "determinize_alg_iadacs.hpp"

#include <climits>
#include <stdexcept>
#include <unordered_set>

using namespace kofola;
using mstate_p = abstract_determinize_alg::mstate_p;
using mstate_col = abstract_determinize_alg::mstate_col;


namespace { // {{{

/// partial macrostate: the (compressed) mapping of slots to runs
class mstate_det_iadacs : public abstract_determinize_alg::mstate
{ // {{{
private: // DATA MEMBERS

  std::vector<unsigned> mapping_;

public: // METHODS

  /// constructor
  explicit mstate_det_iadacs(const std::vector<unsigned>& mapping) : mapping_(mapping) { }

  virtual std::string to_string() const override;
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_det_iadacs() override { }

  friend class kofola::determinize_iadacs;
}; // mstate_det_iadacs }}}


std::string mstate_det_iadacs::to_string() const
{
  return std::string("[IADACs-det: MAPPING=") + std::to_string(this->mapping_) + "]";
}

bool mstate_det_iadacs::eq(const mstate& rhs) const
{
  const mstate_det_iadacs* rhs_iad = dynamic_cast<const mstate_det_iadacs*>(&rhs);
  assert(rhs_iad);
  return this->mapping_ == rhs_iad->mapping_;
}

bool mstate_det_iadacs::lt(const mstate& rhs) const
{
  const mstate_det_iadacs* rhs_iad = dynamic_cast<const mstate_det_iadacs*>(&rhs);
  assert(rhs_iad);
  return this->mapping_ < rhs_iad->mapping_;   // lexicographic
}

/// moves the defined slots of 'mapping' to its front (keeping their order)
std::vector<unsigned> compress(const std::vector<unsigned>& mapping)
{ // {{{
  std::vector<unsigned> compressed(mapping.size(), UINT_MAX);
  unsigned dst = 0;

  for (unsigned run : mapping) {
    if (UINT_MAX != run) { compressed[dst++] = run; }
  }

  return compressed;
} // compress() }}}

} // anonymous namespace }}}


determinize_iadacs::determinize_iadacs(const cmpl_info& info, unsigned part_index, unsigned runs_bound)
  : abstract_determinize_alg(info, part_index),
    runs_bound_(runs_bound),
    acc_cond_(info.aut_->get_acceptance(), runs_bound)
{
  DEBUG_PRINT_LN("IADACs determinization of partition " + std::to_string(part_index) +
                 ": runs bound = " + std::to_string(runs_bound) +
                 ", acc = " + std::to_string(this->acc_cond_.get_acc_cond().get_acceptance()));
}

determinize_iadacs::~determinize_iadacs()
{ }

mstate_p determinize_iadacs::get_init()
{ // {{{
  std::vector<unsigned> init_mapping(this->runs_bound_, UINT_MAX);

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  auto it = this->info_.st_to_part_map_.find(orig_init);
  if (this->info_.st_to_part_map_.end() != it &&
      it->second == static_cast<int>(this->part_index_)) {
    assert(0 < this->runs_bound_);
    init_mapping[0] = orig_init;
  }

  return std::make_shared<mstate_det_iadacs>(init_mapping);
} // get_init() }}}

mstate_col determinize_iadacs::get_succ(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol)
{ // {{{
  const mstate_det_iadacs* src_iad = dynamic_cast<const mstate_det_iadacs*>(src);
  assert(src_iad);

  DEBUG_PRINT_LN("IADACs determinization successor");
  DEBUG_PRINT_LN("glob_reached = " + std::to_string(glob_reached));
  DEBUG_PRINT_LN("src = " + std::to_string(*src_iad));
  DEBUG_PRINT_LN("symbol = " + std::to_string(symbol));

  const std::vector<unsigned>& src_mapping = src_iad->mapping_;
  const unsigned bound = this->runs_bound_;

  // The uncompressed successor mapping: every run keeps its slot, a run that
  // leaves its SCC or merges into a run with a smaller label leaves a hole.
  std::vector<unsigned> uncompressed(bound, UINT_MAX);
  std::unordered_set<unsigned> defined_runs;
  unsigned i = 0;

  for (unsigned run : src_mapping) {
    if (UINT_MAX == run) { break; }   // end of the defined runs

    const unsigned run_scc = this->info_.scc_info_.scc_of(run);
    for (const auto& t : this->info_.aut_->out(run)) {
      if (!bdd_implies(symbol, t.cond)) { continue; }
      if (this->info_.scc_info_.scc_of(t.dst) != run_scc) { continue; }
      if (defined_runs.end() != defined_runs.find(t.dst)) { continue; }   // merged

      uncompressed[i] = t.dst;
      defined_runs.insert(t.dst);
      break;   // the SCC is deterministic inside
    }
    ++i;   // keep the slot even if the run is gone
  }

  // runs entering the partition, in the order of the state numbers
  const auto& sccs_in_part = this->info_.part_to_scc_map_.at(this->part_index_);
  for (unsigned run : glob_reached) {
    if (sccs_in_part.end() == sccs_in_part.find(this->info_.scc_info_.scc_of(run))) { continue; }
    if (defined_runs.end() != defined_runs.find(run)) { continue; }

    // The complementation writes past the mapping here (see the open issue in
    // complement_init_almost_det::get_succ_track()); fail loudly instead,
    // since dropping the run would make the result unsound.
    if (i >= bound) {
      throw std::runtime_error("IADACs determinization: more than " +
        std::to_string(bound) + " simultaneous runs in partition " +
        std::to_string(this->part_index_) + " (runs bound exceeded)");
    }

    uncompressed[i++] = run;
    defined_runs.insert(run);
  }

  std::vector<unsigned> succ_mapping = compress(uncompressed);

  // The colours of the labels that keep carrying the same run: the longest
  // prefix of slots that were not shifted by the compression and whose run
  // existed already in 'src'.
  std::set<unsigned> cols;
  unsigned matching_prefix_len = 0;
  for (; matching_prefix_len < bound; ++matching_prefix_len) {
    const unsigned succ_run = succ_mapping[matching_prefix_len];
    const unsigned src_run = src_mapping[matching_prefix_len];

    if (succ_run != uncompressed[matching_prefix_len] || UINT_MAX == succ_run) { break; }
    if (UINT_MAX == src_run) { break; }

    for (const auto& t : this->info_.aut_->out(src_run)) {
      if (!bdd_implies(symbol, t.cond)) { continue; }
      if (t.dst != succ_run) { continue; }

      for (unsigned col : t.acc.sets()) {
        cols.insert(this->acc_cond_.map_colour(col, matching_prefix_len));
      }
      break;
    }
  }

  // all the other labels are discontinued
  for (unsigned j = matching_prefix_len; j < bound; ++j) {
    cols.insert(this->acc_cond_.get_fin_mark(j));
  }

  return {std::make_shared<mstate_det_iadacs>(succ_mapping), cols};
} // get_succ() }}}
