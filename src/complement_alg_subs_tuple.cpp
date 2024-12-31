#include "complement_alg_subs_tuple.hpp"
#include <algorithm>

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


namespace { // anonymous namespace {{{

/// partial macrostate for the given component
    class mstate_subs_tuple : public abstract_complement_alg::mstate { // {{{
    private: // DATA MEMBERS

        std::deque<std::pair<std::set<unsigned>, int> > greedy_tuple_;
        bool active_;                    // true = active ; false = track
        bool breakpoint_;
        std::set<unsigned> empty_breakpoint {};

    public: // METHODS

        /// constructor
        mstate_subs_tuple(
                const std::deque<std::pair<std::set<unsigned>, int> > &greedy_tuple,
                bool active,
                bool breakpoint
        ) : greedy_tuple_(greedy_tuple),
            active_(active),
            breakpoint_(breakpoint){}

        virtual std::string to_string() const override;

        virtual bool is_active() const override { return this->active_; }

        virtual bool eq(const mstate &rhs) const override;

        virtual bool lt(const mstate &rhs) const override;

        virtual ~mstate_subs_tuple() override {}

        virtual const std::set<unsigned> &get_breakpoint() const override { return empty_breakpoint; }

        virtual void set_breakpoint(const std::set<unsigned> &breakpoint) override {
            (void)breakpoint;
            /*this->breakpoint_ = get_set_intersection(breakpoint, this->check_)*/;
        }

        bool is_upper_part() const;

        friend class kofola::complement_subs_tuple;
    }; // mstate_subs_tuple }}}

    bool is_acc_transition(
            const spot::const_twa_graph_ptr&    aut,
            const spot::scc_info&               scc_info,
            unsigned                     src,
            unsigned                     dst,
            const bdd&                          symbol) { // {{{
        for (const auto &t: aut->out(src)) {
            if (scc_info.scc_of(src) == scc_info.scc_of(t.dst) && dst == t.dst && bdd_implies(symbol, t.cond)) {
                if (t.acc) { return true; }
            }
        }

        return false;
    }
}

std::string mstate_subs_tuple::to_string() const
{
    std::string res = std::string("[subs_tuple(") + ((this->active_)? "A" : "T") + "): ";

    // Create a stringstream object
    std::stringstream ss;

    // Write the input deque to the stringstream
    ss << "greedy_tuple = (";
    for (const auto& p : this->greedy_tuple_) {
        ss << "(";
        ss << "{";
        for (const auto& elem : p.first) {
            ss << elem << ",";
        }
        ss << "},";
        ss << p.second << "),";
    }
    ss << ")";

    // Get the string representation of the stringstream
    std::string str = ss.str();
    res += "C=" + str;

    if (this->active_) {
        res += ", B=" + std::to_string(this->breakpoint_);
    }
    res += "]";
    return res;
}

bool mstate_subs_tuple::is_upper_part() const
{
    std::deque<std::pair<std::set<unsigned>, int> > tuple = greedy_tuple_;

    for(const auto& component : tuple) {
        int color = component.second;
        if(color != -1 && color != -2) {
            return false;
        }
    }

    return true;
}

bool mstate_subs_tuple::eq(const mstate& rhs) const
{
    const mstate_subs_tuple* rhs_subs_tuple = dynamic_cast<const mstate_subs_tuple*>(&rhs);
    assert(rhs_subs_tuple);
    return (this->active_ == rhs_subs_tuple->active_) &&
           (this->greedy_tuple_ == rhs_subs_tuple->greedy_tuple_) &&
            (this->breakpoint_ == rhs_subs_tuple->breakpoint_);
}

bool mstate_subs_tuple::lt(const mstate& rhs) const
{ // {{{
    const mstate_subs_tuple* rhs_subs_tuple = dynamic_cast<const mstate_subs_tuple*>(&rhs);
    assert(rhs_subs_tuple);

    if (this->active_ != rhs_subs_tuple->active_) { return this->active_ < rhs_subs_tuple->active_; }
    if (this->greedy_tuple_ != rhs_subs_tuple->greedy_tuple_) { return this->greedy_tuple_ < rhs_subs_tuple->greedy_tuple_; }
    if(this->breakpoint_ != rhs_subs_tuple->breakpoint_) { return this->breakpoint_ < rhs_subs_tuple->breakpoint_; }
    return false;   // if all are equal
} // lt() }}}


complement_subs_tuple::complement_subs_tuple(const cmpl_info& info, unsigned part_index)
    : abstract_complement_alg(info, part_index)
{ }

mstate_set complement_subs_tuple::get_init()
{ // {{{
    DEBUG_PRINT_LN("init subs_tuple for partition " + std::to_string(this->part_index_));
    std::deque<std::pair<std::set<unsigned>, int> > init_state;

    unsigned orig_init = this->info_.aut_->get_init_state_number();
    if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
        init_state.emplace_back(std::set<unsigned>{orig_init}, -1);
        //get_all_successors_in_scc(this->info_.aut_,);
        //this->info_.aut_->edge_storage() state_is_accepting()
    }

    std::shared_ptr<mstate> ms(new mstate_subs_tuple(init_state, false, true));
    mstate_set result = {ms};

    return result;
} // get_init() }}}

bool is_successor_of(unsigned state, unsigned succ, const cmpl_info& info)
{
    return kofola::is_in(succ, info.reachable_vector_[state]);
}

mstate_col_set complement_subs_tuple::upper_succ(
        const std::set<unsigned>&  glob_reached,
        const mstate*              src,
        const bdd&                 symbol,
        bool track)
{
    const mstate_subs_tuple* src_subs_tuple = dynamic_cast<const mstate_subs_tuple*>(src);

    std::set<unsigned> visited;
    std::deque<std::pair<std::set<unsigned>, int> > succ_tuple;
    std::deque<std::pair<std::set<unsigned>, int> > lower_part;
    bool breakp = true;
    int last_color = -10; // no color

    for(auto component = src_subs_tuple->greedy_tuple_.rbegin(); component != src_subs_tuple->greedy_tuple_.rend(); component++) {
        int color = component->second;

        std::set<unsigned> states = component->first;

        std::set<unsigned> acc_states;
        std::set<unsigned> non_acc_states;

        for(auto state : states) {
            std::set<unsigned> state_succ = kofola::get_all_successors_in_scc(
                    this->info_.aut_, this->info_.scc_info_, std::set<unsigned>{state}, symbol);

            for(auto succ : state_succ) {
                if(!visited.count(succ)) {
                    if(is_acc_transition(
                            this->info_.aut_,
                            this->info_.scc_info_,
                            state, succ, symbol)) {
                        acc_states.insert(succ);
                        non_acc_states.erase(succ);
                    }
                    else {
                        non_acc_states.insert(succ);
                    }
                }
            }
        }
        // has to be done afterwards since we may encounter the non-acc first and to acc after that
        std::set_union(acc_states.begin(), acc_states.end(),
                       non_acc_states.begin(), non_acc_states.end(),
                       std::inserter(visited, visited.begin()));

        if(!acc_states.empty()) {
            succ_tuple.emplace_front(acc_states, -2);
            last_color = -2;

            if(!lower_part.empty() && (lower_part.front().second == 1)) // merging
                lower_part[0].first.insert(acc_states.begin(), acc_states.end());
            else
                lower_part.emplace_front(acc_states, 1);
        }
        if(!non_acc_states.empty()) {
            succ_tuple.emplace_front(non_acc_states, color);
            last_color = color;

            if(color == -2)
                lower_part.emplace_front(non_acc_states, -3);
            else
                lower_part.emplace_front(non_acc_states, 0);
        }
    }
    // fix: only states from our SCC can be there, now takes all
    std::set<unsigned> new_runs;
    std::set_difference(glob_reached.begin(), glob_reached.end(), visited.begin(), visited.end(),
                        std::inserter(new_runs, new_runs.end()));

    auto it = new_runs.begin();
    while (it != new_runs.end()) {
        if (this->info_.st_to_part_map_.at(*it) != static_cast<int>(this->part_index_)) {
            it = new_runs.erase(it);
        } else {
            ++it;
        }
    }

    if(!new_runs.empty()) {
        if(last_color != -1) { // also when there is -10 -> nothing new, so also in this case
            succ_tuple.emplace_front(new_runs, -1);
            lower_part.emplace_front(new_runs, 0);
        }
        else{
            succ_tuple[0].first.insert(new_runs.begin(), new_runs.end()); // union of the first set and the incoming runs
            lower_part[0].first.insert(new_runs.begin(), new_runs.end());
        }
    }

    // we have to set one subset to 3 (each is 1 or 0/-3 now)
    for(auto& pair: lower_part) {
        if(pair.second == 1) {
            pair.second = 3;
            break;
        }
    }

    // only if none of them is empty
    mstate_col_set result;
    if(track)
    {
        std::shared_ptr<mstate> ms(new mstate_subs_tuple(succ_tuple, false, true));
        result = {{ms, {}}};
        return result;
    }
    else
    {
        std::shared_ptr<mstate> ms(new mstate_subs_tuple(succ_tuple, true, true));
        std::shared_ptr<mstate> ms_lower(new mstate_subs_tuple(lower_part, true, breakp));

        if(breakp)
            result = {{ms, {}}, {ms_lower, {0}}};
        else
            result = {{ms, {}}, {ms_lower, {}}};

        return result;
    }
}

mstate_col_set complement_subs_tuple::lower_succ(
        const std::set<unsigned>&  glob_reached,
        const mstate*              src,
        const bdd&                 symbol,
        bool track)
{
    const mstate_subs_tuple* src_subs_tuple = dynamic_cast<const mstate_subs_tuple*>(src);

    std::set<unsigned> visited;
    std::deque<std::pair<std::set<unsigned>, int> > succ_tuple;
    int last_color = -10; // no color
    bool breakp = true;
    bool discont_2 = false;
    bool found_right = false;

    for(auto component = src_subs_tuple->greedy_tuple_.rbegin(); component != src_subs_tuple->greedy_tuple_.rend(); component++) {
        int color = component->second;

        std::set<unsigned> states = component->first;

        std::set<unsigned> acc_states;
        std::set<unsigned> non_acc_states;

        for(auto state : states) {
            std::set<unsigned> state_succ = kofola::get_all_successors_in_scc(
                    this->info_.aut_, this->info_.scc_info_, std::set<unsigned>{state}, symbol);

            for(auto succ : state_succ) {
                if(!visited.count(succ)) {
                    if(is_acc_transition(
                            this->info_.aut_,
                            this->info_.scc_info_,
                            state, succ, symbol)) {
                        acc_states.insert(succ);
                        non_acc_states.erase(succ);
                    }
                    else {
                        non_acc_states.insert(succ);
                    }
                }
            }
        }
        // has to be done afterwards since we may encounter the non-acc first and to acc after that
        std::set_union(acc_states.begin(), acc_states.end(),
                       non_acc_states.begin(), non_acc_states.end(),
                       std::inserter(visited, visited.begin()));

        if(!acc_states.empty()) {
            if((color == 0 || color == -3)) {
                if (!succ_tuple.empty() && (succ_tuple.front().second == 1)) {
                    succ_tuple[0].first.insert(acc_states.begin(), acc_states.end());
                    succ_tuple[0].second = 1;
                }
                else
                    succ_tuple.emplace_front(acc_states, 1);
            }
            else if(color == 1) {
                if(!succ_tuple.empty() && (succ_tuple.front().second == 1)) {
                    succ_tuple[0].first.insert(acc_states.begin(), acc_states.end());
                }
                else
                    succ_tuple.emplace_front(acc_states, 1);
            }
            else if(color == 2) {
                if(!succ_tuple.empty() && (succ_tuple.front().second == 2 || succ_tuple.front().second == 1)) {
                    succ_tuple[0].first.insert(acc_states.begin(), acc_states.end());
                    succ_tuple[0].second = 2;
                }
                else
                    succ_tuple.emplace_front(acc_states, 2);
            }
            else if(color == 3) {
                if (!succ_tuple.empty() && (succ_tuple.front().second == 2 || succ_tuple.front().second == 1)) {
                    succ_tuple[0].first.insert(acc_states.begin(), acc_states.end());
                    succ_tuple[0].second = 2;
                }
                else
                    succ_tuple.emplace_front(acc_states, 2);
            }
            else {
                succ_tuple.emplace_front(acc_states, color); // colors 0 and -3
            }

            last_color = succ_tuple[0].second;
            if(last_color == 2)
                breakp = false;
        }

        if(!non_acc_states.empty()) {
            if(color == 3) {
                if(!succ_tuple.empty() && (succ_tuple.front().second == 2 || succ_tuple.front().second == 1)) {
                    succ_tuple[0].first.insert(non_acc_states.begin(), non_acc_states.end());
                    succ_tuple[0].second = 2;
                }
                else
                    succ_tuple.emplace_front(non_acc_states, 2);
            }
            else if(color == 1) {
                if(!succ_tuple.empty() && (succ_tuple.front().second == 1)) {
                    succ_tuple[0].first.insert(non_acc_states.begin(), non_acc_states.end());
                }
                else
                    succ_tuple.emplace_front(non_acc_states, 1);
            }
            else if(color == 2) {
                if(!succ_tuple.empty() && (succ_tuple.front().second == 2 || succ_tuple.front().second == 1)) {
                    succ_tuple[0].first.insert(non_acc_states.begin(), non_acc_states.end());
                    succ_tuple[0].second = 2;
                }
                else
                    succ_tuple.emplace_front(non_acc_states, 2);
            }
            else
                succ_tuple.emplace_front(non_acc_states, color);

            last_color = succ_tuple[0].second;
            if(last_color == 2)
                breakp = false;
        }

        if(acc_states.empty() && non_acc_states.empty() && (color == 2 || color == 3)) {
            discont_2 = true;
            for(auto& pair: succ_tuple) {
                if(pair.second == 0) {
                    pair.second = 10; // to mark it as the needed 0
                    found_right = true;
                    break;
                }
                if(pair.second == -3) {
                    pair.second = 11; // to mark it as the needed 0
                    found_right = true;
                    break;
                }
            }
        }
    }

    std::set<unsigned> new_runs;
    std::set_difference(glob_reached.begin(), glob_reached.end(), visited.begin(), visited.end(),
                        std::inserter(new_runs, new_runs.end()));

    auto it = new_runs.begin();
    while (it != new_runs.end()) {
        if (this->info_.st_to_part_map_.at(*it) != static_cast<int>(this->part_index_)) {
            it = new_runs.erase(it);
        } else {
            ++it;
        }
    }

    if(!new_runs.empty()) {
        if(last_color == -3 || last_color == -10 || last_color == 1 || last_color == 2 || last_color == 11) { // also when there is -10 -> nothing new, so also in this case
            succ_tuple.emplace_front(new_runs, 0);
        }
        else{
            succ_tuple[0].first.insert(new_runs.begin(), new_runs.end()); // union of the first set and the incoming runs
        }
    }

    // setting color 3
    if(breakp && discont_2) {
        bool was_zero = false;
        bool set_3 = false;

        if(found_right) {
            for (auto &pair: succ_tuple) {
                if (pair.second == 10) {
                    pair.second = 0; // to mark it as the needed 0
                    was_zero = true;
                } else if (pair.second == 11) {
                    pair.second = -3; // to mark it as the needed 0
                    was_zero = true;
                }
                if (was_zero && pair.second == 1) {
                    pair.second = 3;
                    set_3 = true;
                    break;
                }
            }
        }
        if(!set_3) {
            for(auto& pair: succ_tuple) {
                if(pair.second == 1) {
                    pair.second = 3;
                    set_3 = true;
                    break;
                }
            }
        }
    }
    else if(breakp) {
        for(auto& pair: succ_tuple) {
            if(pair.second == 1) {
                pair.second = 3;
                break;
            }
        }
    }
    else {
        for(auto& pair: succ_tuple) {
            if(pair.second == 11) {
                pair.second = -3;
            }
            if(pair.second == 10) {
                pair.second = 0;
            }
        }
    }
    // end of setting color 3

    std::shared_ptr<mstate> ms;
    if(track) {
        ms = std::shared_ptr<mstate>(new mstate_subs_tuple(succ_tuple, false, breakp));
    }
    else {
        ms = std::shared_ptr<mstate>(new mstate_subs_tuple(succ_tuple, true, breakp));
    }

    if(breakp)
    { mstate_col_set result = {{ms, {0}}}; return result; }
    else
    { mstate_col_set result = {{ms, {}}}; return result; }
}

mstate_col_set complement_subs_tuple::get_succ_track(
        const std::set<unsigned>&  glob_reached,
        const mstate*              src,
        const bdd&                 symbol)
{
    const mstate_subs_tuple* src_subs_tuple = dynamic_cast<const mstate_subs_tuple*>(src);
    assert(src_subs_tuple);
    assert(!src_subs_tuple->active_);

    if(src_subs_tuple->is_upper_part()) {
        return upper_succ(glob_reached, src, symbol, true);
    }
    else {
        return lower_succ(glob_reached, src, symbol, true);
    }

} // get_succ_track() }}}

mstate_set complement_subs_tuple::lift_track_to_active(const mstate* src)
{ // {{{
    const mstate_subs_tuple* src_acitve = dynamic_cast<const mstate_subs_tuple*>(src);
    assert(src_acitve);
    assert(!src_acitve->active_);

    std::shared_ptr<mstate> ms(new mstate_subs_tuple(src_acitve->greedy_tuple_, true, src_acitve->breakpoint_));
    return {ms};
} // lift_track_to_active() }}}

mstate_col_set complement_subs_tuple::get_succ_active(
        const std::set<unsigned>&  glob_reached,
        const mstate*              src,
        const bdd&                 symbol,
        bool resample)
{
    (void)resample;
    const mstate_subs_tuple* src_subs_tuple = dynamic_cast<const mstate_subs_tuple*>(src);
    assert(src_subs_tuple);
    assert(src_subs_tuple->active_);

    if(src_subs_tuple->is_upper_part()) {
        return upper_succ(glob_reached, src, symbol, false);
    }
    else {
        return lower_succ(glob_reached, src, symbol, false);
    }
}

complement_subs_tuple::~complement_subs_tuple()
{
}
