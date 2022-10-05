// implementation of Safra trees from Yong Li

#pragma once

namespace kofola { // {{{
namespace safra {

class mstate {
public:
  mstate() {}
  virtual ~mstate() {}
  virtual size_t hash() const { return 0; };
  virtual bool less_than(const mstate *other) const { return false; }
  virtual bool eq(const mstate *other) const { return false; }
  // bool operator != (const mstate& other) {
  //   return ! (*this == &other);
  // }
  virtual mstate *get_ptr() { return this; }
  virtual std::set<unsigned> get_states() const { return std::set<unsigned>(); }
  virtual std::string to_string() const { return ".."; }
};

typedef std::pair<unsigned, int> label;

class order_vec : public mstate {
public:
  // a map from a state to its label
  std::vector<label> labels_;

  order_vec() {}

  order_vec(const order_vec &other) { this->labels_ = other.labels_; }

  order_vec& operator=(const order_vec &other) {
    this->labels_ = other.labels_;
    return *this;
  }

  mstate *get_ptr() { return this; }

  bool operator<(const order_vec &other) const {
    return labels_ < other.labels_;
  }

  bool operator==(const order_vec &other) const {
    return labels_ == other.labels_;
  }

  bool less_than(const mstate *other) const {
    const order_vec *po = static_cast<const order_vec *>(other);
    return *this < *po;
  }
  bool eq(const mstate *other) const {
    const order_vec *po = static_cast<const order_vec *>(other);
    return *this == *po;
  }

  size_t hash() const {
    size_t res = 0;
    for (auto lab : labels_) {
      res = (res << 3) ^ lab.first;
      res = (res << 3) ^ lab.second;
    }
    return res;
  }

  std::set<unsigned> get_states() const {
    std::set<unsigned> res;
    for (const auto &lab : labels_) {
      res.insert(lab.first);
    }
    return res;
  }

  std::string to_string() const {
    std::string res = "C = {";
    res += "}";
    return "";
  }
};

class safra_tree : public order_vec {
public:
  // a map from a node to its parent
  std::vector<int> braces_;

  safra_tree() : order_vec() {}

  safra_tree(const safra_tree &other) : order_vec() {
    this->labels_ = other.labels_;
    this->braces_ = other.braces_;
  }

  bool operator<(const safra_tree &other) const {
    return labels_ < other.labels_ ? true : braces_ < other.braces_;
  }

  bool operator==(const safra_tree &other) const {
    return labels_ == other.labels_ && braces_ == other.braces_;
  }

  safra_tree& operator=(const safra_tree &other) {
    this->labels_ = other.labels_;
    this->braces_ = other.braces_;
    return *this;
  }

  bool less_than(const mstate *other) const {
    const safra_tree *po = static_cast<const safra_tree *>(other);
    return *this < *po;
  }
  bool eq(const mstate *other) const {
    const safra_tree *po = static_cast<const safra_tree *>(other);
    return *this == *po;
  }

  size_t hash() const {
    size_t res = 0;
    for (auto lab : labels_) {
      res = (res << 3) ^ lab.first;
      res = (res << 3) ^ lab.second;
    }
    for (int i : braces_) {
      res ^= (res << 3) ^ i;
    }
    return res;
  }

  std::string to_string() const {
    std::string res = "[";
    for (const auto& pair : labels_) {
      // first state
      unsigned s = pair.first;
      int brace = pair.second;
      // now compute the braces
      std::vector<int> tmp;
      std::string bstr = "";
      while (brace >= 0) {
        // insert in reverse order
        bstr = std::to_string(brace) + ((bstr.length() > 0)? "," : "") + bstr;
        // obtain the i-th braces
        brace = braces_[brace];
      }
      res += ((res.length() > 1)? ", " : "") + std::to_string(s) + "->{" + bstr + "}";
    }
    res += "]";
    return res;
  }
};

}} // namespace kofola::safra }}}
