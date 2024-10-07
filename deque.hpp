#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

template <typename T, typename Allocator = std::allocator<T>>
class Deque {
  const size_t kChunckCap = 10;

 private:
  using alloc_traits = std::allocator_traits<Allocator>;
  using chunck_allocator = typename alloc_traits::template rebind_alloc<T*>;
  using chunck_alloc_traits = std::allocator_traits<chunck_allocator>;
  Allocator allocator_;
  std::vector<T*, chunck_allocator> chuncks_;

  size_t start_chunck_ = 0;
  size_t end_chunck_ = 0;
  size_t start_index_ = 0;
  size_t end_index_ = 0;
  size_t size_ = 0;

  ////////////Allocate Chunck///////////////////////
  T* allocate_chunck() {
    return alloc_traits::allocate(allocator_, kChunckCap);
  }
  ////////////Deallocate Chunck//////////////////////
  void deallocate_chunck(T* chunck) {
    alloc_traits::deallocate(allocator_, chunck, kChunckCap);
  }
  ////////////Construct Element//////////////////////
  void construct_def_element(T* chunck, size_t& count) {
    try {
      for (size_t i = 0; i < kChunckCap; ++i) {
        if (count > 0) {
          try {
            --count;
            alloc_traits::construct(allocator_, chunck + i);
          } catch (...) {
            del_all();
            throw;
          }
          ++size_;
        }
      }
    } catch (...) {
      throw;
    }
  }
  void construct_element(T* elem, const T& value) {
    try {
      alloc_traits::construct(allocator_, elem, value);
    } catch (...) {
      del_all();
      throw;
    }
  }

  void construct_element(T* elem, T&& value) {
    try {
      alloc_traits::construct(allocator_, elem, std::move(value));
    } catch (...) {
      // del_all();
      throw;
    }
  }

  template <class... Args>
  void construct_element(T* elem, Args&&... args) {
    try {
      alloc_traits::construct(allocator_, elem, std::forward<Args>(args)...);
    } catch (...) {
      // del_all();
      throw;
    }
  }

 public:
  /////////////////////////////////////////////
  //////////////CONSTRUCTORS//////////////////
  ///////////////////////////////////////////

  explicit Deque(const Allocator& alloc = Allocator())
      : allocator_(alloc), chuncks_(chunck_allocator(alloc)) {}
  //////////////////////////////
  Deque(Deque&& other) noexcept
      : allocator_(std::move(other.allocator_)),
        chuncks_(std::move(other.chuncks_), chunck_allocator(allocator_)),
        kChunckCap(other.kChunckCap),
        start_chunck_(other.start_chunck_),
        end_chunck_(other.end_chunck_),
        start_index_(other.start_index_),
        end_index_(other.end_index_),
        size_(other.size_) {
    other.start_chunck_ = 0;
    other.end_chunck_ = 0;
    other.start_index_ = 0;
    other.end_index_ = 0;
    other.size_ = 0;
    other.chuncks_.clear();
  }
  //////////////////////////////
  Deque(const Deque& other) try
      : allocator_(std::allocator_traits<Allocator>::
                       select_on_container_copy_construction(other.allocator_)),
        chuncks_(other.chuncks_.size(), nullptr,
                 chunck_allocator(other.allocator_)),
        start_chunck_(other.start_chunck_),
        start_index_(other.start_index_),
        end_chunck_(other.end_chunck_),
        end_index_(other.end_index_) {
    for (size_t i = 0; i < chuncks_.size(); ++i) {
      chuncks_[i] = alloc_traits::allocate(allocator_, kChunckCap);
      for (size_t j = 0; j < kChunckCap; ++j) {
        if ((i == start_chunck_ && j >= start_index_ || i > start_chunck_) &&
            (i == end_chunck_ && j < end_index_ || i < end_chunck_)) {
          try {
            construct_element(chuncks_[i] + j, other.chuncks_[i][j]);
            ++size_;
          } catch (...) {
            del_all();
            throw;
          }
        }
      }
    }
  } catch (...) {
    throw;
  }
  //////////////////////////////
  Deque(size_t count, const Allocator& alloc = Allocator()) try
      : allocator_(alloc),
        chuncks_(count / kChunckCap + 1, nullptr,
                 chunck_allocator(allocator_)) {
    end_chunck_ = count / kChunckCap;
    end_index_ = count % kChunckCap;
    for (auto& chunck : chuncks_) {
      chunck = alloc_traits::allocate(allocator_, kChunckCap);
      construct_def_element(chunck, count);
    }
  } catch (...) {
    throw;
  }

  //////////////////////////////
  Deque(size_t count, const T& value, Allocator alloc = Allocator()) try
      : allocator_(alloc),
        chuncks_(count / kChunckCap + 1, nullptr, chunck_allocator(alloc)) {
    end_chunck_ = count / kChunckCap;
    end_index_ = count % kChunckCap;
    try {
      for (auto& chunck : chuncks_) {
        chunck = alloc_traits::allocate(allocator_, kChunckCap);
        for (size_t j = 0; j < kChunckCap; ++j) {
          if (count != 0) {
            try {
              construct_element(chunck + j, value);
            } catch (...) {
              del_all();
              throw;
            }

            --count;
            ++size_;
          }
        }
      }
    } catch (...) {
      throw 1;
    }
  } catch (...) {
    throw 1;
  }

  Deque(std::initializer_list<T> init, const Allocator& alloc = Allocator()) try
      : allocator_(alloc), chuncks_(chunck_allocator(alloc)) {
    auto iter = init.begin();
    while (iter != init.end()) {
      push_back(*(iter++));
    }
  } catch (...) {
    throw;
  }
  /////////////////////////////////////////////////
  //////////////////DESTRUCTOR/////////////////////
  /////////////////////////////////////////////////
  ~Deque() { del_all(); }

  /////////////////////////////////////////////////
  ///////////////////OPERATORS////////////////////
  ///////////////////////////////////////////////
  Deque& operator=(const Deque& other) {
    if (std::allocator_traits<
            Allocator>::propagate_on_container_copy_assignment::value) {
      allocator_ = other.allocator_;
    }
    Deque copy(other);
    std::swap(chuncks_, copy.chuncks_);
    std::swap(start_chunck_, copy.start_chunck_);
    std::swap(end_chunck_, copy.end_chunck_);
    std::swap(start_index_, copy.start_index_);
    std::swap(end_index_, copy.end_index_);
    std::swap(size_, copy.size_);
    return *this;
  }

  Deque& operator=(Deque&& other) noexcept {
    if (this != &other) {
      chuncks_ = std::move(other.chuncks_);
      start_chunck_ = other.start_chunck_;
      end_chunck_ = other.end_chunck_;
      start_index_ = other.start_index_;
      end_index_ = other.end_index_;
      size_ = other.size_;
      other.chuncks_.clear();
      other.start_chunck_ = 0;
      other.end_chunck_ = 0;
      other.start_index_ = 0;
      other.end_index_ = 0;
      other.size_ = 0;
      if (std::allocator_traits<
              Allocator>::propagate_on_container_copy_assignment::value) {
        allocator_ = std::move(other.allocator_);
      }
    }
    return *this;
  }
  /////////////////////////////////////////////
  T& operator[](int index) {
    if (index < 0) {
      index += size_;
    }
    return chuncks_[start_chunck_ + ((index + start_index_) / kChunckCap)]
                   [(index + start_index_) % kChunckCap];
  }
  const T& operator[](int index) const {
    if (index < 0) {
      index += size_;
    }
    return chuncks_[start_chunck_ + ((index + start_index_) / kChunckCap)]
                   [(index + start_index_) % kChunckCap];
  }

  //////////////////////////////////////////////
  ///////////////PUSHES////////////////////////
  void push_back(const T& val) {
    if (end_chunck_ == chuncks_.size()) {
      resize();
    }
    try {
      construct_element(chuncks_[end_chunck_] + end_index_, val);
    } catch (...) {
      throw;
    }
    size_++;
    end_chunck_ = start_chunck_ + (((size_) + start_index_) / kChunckCap);
    end_index_ = (size_ + start_index_) % kChunckCap;
  }

  void push_back(T&& val) {
    if (end_chunck_ == chuncks_.size()) {
      resize();
    }
    try {
      construct_element(chuncks_[end_chunck_] + end_index_, std::move(val));
    } catch (...) {
      throw;
    }

    size_++;
    end_chunck_ = start_chunck_ + (((size_) + start_index_) / kChunckCap);
    end_index_ = (size_ + start_index_) % kChunckCap;
  }

  void push_front(const T& val) {
    if (start_chunck_ == 0 && start_index_ == 0) {
      resize();
    }
    if (start_index_ == 0) {
      --start_chunck_;
      start_index_ = kChunckCap - 1;
    } else {
      --start_index_;
    }
    construct_element(chuncks_[start_chunck_] + start_index_, val);
    size_++;
  }

  void push_front(T&& val) {
    if (start_chunck_ == 0 && start_index_ == 0) {
      resize();
    }
    if (start_index_ == 0) {
      --start_chunck_;
      start_index_ = kChunckCap - 1;
    } else {
      --start_index_;
    }
    construct_element(chuncks_[start_chunck_] + start_index_, std::move(val));
    size_++;
  }
  /////////////////////////////////////////////
  ///////////////POPES////////////////////////
  void pop_back() {
    if (empty()) {
      return;
    }
    if (end_index_ > 0) {
      --end_index_;
    } else {
      --end_chunck_;
      end_index_ = kChunckCap - 1;
    }
    alloc_traits::destroy(allocator_, chuncks_[end_chunck_] + end_index_);
    --size_;
  }

  void pop_front() {
    if (empty()) {
      return;
    }
    alloc_traits::destroy(allocator_, chuncks_[start_chunck_] + start_index_);

    if (start_index_ < kChunckCap - 1) {
      ++start_index_;
    } else {
      ++start_chunck_;
      start_index_ = 0;
    }
    --size_;
  }

  //////////////////////////////////////////////
  //////////////////METHODS/////////////////////
  /////////////////////////////////////////////
  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  /////////////////////////////////////////////
  //////////////////ATS////////////////////////
  T& at(size_t index) {
    if (index >= size_) {
      throw std::out_of_range("ERROR: OUT OF RANGE");
    }
    return chuncks_[start_chunck_ + ((index + start_index_) / kChunckCap)]
                   [(index + start_index_) % kChunckCap];
  }

  const T& at(size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("ERROR: OUT OF RANGE");
    }
    return chuncks_[start_chunck_ + ((index + start_index_) / kChunckCap)]
                   [(index + start_index_) % kChunckCap];
  }
  //////////////////////////////////////////////

  /////////////////////////////////////////////
 private:
  ////////////METHODS////////////////////
  void del_all() {
    while (!empty()) {
      pop_back();
    }
    for (auto chunck : chuncks_) {
      deallocate_chunck(chunck);
    }
    chuncks_.clear();
  }

  ////////////////////////////////////////
  void resize() {
    if (empty()) {
      chuncks_.resize(1, allocate_chunck());
      return;
    }
    size_t old_size = chuncks_.size();
    chuncks_.resize(old_size * 3, nullptr);
    for (size_t i = 0; i < old_size; ++i) {
      chuncks_[i + old_size] = chuncks_[i];
      chuncks_[i] = alloc_traits::allocate(allocator_, kChunckCap);
      chuncks_[i + 2 * old_size] =
          alloc_traits::allocate(allocator_, kChunckCap);
    }
    start_chunck_ += old_size;
    end_chunck_ += old_size;
  }
  ///////////////////////////////////////

 public:
  /////////////////////////////////////////
  ////////////////////////////////////////
  ////////////ITERATOR///////////////////
  ///////////////////////////////////////
  ///////////////////////////////////////
  template <bool IsConst>
  class DequeIterator
      : public std::iterator<
            std::random_access_iterator_tag,
            typename std::conditional<IsConst, const T, T>::type> {
   private:
    Deque<T, Allocator>* deque_ = nullptr;
    size_t index_ = 0;

   public:
    typedef typename std::conditional<IsConst, const T, T>::type Ttype;
    using value_type = typename std::iterator<std::random_access_iterator_tag,
                                              Ttype>::value_type;
    using pointer =
        typename std::iterator<std::random_access_iterator_tag, Ttype>::pointer;
    using difference_type =
        typename std::iterator<std::random_access_iterator_tag,
                               Ttype>::difference_type;
    using reference = typename std::iterator<std::random_access_iterator_tag,
                                             Ttype>::reference;
    /////DequeIterator constructors/////////////

    DequeIterator() = default;

    DequeIterator(Deque<T, Allocator>* deque, int index)
        : deque_(deque), index_(index) {}

    DequeIterator(const DequeIterator<IsConst>& other)
        : deque_(other.deque_), index_(other.index_) {}

    ///////DequeIterator operators////////////////
    void operator=(const DequeIterator& other) {
      index_ = other.index_;
      deque_ = other.deque_;
    }
    /////////////////инкремент/////////////////
    DequeIterator<IsConst>& operator++() {
      ++index_;
      return *this;
    }

    DequeIterator<IsConst> operator++(int) {
      DequeIterator<IsConst> temp(*this);
      ++index_;
      return temp;
    }
    ////////////////декремент///////////////
    DequeIterator<IsConst>& operator--() {
      --index_;
      return *this;
    }

    DequeIterator<IsConst> operator--(int) {
      DequeIterator<IsConst> temp(*this);
      --index_;
      return temp;
    }
    /////////////////////////////////////////
    DequeIterator<IsConst>& operator+=(difference_type diff) {
      index_ += diff;
      return *this;
    }

    DequeIterator<IsConst>& operator-=(difference_type diff) {
      index_ -= diff;
      return *this;
    }
    //////////////////////////////////////////

    DequeIterator<IsConst> operator+(difference_type diff) const {
      DequeIterator<IsConst> temp(*this);
      temp.index_ += diff;
      return temp;
    }

    DequeIterator<IsConst> operator-(difference_type diff) const {
      DequeIterator<IsConst> temp(*this);
      temp.index_ -= diff;
      return temp;
    }
    ////////////////////////////////////////////

    difference_type operator-(const DequeIterator<IsConst>& iter) const {
      return index_ - iter.index_;
    }

    ////////////////////////////////////////////

    reference operator*() const { return (*deque_)[index_]; }
    pointer operator->() const { return &((*deque_)[index_]); }

    ////////////////////////////////////////////

    friend bool operator==(const DequeIterator<IsConst>& iter1,
                           const DequeIterator<IsConst>& iter2) {
      return iter1.deque_ == iter2.deque_ && iter1.index_ == iter2.index_;
    }

    friend bool operator>(const DequeIterator<IsConst>& iter1,
                          const DequeIterator<IsConst>& iter2) {
      return iter1.index_ > iter2.index_;
    }

    friend bool operator!=(const DequeIterator<IsConst>& iter1,
                           const DequeIterator<IsConst>& iter2) {
      return !(iter1 == iter2);
    }

    friend bool operator<(const DequeIterator<IsConst>& iter1,
                          const DequeIterator<IsConst>& iter2) {
      return (iter1 != iter2 && !(iter1 > iter2));
    }

    friend bool operator<=(const DequeIterator<IsConst>& iter1,
                           const DequeIterator<IsConst>& iter2) {
      return !(iter1 > iter2);
    }

    friend bool operator>=(const DequeIterator<IsConst>& iter1,
                           const DequeIterator<IsConst>& iter2) {
      return !(iter1 < iter2);
    }
  };

  ///////////ITERATOR USINGS////////////////////////////////////////////
  using iterator = DequeIterator<false>;
  using const_iterator = DequeIterator<true>;
  using reverse_iterator = std::reverse_iterator<DequeIterator<false>>;
  using const_reverse_iterator = std::reverse_iterator<DequeIterator<true>>;

  ///////////////////////////////////////////////////////
  ///////////////////BEGINS & ENDS///////////////////////

  iterator begin() { return iterator(this, 0); }

  iterator end() { return iterator(this, size_); }

  const_iterator cbegin() { return const_iterator(this, 0); }

  const_iterator cend() { return const_iterator(this, size_); }

  reverse_iterator rbegin() noexcept {
    return reverse_iterator(iterator(this, size_));
  }

  reverse_iterator rend() noexcept {
    return reverse_iterator(iterator(this, 0));
  }

  const_reverse_iterator crbegin() const noexcept {
    return const_reverse_iterator(const_iterator(this, size_));
  }
  const_reverse_iterator crend() const noexcept {
    return const_reverse_iterator(const_iterator(this, 0));
  }

  ///////////////////////////////////////////////////////
  ///////////////////INSERT & ERASE//////////////////////

  void insert(iterator place, const T& val) {
    size_t index = place - begin();
    push_back(val);
    for (size_t i = size() - 1; i > index; --i) {
      std::swap((*this)[i], (*this)[i - 1]);
    }
  }

  void erase(iterator place) {
    size_t index = place - begin();
    for (size_t i = index; i < size() - 1; ++i) {
      std::swap((*this)[i], (*this)[i + 1]);
    }
    pop_back();
  }
  template <typename... Args>
  void emplace_back(Args&&... args) {
    if (end_chunck_ == chuncks_.size()) {
      resize();
    }
    construct_element(chuncks_[end_chunck_] + end_index_,
                      std::forward<Args>(args)...);
    if (end_index_ == kChunckCap) {
      ++end_chunck_;
      end_index_ = 0;
    } else {
      ++end_index_;
    }

    ++size_;
  }

  template <typename... Args>
  void emplace_front(Args&&... args) {
    if (start_index_ == 0) {
      chuncks_.insert(chuncks_.begin(), allocate_chunck());
      --start_chunck_;
      start_index_ = kChunckCap;
    }
    --start_index_;
    construct_element(chuncks_[start_chunck_] + start_index_,
                      std::forward<Args>(args)...);
    ++size_;
  }

  template <typename... Args>
  void emplace(iterator pos, Args&&... args) {
    if (pos == begin()) {
      emplace_front(std::forward<Args>(args)...);
      return;
    }
    if (pos == end()) {
      emplace_back(std::forward<Args>(args)...);
      return;
    }
    size_t index = pos - begin();
    chuncks_.emplace_back(std::forward<Args>(args)...);
    for (size_t i = index; i < size() - 1; ++i) {
      std::swap((*this)[i], (*this)[i + 1]);
    }
  }
  Allocator get_allocator() { return allocator_; }
};
