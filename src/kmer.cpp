#include "nthash/kmer_functions.hpp"

namespace nthash {

NtHash::NtHash(const char* seq,
               size_t seq_len,
               typedefs::NUM_HASHES_TYPE num_hashes,
               typedefs::K_TYPE k,
               size_t pos)
  : seq(seq, seq_len)
  , num_hashes(num_hashes)
  , k(k)
  , pos(pos)
  , initialized(false)
  , hash_arr(new uint64_t[num_hashes])
{
  if (k == 0) {
    raise_error("NtHash", "k must be greater than 0");
  }
  if (this->seq.size() < k) {
    raise_error("NtHash",
                "sequence length (" + std::to_string(this->seq.size()) +
                  ") is smaller than k (" + std::to_string(k) + ")");
  }
  if (pos > this->seq.size() - k) {
    raise_error("NtHash",
                "passed position (" + std::to_string(pos) +
                  ") is larger than sequence length (" +
                  std::to_string(this->seq.size()) + ")");
  }
}

bool
NtHash::init()
{
  size_t pos_n = 0;
  while (pos <= seq.size() - k + 1 &&
         is_invalid_kmer(seq.data() + pos, k, pos_n)) {
    pos += pos_n + 1;
  }
  if (pos > seq.size() - k) {
    return false;
  }
  fwd_hash = base_forward_hash(seq.data() + pos, k);
  rev_hash = base_reverse_hash(seq.data() + pos, k);
  extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  initialized = true;
  return true;
}

bool
NtHash::roll()
{
  if (!initialized) {
    return init();
  }
  if (pos >= seq.size() - k) {
    return false;
  }
  if (SEED_TAB[(unsigned char)seq[pos + k]] == SEED_N) {
    pos += k;
    return init();
  }
  fwd_hash = next_forward_hash(fwd_hash, k, seq[pos], seq[pos + k]);
  rev_hash = next_reverse_hash(rev_hash, k, seq[pos], seq[pos + k]);
  extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  ++pos;
  return true;
}

bool
NtHash::roll_back()
{
  if (!initialized) {
    return init();
  }
  if (pos == 0) {
    return false;
  }
  if (SEED_TAB[(unsigned char)seq[pos - 1]] == SEED_N && pos >= k) {
    pos -= k;
    return init();
  }
  if (SEED_TAB[(unsigned char)seq[pos - 1]] == SEED_N) {
    return false;
  }
  fwd_hash = prev_forward_hash(fwd_hash, k, seq[pos + k - 1], seq[pos - 1]);
  rev_hash = prev_reverse_hash(rev_hash, k, seq[pos + k - 1], seq[pos - 1]);
  extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  --pos;
  return true;
}

bool
NtHash::peek()
{
  if (pos >= seq.size() - k) {
    return false;
  }
  return peek(seq[pos + k]);
}

bool
NtHash::peek(char char_in)
{
  if (!initialized) {
    return init();
  }
  if (SEED_TAB[(unsigned char)char_in] == SEED_N) {
    return false;
  }
  const uint64_t fwd = next_forward_hash(fwd_hash, k, seq[pos], char_in);
  const uint64_t rev = next_reverse_hash(rev_hash, k, seq[pos], char_in);
  extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
  return true;
}

bool
NtHash::peek_back()
{
  if (pos == 0) {
    return false;
  }
  return peek_back(seq[pos - 1]);
}

bool
NtHash::peek_back(char char_in)
{
  if (!initialized) {
    return init();
  }
  if (SEED_TAB[(unsigned char)char_in] == SEED_N) {
    return false;
  }
  const unsigned char char_out = seq[pos + k - 1];
  const uint64_t fwd = prev_forward_hash(fwd_hash, k, char_out, char_in);
  const uint64_t rev = prev_reverse_hash(rev_hash, k, char_out, char_in);
  extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
  return true;
}

BlindNtHash::BlindNtHash(const char* seq,
                         typedefs::NUM_HASHES_TYPE num_hashes,
                         typedefs::K_TYPE k,
                         ssize_t pos)
  : seq(seq + pos, seq + pos + k)
  , num_hashes(num_hashes)
  , pos(pos)
  , hash_arr(new uint64_t[num_hashes])
{
  if (k == 0) {
    raise_error("BlindNtHash", "k must be greater than 0");
  }
  fwd_hash = base_forward_hash(seq, k);
  rev_hash = base_reverse_hash(seq, k);
  extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
}

void
BlindNtHash::roll(char char_in)
{
  fwd_hash = next_forward_hash(fwd_hash, seq.size(), seq.front(), char_in);
  rev_hash = next_reverse_hash(rev_hash, seq.size(), seq.front(), char_in);
  extend_hashes(fwd_hash, rev_hash, seq.size(), num_hashes, hash_arr.get());
  seq.pop_front();
  seq.push_back(char_in);
  ++pos;
}

void
BlindNtHash::roll_back(char char_in)
{
  fwd_hash = prev_forward_hash(fwd_hash, seq.size(), seq.back(), char_in);
  rev_hash = prev_reverse_hash(rev_hash, seq.size(), seq.back(), char_in);
  extend_hashes(fwd_hash, rev_hash, seq.size(), num_hashes, hash_arr.get());
  seq.pop_back();
  seq.push_front(char_in);
  --pos;
}

void
BlindNtHash::peek(char char_in)
{
  const typedefs::K_TYPE k = seq.size();
  const uint64_t fwd = next_forward_hash(fwd_hash, k, seq.front(), char_in);
  const uint64_t rev = next_reverse_hash(rev_hash, k, seq.front(), char_in);
  extend_hashes(fwd, rev, seq.size(), num_hashes, hash_arr.get());
}

void
BlindNtHash::peek_back(char char_in)
{
  const typedefs::K_TYPE k = seq.size();
  const uint64_t fwd = prev_forward_hash(fwd_hash, k, seq.back(), char_in);
  const uint64_t rev = prev_reverse_hash(rev_hash, k, seq.back(), char_in);
  extend_hashes(fwd, rev, seq.size(), num_hashes, hash_arr.get());
}

} // namespace nthash