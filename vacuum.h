#ifndef CUCKOO_H
#define CUCKOO_H

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include "hashutil.h"

#define memcle(a) memset(a, 0, sizeof(a))
#define sqr(a) ((a) * (a))
#define debug(a) cerr << #a << " = " << a << ' '
#define deln(a) cerr << #a << " = " << a << endl
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ROUNDDOWN(a, b) ((a) - ((a) % (b)))
#define ROUNDUP(a, b) ROUNDDOWN((a) + (b - 1), b)


template <typename fp_t, int fp_len>
class Filter {
public:
    long long n;  // number of buckets
    int m;        // number of slots per bucket
    uint64_t memory_consumption;
    virtual void init(int _n, int _m, int _max_kick_steps) = 0;
    virtual void clear() = 0;
    virtual bool insert(uint64_t ele) = 0;
    virtual bool lookup(uint64_t ele) = 0;
    virtual bool del(uint64_t ele) = 0;
    uint64_t position_hash(long long ele);  // hash to range [0, n - 1]
    virtual double get_load_factor() { return 0; }
    virtual double get_full_bucket_factor() { return 0; }
    virtual void debug_test() {}
};

template <typename fp_t, int fp_len>
uint64_t Filter<fp_t, fp_len>::position_hash(long long ele) {
    return (ele % n + n) % n;
}

template <typename fp_t, int fp_len>
class SemiSortCuckooFilter : public Filter<fp_t, fp_len> {
public:
    int max_2_power;
    int big_seg;
    int len[4];
    virtual void init(int _n, int _m, int _max_kick_steps);
    void clear();
    virtual bool insert(uint64_t ele);
    bool lookup(uint64_t ele);
    virtual bool del(uint64_t ele);
    double get_load_factor();
    double get_full_bucket_factor();
    double get_bits_per_item();

    bool debug_flag = false;
    bool balance = true;
    uint32_t* T;
    uint32_t encode_table[1 << 16];
    uint32_t decode_table[1 << 16];

    ~SemiSortCuckooFilter() { free(T); }

    int filled_cell;
    int full_bucket;
    int max_kick_steps;

    fp_t fingerprint(uint64_t ele);  // 32-bit to 'fp_len'-bit fingerprint

    // interface for semi-sorted bucket
    void get_bucket(int pos, fp_t* store);
    void set_bucket(int pos, fp_t* sotre);
    void test_bucket();
    void make_balance();
    inline int high_bit(fp_t fp);
    inline int low_bit(fp_t fp);
    inline void sort_pair(fp_t& a, fp_t& b);

    virtual int alternate(int pos, fp_t fp) = 0;  // get alternate position
    virtual int insert_to_bucket(
        fp_t* store, fp_t fp);  // insert one fingerprint to bucket [pos]
    virtual int lookup_in_bucket(
        int pos, fp_t fp);  // lookup one fingerprint in  bucket [pos]
    virtual int del_in_bucket(
        int pos, fp_t fp);  // lookup one fingerprint in  bucket [pos]
};

int upperpower2(int x) {
    int ret = 1;
    for (; ret < x;) ret <<= 1;
    return ret;
}

// solve equation : 1 + x(logc - logx + 1) - c = 0
double F_d(double x, double c) { return log(c) - log(x); }
double F(double x, double c) { return 1 + x * (log(c) - log(x) + 1) - c; }
double solve_equation(double c) {
    double x = c + 0.1;
    while (abs(F(x, c)) > 0.001) x -= F(x, c) / F_d(x, c);
    return x;
}
double balls_in_bins_max_load(double balls, double bins) {
    double m = balls;
    double n = bins;
    if (n == 1) return m;

    double c = m / (n * log(n));
    // A more accurate bound..
    if (c < 5)
    {
        double dc = solve_equation(c);
        double ret = (dc - 1 + 2) * log(n);
        return ret;
    }

    double ret = (m / n) + 1.5 * sqrt(2 * m / n * log(n));
    return ret;
}

int proper_alt_range(int M, int i, int* len) {
    double b = 4;      // slots per bucket
    double lf = 0.95;  // target load factor
    int alt_range = 8;
    for (; alt_range < M;) {
        double f = (4 - i) * 0.25;
        if (balls_in_bins_max_load(f * b * lf * M, M * 1.0 / alt_range) <
            0.97 * b * alt_range) {
            return alt_range;
        }
        alt_range <<= 1;
    }
    return alt_range;
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::init(int max_item, int _m, int _step) {
    int _n = (max_item / 0.96 / 4);


    if (_n < 10000) {
        if (_n < 256)
            big_seg = (upperpower2(_n));
        else
            big_seg = (upperpower2(_n / 4));
        _n = ROUNDUP(_n, big_seg);
        len[0] = big_seg - 1;
        len[1] = big_seg - 1;
        len[2] = big_seg - 1;
        len[3] = big_seg - 1;
    } else {
        big_seg = 0;
        big_seg = max(big_seg, proper_alt_range(_n, 0, len));
        int new_n = ROUNDUP(_n, big_seg);
        _n = new_n;

        big_seg--;
        len[0] = big_seg;
        for (int i = 1; i < 4; i++) len[i] = proper_alt_range(_n, i, len) - 1;
        len[0] = max(len[0], 1024);
        len[3] = (len[3] + 1) * 2 - 1;
    }

    this->n = _n;
    this->m = _m;
    this->max_kick_steps = _step;
    this->filled_cell = 0;
    this->full_bucket = 0;

    uint64_t how_many_bit = (uint64_t)this->n * this->m * (fp_len - 1);
    this->memory_consumption =
        ROUNDUP(how_many_bit + 64, 8) / 8 + 8;  // how many bytes !

    max_2_power = 1;
    for (; max_2_power * 2 < _n;) max_2_power <<= 1;
    this->T = (uint32_t*)calloc(this->memory_consumption, sizeof(char));

    int index = 0;
    for (int i = 0; i < 16; i++)
        for (int j = 0; j < ((i == 0) ? 1 : i + 1); j++)
            for (int k = 0; k < ((j == 0) ? 1 : j + 1); k++)
                for (int l = 0; l < ((k == 0) ? 1 : k + 1); l++) {
                    int plain_bit = (i << 12) + (j << 8) + (k << 4) + l;
                    encode_table[plain_bit] = index;
                    decode_table[index] = plain_bit;
                    ++index;
                }
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::clear() {
    this->filled_cell = 0;
    memset(this->T, 0, this->memory_consumption);
}

template <typename fp_t, int fp_len>
fp_t SemiSortCuckooFilter<fp_t, fp_len>::fingerprint(uint64_t ele) {
    fp_t h =
        HashUtil::MurmurHash64(ele ^ 0x192837319273LL) % ((1ull << fp_len) - 1) +
        1;
    return h;
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::get_bucket(int pos, fp_t* store) {
    // Default :
    //
    // Little Endian Store
    // Store by uint32_t
    // store[this -> m] = bucket number

    // 1. read the endcoded bits from memory

    int bucket_length = (fp_len - 1) * 4;
    uint64_t start_bit_pos = (uint64_t)pos * bucket_length;
    uint64_t end_bit_pos = start_bit_pos + bucket_length - 1;
    uint64_t result = 0;

    if (ROUNDDOWN(start_bit_pos, 64) == ROUNDDOWN(end_bit_pos, 64)) {
        uint64_t unit = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64];
        int reading_lower_bound = start_bit_pos & 63;
        int reading_upper_bound = end_bit_pos & 63;

        result = ((uint64_t)unit & ((-1ULL) >> (63 - reading_upper_bound))) >>
                 reading_lower_bound;
    } else {
        uint64_t unit1 = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64];
        uint64_t unit2 = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64 + 1];

        int reading_lower_bound = start_bit_pos & 63;
        int reading_upper_bound = end_bit_pos & 63;

        uint64_t t1 = unit1 >> reading_lower_bound;
        uint64_t t2 = (unit2 & ((1ULL << (reading_upper_bound + 1)) - 1))
                      << (64 - reading_lower_bound);
        result = t1 + t2;
    }

    // 2. read the 4 elements from the encoded bits
    // We use 12 bits to store the 16 most significant bits for the items in
    // bucket, 4 bits per item the low bits are stored in the remaining bits
    //
    // For example, 8 bits per item , require 28 bits to store:
    //
    // Original :
    //
    // hhhh llll
    // hhhh llll
    // hhhh llll
    // hhhh llll
    //
    // encoded :
    //
    //
    // 0 - 11                       12 - 15    16 - 19  20-23   24 - 27
    // HHHHHHHHHHHH                 llll       llll     llll    llll
    //  encoded high bit(12 bits)   item 0     item 1   item 2  item 3
    //
    int decode_result = decode_table[result >> (4 * (fp_len - 4))];

    store[3] = (result & ((1 << (fp_len - 4)) - 1)) +
               ((decode_result & ((1 << 4) - 1)) << (fp_len - 4));
    store[2] = ((result >> (1 * (fp_len - 4))) & ((1 << (fp_len - 4)) - 1)) +
               (((decode_result >> 4) & ((1 << 4) - 1)) << (fp_len - 4));
    store[1] = ((result >> (2 * (fp_len - 4))) & ((1 << (fp_len - 4)) - 1)) +
               (((decode_result >> 8) & ((1 << 4) - 1)) << (fp_len - 4));
    store[0] = ((result >> (3 * (fp_len - 4))) & ((1 << (fp_len - 4)) - 1)) +
               (((decode_result >> 12) & ((1 << 4) - 1)) << (fp_len - 4));

    store[4] = 0;
    store[4] += store[0] != 0;
    store[4] += store[1] != 0;
    store[4] += store[2] != 0;
    store[4] += store[3] != 0;
}

template <typename fp_t, int fp_len>
inline void SemiSortCuckooFilter<fp_t, fp_len>::sort_pair(fp_t& a, fp_t& b) {
    if ((a) < (b)) swap(a, b);
}

template <typename fp_t, int fp_len>
void SemiSortCuckooFilter<fp_t, fp_len>::set_bucket(int pos, fp_t* store) {
    // 0. sort store ! descendant order >>>>>>

    sort_pair(store[0], store[2]);
    sort_pair(store[1], store[3]);
    sort_pair(store[0], store[1]);
    sort_pair(store[2], store[3]);
    sort_pair(store[1], store[2]);

    // 1. compute the encode

    uint64_t high_bit = 0;
    uint64_t low_bit = 0;

    low_bit =
        (store[3] & ((1 << (fp_len - 4)) - 1)) +
        ((store[2] & ((1 << (fp_len - 4)) - 1)) << (1 * (fp_len - 4))) +
        (((uint64_t)store[1] & ((1 << (fp_len - 4)) - 1)) << (2 * (fp_len - 4))) +
        (((uint64_t)store[0] & ((1 << (fp_len - 4)) - 1)) << (3 * (fp_len - 4)));

    high_bit = ((store[3] >> (fp_len - 4)) & ((1 << 4) - 1)) +
               (((store[2] >> (fp_len - 4)) & ((1 << 4) - 1)) << 4) +
               (((store[1] >> (fp_len - 4)) & ((1 << 4) - 1)) << 8) +
               (((store[0] >> (fp_len - 4)) & ((1 << 4) - 1)) << 12);


    // 2. store into memory
    uint64_t high_encode = encode_table[high_bit];
    uint64_t all_encode = (high_encode << (4 * (fp_len - 4))) + low_bit;

    int bucket_length = (fp_len - 1) * 4;
    uint64_t start_bit_pos = (uint64_t)pos * bucket_length;
    uint64_t end_bit_pos = start_bit_pos + bucket_length - 1;

    if (ROUNDDOWN(start_bit_pos, 64) == ROUNDDOWN(end_bit_pos, 64)) {
        uint64_t unit = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64];
        int writing_lower_bound = start_bit_pos & 63;
        int writing_upper_bound = end_bit_pos & 63;

        ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64] =
            (unit & (((1ULL << writing_lower_bound) - 1) +
                     ((-1ULL) - ((-1ULL) >> (63 - writing_upper_bound))))) +
            ((all_encode &
              ((1ULL << (writing_upper_bound - writing_lower_bound + 1)) - 1))
             << writing_lower_bound);
    } else {
        uint64_t unit1 = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64];
        uint64_t unit2 = ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64 + 1];
        int writing_lower_bound = start_bit_pos & 63;
        int writing_upper_bound = end_bit_pos & 63;
        uint64_t lower_part =
            all_encode & ((1LL << (64 - writing_lower_bound)) - 1);
        uint64_t higher_part = all_encode >> (64 - writing_lower_bound);
        ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64] =
            (unit1 & ((1LL << writing_lower_bound) - 1)) +
            (lower_part << writing_lower_bound);
        ((uint64_t*)T)[ROUNDDOWN(start_bit_pos, 64) / 64 + 1] =
            ((unit2 >> (writing_upper_bound + 1)) << (writing_upper_bound + 1)) +
            higher_part;
    }

}


template <typename fp_t, int fp_len>
inline int SemiSortCuckooFilter<fp_t, fp_len>::high_bit(fp_t fp) {
    return (fp >> (fp_len - 4)) & ((1 << 4) - 1);
}

template <typename fp_t, int fp_len>
inline int SemiSortCuckooFilter<fp_t, fp_len>::low_bit(fp_t fp) {
    return fp & ((1 << (fp_len - 4)) - 1);
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::insert_to_bucket(fp_t* store, fp_t fp) {
    // if success return 0
    // if find collision : return 1 + position
    // if full : return 1 + 4

    if (store[this->m] == this->m)
        return 1 + 4;
    else {
        store[3] = fp; // sorted -- store[3] must be empty !
        return 0;
    }
}

template <typename fp_t, int fp_len>
bool SemiSortCuckooFilter<fp_t, fp_len>::insert(uint64_t ele) {
    return false;
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::lookup_in_bucket(int pos, fp_t fp) {
    // If lookup success return 1
    // If lookup fail and the bucket is full return 2
    // If lookup fail and the bucket is not full return 3

    fp_t store[8];
    get_bucket(pos, store);

    int isFull = 1;
    for (int i = 0; i < this->m; i++) {
        fp_t t = store[i];
        if (t == fp) return 1;
        isFull &= (t != 0);
    }
    return (isFull) ? 2 : 3;
}

template <typename fp_t, int fp_len>
bool SemiSortCuckooFilter<fp_t, fp_len>::lookup(uint64_t ele) {
    // If ele is positive, return true
    // negative -- return false

    ele = HashUtil::MurmurHash64(ele ^ 0x12891927);

    fp_t fp = fingerprint(ele);
    int pos1 = this->position_hash(ele);

    int ok1 = lookup_in_bucket(pos1, fp);

    if (ok1 == 1) return true;
    // if (ok1 == 3) return false;

    int pos2 = alternate(pos1, fp);
    assert(pos1 == alternate(pos2, fp));
    int ok2 = lookup_in_bucket(pos2, fp);

    return ok2 == 1;
}

template <typename fp_t, int fp_len>
int SemiSortCuckooFilter<fp_t, fp_len>::del_in_bucket(int pos, fp_t fp) {
    fp_t store[8];
    get_bucket(pos, store);

    for (int i = 0; i < this->m; i++) {
        fp_t t = store[i];
        if (t == fp) {
            store[i] = 0;
            --this->filled_cell;
            set_bucket(pos, store);
            return 1;
        }
    }
    return 0;
}

template <typename fp_t, int fp_len>
bool SemiSortCuckooFilter<fp_t, fp_len>::del(uint64_t ele) {
    // If ele is positive, return true
    // negative -- return false

    ele = HashUtil::MurmurHash64(ele ^ 0x12891927);
    fp_t fp = fingerprint(ele);
    int pos1 = this->position_hash(ele);

    int ok1 = del_in_bucket(pos1, fp);

    if (ok1 == 1) return true;
    // if (ok1 == 3) return false;

    int pos2 = alternate(pos1, fp);
    int ok2 = del_in_bucket(pos2, fp);

    return ok2 == 1;
}
template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_load_factor() {
    return filled_cell * 1.0 / this->n / this->m;
}


template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_bits_per_item() {
    return double(this->memory_consumption) * 8 / filled_cell;
}

template <typename fp_t, int fp_len>
double SemiSortCuckooFilter<fp_t, fp_len>::get_full_bucket_factor() {
    return full_bucket * 1.0 / this->n;
}

template <typename fp_t, int fp_len>
class VacuumFilter : public SemiSortCuckooFilter<fp_t, fp_len> {
private:
    virtual int alternate(int pos, fp_t fp)  // get alternate position
    {
        uint32_t fp_hash = fp * 0x5bd1e995;
        int seg = this->len[fp & 3];
        return pos ^ (fp_hash & seg);
    }

public:
    bool insert(uint64_t ele) {
        // If insert success return true
        // If insert fail return false

        ele = HashUtil::MurmurHash64(ele ^ 0x12891927);

        fp_t fp = this->fingerprint(ele);
        int cur1 = this->position_hash(ele);
        int cur2 = alternate(cur1, fp);

        fp_t store1[8];
        fp_t store2[8];

        this->get_bucket(cur1, store1);
        this->get_bucket(cur2, store2);

        if (store1[this->m] <= store2[this->m]) {
            if (this->insert_to_bucket(store1, fp) == 0) {
                this->filled_cell++;
                this->set_bucket(cur1, store1);
                return true;
            }
        } else {
            if (this->insert_to_bucket(store2, fp) == 0) {
                this->filled_cell++;
                this->set_bucket(cur2, store2);
                return true;
            }
        }

        // randomly choose one bucket's elements to kick
        int rk = rand() % this->m;

        // get those item
        int cur;
        fp_t* cur_store;

        if (rand() & 1)
            cur = cur1, cur_store = store1;
        else
            cur = cur2, cur_store = store2;

        fp_t tmp_fp = cur_store[rk];
        cur_store[rk] = fp;
        this->set_bucket(cur, cur_store);

        int alt = alternate(cur, tmp_fp);

        for (int i = 0; i < this->max_kick_steps; i++) {
            memset(store1, 0, sizeof(store1));
            this->get_bucket(alt, store1);
            if (store1[this->m] == this->m) {
                for (int j = 0; j < this->m; j++) {
                    int nex = alternate(alt, store1[j]);
                    this->get_bucket(nex, store2);
                    if (store2[this->m] < this->m) {
                        store2[this->m - 1] = store1[j];
                        store1[j] = tmp_fp;
                        this->filled_cell++;
                        this->set_bucket(nex, store2);
                        this->set_bucket(alt, store1);
                        return true;
                    }
                }

                rk = rand() % this->m;
                fp = store1[rk];
                store1[rk] = tmp_fp;
                this->set_bucket(alt, store1);

                tmp_fp = fp;
                alt = alternate(alt, tmp_fp);
            } else {
                store1[this->m - 1] = tmp_fp;
                this->filled_cell++;
                this->set_bucket(alt, store1);
                return true;
            }
        }
        return false;
    }
};

#endif
