#include <bits/stdc++.h>
#include <time.h>
#include <unistd.h>
#include <chrono>
#include <random>
#include <ratio>
#include "./ModifiedCuckooFilter/src/cuckoofilter.h"
#include "cuckoo.h"
#include "hashutil.h"

using namespace std;

// generate n 64-bit random numbers
void random_gen(int n, vector<uint64_t>& store, mt19937& rd) {
    store.resize(n);
    for (int i = 0; i < n; i++)
        store[i] = (uint64_t(rd()) << 32) + rd();
}

// generate n 64-bit random numbers
void random_gen_1(int n, uint64_t** store, mt19937& rd) {
    *store = new uint64_t[n + 128];
    for (int i = 0; i < n; i++)
        (*store)[i] = (uint64_t(rd()) << 32) + rd();
}
void test_vf_no_padding() {

    /*
        We implemented VF_no_padding from scratch.
        It supports fingerprint length from 4 to 16 bits, but we recommend to use fingerprint longer than 8 bits.
        This version aims at flexibility, so it is slower than VF_with_padding.
    */

    cout << "Testing vacuum filter(no padding)..." << endl;

    int n = 1 << 25; // number of inserted keys
    int q = 10000000; // number of queries

    cout << "Keys number = " << n << endl;
    cout << "Queries number = " << q << endl;

    mt19937 rd(12821);
    vector<uint64_t> insKey;
    vector<uint64_t> alienKey;
    random_gen(n, insKey, rd);
    random_gen(q, alienKey, rd);

    VacuumFilter<uint16_t, 16> vf;

    // vf.init(max_item_numbers, slots per bucket, max_kick_steps)
    vf.init(n, 4, 400);

    for (int i = 0; i < n; i++)
        if (vf.insert(insKey[i]) == false)
            cout << "Insertion fails when inserting " << i << "th key: " << insKey[i] << endl;

    cout << "Load factor = " << vf.get_load_factor() << endl;

    for (int i = 0; i < n; i++)
        if (vf.lookup(insKey[i]) == false)
            cout << "False negative happens at " << i << "th key: " << insKey[i] << endl;
    
    int false_positive_cnt = 0;

    for (int i = 0; i < q; i++)
        if (vf.lookup(alienKey[i]) == true)
            false_positive_cnt++;

    cout << "False positive rate = " << double(false_positive_cnt) / q << endl;
    cout << "Bits per key = " << vf.get_bits_per_item() << endl;

    for (int i = 0; i < n; i++)
        if (vf.del(insKey[i]) == false)
            cout << "Deletion fails when inserting " << i << "th key: " << insKey[i] << endl;

    cout << endl;
}

void test_vf_with_padding() {

    /*
    We also implemented VF_with_padding based on Bin Fan's cuckoo filter.
    It support fingerprint length 4, 8, 12 and 16 (if you enable semi-sorting, it should be 5, 9, 13, 17).
    */

    cout << "Testing vacuum filter(with padding)..." << endl;

    int n = 1 << 25; // number of inserted keys
    int q = 10000000; // number of queries

    cout << "Keys number = " << n << endl;
    cout << "Queries number = " << q << endl;

    mt19937 rd(12821);
    vector<uint64_t> insKey;
    vector<uint64_t> alienKey;
    random_gen(n, insKey, rd);
    random_gen(q, alienKey, rd);

    cuckoofilter::VacuumFilter<size_t, 16> vf(n);

    // If you want to enable semi-sorting to save memory and allow some loss on throughput, use
    // cuckoofilter::VacuumFilter<size_t, 17, cuckoofilter::PackedTable> vf(n);

    for (int i = 0; i < n; i++)
        if (vf.Add(insKey[i]) != cuckoofilter::Ok) {
            cout << "Load factor = " << vf.LoadFactor() << endl;
            cout << "Insertion fails when inserting " << i << "th key: " << insKey[i] << endl;
        }

    cout << "Load factor = " << vf.LoadFactor() << endl;

    for (int i = 0; i < n; i++)
        if (vf.Contain(insKey[i]) != cuckoofilter::Ok)
            cout << "False negative happens at " << i << "th key: " << insKey[i] << endl;
    
    int false_positive_cnt = 0;

    for (int i = 0; i < q; i++)
        if (vf.Contain(alienKey[i]) == cuckoofilter::Ok)
            false_positive_cnt++;

    cout << "False positive rate = " << double(false_positive_cnt) / q << endl;
    cout << "Bits per key = " << vf.BitsPerItem() << endl;

    for (int i = 0; i < n; i++)
        if (vf.Delete(insKey[i]) != cuckoofilter::Ok)
            cout << "Deletion fails when inserting " << i << "th key: " << insKey[i] << endl;

    cout << endl;
}

void test_batch() {

    /*
    We also implemented batching mode for VF.
    Given an array of keys, VF slices the array into multiple batches. Each batch contains 128 keys.
    Then VF performs insertion/deletion/lookup operation for those batches.
    */

    cout << "Testing VF in batching mode..." << endl;

    int n = (1 << 25);
    int q = 10000000;
    cout << "Keys number = " << n << endl;
    cout << "Queries number = " << q << endl;

    mt19937 rd(112983);
    uint64_t* insKey;
    uint64_t* alienKey;
    bool* res;
    res = new bool[max(n, q)];

    random_gen_1(n, &insKey, rd);
    random_gen_1(q, &alienKey, rd);

    cuckoofilter::VacuumFilter<size_t, 16> vf(n);

    // If you want to enable semi-sorting to save memory and allow some loss on throughput, use
    // cuckoofilter::VacuumFilter<size_t, 17, cuckoofilter::PackedTable> vf(n);

    vf.Add_many(insKey, res, n);
    for (int i = 0; i < n; i++)
        if (res[i] == false) {
            cout << "Insertion fails when inserting " << i << "th key: " << insKey[i] << endl;
            break;
        }

    cout << "Load factor = " << vf.LoadFactor() << endl;

    vf.Contain_many(insKey, res, n);
    for (int i = 0; i < n; i++)
        if (res[i] == false) {
            cout << "False negative happens at " << i << "th key: " << insKey[i] << endl;
            break;
        }

    int cnt = 0;
    vf.Contain_many(alienKey, res, q);
    for (int i = 0; i < q; i++) if (res[i] == true) cnt++;

    cout << "False positive rate = " << double(cnt) / q << endl;
    cout << "Bits per key = " << vf.BitsPerItem() << endl;

    vf.Delete_many(insKey, res, n);
    for (int i = 0; i < n; i++)
        if (res[i] == false) {
            cout << "Deletion fails when inserting " << i << "th key: " << insKey[i] << endl;
            break;
        }

    delete insKey;
    delete alienKey;
    delete res;
}

int main() {
    test_vf_no_padding();
    test_vf_with_padding();
    test_batch();
    return 0;
}
