#include <cuckoofilter.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <iostream>
using namespace std;
#define deln(a) cerr << #a << " " << a << endl

using cuckoofilter::CuckooFilter;
using cuckoofilter::VacuumFilter;

int main(int argc, char **argv) {
  size_t total_items = int(1 << 26) * 0.95;

  // Create a cuckoo filter where each item is of type size_t and
  // use 12 bits for each item:
  //    CuckooFilter<size_t, 12> filter(total_items);
  // To enable semi-sorting, define the storage of cuckoo filter to be
  // PackedTable, accepting keys of size_t type and making 13 bits
  // for each key:
  //   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
  deln(100);
  VacuumFilter<size_t, 12> filter(total_items);

  // Insert items to this cuckoo filter
  size_t num_inserted = 0;
  for (size_t i = 0; i < total_items; i++, num_inserted++) {
    if (filter.Add(i) != cuckoofilter::Ok) {
        deln(i);
      break;
    }
  }

  // Check if previously inserted items are in the filter, expected
  // true for all items
  for (size_t i = 0; i < num_inserted; i++) {
    if (filter.Contain(i) != cuckoofilter::Ok)
    {
        printf("fail\n");
        deln(i);
        filter.Contain(i);
        break;
    }
  }

  // Check non-existing items, a few false positives expected
  size_t total_queries = 0;
  size_t false_queries = 0;
  for (size_t i = total_items; i < 2 * total_items; i++) {
    if (filter.Contain(i) == cuckoofilter::Ok) {
      false_queries++;
    }
    total_queries++;
  }

  // Output the measured false positive rate
  std::cout << "false positive rate is "
            << 100.0 * false_queries / total_queries << "%\n";

  return 0;
}
