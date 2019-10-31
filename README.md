# Vacuum Filter

## Overview

**Vacuum filter** is an approximate set-membership query data structure based on **cuckoo filter**. Vacuum filters costs the smallest space among all known AMQ data structure (when enabling semi-sorting optimization) and provides higher insertion and lookup throughput in most situations. It also support batch-mode operations that significantly increase the throughput.

Our paper will appear in VLDB 2020.  

**Paper Link:** http://www.vldb.org/pvldb/vol13/p197-wang.pdf

## Repository Structure

We implemented two versions of vacuum filters. One is implemented from scratch, which supports arbitrary fingerprint length from 6 to 16. Another one is implemented based on Bin Fan's cuckoo filter(https://github.com/efficient/cuckoofilter), which provides higher throughputs.

1. `vacuum.h` : The VF implemented from scratch (more flexible).
2. `ModifiedCuckooFilter` : The VF implemented based on cuckoo filter (faster).

## Demo

```
make test
./test
```

## API

Here is a small example of VF's APIs.

```c++
VacuumFilter<uint16_t, 16> vf; // vacuum filter with 16-bit fingerprint
vf.init(n, 4, 400); // vf.init(max_item_numbers, slots per bucket, max_kick_steps)
vf.insert(999); // vf.insert(item), item is a 64-bit integer
vf.lookup(999); // true
vf.del(999); // vf.del(item), item should exist in the filter
vf.lookup(999); // false
```

For the details, check out `test.cpp`. It contains two detailed examples of those two versions of VFs and one example of the batch-mode VF.

## Authors

- Mingxun Zhou(zhoumingxun@pku.edu.cn)
- Minmei Wang(mwang107@ucsc.edu)
- Chen Qian(cqian12@ucsc.edu)
- Shouqian Shi(sshi27@ucsc.edu)



