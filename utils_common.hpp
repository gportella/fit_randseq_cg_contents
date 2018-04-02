
#ifndef PERIOD_ELASTIC_H
#define PERIOD_ELASTIC_H

#include <functional>
#include <iostream>
#include <math.h>
#include <seqan/arg_parse.h>
#include <seqan/bed_io.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <stdlib.h>
#include <string>
#include <vector>

#define NUC_LEN 147
#define NUC_CORE 74

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// exponentiate, any faster than loops? hard to believe
// from http://stackoverflow.com/a/40642375
struct op {
  double operator()(double d) const { return std::exp(d); }
};
std::vector<double> exp_vect(const std::vector<double> &x) {
  const int n = x.size();
  std::vector<double> y(n);
  std::transform(x.begin(), x.end(), y.begin(), op());
  return y;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename tt1> bool checkArguments(tt1 cond) {
  if (cond.cutoff > 147) {
    std::cerr << "A cutoff longer than the length of a nucleosome does not "
                 "make sense. "
              << std::endl;
    exit(0);
  } else {
    return true;
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// check if char is not in Dna5String
template <typename tt1> bool notNInside(tt1 seq) {
  Iterator<Dna5String>::Type it = seqan::begin(seq);
  Iterator<Dna5String>::Type itEnd = seqan::end(seq);
  for (; it != itEnd; goNext(it)) {
    if (getValue(it) == 'N') {
      return false;
    }
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Functions to get random sequence

typedef std::vector<char> char_array;
char_array charset() {
  // Change this to suit
  return char_array({'A', 'T', 'C', 'G'});
};

// given a function that generates a random character,
// return a string of the requested length
std::string random_string(size_t length, std::function<char(void)> rand_char) {
  std::string str(length, 0);
  std::generate_n(str.begin(), length, rand_char);
  return str;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Dna5String genRandSeq(int lbp) {
  // 0) create the character set.
  //   yes, you can use an array here,
  //   but a function is cleaner and more flexible
  const auto ch_set = charset();

  // 1) create a non-deterministic random number generator
  std::default_random_engine rng(std::random_device{}());
  // 2) create a random number "shaper" that will give
  //   us uniformly distributed indices into the character set
  std::uniform_int_distribution<> dist(0, ch_set.size() - 1);
  // 3) create a function that ties them together, to get:
  //   a non-deterministic uniform distribution from the
  //   character set of your choice.
  auto randchar = [ch_set, &dist, &rng]() { return ch_set[dist(rng)]; };

  Dna5String seqout = random_string(lbp, randchar);
  return seqout;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Adds two vectors, first checks if the are of the same length
template <typename T>
std::vector<T> safe_add_vector(std::vector<T> const &v1,
                               std::vector<T> const &v2) {

  std::vector<T> result(v1.size());

  // better safe than sorry
  if (v1.size() == v2.size()) {
    for (unsigned i = 0; i < v1.size(); ++i) {
      result[i] = v1[i] + v2[i];
    }

  } else {
    std::cerr << "Trying to add two vectors of different length." << std::endl;
    std::cerr << "v1: " << v1.size() << " vs "
              << "v2: " << v2.size() << std::endl;
    exit(1);
  }
  return result;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Adds two vectors, ASSUMES they are both of the same size. Use only if
// you are completely sure that both have the same length or expect undefined
// behaviour
template <typename T>
std::vector<T> brave_add_vector(std::vector<T> const &v1,
                                std::vector<T> const &v2) {

  std::vector<T> result(v1.size());

  for (unsigned i = 0; i < v1.size(); ++i) {
    result[i] = v1[i] + v2[i];
  }

  return result;
}


///////////////////////////////////////////////////////////////////////////////
// Return a std vect of bed3 records from a bed3 file name
///////////////////////////////////////////////////////////////////////////////
template <typename tt1> std::vector<BedRecord<Bed3>> read_in_bed3(tt1 fn) {
  // Open input bed file.
  BedFileIn bedIn;
  if (!open(bedIn, toCString(fn))) {
    std::cerr << "ERROR: Could not open bed file " << toCString(fn)
              << std::endl;
    exit(1);
  }
  BedRecord<Bed3> inBedrecord;
  std::vector<BedRecord<Bed3>> vect_bed3_records;
  try {
    while (!atEnd(bedIn)) {
      readRecord(inBedrecord, bedIn);
      vect_bed3_records.push_back(inBedrecord);
    }
  } catch (Exception const &e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    exit(1);
  }
  return vect_bed3_records;
}

///////////////////////////////////////////////////////////////////////////////
// Return a StringSet<Dna5String> from a vector of bed3 records
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2>
StringSet<Dna5String> get_seqs_from_beds(tt1 bedRecords, tt2 genomeFileName) {

  // Translate sequence name to index.
  // Try to load index and create on the fly if necessary.
  FaiIndex faiIndex;
  if (!open(faiIndex, toCString(genomeFileName))) {
    if (!build(faiIndex, toCString(genomeFileName))) {
      std::cerr << "ERROR: Index could not be loaded or built.\n";
      exit(1);
    }
    if (!save(faiIndex)) // Name is stored from when reading.
    {
      std::cerr << "ERROR: Index could not be written do disk.\n";
      exit(1);
    }
  }

  StringSet<Dna5String> seqs;
  for (auto const &rec : bedRecords) {
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, rec.ref)) {
      std::cerr << "ERROR: Index does not know about sequence " << rec.ref
                << "\n";
      exit(1);
    }

    // Convert positions into integers.
    unsigned beginPos = rec.beginPos, endPos = rec.endPos;

    // Make sure begin and end pos are on the sequence and begin <= end.
    if (beginPos > sequenceLength(faiIndex, idx))
      beginPos = sequenceLength(faiIndex, idx);
    if (endPos > sequenceLength(faiIndex, idx))
      endPos = sequenceLength(faiIndex, idx);
    if (beginPos > endPos)
      endPos = beginPos;

    // Finally, get infix of sequence.
    Dna5String sequenceInfix;
    readRegion(sequenceInfix, faiIndex, idx, beginPos, endPos);
    appendValue(seqs, sequenceInfix);
  }
  return seqs;
}

///////////////////////////////////////////////////////////////////////////////
// write out x y
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2> void write_xy(tt1 fout, tt2 x, tt2 y) {
  std::ofstream file;
  file.open(toCString(fout));
  for (unsigned i = 0; i < x.size(); ++i) {
    file << x[i] << " " << y[i] << std::endl;
  }
  file.close();
}

///////////////////////////////////////////////////////////////////////////////
// write out x
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2> void write_x(tt1 fout, tt2 x) {
  std::ofstream file;
  file.open(toCString(fout));
  for (unsigned i = 0; i < x.size(); ++i) {
    file << x[i] << std::endl;
  }
  file.close();
}

#endif /* end protective inclusion */
