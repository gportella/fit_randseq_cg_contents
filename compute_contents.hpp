

#ifndef CM_CONTENTS
#define CM_CONTENTS
#include "utils_common.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <math.h>
#include <seqan/arg_parse.h>
#include <seqan/bed_io.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace seqan;

// for the conditions
typedef struct my_results {
  std::vector<double> reference_cg;
  std::vector<double> random_cg;
  std::vector<BedRecord<Bed3>> random_bed;
} mimic_results;

///////////////////////////////////////////////////////////////////////////////
// Calculate weighted error between two curves
///////////////////////////////////////////////////////////////////////////////
template <typename tt1> double weighted_error(tt1 v1, tt1 v2) {

  double result = 0.0;
  // better safe than sorry
  if (v1.size() == v2.size() and v1.size() > 0) {
    for (unsigned i = 0; i < v1.size(); ++i) {
      result += ((v1[i] - v2[i]) * (v1[i] - v2[i]));
    }

  } else {
    std::cerr << "Trying to divide two vectors of different length or empty"
              << std::endl;
    std::cerr << "v1: " << v1.size() << " vs "
              << "v2: " << v2.size() << std::endl;
    exit(1);
  }
  return sqrt(result / double(v1.size()));
}

///////////////////////////////////////////////////////////////////////////////
// Compute the fraction of CG for a given sequence
///////////////////////////////////////////////////////////////////////////////
template <typename tt1>
std::vector<double> find_cg_composition(tt1 seq, int window, float cg_content_bkground) {

  std::vector<double> cg_comp;

  for (unsigned i = 0; i < length(seq) - window + 1; ++i) {
    Infix<Dna5String>::Type genomeFragment = infix(seq, i, i + window);

    Iterator<Dna5String>::Type it = begin(genomeFragment);
    Iterator<Dna5String>::Type itEnd = end(genomeFragment);

    float n_count = 0, a_count = 0, t_count = 0, c_count = 0, g_count = 0;

    for (; it != itEnd; goNext(it)) {
      // we can't use switch because of SeqAn
      char candidate = getValue(it);
      switch (candidate) {
      case 'A':
        ++a_count;
        break;
      case 'T':
        ++t_count;
        break;
      case 'C':
        ++c_count;
        break;
      case 'G':
        ++g_count;
        break;
      case 'N':
        ++n_count;
        break;
      }
    }
    double gc_content = ((g_count + c_count) /
                         (a_count + t_count + c_count + g_count + n_count)) -
                        cg_content_bkground;

    cg_comp.push_back(gc_content);
  }
  return cg_comp;
}

///////////////////////////////////////////////////////////////////////////////
// Compute the average fraction of CG for a bunch of sequences
// The sequences should all have the same length.
// TODO:  check that the seqs are actually of the same length
//        alternatively use the safe_add_vector function instead of the
//        brave_add_vector one
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2, typename tt3>
std::vector<double> compute_cg_contents(tt1 seqs, tt2 kmer_window, tt3 cg_content_bkground) {

  std::vector<double> av_cg(length(seqs[0]) - kmer_window, 0.0);
  std::vector<double> cg;
  for (unsigned i = 0; i < length(seqs); ++i) {
    cg = find_cg_composition(seqs[i], kmer_window, cg_content_bkground);
    av_cg = brave_add_vector(av_cg, cg);
  }
  // average
  std::transform(av_cg.begin(), av_cg.end(), av_cg.begin(),
                 std::bind2nd(std::divides<double>(), length(seqs)));
  return av_cg;
}

///////////////////////////////////////////////////////////////////////////////
// random vector
// generates numWantedNumbers indices from a maximum of maxSize
// from https://stackoverflow.com/a/12182599
///////////////////////////////////////////////////////////////////////////////
template <typename tt1>
std::vector<int> random_int_vector(tt1 maxSize, tt1 numWantedNumbers) {

  // Generate a vector with all possible numbers and shuffle it.
  std::vector<int> result;
  if (numWantedNumbers < maxSize) {
    std::random_device device;
    std::mt19937 g(device());
    for (unsigned long i = 0; i < maxSize; ++i)
      result.push_back(i);

    std::shuffle(result.begin(), result.end(), g);
    result.resize(numWantedNumbers);
  } else {
    std::cerr << "In random_int_vector, you ask for more elemnts that allowed."
              << std::endl;
    std::cerr << "maxSize: " << maxSize << " vs "
              << "numWantedNumbers: " << numWantedNumbers << std::endl;
    exit(1);
  }
  // sort the vector in descending order
  std::sort(std::begin(result), std::end(result), std::greater<int>());
  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Given a list of indices, return a string of sequences from a lager array
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2>
StringSet<Dna5String> get_seq_from_indices(tt2 randvect, tt1 seqs) {

  StringSet<Dna5String> rand_seqs;
  for (auto const &i : randvect) {
    appendValue(rand_seqs, seqs[i]);
  }

  return rand_seqs;
}

///////////////////////////////////////////////////////////////////////////////
// Given a list of indices, return a bed record from a lager array
///////////////////////////////////////////////////////////////////////////////
template <typename tt1, typename tt2>
std::vector<BedRecord<Bed3>> get_bed_from_indices(tt2 randvect, tt1 v_beds) {

  std::vector<BedRecord<Bed3>> rand_beds;
  for (auto const &i : randvect) {
    appendValue(rand_beds, v_beds[i]);
  }

  return rand_beds;
}

template <typename tt0, typename tt1, typename tt2>
mimic_results do_mimic_profile(tt0 parseOptions, tt1 dseqs, tt1 dseqs_rand,
                               tt2 v_randBedrecords) {

  // obtain the reference CG content curve
  std::vector<double> reference_cg =
      compute_cg_contents(dseqs, parseOptions.kmer_window, parseOptions.cg_background);
  // from the random sequences, extract random indices, as many
  // as sequences we have in the reference, in descending order.
  std::vector<int> init_proposed_ind =
      random_int_vector(length(dseqs_rand), length(dseqs));
  // using the random indices just generated, produce a vector
  // of sequences for those indices
  StringSet<Dna5String> init_proposed_seq =
      get_seq_from_indices(init_proposed_ind, dseqs_rand);
  std::vector<BedRecord<Bed3>> init_proposed_bed =
      get_bed_from_indices(init_proposed_ind, v_randBedrecords);

  // remove init_proposed_ind from dseqs_rand
  // WARNING --> this method requires the indices to
  // be ordered from larger to smaller. This is now the case because
  // random_int_vector returns such a vector
  // doing it like this reorders the vector, but I don't mind
  // other alternatives here https://stackoverflow.com/a/3487736
  for (auto const &i : init_proposed_ind) {
    dseqs_rand[i] = back(dseqs_rand);
    eraseBack(dseqs_rand);
    // now the bedrecords
    v_randBedrecords[i] = v_randBedrecords.back();
    v_randBedrecords.pop_back();
  }

  // the CG content of the initial trial bunch of sequences
  std::vector<double> rand_cg =
      compute_cg_contents(init_proposed_seq, parseOptions.kmer_window, parseOptions.cg_background);
  // how far are we from the target curve?
  double prev_diff = weighted_error(reference_cg, rand_cg);
  if (parseOptions.b_verbose) {
    std::cout << "Initial weighted error is " << prev_diff << std::endl;
  }

  // initialize random generator
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());

  std::vector<double> test_cg(length(dseqs), 0);

  unsigned int iterations = 0;

  while (prev_diff > parseOptions.disimilarity_cutoff and
         iterations < parseOptions.max_iter) {

    ++iterations;

    // uniform_int_distribution goes from [a, b], so you need the -1
    std::uniform_int_distribution<int> distr(0, length(init_proposed_seq) - 1);
    std::uniform_int_distribution<int> distr2(0, length(dseqs_rand) - 1);
    int rand_trial_ref = distr(generator);
    int rand_trial_rand = distr2(generator);
    std::cout << "1 " << rand_trial_ref << std::endl;

    StringSet<Dna5String> proposed_seq = init_proposed_seq;
    std::vector<BedRecord<Bed3>> proposed_bed = init_proposed_bed;
    // remove random from proposed, both seq and bed
    erase(proposed_seq, rand_trial_ref);
    proposed_bed.erase(proposed_bed.begin() + rand_trial_ref);
    // add one random from random source of seq, and the equivalent from
    // bed
    appendValue(proposed_seq, dseqs_rand[rand_trial_rand]);
    appendValue(proposed_bed, v_randBedrecords[rand_trial_rand]);

    test_cg = compute_cg_contents(proposed_seq, parseOptions.kmer_window, parseOptions.cg_background);

    double diff = weighted_error(reference_cg, test_cg);

    if (diff < prev_diff) {
      // remove proposed addition form bag of rand_seqs
      erase(dseqs_rand, rand_trial_rand);
      // remove proposed addition from the bag of bed records
      v_randBedrecords.erase(v_randBedrecords.begin() + rand_trial_rand);
      prev_diff = diff;
      init_proposed_seq = proposed_seq;
      init_proposed_bed = proposed_bed;
      if (parseOptions.b_verbose) {
        std::cout << "Accepted; current diff is " << prev_diff << std::endl;
      }
    }
  }

  mimic_results results_back = {reference_cg, test_cg, init_proposed_bed};
  return results_back;
}

#endif /* end of include guard:  */
