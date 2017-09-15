
#ifndef CG_FIT_MAIN
#define CG_FIT_MAIN
#include "commandline_parse.hpp"
#include "compute_contents.hpp"
#include "utils_common.hpp"
#include <Eigen/Dense>
#include <functional>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <seqan/arg_parse.h>
#include <seqan/bed_io.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

int main(int argc, char const **argv) {

  // OpenMP setting the maximum number of threads possible
  int nProcessors = omp_get_max_threads();
  omp_set_num_threads(nProcessors);

  // Parse the command line
  Options parseOptions;
  ArgumentParser::ParseResult res = parseCommandLine(parseOptions, argc, argv);

  // If parsing did not work, then exit with code 1
  // Otherwise exit with code 0
  if (res != ArgumentParser::PARSE_OK)
    return res == ArgumentParser::PARSE_ERROR;

  // input files
  CharString genomeFileName = parseOptions.genomeFileName;
  CharString refbedFileName = parseOptions.refbedFileName;
  CharString randbedFileName = parseOptions.randbedFileName;
  // declarations for fasta inputs

  /////////////////////////////////////////////////////////////////////////////
  // Reading inputs
  /////////////////////////////////////////////////////////////////////////////

  std::vector<BedRecord<Bed3>> refBedRecords = read_in_bed3(refbedFileName);
  std::vector<BedRecord<Bed3>> randBedRecords = read_in_bed3(randbedFileName);

  // get ref_seqs
  // Try to load index and create on the fly if necessary.
  StringSet<Dna5String> ref_seqs =
      get_seqs_from_beds(refBedRecords, genomeFileName);

  // get rand_seqs
  // Try to load index and create on the fly if necessary.
  StringSet<Dna5String> rand_seqs =
      get_seqs_from_beds(randBedRecords, genomeFileName);

  /////////////////////////////////////////////////////////////////////////////
  // Carry out the calculations
  /////////////////////////////////////////////////////////////////////////////
  mimic_results results_sampling =
      do_mimic_profile(parseOptions, ref_seqs, rand_seqs, randBedRecords);

  /////////////////////////////////////////////////////////////////////////////
  // Write out the results
  /////////////////////////////////////////////////////////////////////////////
  write_x("reference.txt", results_sampling.reference_cg);
  write_x("random.txt", results_sampling.random_cg);
  BedFileOut bedFileOut;

  if (!open(bedFileOut, toCString(parseOptions.outbedFileName))) {
    std::cerr << "ERROR: Cound not open output bed.\n";
    return 1;
  }
  try {
    for (auto record : results_sampling.random_bed) {
      writeRecord(bedFileOut, record);
    }
  } catch (Exception const &e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return 1;
  }

  // don't forget to close the output bed file!!
  close(bedFileOut);

  exit(0);
}

#endif /* end protective inclusion */
