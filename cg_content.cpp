
#ifndef CG_FIT_MAIN
#define CG_FIT_MAIN
#include "commandline_parse.hpp"
#include "compute_contents.hpp"
#include "utils_common.hpp"
#include <functional>
#include <iostream>
#include <math.h>
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
  CharString outCGrefFileName = parseOptions.outCGrefFileName;
  CharString outCGrandFileName = parseOptions.outCGrandFileName;

  bool b_onlyprof = parseOptions.b_onlyprof;

  /////////////////////////////////////////////////////////////////////////////
  // Reading inputs
  /////////////////////////////////////////////////////////////////////////////

  std::vector<BedRecord<Bed3>> refBedRecords = read_in_bed3(refbedFileName);

  // get ref_seqs
  // Try to load index and create on the fly if necessary.
  StringSet<Dna5String> ref_seqs =
      get_seqs_from_beds(refBedRecords, genomeFileName);

  // get rand_seqs
  // Try to load index and create on the fly if necessary.
  if (!b_onlyprof) {
    std::vector<BedRecord<Bed3>> randBedRecords = read_in_bed3(randbedFileName);
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
    write_x(outCGrefFileName, results_sampling.reference_cg);
    write_x(outCGrandFileName, results_sampling.reference_cg);
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
    close(bedFileOut);
  } else {
    std::vector<double> reference_cg =
        compute_cg_contents(ref_seqs, parseOptions.kmer_window, parseOptions.cg_background);
    write_x(outCGrefFileName, reference_cg);
  }

  // don't forget to close the output bed file!!

  exit(0);
}

#endif /* end protective inclusion */
