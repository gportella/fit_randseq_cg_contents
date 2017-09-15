
#ifndef CL_PARSER_H
#define CL_PARSER_H
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
// using namespace std;

struct Options {
  bool b_verbose;
  unsigned int max_iter;
  double disimilarity_cutoff;
  unsigned int kmer_window;
  CharString genomeFileName;
  CharString refbedFileName;
  CharString randbedFileName;
  CharString outbedFileName;

  // I guess this is how to initialize
  Options() : b_verbose(false) {}
};

ArgumentParser::ParseResult parseCommandLine(Options &parseOptions, int argc,
                                             char const **argv) {
  // Setup ArgumentParser.
  ArgumentParser parser("cg_contents");
  setShortDescription(parser, "Computes CG contents ");
  addDescription(parser, "Computes CG contents"
                         "and lalala.");
  addUsageLine(parser, "[\\fIOPTIONS\\fP] ");
  setVersion(parser, "0.1");
  setDate(parser, "September 2017");

  // Define Options
  addOption(parser, ArgParseOption("i", "genome_file", "A fasta input file "
                                                       "with sequence to "
                                                       "search for.",
                                   ArgParseArgument::INPUT_FILE));
  addOption(parser, ArgParseOption(
                        "b", "ref_bed_file",
                        "A bz2 BED file with the reference intervals. Should "
                        "have the same lenght as rand seq file.",
                        ArgParseArgument::INPUT_FILE));
  addOption(parser, ArgParseOption("br", "rand_bed_file",
                                   "A bz2 BED file random intervals. Should "
                                   "have the same lenght as rand seq file.",
                                   ArgParseArgument::INPUT_FILE));
  addOption(parser, ArgParseOption("bo", "out_bed_file", "Output bed file.",
                                   ArgParseArgument::OUTPUT_FILE));

  addOption(parser, ArgParseOption("v", "be_verbose", "Be verbose in "
                                                      "what you do."));
  addOption(parser, ArgParseOption(
                        "max_iter", "max_iter", "Maximum number of iterations",
                        seqan::ArgParseArgument::INTEGER, "INTEGER"));
  addOption(parser, ArgParseOption("cutoff", "disimilarity_cutoff", "The ",
                                   seqan::ArgParseArgument::DOUBLE, "FLOAT"));
  addOption(parser,
            ArgParseOption("w", "kmer_window",
                           "Window size for extracting CG composition",
                           seqan::ArgParseArgument::INTEGER, "INTEGER"));

  setDefaultValue(parser, "ref_bed_file", "reference_bed.fa");
  setDefaultValue(parser, "rand_bed_file", "random_bed.fa");
  setDefaultValue(parser, "out_bed_file", "random_bed_fit_CG.bed");
  setDefaultValue(parser, "max_iter", "10000");
  setDefaultValue(parser, "disimilarity_cutoff", "0.0001");
  setDefaultValue(parser, "kmer_window", "40");
  setValidValues(parser, "genome_file",
                 "fna fa fna.gz fasta fa.gz fasta.bz2 fa.bz2");
  setValidValues(parser, "ref_bed_file", "bed bed.bz2 bed.gz");
  setValidValues(parser, "rand_bed_file", "bed bed.bz2 bed.gz");

  // Parse command line.
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // Only extract options if the
  // program continues after
  // parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

  // Extract option values.
  getOptionValue(parseOptions.genomeFileName, parser, "genome_file");
  getOptionValue(parseOptions.refbedFileName, parser, "ref_bed_file");
  getOptionValue(parseOptions.randbedFileName, parser, "rand_bed_file");
  getOptionValue(parseOptions.outbedFileName, parser, "out_bed_file");
  getOptionValue(parseOptions.max_iter, parser, "max_iter");
  getOptionValue(parseOptions.disimilarity_cutoff, parser,
                 "disimilarity_cutoff");
  getOptionValue(parseOptions.kmer_window, parser, "kmer_window");
  parseOptions.b_verbose = isSet(parser, "be_verbose");

  return ArgumentParser::PARSE_OK;
}
#endif /* end of protective declarations */
