
#ifndef CL_PARSER_H
#define CL_PARSER_H
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
// using namespace std;

struct Options {
  bool b_verbose;
  bool b_onlyprof;
  unsigned int max_iter;
  double disimilarity_cutoff;
  double cg_background;
  unsigned int kmer_window;
  CharString genomeFileName;
  CharString refbedFileName;
  CharString randbedFileName;
  CharString outbedFileName;
  CharString outCGrefFileName;
  CharString outCGrandFileName;

  // I guess this is how to initialize
  Options() : b_verbose(false) {}
};

ArgumentParser::ParseResult parseCommandLine(Options &parseOptions, int argc,
                                             char const **argv) {
  // Setup ArgumentParser.
  ArgumentParser parser("cg_contents");
  setShortDescription(parser, "Computes CG contents ");
  addDescription(parser, "Computes CG contents "
                         "and generate random regions with similar CG profile.");
  addUsageLine(parser, "[\\fIOPTIONS\\fP] ");
  setVersion(parser, "0.4");
  setDate(parser, "November 2017");

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
  addOption(parser,
            ArgParseOption("cgref", "cg_ref_file", "Output CG reference.",
                           ArgParseArgument::OUTPUT_FILE));
  addOption(parser,
            ArgParseOption("cgrand", "cg_rand_file", "Output CG random.",
                           ArgParseArgument::OUTPUT_FILE));

  addOption(parser, ArgParseOption("v", "be_verbose", "Be verbose in "
                                                      "what you do."));
  addOption(parser,
            ArgParseOption("pf", "only_profile", "Compute only the profile "
                                                 "of the reference bed file."));
  addOption(parser, ArgParseOption(
                        "max_iter", "max_iter", "Maximum number of iterations",
                        seqan::ArgParseArgument::INTEGER, "INTEGER"));
  addOption(parser, ArgParseOption("cutoff", "disimilarity_cutoff", "The RMSD cutoff, smaller means that profiles are more similar.",
                                   seqan::ArgParseArgument::DOUBLE, "FLOAT"));
  addOption(parser, ArgParseOption("cg_bkg", "cg_background", "The average CG contents in your reference genome, used to "
                                    "calculate enrichment.",
                                   seqan::ArgParseArgument::DOUBLE, "FLOAT"));
  addOption(parser,
            ArgParseOption("w", "kmer_window",
                           "Window size for extracting CG composition",
                           seqan::ArgParseArgument::INTEGER, "INTEGER"));

  setDefaultValue(parser, "ref_bed_file", "reference_bed.fa");
  setDefaultValue(parser, "rand_bed_file", "random_bed.fa");
  setDefaultValue(parser, "out_bed_file", "random_bed_fit_CG.bed");
  setDefaultValue(parser, "cg_ref_file", "cg_contents_reference.dat");
  setDefaultValue(parser, "cg_rand_file", "cg_contents_random.dat");
  setDefaultValue(parser, "max_iter", "10000");
  setDefaultValue(parser, "disimilarity_cutoff", "0.0001");
  setDefaultValue(parser, "cg_background", "0.42");
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
  getOptionValue(parseOptions.outCGrefFileName, parser, "cg_ref_file");
  getOptionValue(parseOptions.outCGrandFileName, parser, "cg_rand_file");
  getOptionValue(parseOptions.max_iter, parser, "max_iter");
  getOptionValue(parseOptions.disimilarity_cutoff, parser,
                 "disimilarity_cutoff");
  getOptionValue(parseOptions.disimilarity_cutoff, parser,
                 "cg_background");
  getOptionValue(parseOptions.kmer_window, parser, "kmer_window");
  parseOptions.b_verbose = isSet(parser, "be_verbose");
  parseOptions.b_onlyprof = isSet(parser, "only_profile");

  return ArgumentParser::PARSE_OK;
}
#endif /* end of protective declarations */
