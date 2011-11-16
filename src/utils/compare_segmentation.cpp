/* compare_segmentation: segmentation schemes evaluating program
 * Song Qiang <qiang.song@usc.edu> 2009
 *
 * This program expects two segmentation schemes and yields
 * various measurements of their similarity and their own 
 * properties
 *
 */ 

#include <string>
#include <vector>
#include <iostream>

using namespace std;

#include "GenomicRegion.hpp"
#include "GenomicRegionTool.hpp"
#include "OptionParser.hpp"

void 
compute_mutual_information(const vector<SimpleGenomicRegion> &genomic_region_set_1,
						   const vector<SimpleGenomicRegion> &genomic_region_set_2)
{
		vector< vector<SimpleGenomicRegion> > region_set_by_chrom_1;
		vector< vector<SimpleGenomicRegion> > region_set_by_chrom_2;
		GenomicRegionTool::seperate_chromosomes(genomic_region_set_1,
												region_set_by_chrom_1);
		GenomicRegionTool::seperate_chromosomes(genomic_region_set_2,
												region_set_by_chrom_2);
		
		size_t lhs_counter = 0, rhs_counter = 0;
		string lhs_chrom_name, rhs_chrom_name;
		while(lhs_counter < region_set_by_chrom_1.size()
			  && rhs_counter < region_set_by_chrom_2.size())
		{
				lhs_chrom_name = region_set_by_chrom_1[lhs_counter][0].get_chrom();
				rhs_chrom_name = region_set_by_chrom_2[rhs_counter][0].get_chrom();
						
				if(lhs_chrom_name == rhs_chrom_name)
				{
						SimpleGenomicRegion chrom(lhs_chrom_name,
												  0,
												  std::max(region_set_by_chrom_1[lhs_counter].back().get_end(),
														   region_set_by_chrom_2[rhs_counter].back().get_end()));
						double mutual_infor =
								GenomicRegionTool::mutual_information(
										  region_set_by_chrom_1[lhs_counter],
										  region_set_by_chrom_2[rhs_counter],
										  chrom);
						double entropy =
								GenomicRegionTool::entropy(region_set_by_chrom_1[lhs_counter],
														   chrom)
								+
								GenomicRegionTool::entropy(region_set_by_chrom_2[rhs_counter],
															 chrom);

						std::cout << lhs_chrom_name << "\t"
								  << 2 * mutual_infor / entropy
								  << std::endl;
						lhs_counter++;
						rhs_counter++;
						continue;
				}

				if(lhs_chrom_name < rhs_chrom_name)
				{				
						std::cout << lhs_chrom_name
								  << "\t" << 0.0 << std::endl;
						lhs_counter++;
				}
				else
				{
						std::cout << rhs_chrom_name
								  << "\t" << 0.0 << std::endl;
						rhs_counter++;
				}
		}

		while(lhs_counter < region_set_by_chrom_1.size())
		{
				std::cout << region_set_by_chrom_1[lhs_counter].front().get_chrom()
						  << "\t" << 0.0 << std::endl;
				lhs_counter++;
		}

		while(rhs_counter < region_set_by_chrom_2.size())
		{
				std::cout << region_set_by_chrom_2[rhs_counter].front().get_chrom()
						  << "\t" << 0.0 << std::endl;
				rhs_counter++;
		}

}

void 
compute_distance(const vector<SimpleGenomicRegion> &genomic_region_set_1,
				 const vector<SimpleGenomicRegion> &genomic_region_set_2)
{
		vector< vector<SimpleGenomicRegion> > region_set_by_chrom_1;
		vector< vector<SimpleGenomicRegion> > region_set_by_chrom_2;
		GenomicRegionTool::seperate_chromosomes(genomic_region_set_1,
												region_set_by_chrom_1);
		GenomicRegionTool::seperate_chromosomes(genomic_region_set_2,
												region_set_by_chrom_2);
		
		size_t lhs_counter = 0, rhs_counter = 0;
		string lhs_chrom_name, rhs_chrom_name;
		while(lhs_counter < region_set_by_chrom_1.size()
			  && rhs_counter < region_set_by_chrom_2.size())
		{
				lhs_chrom_name = region_set_by_chrom_1[lhs_counter][0].get_chrom();
				rhs_chrom_name = region_set_by_chrom_2[rhs_counter][0].get_chrom();
						
				if(lhs_chrom_name == rhs_chrom_name)
				{
						SimpleGenomicRegion chrom(lhs_chrom_name,
												  0,
												  std::max(region_set_by_chrom_1[lhs_counter].back().get_end(),
														   region_set_by_chrom_2[rhs_counter].back().get_end()));
						std::cout << lhs_chrom_name << "\t"
								  << GenomicRegionTool::distance(
										  region_set_by_chrom_1[lhs_counter],
										  region_set_by_chrom_2[rhs_counter],
										  chrom)
								  << std::endl;
						lhs_counter++;
						rhs_counter++;
						continue;
				}

				if(lhs_chrom_name < rhs_chrom_name)
				{				
						std::cout << lhs_chrom_name
								  << "\t" << 1.0 << std::endl;
						lhs_counter++;
				}
				else
				{
						std::cout << rhs_chrom_name
								  << "\t" << 1.0 << std::endl;
						rhs_counter++;
				}
		}

		while(lhs_counter < region_set_by_chrom_1.size())
		{
				std::cout << region_set_by_chrom_1[lhs_counter].front().get_chrom()
						  << "\t" << 1.0 << std::endl;
				lhs_counter++;
		}

		while(rhs_counter < region_set_by_chrom_2.size())
		{
				std::cout << region_set_by_chrom_2[rhs_counter].front().get_chrom()
						  << "\t" << 1.0 << std::endl;
				rhs_counter++;
		}

}

int
main(int argc, const char **argv)
{
		const int Exit_Success = 0;
		const int Exit_Failure = -1;
		const bool Option_Not_Required = false;

		bool BOOL_normalized_intersection_score = false;
		bool BOOL_mutual_information_score = false;
		bool BOOL_mutual_information_distance = false;
		bool BOOL_summary = false;
		bool VERBOSE = false;

		string file_name_1, file_name_2;

		// Parsing Command line options 
		try
		{
				OptionParser opt_parse("compare_segmentation",
									   "Compare two segmentation schemes and get their info",
									   "file_name_1 file_name_2");
				opt_parse.add_opt("intersection", 'i',
								  "Output normalized intersection score",
								  Option_Not_Required,
								  BOOL_normalized_intersection_score);
				opt_parse.add_opt("mutual_information", 'm',
								  "Output mutual information score",
								  Option_Not_Required,
								  BOOL_mutual_information_score);
				opt_parse.add_opt("summary", 's',
								  "Output summary information for each segmentation",
								  Option_Not_Required,
								  BOOL_summary);
				opt_parse.add_opt("distance", 'd',
								  "Output distance based on mutual infomation",
								  Option_Not_Required,
								  BOOL_mutual_information_distance);
				opt_parse.add_opt("verbose", 'v',
								  "Output verbose information",
								  Option_Not_Required,
								  VERBOSE);
				
				vector<string> leftover_args;
				opt_parse.parse(argc, argv, leftover_args);
						
				if (argc == 1 || opt_parse.help_requested())
				{
						cerr << opt_parse.help_message() << endl;
						return Exit_Success;
				}
				if (opt_parse.about_requested())
				{
						cerr << opt_parse.about_message() << endl;
						return Exit_Success;
				}
				if (opt_parse.option_missing())
				{
						cerr << opt_parse.option_missing_message() << endl;
						return Exit_Success;
				}
				if (leftover_args.size() < 2) 
				{
						cerr << opt_parse.help_message() << endl;
						return Exit_Success;
				}
		
				// reads in segmentation files
				if(VERBOSE)
						std::cout << "Reading segmentation files" << std::endl;
		
				file_name_1 = leftover_args[0];
				file_name_2 = leftover_args[1];
		}
		catch(RMAPOptionException e)
		{
				std::cerr << e.what() << std::endl;
				return Exit_Failure;
		}
		
		
		// reads in segmentation files
		if(VERBOSE)
				std::cout << "Reading segmentation files" << std::endl;
		
		vector<SimpleGenomicRegion> genomic_region_set_1, genomic_region_set_2;
		try
		{
				GenomicRegionTool::ReadBEDFile(file_name_1,
											   genomic_region_set_1);
				GenomicRegionTool::ReadBEDFile(file_name_2,
											   genomic_region_set_2);
		}
		catch(...)
		{
				std::cerr << "Error: can not open segmentation files" << std::endl;
				return Exit_Failure;
		}
		
		// Output summary information for each segmentation
		if(BOOL_summary)
		{
				std::cout << "Summary for segmenttion " + file_name_1
						  << std::endl;
				GenomicRegionTool::summary(genomic_region_set_1);
				std::cout << "Summary for segmenttion " + file_name_2
						  << std::endl;
				GenomicRegionTool::summary(genomic_region_set_2);
		}

		// Output mutual_information score 
		if(BOOL_mutual_information_score)
		{
				std::cout << "calculating normalized mutual information (symmetric uncertainty)" << std::endl;
				compute_mutual_information(genomic_region_set_1,
										   genomic_region_set_2);
		}

		// Output mutual_information based distance 
		if(BOOL_mutual_information_distance)
		{
				std::cout << "calculating mutual information based distance" << std::endl;
				compute_distance(genomic_region_set_1,
								 genomic_region_set_2);
		}

		// Output normalized intersection score
		if(BOOL_normalized_intersection_score)
		{
				std::cout << "calculating normalized intersection score" << std::endl;
				std::cout << GenomicRegionTool::intersection2union_ratio(genomic_region_set_1,
																		 genomic_region_set_2) 
						  << std::endl;
		}
		
		return Exit_Success;
}

