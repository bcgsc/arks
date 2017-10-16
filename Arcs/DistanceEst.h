#ifndef _DISTANCE_EST_H_
#define _DISTANCE_EST_H_ 1

#include "Arcs/Arcs.h"
#include "Common/StatUtil.h"

/**
 * Calculate stats regarding number of shared barcodes
 * vs. distance.  These stats are based on examining
 * the number of shared barcodes between the ends of
 * the same contig.
 */
void calcBarcodeToDistStats(const ARCS::IndexMap& imap,
	const ARCS::ContigToLength& contigToLength,
	const std::unordered_map<std::string, int>& indexMultMap,
	ARCS::BarcodeToDistStats& barcodeToDistStats,
	ARCS::ArcsParams& params)
{
	typedef std::string ContigID;
	typedef std::unordered_map<ContigID, unsigned> ContigToNumBarcodes;
	typedef typename ContigToNumBarcodes::const_iterator ContigToNumBarcodesIt;

	/* contig ID => number of barcodes shared between head and tail */
	ContigToNumBarcodes numSharedBarcodes;

	/* tally number of barcodes shared between heads/tails of same contig */

	/* for each chromium barcode */
	for (auto it = imap.begin(); it != imap.end(); ++it) {

		/* skip barcodes outside of min/max multiplicity range */
		std::string index = it->first;
		int indexMult = indexMultMap.at(index);
		if (indexMult < params.min_mult || indexMult > params.max_mult)
			continue;

		/* contig head/tail => number of read pairs from current barcode */
		const ARCS::ScafMap& contigEndToPairs = it->second;

		for (auto o = contigEndToPairs.begin();
			 o != contigEndToPairs.end(); ++o) {

			for (auto p = contigEndToPairs.begin();
				 p != contigEndToPairs.end(); ++p) {

				std::string scafA, scafB;
				bool scafAflag, scafBflag;
				std::tie(scafA, scafAflag) = o->first;
				std::tie(scafB, scafBflag) = p->first;

				/*
				 * skip over pairs that are not head and
				 * tail of same contig
				 */
				if (scafA != scafB || scafAflag == scafBflag)
					continue;

				std::string contigID = scafA;

                ARCS::ContigToLengthIt lengthIt = contigToLength.find(contigID);
				assert(lengthIt != contigToLength.end());
				unsigned l = lengthIt->second;

				/*
				 * skip contigs shorter than 2 times the contig
				 * end length, because we want our distance samples
				 * to be based on a uniform head/tail length
				 */
				if ((int)l < 2 * params.end_length)
					continue;

				numSharedBarcodes[contigID]++;
			}
		}
	}

	/* build map of num shared barcodes => distance between head/tail */

	ofstream samplesOut;
	if (!params.dist_samples_tsv.empty()) {
		samplesOut.open(params.dist_samples_tsv.c_str());
		assert(samplesOut);
		samplesOut << "contig_id" << '\t'
			<< "shared_barcodes" << '\t'
			<< "distance" << '\n';
		assert(samplesOut);
	}

	ARCS::BarcodeToDist barcodeToDist;

	for (ContigToNumBarcodesIt it = numSharedBarcodes.begin();
		 it != numSharedBarcodes.end(); ++it)
	{
		std::string contigID = it->first;
        unsigned numBarcodes = it->second;

		assert(numBarcodes > 0);
		ARCS::BarcodesBinIndex binIndex =
			(numBarcodes - 1) / params.barcodes_bin_size;

		/* get length of contig */

        ARCS::ContigToLengthIt lengthIt = contigToLength.find(contigID);
		assert(lengthIt != contigToLength.end());
		unsigned l = lengthIt->second;
		assert((int)l >= 2 * params.end_length);

		/* record distance sample */

	    unsigned d = l - 2 * params.end_length;
		barcodeToDist[binIndex].push_back(d);

		if (!params.dist_stats_tsv.empty()) {
			samplesOut << contigID << '\t'
				<< numBarcodes << '\t'
				<< d << '\n';
			assert(samplesOut);
		}
	}

	if (!params.dist_samples_tsv.empty()) {
		samplesOut.close();
	}

	/* calc distance median and IQR for each barcode bin */

	ofstream statsOut;
	if (!params.dist_stats_tsv.empty()) {
		statsOut.open(params.dist_stats_tsv.c_str());
		assert(statsOut);
		statsOut << "barcodes" << '\t'
			<< "q1" << '\t'
			<< "q2" << '\t'
			<< "q3" << '\t'
			<< "n" << '\n';
		assert(statsOut);
	}

	for (ARCS::BarcodeToDistIt it = barcodeToDist.begin();
		it != barcodeToDist.end(); ++it)
	{
		unsigned binIndex = it->first;
		ARCS::DistanceSamples& distSamples = it->second;

		ARCS::DistStats stats;
 		stats.q2 = median(distSamples.begin(), distSamples.end());
		boost::tie(stats.q1, stats.q3) =
			IQR(distSamples.begin(), distSamples.end());
		stats.n = distSamples.size();

		barcodeToDistStats[binIndex] = stats;

		if (!params.dist_stats_tsv.empty()) {
			statsOut << '['
				<< binIndex * params.barcodes_bin_size + 1
				<< ','
				<< (binIndex + 1) * params.barcodes_bin_size
				<< ']' << '\t'
				<< stats.q1 << '\t'
				<< stats.q2 << '\t'
				<< stats.q3 << '\t'
				<< stats.n << '\n';
			assert(statsOut);
		}
	}

	if (!params.dist_stats_tsv.empty()) {
		statsOut.close();
	}
}


#endif
