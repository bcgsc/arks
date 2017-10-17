#ifndef _DISTANCE_EST_H_
#define _DISTANCE_EST_H_ 1

#include "Arcs/Arcs.h"
#include "Common/StatUtil.h"
#include <limits>
#include <iostream>

/**
 * Records the distance between the head/tail regions of the same
 * contig vs. barcode union size, barcode intersection size,
 * and number of distinct barcodes mapped to each end.
 */
struct DistSample
{
	unsigned distance;
	unsigned barcodesHead;
	unsigned barcodesTail;
	unsigned barcodesUnion;
	unsigned barcodesIntersect;

	DistSample() :
		distance(std::numeric_limits<unsigned>::max()),
		barcodesHead(0),
		barcodesTail(0),
		barcodesUnion(0),
		barcodesIntersect(0)
	{}
};

typedef std::unordered_map<std::string, DistSample> DistSampleMap;
typedef typename DistSampleMap::const_iterator DistSampleConstIt;

typedef std::map<double, DistSample> JaccardToDistSample;
typedef typename JaccardToDistSample::const_iterator JaccardToDistSampleConstIt;

/**
 * Measure distance between contig ends vs.
 * barcode intersection size and barcode union size.
 */
void calcDistSamples(const ARCS::IndexMap& imap,
	const ARCS::ContigToLength& contigToLength,
	const std::unordered_map<std::string, int>& indexMultMap,
	const ARCS::ArcsParams& params,
	DistSampleMap& distSamples)
{
	/* for each chromium barcode */
	for (auto barcodeIt = imap.begin(); barcodeIt != imap.end();
		++barcodeIt)
	{
		/* skip barcodes outside of min/max multiplicity range */
		std::string index = barcodeIt->first;
		int indexMult = indexMultMap.at(index);
		if (indexMult < params.min_mult || indexMult > params.max_mult)
			continue;

		/* contig head/tail => number of mapped read pairs */
		const ARCS::ScafMap& contigToCount = barcodeIt->second;

		for (auto contigIt = contigToCount.begin();
			contigIt != contigToCount.end(); ++contigIt)
		{
			std::string contigID;
			bool isHead;
			std::tie(contigID, isHead) = contigIt->first;
			int readPairs = contigIt->second;

			/*
			 * skip contigs with less than required number of
			 * mapped read pairs (-c option)
			 */
			if (readPairs < params.min_reads)
				continue;

			/*
			 * skip contigs shorter than 2 times the contig
			 * end length, because we want our distance samples
			 * to be based on a uniform head/tail length
			 */

			unsigned l = contigToLength.at(contigID);
			if (l < (unsigned) 2 * params.end_length)
				continue;

			DistSample& distSample = distSamples[contigID];
			distSample.distance = l - 2 * params.end_length;

			if (isHead)
				distSample.barcodesHead++;
			else
				distSample.barcodesTail++;

			/*
			 * Check if barcode also maps to other end of contig
			 * with sufficient number of read pairs.
			 *
			 * The `isHead` part of the `if` condition prevents
			 * double-counting when a barcode maps to both
			 * ends of a contig.
			 */

			ARCS::CI otherEnd(contigID, !isHead);
			ARCS::ScafMapConstIt otherIt = contigToCount.find(otherEnd);
			bool foundOther = otherIt != contigToCount.end()
				&& otherIt->second >= params.min_reads;

			if (foundOther && isHead) {
				distSample.barcodesIntersect++;
				distSample.barcodesUnion++;
			} else if (!foundOther) {
				distSample.barcodesUnion++;
			}
		}
	}
}

/**
 * Build a ordered map from barcode Jaccard index to
 * distance sample. Each distance sample comes from
 * measuring the distance between the head/tail of the
 * same contig, along with associated head/tail barcode
 * counts.
 */
static inline void buildJaccardToDistSamples(
	const DistSampleMap& distSamples,
	JaccardToDistSample& jaccardToSamples)
{
	for (DistSampleConstIt it = distSamples.begin();
		it != distSamples.end(); ++it)
	{
		const DistSample& sample = it->second;
		double jaccard = double(sample.barcodesIntersect)
			/ sample.barcodesUnion;
		jaccardToSamples.insert(
			JaccardToDistSample::value_type(jaccard, sample));
	}
}

/**
 * Write distance samples to an output stream.  The distance
 * samples record the distance between the head and tail regions
 * of the same contig with associated barcode stats (e.g.
 * barcode intersection size).
 */
static inline std::ostream& writeDistSamples(std::ostream& out,
	const DistSampleMap& distSamples)
{
	out << "contig_id" << '\t'
		<< "distance" << '\t'
		<< "barcodes_head" << '\t'
		<< "barcodes_tail" << '\t'
		<< "barcodes_union" << '\t'
		<< "barcodes_intersect" << '\n';

	for (DistSampleConstIt it = distSamples.begin();
		it != distSamples.end(); ++it)
	{
		const std::string& contigID = it->first;
		const DistSample& sample = it->second;

		out << contigID << '\t'
			<< sample.distance << '\t'
			<< sample.barcodesHead << '\t'
			<< sample.barcodesTail << '\t'
			<< sample.barcodesUnion << '\t'
			<< sample.barcodesIntersect << '\n';
	}

	return out;
}

#endif
