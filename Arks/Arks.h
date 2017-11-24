#ifndef ARCS_H
#define ARCS_H 1

#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <iomanip>
#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <limits>
#include <utility> 
#include <vector>
#include <iterator>
#include <time.h> 
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/graphviz.hpp>
//#include "Common/Uncompress.h"
#include "DataLayer/FastaReader.h"
#include "DataLayer/FastaReader.cpp"
#include "Common/ReadsProcessor.h"
// using sparse hash maps for k-merization
#include <google/sparse_hash_map>
#include "city.h"

namespace ARCS {

/**
 * Parameters controlling ARCS run
 */
struct ArcsParams {

	std::string program; 
	std::string file;
	std::string multfile;
	std::string conrecfile; 
	std::string kmapfile; 
	std::string imapfile; 
	int checkpoint_outs;
	int min_reads;
	int k_value;
	int k_shift;
	double j_index;
	int min_links;
	unsigned min_size;
	std::string base_name;
	int min_mult;
	int max_mult;
	int max_degree;
	int end_length;
	float error_percent;
	int verbose;
	unsigned threads;
	bool distance_est;
	std::string intra_contig_tsv;
	std::string inter_contig_tsv;
	unsigned dist_bin_size;

	ArcsParams() :
			program(), file(), multfile(), conrecfile(), kmapfile(), imapfile(), checkpoint_outs(0), min_reads(5), k_value(
					30), k_shift(1), j_index(0.55), min_links(0), min_size(500), base_name(
					""), min_mult(50), max_mult(10000), max_degree(0), end_length(
					0), error_percent(0.05), verbose(0), threads(1), distance_est(false), dist_bin_size(20) {
	}

};

/* SIMPLIFYING VARIABLES: */
typedef const char* Kmer;
typedef std::pair<std::string, bool> CI;

/* MAP DATA STRUCTURES: */

/* ContigKMap: <k-mer, pair(contig id, bool), hash<k-mer>, eqstr>
 * 	k-mer = string sequence
 *  contig id = string
 *  bool = True for Head; False for Tail
 *  eqstr = equal key
 */

// simple hash adapter for types without pointers
template<typename T>
struct CityHasher {
	size_t operator()(const T& t) const {
		return CityHash64(t, sizeof(t));
	}
};

// specialization for strings
template<>
struct CityHasher<std::string> {
	size_t operator()(const string t) const {
		return CityHash64(t.c_str(), t.size());
	}
};

struct eqstr {
	bool operator()(std::string s1, std::string s2) const {
		return (s1 == s2);
	}
};

typedef google::sparse_hash_map<std::string, int, CityHasher<std::string>, eqstr> ContigKMap;

/* ScafMap: <pair(scaffold id, bool), count>, cout =  # times index maps to scaffold (c), bool = true-head, false-tail*/
typedef std::map<CI, int> ScafMap;
typedef typename ScafMap::const_iterator ScafMapConstIt;

/* IndexMap: key = index sequence, value = ScafMap */
typedef std::unordered_map<std::string, ScafMap> IndexMap;

/** a pair of contig IDs */
typedef std::pair<std::string, std::string> ContigPair;

/* PairMap: key = pair of scaf sequence id, value = num links*/
typedef std::map<std::pair<std::string, std::string>, std::vector<int>> PairMap;
typedef typename PairMap::iterator PairMapIt;

/** maps contig FASTA ID to contig length (bp) */
typedef std::unordered_map<std::string, int> ContigToLength;
typedef typename ContigToLength::const_iterator ContigToLengthIt;

/* GRAPH DATA STRUCTURES: */

struct VertexProperties {
	std::string id;
};

/* Orientation: 0-HH, 1-HT, 2-TH, 3-TT */
struct EdgeProperties {
	int orientation;
	int weight;
	int minDist;
	int maxDist;
	float jaccard;
	EdgeProperties() :
		orientation(0), weight(0),
		minDist(std::numeric_limits<int>::min()),
		maxDist(std::numeric_limits<int>::max()),
		jaccard(-1.0f)
	{}
};

template <class GraphT>
struct EdgePropertyWriter
{
	typedef typename boost::graph_traits<GraphT>::edge_descriptor E;
	typedef typename boost::edge_property<GraphT>::type EP;

	GraphT& m_g;

	EdgePropertyWriter(GraphT& g) : m_g(g) {}

	void operator()(std::ostream& out, const E& e) const
	{
		EP ep = m_g[e];
		out << '['
			<< "label=" << ep.orientation << ','
			<< "weight=" << ep.weight;

		if (ep.minDist != std::numeric_limits<int>::min()) {
			assert(ep.maxDist != std::numeric_limits<int>::max());
			assert(ep.jaccard >= 0.0f);
			out << ','
				<< "mind=" << ep.minDist << ','
				<< "maxd=" << ep.maxDist << ','
				<< std::fixed << std::setprecision(2)
				<< "j=" << ep.jaccard << '\n';
		}
		out << ']';
	}
};

template <class GraphT>
struct VertexPropertyWriter
{
	typedef typename boost::graph_traits<GraphT>::vertex_descriptor V;
	typedef typename boost::vertex_property<GraphT>::type VP;

	GraphT& m_g;

	VertexPropertyWriter(GraphT& g) : m_g(g) {}
	void operator()(std::ostream& out, const V& v) const
	{
		out << "[id=" << m_g[v].id << "]";
	}
};

typedef boost::undirected_graph<VertexProperties, EdgeProperties> Graph;
typedef std::unordered_map<std::string, Graph::vertex_descriptor> VidVdesMap;
typedef boost::graph_traits<ARCS::Graph>::vertex_descriptor VertexDes;
}

#endif
