#ifndef _STAT_UTIL_H_
#define _STAT_UTIL_H_ 1

#include <algorithm>

/** Calculate the median of an unsorted container */
template <typename IteratorT>
static inline float median(const IteratorT& it1,
	const IteratorT& it2)
{
	size_t size = it2 - it1;
	if (size == 0)
		return 0.0f;

	if (size % 2 == 0) {
		/* return average of middle two elements */
		IteratorT mid1 = it1 + size / 2 - 1;
		IteratorT mid2 = it1 + size / 2;
		std::nth_element(it1, mid1, it2);
		std::nth_element(it1, mid2, it2);
		return float(*mid1 + *mid2) / 2.0f;
	} else {
		/* return middle element */
		IteratorT mid = it1 + size / 2;
		std::nth_element(it1, mid, it2);
		return float(*mid);
	}
}

/** Calculate the inter-quartile range of an unsorted container */
template <typename IteratorT>
static inline
std::pair<typename IteratorT::value_type, typename IteratorT::value_type>
IQR(const IteratorT& it1, const IteratorT& it2)
{
	typedef typename IteratorT::value_type value_type;

	size_t size = it2 - it1;
	if (size == 0)
		return std::make_pair(0, 0);
	if (size == 1)
		return std::make_pair(*it1, *it1);

	value_type q1, q3;
	IteratorT mid = it1 + size / 2;

	if (size % 2 == 0) {
		q1 = median(it1, mid);
		q3 = median(mid, it2);
	} else {
		q1 = median(it1, mid + 1);
		q3 = median(mid, it2);
	}

	return std::make_pair(q1, q3);
}

#endif
