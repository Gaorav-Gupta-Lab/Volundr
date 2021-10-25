#!python
#cython: language_level=3

import collections
from Levenshtein import distance
from natsort import natsort

cdef dict temp_dict = {}
cdef str fq1_name, fq1_seq, fq1_qual, marker_read, umi_fq, umi_anchor
# data_dict = {}
data_dict = collections.defaultdict(lambda: collections.defaultdict(list))

cpdef object file_writer(object self, fq1_batch):

    for fq1_read in fq1_batch:
        fq1_name = fq1_read[0]
        fq1_seq = fq1_read[1]
        fq1_qual = fq1_read[3]

        # Apply Filters
        min_length = self.args.MinimumReadLength

        sample_index, index_mismatch = index_search(self.master_index_dict, fq1_name, temp_dict)
        marker_read = ""
        if not self.sample_manifest_dictionary[sample_index] and sample_index != "Unknown":
            sample_index = "GhostIndex"

        if sample_index is not "Unknown" and sample_index is not "GhostIndex":
            marker_read = self.sample_manifest_dictionary[sample_index][1]

        # Filter reads based on length and number of N's.
        filtered = False
        if len(fq1_seq) < min_length or fq1_seq.count("N") / len(fq1_seq) >= self.args.N_Limit:
            data_dict[sample_index]["QC"].append("Filtered")
            filtered = True

        if not filtered:
            if sample_index in data_dict:
                data_dict[sample_index]["R1"].append("@{}\n{}\n+\n{}\n".format(fq1_name, fq1_seq, fq1_qual))
                if sample_index is not "Unknown":
                    data_dict[sample_index]["QC"][index_mismatch] += 1
                # data_dict[sample_index]["R1"].append("@{}\n{}\n+\n{}\n".format(fq1_name, fq1_seq, fq1_qual))

            else:
                data_dict[sample_index]["R1"] = ["@{}\n{}\n+\n{}\n".format(fq1_name, fq1_seq, fq1_qual)]
                data_dict[sample_index]["QC"] = [0, 0, 0]
                if sample_index is not "Unknown":
                    data_dict[sample_index]["QC"][index_mismatch] += 1

            # data_dict[sample_index]["R1"] = ["@{}\n{}\n+\n{}\n".format(fq1_name, fq1_seq, fq1_qual)]

    return data_dict


cdef umi_search(umi_anchor, seq):
    rt_pos = len(umi_anchor)
    lft_pos = 0
    umi = "No_UMI"

    while rt_pos < 26:
        query = seq[lft_pos:rt_pos]
        query_mismatch = distance(query, umi_anchor)
        tmp_mismatch = ""
        if query_mismatch <= 1:
            umi = seq[lft_pos+2:rt_pos+2]
            if query_mismatch == 0:
                return umi, rt_pos+2
            else:
                tmp_mismatch = query_mismatch
        lft_pos += 1
        rt_pos += 1

    return umi, rt_pos+1


cdef index_search(master_index_dict, fastq_name, temp_dict):
    def match_maker(query, unknown):
        """
        This little ditty gives us some wiggle room in identifying our indices and any other small targets.
        :param query
        :param unknown
        :return:
        """
        query_mismatch = distance(query, unknown)

        # Unknown length can be longer than target length.  Need to adjust mismatch index to reflect this.
        adjusted_query_mismatch = query_mismatch - (len(unknown) - len(query))

        return adjusted_query_mismatch

    left_query = fastq_name.split(":")[-1].split("+")[0]
    right_query = fastq_name.split(":")[-1].split("+")[1]

    # An iSeq100 can give this error.  Speeds up processing.
    if left_query == "TTTTTTTT" or right_query == "TTTTTTTT":
        return "Unknown", ""

    for index_key in master_index_dict:
        left_index = master_index_dict[index_key][0]
        right_index = master_index_dict[index_key][1]

        left_match = match_maker(left_index, left_query)
        right_match = match_maker(right_index, right_query)

        if left_match == 0 == right_match:
            return index_key, 0
        elif left_match <= 1 and right_match <= 1:
            temp_dict[left_match, right_match] = index_key

    if temp_dict:
        natsort.natsorted(temp_dict)
        return next(iter(temp_dict.values())), 1

    else:
        return "Unknown", ""