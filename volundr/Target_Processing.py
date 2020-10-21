import csv
import pathlib
import natsort
import collections
from contextlib import suppress

__version__ = "0.1.0"


def process_target(o):
    """
    This was written to trim the PAM sequence from the sgRNA in a target file.
    :param o:
    :return:
    """
    target_data_file = list(csv.reader(open(o.target_data_file), delimiter="\t"))
    data_dict = collections.defaultdict(lambda: collections.defaultdict(float))
    line_count = 0
    header = ""

    for line in target_data_file:
        if line_count > 0:
            gene_name = line[0].split("_")[0]

            if gene_name in data_dict:
                counter = 0
                for v in line:
                    if counter > 0:
                        with suppress(ValueError):
                            with suppress(IndexError):
                                data_dict[gene_name][counter-1] += float(v)

                    counter += 1
            else:
                line_list = []
                for i in range(len(line)):
                    if i > 0:
                        try:
                            line_list.append(float(line[i]))
                        except ValueError:
                            pass
                data_dict[gene_name] = line_list

        elif line_count == 0:
            for i in range(len(line)):
                if i+1 < len(line):
                    header += line[i]+"\t"
                else:
                    header += line[i]+"\n"
            line_count += 1

    filepath = str(pathlib.PurePath(pathlib.PurePath(o.target_data_file).parent, ""))
    out_file = open(filepath+"outfile.txt", "w")
    out_file.write(header)

    for gene in natsort.natsorted(data_dict):
        temp_row = ""
        temp_row += gene + "\t"
        for i in range(len(data_dict[gene])):
            if i+1 < len(data_dict[gene]):
                temp_row += str(data_dict[gene][i]) + "\t"
            else:
                temp_row += str(data_dict[gene][i]) + "\n"
        out_file.write(temp_row)
