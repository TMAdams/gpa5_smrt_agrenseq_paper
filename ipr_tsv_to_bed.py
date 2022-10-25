#!/mnt/shared/scratch/tadams/apps/conda/bin/python

# Script to create bed file from interpro tsv

# Import python modules

import argparse
from collections import defaultdict

# Prepare function to parse CLI arguments


def parse_args():
    parser = argparse.ArgumentParser(description='Convert IPR tsv to BED for \
    specific feature predictions')
    parser.add_argument('--tsv', required=True, help='interpro tsv file')
    parser.add_argument('--feature', required=True, help='Interproscan feature \
    to be extracted to BED')
    parser.add_argument('--bed_out', required=True, help='Location to write \
    bed file of feature locations to')
    return parser.parse_args()

# Prepare function to load data into a dictionary


def parse_data(tsv, feature):
    location_dict = defaultdict(list)
    previous_gene = ''
    count = 0
    for line in tsv:
        line = line.rstrip()
        split_line = line.split('\t')
        ipr_feat = split_line[11]
        if ipr_feat == feature:
            gene = split_line[0]
            if gene == previous_gene:
                count += 1
                key = gene + '_' + str(count)
            else:
                key = gene
                count = 0
            start = int(split_line[6]) - 1
            stop = split_line[7]
            previous_gene = gene
            location_dict[key] = [gene, str(start), stop]
    return location_dict

# Prepare function to write output file


def write_file(tsv, feature, bed_out):
    location_dict = parse_data(tsv, feature)
    with open(bed_out, 'w') as o:
        for key in location_dict.keys():
            list_to_write = location_dict[key]
            string_to_write = '\t'.join(list_to_write)
            o.write(string_to_write)
            o.write('\n')
    o.close()

# Prepare main function


def main():
    args = parse_args()
    tsv = open(args.tsv).readlines()
    feature = args.feature
    bed_out = args.bed_out
    write_file(tsv, feature, bed_out)


if __name__ == '__main__':
    main()
