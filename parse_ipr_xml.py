#!/mnt/shared/scratch/tadams/apps/conda/bin/python

# Script to parse output of nucleotide interproscan xml and extract an AA fasta
# and a bed file of locations of requested features
# WARNING!!!! This parses an XML file produced by interproscan. This assumes
# interproscan version 5.x-x.x is being used with XML schema version 4.5.
# This was the most recent at the time, but if your version number is HIGHER
# than 5.54-87.0 you'll need to check if the namespace heading is different and
# if the schema has changed. This only looks for hits found by hmme3, if you
# need ones called by a different analysis you will need to change the script.
# This looks to be very memory inefficient, possible TODO: improve this.

# Import python  modules

from collections import defaultdict
import xml.etree.ElementTree as ET
import argparse

# Prepare function to load CLI arguments


def parse_args():
    parser = argparse.ArgumentParser(description='Extract AA fasta and BED of \
    requested features')
    parser.add_argument('--xml', required=True, help='XML interproscan output \
    to be parsed')
    parser.add_argument('--feature', required=True, help='Interproscan feature \
    to be extracted to BED')
    parser.add_argument('--pep_out', required=True, help='Location to write \
    pep fasta file to')
    parser.add_argument('--bed_out', required=True, help='Location to write \
    bed file of feature locations to')
    return parser.parse_args()

# Prepare function to create data structures


def create_structures(xml, feature):
    AA_dict = defaultdict(list)
    Bed_dict = defaultdict(list)
    tree = ET.parse(xml)  # Get tree structure of XML file
    root = tree.getroot()  # Get data from XML tree root level
    ns = {'xmlns':
          'http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5'}
# This aids with parsing entries downstream
# There might be a way to use full paths, but I couldn't get it to work
# Possible TODO: Clean up searches with XPATH commands
    for query in root.findall('xmlns:nucleotide-sequence', ns):
        for orf in query.findall('xmlns:orf', ns):
            protein = orf.find('xmlns:protein', ns)
            sequence = protein.find('xmlns:sequence', ns)
            sequence_string = sequence.text
            name_list = protein.findall('xmlns:xref', ns)
            for name_entry in name_list:
                name_dict = name_entry.attrib
                orf_id = name_dict['id']
                query_name = name_dict['name'].split()[1].split('=')[1]
                orf_name = query_name + '_' + orf_id
                AA_dict[orf_name] = [orf_name, sequence_string]
                matches = protein.find('xmlns:matches', ns)
                for hmmer_hit in matches.findall('xmlns:hmmer3-match', ns):
                    signature = hmmer_hit.find('xmlns:signature', ns)
                    entry = signature.find('xmlns:entry', ns)
                    if entry is not None:
                        entry_dict = entry.attrib
                        entry_accession = entry_dict['ac']
                        if entry_accession == feature:
                            location = hmmer_hit.find('xmlns:locations', ns)
                            hmm = location.findall('xmlns:hmmer3-location', ns)
                            count = 0
                            for hit in hmm:
                                location_dict = hit.attrib
                                start = int(location_dict['start']) - 1
                                end = location_dict['end']
                                entry = orf_name + '_' + str(count)
                                count += 1
                                Bed_dict[entry] = [orf_name, str(start), end]
    return(AA_dict, Bed_dict)

# Prepare function to write output files


def write_files(xml, feature, pep_out, bed_out):
    dictionary_list = create_structures(xml, feature)
    AA_dict = dictionary_list[0]
    Bed_dict = dictionary_list[1]
    with open(pep_out, 'w') as po:
        for entry in AA_dict.keys():
            orf_name = AA_dict[entry][0]
            fasta_header = '>' + orf_name
            sequence = AA_dict[entry][1]
            po.write(fasta_header)
            po.write('\n')
            po.write(sequence)
            po.write('\n')
        po.close()
    with open(bed_out, 'w') as bo:
        for dict_key in Bed_dict.keys():
            list_to_write = Bed_dict[dict_key]
            string_to_write = '\t'.join(list_to_write)
            bo.write(string_to_write)
            bo.write('\n')
        bo.close()

# Preapre main function


def main():
    args = parse_args()
    xml = args.xml
    feature = args.feature
    pep_out = args.pep_out
    bed_out = args.bed_out
    write_files(xml, feature, pep_out, bed_out)


if __name__ == '__main__':
    main()
