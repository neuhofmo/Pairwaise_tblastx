#! /usr/bin/env python2

import argparse
import os
import subprocess
from time import strftime, sleep

from Bio import SeqIO

from RecBlastUtils import debug_s, create_folder_if_needed, write_run_script, exists_not_empty

# THIS PROGRAM RECEIVES A FOLDER CONTAINING FASTA FILES,
# PERFORMS TBLASTX BETWEEN DIFFERENT COMBINATIONS OF THE FILES (NON-REDUNDANT)
# THEN FILTERS THE RESULTS ACCORDING TO SPECIFIC PARAMETERS
# AND COMBINES THE HOMOLOGUE SEQUENCES

__version__ = '0.1'

# for efficiency
strip = str.strip
split = str.split
join_folder = os.path.join


def write_and_run_script(command, script_name, check_call=False):
    """"Write the script to the script_name path and run it using subprocess."""
    try:
        script_path = write_run_script(command, script_name)
        if check_call:
            return subprocess.check_call(script_path)
        else:
            return strip(subprocess.check_output(script_path, universal_newlines=True))  #
        # subprocess.check_call(script_path)   # switched to check_output
    except subprocess.CalledProcessError, e:  # restarting the process (with a little sleep period)
        debug("Had a little problem with running this command: "
              "{}\nSo we are running it again.\nThis was the error:".format(command, e))
        sleep(10)
        script_path = write_run_script(command, script_name)
        sleep(20)
        # subprocess.check_call(script_path)  # switched to check_output
        if check_call:
            return subprocess.check_call(script_path)
        else:
            return strip(subprocess.check_output(script_path, universal_newlines=True))  #


def run_tblastx(query_path, db_path, outfmt, e_value_thresh, coverage_threshold, cpu,
                blast_output_file):
    """
    This function runs tblastx.
    :param query_path: a fasta query file path
    :param db_path: ncbi blast+ db file path
    :param outfmt: out format (outfmt) for the blast results
    :param e_value_thresh: e_value threshold
    :param coverage_threshold: coverage_threshold
    :param cpu: number CPUs to run blast on
    :param blast_output_file: output file path
    :return:
    """
    command_line = "{0} -query {1} -db {2} -outfmt '{3}' -evalue {4} -max_hsps 1 -qcov_hsp_perc {5} -num_threads {6} " \
                   "-out {7}\n".format(TBLASTX_PATH, query_path, db_path, outfmt, e_value_thresh, coverage_threshold,
                                       cpu, blast_output_file)
    debug("Running the following line:\n{}".format(command_line))
    # writing the run script
    try:
        script_path = write_run_script(command_line)
        subprocess.check_call(script_path)
    except subprocess.CalledProcessError, e:  # restarting the process (with a little sleep period)
        debug("Had a little problem with running this command: "
              "{}\nSo we are running it again.\nThis was the error:".format(command_line, e))
        sleep(10)
        script_path = write_run_script(command_line)
        sleep(20)
        subprocess.check_call(script_path)
    return True


def merge_dicts_by_values(target_dict):
    """
    Receives a dict and merges values and keys.
    Then returns a dictionary with a uniq list of values.
    :param target_dict:
    :return:
    """
    # TODO: make more efficient if possible
    # if key(1) is in values of another dict(2),
    # add the values of key(1) to values of dict(2)

    all_keys = target_dict.keys()
    debug("Running on {} sequences:".format(len(all_keys)))
    for key in all_keys:  # going over every key in the dictionary
        for value_tup in target_dict[key]:  # in every key, going over the tuples in the value list
            if value_tup in all_keys:  # if the value tuple is also a key in the same dict
                debug("{} is already in the dictionary, merging.".format(value_tup))
                # add the value list of the key to the original
                target_dict[value_tup].extend(target_dict[key].remove(value_tup))

    debug("Ran and merged all keys.")
    # turning the values in the dictionary to unique
    for key in target_dict:
        target_dict[key] = list(set(target_dict[key]))

    return target_dict


def make_blast_db(fasta_path, makeblastdb_path='/usr/local/bin/makeblastdb'):
    """
    This function receives a path to a fasta file, and makes nucleotide blast db out of it.
    It fetches the taxid and organism name from the title (assuming all files have it,
    which is what we are trying to assume here.
    :param fasta_path:
    :param makeblastdb_path:
    :return:
    """
    # get taxid:
    taxid_script = '/tmp/get_taxid.sh'
    taxid_command = "cat {0} | grep ORGANISM | cut -d'=' -f3 | cut -d' ' -f1 | sort | uniq\n".format(
        fasta_path)  # according to our headers
    taxid = write_and_run_script(taxid_command, taxid_script)
    debug("Got taxid: {}".format(taxid))
    # get organism name:
    organism_script = '/tmp/get_organism.sh'
    organism_command = "cat {0} | grep ORGANISM | cut -d'=' -f4 | cut -d'\"' -f2 | sort | uniq\n".format(
        fasta_path)  # according to our headers
    organism = write_and_run_script(organism_command, organism_script)
    debug("Got organism name: {}".format(organism))
    # run makeblastdb command
    makeblastdb_command = '{0} -in {1} -dbtype nucl -taxid {2}\n'.format(makeblastdb_path, fasta_path, taxid)
    makeblastdb_script = '/tmp/run_makeblastdb.sh'
    if write_and_run_script(makeblastdb_command, makeblastdb_script, check_call=True):  # the building run properly
        return taxid, organism


######################
# Default parameters #
######################

E_VALUE_THRESH = 1e-4
COVERAGE_THRESHOLD = 10  # TODO change
IDENTITY_THRESHOLD = 15  # TODO: DO we need that?
MATCH_LENGTH = 40  # TODO change and check (it's 40aa, but not sure how it looks like in tbalstx)
MAX_TARGET_SEQS = '50'
# OUTFMT = '6 staxids sseqid pident qcovs evalue length sscinames sblastnames'  # no need for names
OUTFMT = '6 sseqid sseqid pident qcovs evalue length'
BLAST_SUF = ['nhr', 'nin', 'nsq']

####################
#  parse arguments #
####################

parser = argparse.ArgumentParser()
parser.add_argument("fasta_folder", help="Folder containing the fasta nt files you want to blast.")
parser.add_argument("--num_threads", help="The number of threads (CPU) dedicated for parallel blast run.",
                    default=1, type=int)
parser.add_argument("--evalue", help="The e-value threshold for matches of the first blast.", default=E_VALUE_THRESH)
parser.add_argument("--identity", help="The minimum identity required for blast matches.",
                    default=IDENTITY_THRESHOLD)
parser.add_argument("--coverage", help="The minimum query and hit coverage required for blast matches.",
                    default=COVERAGE_THRESHOLD)
parser.add_argument("--max_seqs", help="The maximum number of sequences reported by blast.", default=MAX_TARGET_SEQS)
parser.add_argument("--keep_files", help="Keeps intermediate files after completion", action="store_true",
                    default=False)
parser.add_argument("--skip_blast", help="Don't run blast, in case you want to continue an old run",
                    action="store_true", default=False)
parser.add_argument("--debug", help="Adds debug prints in various stages of the run.", action="store_true",
                    default=False)
parser.add_argument("-v", "--version", help="Prints version in formation.", action="store_true",
                    default=False)

args = parser.parse_args()

# DEBUG flags
DEBUG = args.debug

if args.version:
    print __version__


def debug(s):
    """Wrapper for the debug function. Debug prints if DEBUG is True."""
    return debug_s(s, DEBUG)

fasta_folder = args.fasta_folder  # input folder

# Parsing input arguments and local parameters:
# locating TBLASTX path on your system
TBLASTX_PATH = "Not valid"
try:
    TBLASTX_PATH = strip(subprocess.check_output(["which", "tblastx"], universal_newlines=True))
    debug("TBLASTX found in {}".format(TBLASTX_PATH))
except subprocess.CalledProcessError:
    print("No TBLASTX found. Please check install blast properly or make sure it's in $PATH. Aborting.")
    exit(1)

CPU = args.num_threads

# TODO: use argparse arguments to change input fixed parameters

#######################
#  Preparing folders  #
#######################

# creating output folder:
output_fasta_folder = join_folder(fasta_folder, 'fasta_output')
create_folder_if_needed(output_fasta_folder)

# go over fasta files in folder
only_fasta_files = [f for f in os.listdir(fasta_folder) if os.path.isfile(join_folder(fasta_folder, f)) and
                    os.path.splitext(f)[1] in ['.fna', '.fasta', '.fa']]
# debug(only_fasta_files)  # DEBUG
# exit(1)  # DEBUG

# making blastdb if doesn't exist already:
for fasta_file in only_fasta_files:
    full_path = join_folder(fasta_folder, fasta_file)
    if all(os.path.exists(".".join([full_path, x])) for x in BLAST_SUF):
        debug("Blast DB for the file {} already exists.".format(full_path))
    else:  # if db file doesnt exist..:
        debug("Blast DB for file {} does not exist, creating it.".format(full_path))
        if make_blast_db(full_path):
            debug("Blast DB made for file {}".format(full_path))
debug("Made BLASTDBs where needed.")

blast_out_files = {}  # initialize

# blast a-b, a-c... b-c, b-d... c-d, etc. Use standard parameters.
num_fastas = len(only_fasta_files)
debug("Working on {} files:".format(num_fastas))

for i in range(num_fastas):  # first loop - query file
    for j in range(i + 1, num_fastas):  # 2nd loop - db file
        query_db_tuple = (only_fasta_files[i], only_fasta_files[j])  # files involved
        query_base = only_fasta_files[i]
        db_base = only_fasta_files[j]
        # blast output name for this run
        blast_out_fname = join_folder(output_fasta_folder, "{0}_VS_{1}.blast.out".format(query_base, db_base))

        # connect blast_out_name to the files it involves
        blast_out_files[blast_out_fname] = query_db_tuple

        # run tblastx
        full_path_query = join_folder(fasta_folder, query_db_tuple[0])
        full_path_db = join_folder(fasta_folder, query_db_tuple[1])
        if args.skip_blast:
            print "Skipping blast run. Using output from {}".format(blast_out_fname)
        elif exists_not_empty(blast_out_fname):
            print "Output file exists and not empty. Skipping blast run. Using output from {}".format(blast_out_fname)
        else:
            if run_tblastx(full_path_query, full_path_db, OUTFMT, E_VALUE_THRESH, COVERAGE_THRESHOLD, CPU,
                           blast_out_fname):
                debug("Ran tblastx on query {0} and db {1}.\nOutput in {2}.".format(query_base, db_base,
                                                                                    blast_out_fname))

# At this point we should have ran all non redundant combinations.
# when it's done:

###########################
#  Analyze blast results  #
###########################

homolog_dict = {}

# go over output files
debug("Analyzing output files:")
for blast_file in blast_out_files:
    debug("Working on file {}".format(blast_file))
    # filter results - above evalue threshold, length...
    try:
        with open(blast_file, 'r') as infile:  # open file for reading
            for line in infile:  # read every line
                # parse every line and take relevant information
                query_seq_id, target_seq_id, identity, coverage, evalue, length = split(strip(line))
                debug("query_seq_id: {0}; target_seq_id: {1}; identity: {2}; "
                      "coverage: {3}; evalue: {4}; length: {5}".format(query_seq_id, target_seq_id, identity, coverage,
                                                                       evalue, length))
                query_file = blast_out_files[blast_file][0]
                target_file = blast_out_files[blast_file][1]
                # if sequence is a good match:
                if float(evalue) < E_VALUE_THRESH \
                        and float(identity) >= IDENTITY_THRESHOLD \
                        and float(coverage) >= COVERAGE_THRESHOLD \
                        and int(length) >= MATCH_LENGTH:
                    debug("**PASSED**")
                    try:
                        homolog_dict[(query_file, query_seq_id)].append((target_file, target_seq_id))  # adding the file
                    except KeyError:  # the key doesn't exist
                        homolog_dict[(query_file, query_seq_id)] = [(target_file, target_seq_id)]  # start a new list

    except IOError:
        debug("No results at all for file {}".format(blast_file))

# merging values and keys in the dict to remove redundancy
merged_homolog_dict = merge_dicts_by_values(homolog_dict)
debug("Merged dictionaries.")

debug("Now writing output fastas to {}".format(output_fasta_folder))
for seq in merged_homolog_dict:
    # open a fasta file for that seq
    output_fasta_file = join_folder(output_fasta_folder, "{0}__{1}_homologues.fasta".format(seq[0], seq[1]))

    # a list of all the files that contain sequences we need
    get_files = seq[0] + list(set([x[0] for x in merged_homolog_dict[seq]]))  # all files

    # fetching the file from
    with open(output_fasta_file, 'wa') as outfile:
        for filename in get_files:
            seqs_to_fetch_by_filename = [x[1] for x in merged_homolog_dict[seq] if x[0] == filename]
            seq_iterator = SeqIO.parse(open(join_folder(fasta_folder, filename), 'r'), 'fasta')
            SeqIO.write((seq for seq in seq_iterator if seq.id in seqs_to_fetch_by_filename), outfile, "fasta")
        debug("Done writing sequences to output file {}".format(output_fasta_file))
    # sorted(merged_homolog_dict[seq], key=lambda x: x[0])  # no need

print "Done writing all sequences."
print "Program done at {}".format(strftime('%H:%M:%S'))
exit(0)

# that's it
