import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import sys, os
import gzip
import random
from multiprocessing import Pool
import argparse
import time
import subprocess
import pkg_resources
import tempfile
import pysam
from scipy import stats
from canvasPainter import paint_canvas
from pysamReader import pysamReader
from repeatAnnotation import repeats_annotation
from gnomadAnnotation import gnomad_annotation
from homopolymerDinucleotideAnnotation import homopolymer_dinucleotide_annotation

MAX_DP = 500
WIDTH = 300


def list_to_string(info_list):
    string = ""
    for item in info_list:
        string += item[1] + ":" + str(item[0]) + " "
    return string


def wilson_binom_interval(success, total, alpha = 0.05):
    q_ = success / total
    crit = stats.norm.isf( alpha / 2.)
    crit2 = crit**2
    denom = 1 + crit2 / total
    center = (q_ + crit2 / (2 * total)) / denom
    dist = crit * np.sqrt(q_ * (1. - q_) / total + crit2 / (4. * total**2))
    dist /= denom
    ci_low = center - dist
    ci_upp = center + dist
    return ci_low, ci_upp


def check_x_region(position):
    position = int(position)
    in_par1 = (position >= x_par1_region[0]) and (position <= x_par1_region[1])
    in_par2 = (position >= x_par2_region[0]) and (position <= x_par2_region[1])
    return (not in_par1) and (not in_par2)


def check_y_region(position):
    position = int(position)
    in_par1 = (position >= y_par1_region[0]) and (position <= y_par1_region[1])
    in_par2 = (position >= y_par2_region[0]) and (position <= y_par2_region[1])
    return (not in_par1) and (not in_par2)


def normalize_fasta_chrom(fasta, chrom):
    references = set(fasta.references)
    if chrom in references:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in references:
        return chrom[3:]
    if ("chr" + chrom) in references:
        return "chr" + chrom
    return chrom


def custom_homopolymer_dinucleotide_annotation(chrom, pos, reference_fasta):
    import re

    def check_if_in_homopolymer(seq_str):
        pattern = re.compile(r'([ACGT])\1{3,}')
        matches = [m.group() for m in re.finditer(pattern, seq_str)]
        if len(matches) > 0:
            is_homopolymer = 1
        else:
            is_homopolymer = 0
        return is_homopolymer

    def check_if_in_dinucleotide_repeat(seq_str):
        pattern = re.compile(r'([ACGT]{2})\1{3,}')
        matches = [m.group() for m in re.finditer(pattern, seq_str)]
        if len(matches) > 0:
            is_dinucleotide = 1
        else:
            is_dinucleotide = 0
        return is_dinucleotide

    pos = int(pos)
    fasta = pysam.FastaFile(reference_fasta)
    fasta_chrom = normalize_fasta_chrom(fasta, str(chrom))

    seq_str_9bp = fasta.fetch(fasta_chrom, max(pos - 5, 0), pos + 4).upper()
    seq_str_17bp = fasta.fetch(fasta_chrom, max(pos - 9, 0), pos + 8).upper()
    fasta.close()

    return check_if_in_homopolymer(seq_str_9bp), check_if_in_dinucleotide_repeat(seq_str_17bp)


def custom_repeats_annotation(all_variants, output_dir, repeat_bed, segdup_bed):
    rp_fd, rp_path = tempfile.mkstemp()
    try:
        with os.fdopen(rp_fd, 'w') as tmp:
            for variant in all_variants:
                sample_name, bam, chrom, pos, ref, alt, depth, sex = variant
                key = "_".join([chrom, pos, ref, alt])
                line = "\t".join(map(str, [chrom, int(pos)-1, int(pos) + len(ref)-2, ref, alt, key])) + "\n"
                tmp.write(line)

        command = "bedtools annotate -i " + rp_path + " -files " + repeat_bed + " " + segdup_bed + " > " + output_dir + "repeats_annotation.bed"
        subprocess.call(command, shell=True)
        os.remove(rp_path)
        df = pd.read_csv(output_dir + "repeats_annotation.bed", header=None, sep="\t", dtype={0: str})
        repeats_dict = dict(zip(df[5], zip(df[6], df[7])))
        return repeats_dict
    except:
        sys.stderr.write("Error with repeat annotation. Check if you have module loaded bedtools.\n")
        sys.exit(2)


def multiprocess_iterator(args):
    line, image_outdir, matrix_outdir, build, reference_fasta, apply_sex_chromosome_rules = args
    sample_name, bam, chrom, pos, ref, alt, sequencing_depth, sex, cram_ref_dir = line

    pysam_reader = pysamReader(bam, chrom, pos, cram_ref_dir, ref, alt)
    pysam_reader.downsample_to_max_depth()
    pysam_reader.build_reads_dict()

    if ref != None and alt != None:
        pysam_reader.rearrange_reads_ref_alt()
    elif ref == None and alt == None:
        pysam_reader.rearrange_reads_no_ref()
    elif ref != None and alt == None:
        pysam_reader.rearrange_reads_no_alt()

    reads, reads_count, base_info = pysam_reader.close()
    del pysam_reader
    canvas = paint_canvas(reads, int(pos))

    #compute depth fraction
    depth_fraction = reads_count/int(sequencing_depth)
    if apply_sex_chromosome_rules and sex == "M" and chrom == "X" and check_x_region(int(pos)):
        depth_fraction = depth_fraction * 2
    if apply_sex_chromosome_rules and sex == "M" and chrom == "Y" and check_y_region(int(pos)):
        depth_fraction = depth_fraction * 2

    #compute maf and binomial CI
    ref_count = base_info[0][0]
    alt_count = base_info[1][0]
    if ref_count + alt_count == 0:
        maf = 0
        lower_CI = 0
        upper_CI = 0
    else:
        maf = alt_count /(ref_count + alt_count)
        lower_CI, upper_CI = wilson_binom_interval(alt_count, ref_count + alt_count)

    #save images
    key =  "_".join(list(map(str,[chrom, pos, ref, alt])))
    filename = sample_name + "-" + key
    matrix_file = matrix_outdir + filename
    image_file = image_outdir + filename + ".jpg"

    #save matrix
    np.save(matrix_file, canvas)

    #save image
    canvas = canvas.astype(int)
    fig1 = plt.figure()
    plt.imshow(canvas)
    plt.title(chrom + ":" +  str(pos) + "  "+ list_to_string(base_info))
    fig1.savefig(image_file)
    plt.close(fig1)

    #check if homopolymer
    if build == "custom":
        is_homopolymer, is_dinucleotide = custom_homopolymer_dinucleotide_annotation(chrom, pos, reference_fasta)
    else:
        is_homopolymer, is_dinucleotide = homopolymer_dinucleotide_annotation(chrom, pos, build)

    return sample_name, sex, key, maf, lower_CI, upper_CI, depth_fraction, is_homopolymer, is_dinucleotide, image_file, matrix_file


def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-i", "--input_file", required=True, help="Input file (input.txt). [bam],[vcf],[sequencing_depth],[sex]")
    parser.add_argument("-f", "--vcf_filters", required=False, default=None, help="Filter the vcf file by INFO column, e.g. PASS. Default is no filtering")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory (output)")
    parser.add_argument("-a", "--annovar_path", required=False, help="Absolute path to the annovar package \
                                                                     (humandb directory should already be specified inside)")
    parser.add_argument("-db", "--dbtype", required=False, default="gnomad_genome", help="db file located in annovar directory,  this feeds directly into the annovar parameter --dbtype, default: gnomad_genome")
    parser.add_argument("-b", "--build", required=False, default="hg19", help="Version of genome build, options: hg19, hg38, custom")
    parser.add_argument("--reference-fasta", required=False, help="Reference FASTA file for custom genome build mode")
    parser.add_argument("--repeat-bed", required=False, help="Repeat annotation BED file for custom genome build mode")
    parser.add_argument("--segdup-bed", required=False, help="Segmental duplication BED file for custom genome build mode")
    parser.add_argument("--skip-annovar", action="store_true", help="Skip ANNOVAR/gnomAD annotation. Useful for custom or non-human assemblies")
    parser.add_argument("-c", "--cram_ref_dir", required=False, help="CRAM file reference path")
    options = parser.parse_args(args)
    return options


def main():
    since = time.time()
    options = getOptions(sys.argv[1:])

    global x_par1_region
    global y_par1_region
    global x_par2_region
    global y_par2_region

    apply_sex_chromosome_rules = True

    if options.build == 'hg19':
        x_par1_region = [60001, 2699520]
        y_par1_region = [10001, 2649520]
        x_par2_region = [154931044, 155260560]
        y_par2_region = [59034050, 59363566]
    elif options.build == 'hg38':
        x_par1_region = [10001, 2781479]
        y_par1_region = [10001, 2781479]
        x_par2_region = [155701383, 156030895]
        y_par2_region = [56887903, 57217415]
    elif options.build == 'custom':
        apply_sex_chromosome_rules = False
        x_par1_region = None
        y_par1_region = None
        x_par2_region = None
        y_par2_region = None
        sys.stderr.write(
            "WARNING: custom genome build mode is experimental. "
            "Please ensure that the BAM, VCF, reference FASTA, repeat BED, "
            "and segdup BED use the same coordinate system. "
            "Note: hg19/hg38-specific male X/Y depth_fraction adjustment will not be applied in custom mode.\n"
        )
    else:
        sys.stderr.write(options.build + " is an invalid genome build, please use hg19, hg38, or custom")
        sys.exit(3)

    input_file = options.input_file
    filters = options.vcf_filters
    output_dir = options.output_dir
    annovar_path = options.annovar_path

    global build
    build = options.build

    dbtype = options.dbtype
    cram_ref_dir = options.cram_ref_dir
    reference_fasta = options.reference_fasta
    repeat_bed = options.repeat_bed
    segdup_bed = options.segdup_bed
    skip_annovar = options.skip_annovar

    #check if input path is valid
    if not os.path.exists(input_file):
        sys.stderr.write("Please provide a valid input file.\n")
        sys.exit(2)

    #check if the build version is valid
    if build not in ["hg19", "hg38", "custom"]:
        sys.stderr.write("Please provide a valid genome build version. Supported options: hg19, hg38, custom.\n")
        sys.exit(2)

    if build == "custom":
        if not reference_fasta:
            sys.stderr.write("--reference-fasta is required when --build custom.\n")
            sys.exit(2)
        if not repeat_bed:
            sys.stderr.write("--repeat-bed is required when --build custom.\n")
            sys.exit(2)
        if not segdup_bed:
            sys.stderr.write("--segdup-bed is required when --build custom.\n")
            sys.exit(2)
        if not os.path.exists(reference_fasta):
            sys.stderr.write("Please provide a valid --reference-fasta file.\n")
            sys.exit(2)
        if not os.path.exists(reference_fasta + ".fai"):
            sys.stderr.write("Please provide an indexed --reference-fasta file. Missing .fai index. You can create it with: samtools faidx <reference.fa>\n")
            sys.exit(2)
        if not os.path.exists(repeat_bed):
            sys.stderr.write("Please provide a valid --repeat-bed file.\n")
            sys.exit(2)
        if not os.path.exists(segdup_bed):
            sys.stderr.write("Please provide a valid --segdup-bed file.\n")
            sys.exit(2)

    annovar = None
    annovar_db = None

    if skip_annovar:
        sys.stderr.write(
            "WARNING: ANNOVAR/gnomAD annotation skipped. "
            "Gene annotation and population-frequency filtering may be unavailable. "
            "Results should be interpreted with caution.\n"
        )
    else:
        if not annovar_path:
            sys.stderr.write("Please provide --annovar_path, or use --skip-annovar.\n")
            sys.exit(2)

        if annovar_path.endswith("/"):
            annovar = annovar_path + "annotate_variation.pl"
            annovar_db = annovar_path + "humandb"
        else:
            annovar = annovar_path + "/annotate_variation.pl"
            annovar_db = annovar_path + "/humandb"

        #check if annovar path is valid
        if not os.path.exists(annovar) or not os.path.exists(annovar_db):
            sys.stderr.write("Please provide a valid annovar program directory.\n")
            sys.exit(2)

        #check if db path is valid
        if not os.path.exists(annovar_db + '/' + build + '_' + dbtype + '.txt'):
            sys.stderr.write("Please provide a valid annovar dbtype. db file should have the location annovar/humandb/<build>_<dbtype>.txt\n")
            sys.exit(2)

    #make dir if output_dir does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not output_dir.endswith("/"):
        output_dir += "/"

    #check if cram path is valid
    if cram_ref_dir != None:
        if not os.path.exists(cram_ref_dir):
            sys.stderr.write("Please provide a valid CRAM reference path.\n")
            sys.exit(2)

    global image_outdir
    image_outdir= output_dir + "images/"

    global matrix_outdir
    matrix_outdir = output_dir + "matrices/"

    wfile = open(output_dir + "features.txt", "w")
    header = ["#sample_name","sex", "chrom", "pos", "ref", "alt", "variant", "maf", "lower_CI", "upper_CI", "variant_type", "gene_id",
              "gnomad", "all_repeat", "segdup", "homopolymer", "dinucluotide", "depth_fraction", "image_filepath", "npy_filepath"]
    wfile.write("\t".join(header) + "\n")

    if not os.path.exists(image_outdir):
        os.makedirs(image_outdir)

    if not os.path.exists(matrix_outdir):
        os.makedirs(matrix_outdir)

    #process input files
    all_variants = []
    if filters != None:
        filters = filters.split(",")

    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            sample_name, bam, vcf, depth, sex = line.rstrip().split("\t")

            if bam.endswith(".cram"):
                sys.stdout.write("NOTICE: Input file is CRAM, checking if -c is used... ")
                sys.stdout.write('\n')
                if cram_ref_dir is None:
                    raise Exception("CRAM input must have a reference path. Use -c and make sure to put in the correct path to the reference file.")
                    sys.exit(2)
                else:
                    sys.stdout.write("CRAM reference path has been given for " + str(sample_name))
                    sys.stdout.write('\n')

            if vcf.endswith(".vcf.gz"):
                vcf_file = gzip.open(vcf, "rt")
            elif vcf.endswith(".vcf"):
                vcf_file = open(vcf, "r")
            else:
                raise Exception("input file must contains valid vcf files ending with '.vcf' or '.vcf.gz'")
                sys.exit(2)

            for vcf_line in vcf_file:
                if vcf_line.startswith("#"):
                   continue

                vcf_line = vcf_line.rstrip().split("\t")
                chrom, pos = vcf_line[:2]
                ref, alt = vcf_line[3:5]

                if len(ref) != 1 or len(alt) != 1:
                    continue

                if filters != None and vcf_line[6] not in filters:
                    continue

                all_variants.append([sample_name, bam, chrom, pos, ref, alt, depth, sex])

            vcf_file.close()

    #annotation repeat and segdup
    if build == "custom":
        repeats_dict = custom_repeats_annotation(all_variants, output_dir, repeat_bed, segdup_bed)
    else:
        repeats_dict = repeats_annotation(all_variants, output_dir, build)

    #annovar annotation for gnomad
    if skip_annovar:
        function_dict = {}
        gnomad_dict = {}
        for variant in all_variants:
            sample_name, bam, chrom, pos, ref, alt, depth, sex = variant
            key = "_".join([chrom, pos, ref, alt])
            function_dict[key] = [".", "."]
    else:
        function_dict, gnomad_dict = gnomad_annotation(all_variants, output_dir, annovar, annovar_db, build, dbtype)

    #Add cram ref path to variants
    for variants in all_variants:
        variants.append(cram_ref_dir)

    #prepare args list
    args_list = [(variant, image_outdir, matrix_outdir, build, reference_fasta, apply_sex_chromosome_rules) for variant in all_variants]

    #draw images
    try:
        pool = Pool(8) # on 8 processors
        results = pool.map(multiprocess_iterator, args_list, chunksize=8)

        for result in results:
            if not result:
                continue

            sample_name, sex, key, maf, lower_CI, upper_CI, depth_fraction, is_homopolymer, is_dinucleotide, image_file, npy_file = result
            npy_file = os.path.abspath(npy_file + ".npy")
            image_file = os.path.abspath(image_file)

            all_repeat, segdup = repeats_dict[key]

            if key in gnomad_dict.keys():
                gnomad = gnomad_dict[key]
            else:
                gnomad = 0

            var_type, gene = function_dict.get(key, [".", "."])
            chrom, pos, ref, alt = key.split("_")

            wfile.write("\t".join(list(map(str,[sample_name, sex, chrom, pos, ref, alt, key, maf, lower_CI, upper_CI,
                                                var_type, gene, gnomad, int(all_repeat), int(segdup),
                                                int(is_homopolymer), int(is_dinucleotide), depth_fraction,
                                                image_file, npy_file]))) + "\n")

        wfile.close()

    #except:
     #   sys.stderr.write("Error during multiprocess feature extraction.\n")
    finally: # To make sure processes are closed in the end, even if errors happen
        pool.close()
        pool.join()

    #merge features together
    time_elapsed = time.time() - since
    sys.stdout.write("complete image recoding in {:.0f}m {:.0f}s".format(time_elapsed // 60, time_elapsed % 60) + "\n")
