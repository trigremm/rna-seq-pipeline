"""
example:
    #  check settings
    nano easy_STAR.py

    # print settings
    python easy_STAR.py

    # generate scripts
    python easy_STAR.py run

    # dangerous
    # to reset to last git version
    # git fetch --all; git reset --hard origin/master
"""


__VERSION__ = "0.1.0"


import os, sys
import pprint
pp = pprint.PrettyPrinter(indent=4)

from collections import defaultdict


def get_settings():
    settings = {
      "number_of_threads": "4",
      "project_root": "/mnt/ds2413p/aigul/KAZ_RNA/KAZ_RNA_STAR",
      "fastq_dirs_list": ["/mnt/ds2413p/aigul/archive", "/mnt/ds2413p/q-symphony/icebox/Archive/KAZ_WT/", ],
      "sample_delimiter": "_",
      "fastq_extension": ".fastq.gz",
      "R1_fastq_extension": ".R1.fastq.gz",
      "R2_fastq_extension": ".R2.fastq.gz",
      "STAR": "/home/adminrig/anaconda3/envs/rna/bin/STAR",
      "ref_dir": "~/PublicData/ensembl_GRCh37_75/",
      "ref": "~/PublicData/ensembl_GRCh37_75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa",
      "ref_gtf": "~/PublicData/ensembl_GRCh37_75/Homo_sapiens.GRCh37.75.gtf",
      "add_tokens": True,
    }
    samples_dict = load_fastq_samples(settings)
    settings.update({
        "samples_dict": {k:samples_dict[k] for k in sorted(samples_dict)},
        "project_scripts_dir": os.path.join(settings["project_root"], "scripts"),
    })
    return settings


def run_pipeline(settings):
    _ = -1
    for _, sample in enumerate(sorted(settings["samples_dict"])):
        sample_settings = get_sample_settings(sample, settings)
        cmd_list = get_cmd_list(sample_settings)
        write_cmd_list_to_file(sample_settings, cmd_list)
    return _


def get_sample_settings(sample, settings):
    sample_dir = os.path.join(settings["project_root"], sample)
    sample_prefix = os.path.join(sample_dir, sample)
    settings.update({
        "sample": sample,
        "sample_dir": sample_dir,
        "read1": settings["samples_dict"][sample]["read1"],
        "read2": settings["samples_dict"][sample]["read2"],
        "star_out_prefix": sample_prefix + '.STAR.',
        "star_out_sam": sample_prefix + '.STAR.Aligned.out.sam',
        "star_out_genecounts":sample_prefix + ".STAR.ReadsPerGene.out.tab",
    })
    return settings


def get_mkdir_cmd(d):
    return "mkdir -p {sample_dir}".format(**d)


def get_cmd_list(sample_settings):
    cmd_list = [
        get_mkdir_cmd(sample_settings),
        get_cmd_star(sample_settings),
        bash_delete_sam(sample_settings)
    ]
    return cmd_list


def write_cmd_list_to_file(sample_settings, cmd_list):
    script_file = os.path.join(
        sample_settings["project_scripts_dir"],
        sample_settings["sample"] + ".ss.sh",
    )
    with open(script_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        for cmd in cmd_list:
            new_line = cmd + "\n\n"
            f.write(new_line)


def load_fastq_samples(settings):
    fastq_dirs_list = settings["fastq_dirs_list"]
    sample_delimiter = settings["sample_delimiter"]
    fastq_extension = settings["fastq_extension"]
    R1_fastq_extension = settings["R1_fastq_extension"]
    R2_fastq_extension = settings["R2_fastq_extension"]

    res = defaultdict(lambda: defaultdict(str))
    for fastq in get_files_generator(fastq_dirs_list, fastq_extension):
        sample = os.path.basename(fastq).split(sample_delimiter)[0]
        if fastq.endswith(R1_fastq_extension):
            res[sample]["read1"] = fastq
        elif fastq.endswith(R2_fastq_extension):
            res[sample]["read2"] = fastq
    return res


def get_files_generator(dirs_list, extension=""):
    for path in dirs_list:
        for data_file in os.listdir(path):
            if data_file:
                data_path = os.path.join(path, data_file)
                if os.path.isfile(data_path) and data_path.endswith(extension):
                    yield data_path
                elif os.path.isdir(data_path):
                    yield from get_files_generator([data_path], extension)


def reduce_spaces_and_newlines(func):
    def wrapper(*args, **kwargs):
        s = func(*args, **kwargs)
        s = s.replace("\n", " ")
        s = " ".join([i for i in s.split(" ") if i])
        return s
    return wrapper


@reduce_spaces_and_newlines
def get_cmd(d):
    d["token"] = "{sample_dir}/token.{sample}.{token_suffix}".format(**d)
    d["flags"] = " && ".join([" [ -f {:s} ] ".format(i) for i in d["files_list"]]) + " && [ ! -f {token} ] ".format(**d)
    if d["add_tokens"]:
        cmd = """
            {flags} &&
            dt1=`date +%y%m%d_%H%M%S` && echo $dt1 {token} &&
            {main_cmd} &&
            dt2=`date +%y%m%d_%H%M%S` &&
            echo $dt1 $dt2 > {token} ||
            echo "TOKEN SKIPPED {token}"
            """.format(**d)
    else:
        cmd = "{main_cmd}".format(**d)
    return cmd


def get_cmd_star(d):
    d["token_suffix"] = "star"
    d["files_list"] = [d["read1"], d["read2"]]
    d["out_file"] = d["star_out_genecounts"]
    d["main_cmd"] = bash_star(d)
    return get_cmd(d)


@reduce_spaces_and_newlines
def bash_star(d):
    return """nice {STAR}
             --runThreadN  8
             --genomeDir    {ref_dir}
             --sjdbGTFfile  {ref_gtf}
             --readFilesIn  {read1}  {read2}
             --outFileNamePrefix  {star_out_prefix}
             --quantMode  TranscriptomeSAM  GeneCounts
             --readFilesCommand gunzip -c""".format(**d)


def bash_delete_sam(d):
    return """rm {star_out_sam}""".format(**d)


if __name__ == '__main__':
    settings = get_settings()
    if 'run' in sys.argv:
        num = run_pipeline(settings)
        print (f"# processed {num} samples")
        print (f"# cd {settings['project_scripts_dir']}")
    else:
        print (__doc__)
        pp.pprint(settings)
