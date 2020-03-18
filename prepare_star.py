"""
example:
    #  check run options
    python prepare_star.py
    # initial config generation
    python prepare_star.py \\
        --project_root ~/projects/KAZ_RNA/KAZ_RNA_STAR \\
        --fastq_dirs_list ~/icebox/fastq_gz/KAZ_RNA/ \\
        --sample_delimiter . \\
        --fastq_extension .fastq.gz \\
        --R1_fastq_extension .R1.fastq.gz \\
        --R2_fastq_extension .R2.fastq.gz \\
        --script_dir_name scripts
    # optional: generate config with --add_tokens option
    # to add tokens to each step, for rerun from last failed step
    # precise config tuning
    # scripts generation
    python prepare_star.py -j ~/projects/KAZ_RNA/KAZ_RNA_STAR/scripts/default_settings.json

"""
__VERSION__ = "0.1.0"

import os
import argparse
import json

from collections import defaultdict


__NOT_READY__ = "NOT_READY"
__READY__ = "READY"
__ALMOST_READY__ = "ALMOST_READY"


def main():
    settings = parse_arguments_to_settings()
    if settings["ready"] == __ALMOST_READY__:
        save_settings(settings)
    elif settings["ready"] == __READY__:
        run_pipeline(settings)
    else:
        print(__doc__)


def parse_arguments_to_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--settings_json", default=None, required=False)
    parser.add_argument("--project_root", default=None, required=False)
    parser.add_argument("--fastq_dirs_list", default=[], required=False, nargs="+")
    parser.add_argument("--sample_delimiter", default="_", required=False)
    parser.add_argument("--fastq_extension", default=".fastq.gz", required=False)
    parser.add_argument("--R1_fastq_extension", default=".R1.fastq.gz", required=False)
    parser.add_argument("--R2_fastq_extension", default=".R2.fastq.gz", required=False)
    parser.add_argument("--script_dir_name", default="scripts", required=False)
    parser.add_argument("--run_annovar", action="store_true")
    parser.add_argument("--add_tokens", action="store_true")
    parser.add_argument("--debug", action="store_true")
    #
    args = parser.parse_args()
    if args.settings_json:
        settings = json.load(open(args.settings_json))   # config exist, import it
        __samples_dict__ = load_fastq_samples(settings)  # find all fastq files
        __samples_list__ = settings["samples_list"]  # get list of target samples from config
        settings["samples_dict"] = {  # filter target samples from all fastq
            list_key: dict_value
            for list_key in __samples_list__
            for dict_key, dict_value in __samples_dict__.items()
            if list_key == dict_key or list_key + "_m" == dict_key
        }
        settings["ready"] = __READY__
    elif args.project_root:
        settings = {
            "settings_json": args.settings_json,
            "project_root": args.project_root,
            "fastq_dirs_list": args.fastq_dirs_list,
            "sample_delimiter": args.sample_delimiter,
            "fastq_extension": args.fastq_extension,
            "R1_fastq_extension": args.R1_fastq_extension,
            "R2_fastq_extension": args.R2_fastq_extension,
            "script_dir_name": args.script_dir_name,
            "run_annovar": args.run_annovar,
            "add_tokens": args.add_tokens,
            "debug": args.debug,
            "ready":__ALMOST_READY__,
        }
    else:
        settings = {
            "ready":__NOT_READY__,
        }
    return settings


def load_fastq_samples(settings):
    fastq_dirs_list = settings["fastq_dirs_list"]
    sample_delimiter = settings["sample_delimiter"]
    fastq_extension = settings["fastq_extension"]
    R1_fastq_extension = settings["R1_fastq_extension"]
    R2_fastq_extension = settings["R2_fastq_extension"]
    #
    res = defaultdict(lambda: defaultdict(str))
    for fastq in get_files_generator(fastq_dirs_list, fastq_extension):
        sample = os.path.basename(fastq).split(sample_delimiter)[0]
        if fastq.endswith(R1_fastq_extension):
            res[sample]["read1"] = fastq
        elif fastq.endswith(R2_fastq_extension):
            res[sample]["read2"] = fastq
    res = {
        key: value
        for key, value in res.items()
        if key + "_m" not in res
    }
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


def save_settings(settings):
    settings["project_script_dir"] = os.path.join(
        settings["project_root"],
        settings["script_dir_name"],
    )
    settings.update(get_default_settings(settings))
    mkdir(settings["project_script_dir"])
    json_file = os.path.join(
        settings["project_script_dir"],
        "default_settings.json",
    )
    with open(json_file, "w") as f:
        default_settings_str = json.dumps(settings, indent=4, sort_keys=True)
        f.write(default_settings_str)
    # debug print
    print(f"# ls   {settings['project_script_dir']}")
    print(f"# more {json_file}")
    print(f"# nano {json_file}")
    print(f"# python prepare_star.py -j {json_file}")


def get_default_settings(d):
    default_settings_dict = {
        "number_of_threads": "8",
        "fastq_dirs_list": d["fastq_dirs_list"],
        "sample_delimiter": d["sample_delimiter"],
        "fastq_extension": d["fastq_extension"],
        "R1_fastq_extension": d["R1_fastq_extension"],
        "R2_fastq_extension": d["R2_fastq_extension"],
        "samples_list": sorted(load_fastq_samples(d)) if  d["fastq_dirs_list"] else [],
        "project_root": d["project_root"],

        "tools": {
            "STAR": "~/anaconda3/envs/rna/bin/STAR",
        },
        "databases": {
            "ref_dir": "~/PublicData/ensembl_GRCh37_75/",
            "ref": "~/PublicData/ensembl_GRCh37_75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa",
            "ref_gtf": "~/PublicData/ensembl_GRCh37_75/Homo_sapiens.GRCh37.75.gtf"
        },
    }
    return default_settings_dict


def mkdir(dir_name):
    os.makedirs(dir_name, exist_ok=True)


def run_pipeline(settings):
    for sample in sorted(settings["samples_dict"]):
        sample_settings = settings
        sample_settings["sample"] = sample
        sample_settings = get_settings_for_RNA_SEQ_PIPELINE(sample_settings)
        cmd_list = get_cmd_list_for_RNA_SEQ_PIPELINE(sample_settings)
        write_cmd_list_to_file(sample_settings, cmd_list)
    # debug print
    print(f"# ls {settings['project_script_dir']}")
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ss.sh ); do echo $i; done" )
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ss.sh ); do qsub $i; done" )


def get_settings_for_RNA_SEQ_PIPELINE(sample_dict):
    sample = sample_dict["sample"]
    sample_dir = os.path.join(sample_dict["project_root"], sample)
    sample_path_prefix = os.path.join(sample_dir, sample)
    _dict = {
        "sample": sample,
        "sample_dir": sample_dir,
        "project_script_dir" : sample_dict["project_script_dir"],
        "number_of_threads": sample_dict["number_of_threads"],

        "STAR": sample_dict["tools"]["STAR"],

        "ref_dir": sample_dict["databases"]["ref_dir"],
        "ref": sample_dict["databases"]["ref"],
        "ref_gtf": sample_dict["databases"]["ref_gtf"],

        "read1": sample_dict["samples_dict"][sample]["read1"],
        "read2": sample_dict["samples_dict"][sample]["read2"],

        "star_out_prefix": sample_path_prefix,

        "add_tokens": sample_dict["add_tokens"],
        "debug": sample_dict["debug"],
    }
    return _dict

###############################################################################
def get_cmd_list_for_RNA_SEQ_PIPELINE(sample_settings):
    if sample_settings["debug"]:
        print("\n\n# DEBUG", sample_settings)
    cmd_list = [
        get_mkdir_cmd(sample_settings),
        get_cmd_star(sample_settings),
    ]
    return cmd_list


###############################################################################
def get_mkdir_cmd(d):
    return "mkdir -p {sample_dir}".format(**d)


###############################################################################
def reduce_spaces_and_newlines(s):
    s = s.replace("\n", " ")
    s = " ".join([i for i in s.split(" ") if i])
    return s

def get_cmd(d):
    d["token"] = "{sample_dir}/token.{sample}.{token_suffix}".format(**d)
    d["flags"] = " && ".join([" [ -f {:s} ] ".format(i) for i in d["files_list"]]) + " && [ ! -f {token} ] ".format(**d)
    if d["add_tokens"]:
        cmd = """
            {flags} &&
            dt1=`date +%y%m%d_%H%M%S` && echo $dt1 {token} &&
            {main_cmd} &&
            du {out_file} > {out_file}.$dt1.du &&
            md5sum {out_file} > {out_file}.$dt1.md5 &&
            dt2=`date +%y%m%d_%H%M%S` &&
            echo $dt1 $dt2 > {token} ||
            echo "TOKEN SKIPPED {token}"
            """.format(**d)
    else:
        cmd = "{main_cmd}".format(**d)
    return reduce_spaces_and_newlines(cmd)


###############################################################################
def get_cmd_star(d):
    d["token_suffix"] = "star"
    d["files_list"] = [d["read1"], d["read2"]]
    d["out_file"] = d["star_out_prefix"]
    d["main_cmd"] = bash_star(d)
    return get_cmd(d)

def bash_star(d):
    return """nice {STAR}
             --runThreadN  8
             --genomeDir    {ref_dir}
             --sjdbGTFfile  {ref_gtf}
             --readFilesIn  {read1}  {read2}
             --outFileNamePrefix  {star_out_prefix}
             --quantMode  TranscriptomeSAM  GeneCounts
             --readFilesCommand gunzip -c""".format(**d)


###############################################################################
def write_cmd_list_to_file(sample_settings, cmd_list):
    script_file = os.path.join(
        sample_settings["project_script_dir"],
        sample_settings["sample"] + ".ss.sh",
    )
    with open(script_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        for cmd in cmd_list:
            new_line = cmd + "\n\n"
            f.write(new_line)


###############################################################################
if __name__ == "__main__":
    main()
