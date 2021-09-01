import fileinput
import logging
import os
import shlex
import shutil
import subprocess
import sys
import uuid
import gzip
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path
from typing import Tuple, List, Dict

import click


class ConcatenationError(Exception):
    pass


class MafftError(Exception):
    pass


def extract_name_from_path(path: Path) -> str:
    return path.name.split("_na_")[0]


def concatenate(infiles: List[Path], outfile: Path):
    with outfile.open("w") as outstream, fileinput.input(files=infiles) as instream:
        for line in map(str.rstrip, instream):
            if line.startswith(">"):
                old_header = line[1:]
                new_header = f">{uuid.uuid4()} {old_header}"
                print(new_header, file=outstream)
            else:
                print(line, file=outstream)


def update_with_new_sequences(msa: Path, new_sequences: List[Path], outdir: Path):
    name = extract_name_from_path(msa)
    logging.debug(f"Updating MSA for {name}...")

    new_sequence_file = outdir / f"tmp.new_sequences.{name}.fa"
    concatenate(new_sequences, new_sequence_file)

    # make paths shell-safe
    new_msa = shlex.quote(str(outdir / f"{name}.fa"))
    existing_msa = shlex.quote(str(msa))
    new_sequence_shstr = shlex.quote(str(new_sequence_file))

    # have run into issues on noah with scratch being full
    mafft_tmpdir = outdir / str(uuid.uuid4())
    if mafft_tmpdir.exists():
        shutil.rmtree(mafft_tmpdir)
    mafft_tmpdir.mkdir()
    env = os.environ
    env["TMPDIR"] = str(mafft_tmpdir)

    args = " ".join(
        [
            "mafft",
            "--auto",
            "--thread",
            "-1",
            "--add",
            new_sequence_shstr,
            existing_msa,
            ">",
            new_msa,
        ]
    )
    logging.debug(f"Running mafft with args: {args}")
    process = subprocess.Popen(
        args,
        stderr=sys.stderr,
        encoding="utf-8",
        shell=True,
        env=env,
    )
    exit_code = process.wait()
    shutil.rmtree(mafft_tmpdir)
    if exit_code != 0:
        raise MafftError(
            f"Failed to execute mafft for {name} - check the logs"
        )
    logging.debug(f"Finished updating MSA for {name}")
    new_sequence_file.unlink(missing_ok=True)


def copy_compressed_files(infile: Path, outfile: Path):
    content = gzip.open(infile, "rt").read()
    with open(outfile, "w") as f:
        f.write(content)


@click.command()
@click.help_option("--help", "-h")
@click.argument(
    "indirs", nargs=-1, type=click.Path(exists=True, file_okay=False, resolve_path=True)
)
@click.option(
    "-M",
    "--msa-dir",
    help="Directory containing the MSAs to update.",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    required=True,
)
@click.option(
    "-o",
    "--outdir",
    help="Directory to write the updated MSAs to.",
    type=click.Path(file_okay=False, writable=True, resolve_path=True),
    required=True,
)
@click.option(
    "-e",
    "--extensions",
    help=(
        "Valid extensions for MSAs and new sequences. Use option multiple times for "
        "more than one extension. Note, there is no need to specify '.gz' if the files "
        "are compressed."
    ),
    default=[".fa", ".fasta"],
    show_default=True,
    multiple=True,
    metavar="EXT",
)
@click.option(
    "-j",
    "-t",
    "--processes",
    help=(
        "Maximum number of processes to use. A value of 0 will use the number of "
        "processes on the machine."
    ),
    default=1,
    show_default=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    indirs: Tuple[str],
    msa_dir: str,
    outdir: str,
    verbose: bool,
    extensions: Tuple[str],
    processes: int,
):
    """A script to update multiple sequence alignments (MSAs) with newly discovered
    sequences (de novo paths).

    The files are matched up by comparing the part of the file name *before* the first
    `.` character.

    INDIRS: Any number of directories containing fasta files with the new sequences to
    add.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    logging.info("Searching for MSAs...")
    msa_lookup: Dict[str, Path] = dict()
    for file in Path(msa_dir).rglob("*"):
        if any(ext in extensions for ext in file.suffixes):
            name = extract_name_from_path(file)
            msa_lookup[name] = file
    logging.info(f"Found {len(msa_lookup)} MSA files.")

    logging.info("Searching for discovered sequence files...")
    denovo_lookup: Dict[str, List[Path]] = defaultdict(list)
    for directory in indirs:
        for file in Path(directory).rglob("*"):
            if any(ext in extensions for ext in file.suffixes):
                name = file.name.split(".")[0]
                denovo_lookup[name].append(file)
    logging.info(f"{len(denovo_lookup)} MSAs have discovered sequence files.")

    no_update_msas = msa_lookup.keys() - denovo_lookup.keys()
    to_update_msas = msa_lookup.keys() - no_update_msas
    logging.info(
        f"{len(no_update_msas)} MSAs have no discovered sequences. Skipping MSA for "
        f"those..."
    )
    if processes == 0:
        logging.info(f"{processes} processed requested. Using all available...")
        processes = None

    jobs = []
    logging.info("Copying files that don't require update to outdir...")
    for name in no_update_msas:
        original_file = msa_lookup[name]
        new_file = outdir / f"{name}.fa"
        if original_file.suffix == ".gz":
            jobs.append((original_file, new_file))
        else:
            new_file = shutil.copy(str(original_file), str(new_file))
            logging.debug(f"{original_file} copied to {new_file}")

    if jobs:
        logging.info("Copying compressed files that don't require update to outdir...")
        with Pool(processes=processes) as pool:
            pool.starmap(copy_compressed_files, jobs)
        logging.info("Successfully copied compressed files")

    jobs = []
    for name in to_update_msas:
        msa = msa_lookup[name]
        new_sequences = denovo_lookup[name]
        jobs.append((msa, new_sequences, outdir))

    with Pool(processes=processes) as pool:
        logging.info(f"Updating {len(jobs)} MSAs...")
        pool.starmap(update_with_new_sequences, jobs)

    logging.info("All MSAs updated!")


if __name__ == "__main__":
    main()
