import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable

# Config Parsers 
# try to import the python2 ConfigParser
# if unable to import, then try to import the python3 configparser
try:
    import ConfigParser as configparser
except ImportError:
    import configparser

config = configparser.ConfigParser()
config_file = os.path.abspath("baqlava/configs/baqlava.cfg")
config.read(config_file)


# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.0.1",                    #Update the version as needed
    description="Download and place full BAQLaVa databases"     #Update the description as needed
    )

###############
# custom args #
###############
workflow.add_argument(
    name = "nucdb",
    desc = "nucleotide database folder to use",
    default = config.get('database','nucleotide_db_full'))

workflow.add_argument(
    name = "protdb",
    desc = "protein database folder to use",
    default = config.get('database','uniref_db_full'))

args = workflow.parse_args()


def main():
    ############################
    # function to run workflow #
    ############################

    workflow.add_task(
    "wget [n_db]",
    targets = ["BAQLaVa.V0.1.nucleotide.tar.gz"],
    n_db = os.path.abspath(args.nucdb))

    workflow.go()


if __name__ == '__main__':
    main()
