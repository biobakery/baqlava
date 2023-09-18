# python anadama_pipeline_v3.py -i ../nucleotide_db/raw_databases/databases/ -o pipeline_v3
# before running: source activate cluster (mmseqs installed)
# hutlab load centos7/python3/anadama2/0.10.0-devel
# module load samtools/1.10-fasrc01
# module load blast/2.6.0+-fasrc01
# after fasrc updates:
# module load centos6/0.0.1-fasrc01
# module load bowtie2/2.3.2-fasrc01
# module load diamond/2.0.4-fasrc01


import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable

# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.0.1",                    #Update the version as needed
    description="Generate Nucleotide Database"     #Update the description as needed
    )

###############
# custom args #
###############

#workflow.add_argument(
#    name = "databases",
#    desc = "directory where all and ONLY raw databses are saved",
#    default = "/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava_v1/primary_files/nucleotide_db/raw_databases/databases/"

#workflow.add_argument(
#    name = "fastANI",
#    desc = "location of fastANI executable",
#    default = "/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava_V0_sandbox/20230113_fastANI/FastANI/fastANI")

args = workflow.parse_args()
  
def main():
    ############################
    # function to run workflow #
    ############################
    final_database_dir = args.output + "/final_databases/"
    mmseqs_dir = args.output + "/mmseqs/" 
    bt2_dir = args.output + "/bowtie2/"
    cluster_dir = args.output + "/cluster/"
    marker_dir = args.output + "/markers/"
    marker_info_dir = args.output + "/marker_information/"
    demo_dir = args.output + "/demo/"

    os.system("mkdir " + args.output)
    os.system("mkdir " + final_database_dir)
    os.system("mkdir " + bt2_dir)
    os.system("mkdir " + mmseqs_dir)
    os.system("mkdir " + cluster_dir)
    os.system("mkdir " + marker_dir)
    os.system("mkdir " + marker_info_dir)
    os.system("mkdir " + demo_dir)
    
    #workflow.add_task("python scripts/d01_format_database.py [args[0]] [args[1]] [targets[0]]",
    #args = [args.input, final_database_dir],
    #targets = [final_database_dir + "BAQLaVa_nucleotidedb.fasta", final_database_dir + "BAQLaVa_nucleotidedb_reference_file.txt"])
        
    #workflow.add_task("mmseqs createdb [depends[0]] [targets[0]]",
    #depends = [final_database_dir + "BAQLaVa_nucleotidedb.fasta"],
    #targets = [mmseqs_dir + "BAQ"])
    
    #workflow.add_task("mmseqs cluster [depends[0]] [args[0]] [args[1]] -c 0.85 --cov-mode 1 --min-seq-id 0.95 --seq-id-mode 0 --threads 10 --split-memory-limit 32G",
    #depends = [mmseqs_dir + "BAQ"],
    #args = [mmseqs_dir + "BAQ_clu", mmseqs_dir + "tmp"],
    #targets = [mmseqs_dir + "BAQ_clu.index"])
    
    #workflow.add_task("mmseqs createtsv [depends[0]] [depends[0]] [args[0]] [targets[0]]",
    #depends = [mmseqs_dir + "BAQ"],
    #args = [mmseqs_dir + "BAQ_clu"],
    #targets = [mmseqs_dir + "BAQ_clu.tsv"])
    
    workflow.add_task("python scripts/d02_make_derep_database.py [depends[0]] [depends[1]] [args[0]]",
    depends = [final_database_dir + "BAQLaVa_nucleotidedb.fasta", mmseqs_dir + "BAQ_clu.tsv"],
    args = [bt2_dir],
    targets = [bt2_dir + "BAQLaVa_nucleotidedb_dereplicated.fasta"])

    #workflow.add_task("python scripts/d03_make_100mers.py [depends[0]] [depends[1]] [args[0]]",
    #depends = [final_database_dir + "BAQLaVa_nucleotidedb.fasta", mmseqs_dir + "BAQ_clu.tsv"],
    #args = [bt2_dir],
    #targets = [bt2_dir + "BAQLaVa_dereplicated_100mers.fasta"])

    #workflow.add_task("bowtie2-build -f [depends[0]] [args[0]]",
    #depends = [final_database_dir + "BAQLaVa_nucleotidedb.fasta"],
    #args = [bt2_dir + "BAQLaVa_bt2"],
    #targets = [bt2_dir + "BAQLaVa_bt2.1.bt2"])    

    #workflow.add_task("bowtie2 -x [args[0]] -U [depends[0]] -f -a --very-sensitive -S [targets[0]] --threads 12",
    #depends = [bt2_dir + "BAQLaVa_dereplicated_100mers.fasta", bt2_dir + "BAQLaVa_bt2.1.bt2"],
    #args = [bt2_dir + "BAQLaVa_bt2"],
    #targets = [bt2_dir + "BAQ_100mers_mapped.sam"])

    #workflow.add_task("samtools sort -o [targets[0]] -n -@ 5 [depends[0]]",
    #depends = [bt2_dir + "BAQ_100mers_mapped.sam"],
    #targets = [bt2_dir + "BAQ_100mers_mapped_sorted.sam"])

    #workflow.add_task("python scripts/d04_loop_through_sam.py [depends[0]] [depends[1]] [args[0]]",
    #depends = [bt2_dir + "BAQ_100mers_mapped_sorted.sam", mmseqs_dir + "BAQ_clu.index"], 
    #args = [cluster_dir],
    #targets = [cluster_dir + "mapped_reads.txt", cluster_dir + "self_mapped.txt", cluster_dir + "not_marker_candidates.txt"])

    #workflow.add_task("python scripts/d05_duplicated_genomes.py [depends[0]] [args[0]]",
    #depends = [cluster_dir + "self_mapped.txt"],
    #args = [cluster_dir],
    #targets = [cluster_dir + "duplicated_genomes.txt"])    

    #workflow.add_task("cut -f1,2,3 [depends[0]] > [targets[0]]",
    #depends = [cluster_dir + "mapped_reads.txt"],
    #targets = [cluster_dir + "mapped_reads_columns123.txt"])

    #workflow.add_task("python scripts/d06_bowtie2_cluster.py [depends[0]] [depends[1]] [depends[2]] [depends[3]] [args[0]]",
    #depends = [cluster_dir + "mapped_reads_columns123.txt", bt2_dir + "BAQLaVa_dereplicated_100mers.fasta", final_database_dir + "BAQLaVa_nucleotidedb_reference_file.txt", mmseqs_dir + "BAQ_clu.tsv"],
    #args = [cluster_dir],
    #targets = [cluster_dir + "bowtie2_clu.txt", cluster_dir + "bowtie2_cc.txt"])

    #workflow.add_task("python scripts/d07_combine_mmseqs_bowtie2_cluster.py [depends[0]] [depends[1]] [depends[2]] [depends[3]] [depends[4]] [args[0]]",
    #depends = [final_database_dir + "BAQLaVa_nucleotidedb_reference_file.txt", mmseqs_dir + "BAQ_clu.tsv", cluster_dir + "bowtie2_clu.txt", cluster_dir + "duplicated_genomes.txt", cluster_dir + "bowtie2_cc.txt"],
    #args = [cluster_dir],
    #targets = [cluster_dir + "final_clusters.txt"])

    #workflow.add_task("python scripts/d08_loop_2nd.py [depends[0]] [depends[1]] [depends[2]] [args[0]]",
    #depends = [cluster_dir + "final_clusters.txt", bt2_dir + "BAQLaVa_dereplicated_100mers.fasta", bt2_dir + "BAQ_100mers_mapped_sorted.sam"],
    #args = [marker_dir],
    #targets = [marker_dir + "cluster_candidate100mers.txt", marker_dir + "supercluster_candidate100mers.txt", marker_dir + "supercluster_alignments.txt"])

    workflow.add_task("python scripts/d09_make_clusterrep_markers.py [depends[0]] [depends[1]] [depends[2]] [args[0]]",
    depends = [marker_dir + "cluster_candidate100mers.txt", bt2_dir + "BAQLaVa_dereplicated_100mers.fasta", cluster_dir + "final_clusters.txt"],
    args = [marker_dir],
    targets = [marker_dir + "clusterrep_markers.fasta", marker_dir + "clusterrep_marker_stats.txt", marker_dir + "clusterrep_marker_lengths.txt"])

    workflow.add_task("python scripts/d10_make_supercluster_markers.py [depends[0]] [depends[1]] [depends[2]] [args[0]]",
    depends = [marker_dir + "supercluster_candidate100mers.txt", cluster_dir + "final_clusters.txt", bt2_dir + "BAQLaVa_dereplicated_100mers.fasta"],
    args = [marker_dir],
    targets= [marker_dir + "supercluster_markers.fasta", marker_dir + "supercluster_marker_stats.txt", marker_dir + "supercluster_marker_lengths.txt"])
    
    workflow.add_task("mmseqs createdb [depends[0]] [targets[0]]",
    depends = [marker_dir + "supercluster_markers.fasta"],
    targets = [marker_dir + "SC"])

    workflow.add_task("mmseqs cluster [depends[0]] [args[0]] [args[1]] -c 0.8 --cov-mode 2 --min-seq-id 0.95 --seq-id-mode 0",
    depends = [marker_dir + "SC"],
    args = [marker_dir + "SC_clu", marker_dir + "tmp"],
    targets = [marker_dir + "SC_clu.index"])

    workflow.add_task("mmseqs createtsv [depends[0]] [depends[0]] [args[0]] [targets[0]]",
    depends = [marker_dir + "SC"],
    args = [marker_dir + "SC_clu"],
    targets = [marker_dir + "SC_clu.tsv"])

    workflow.add_task("python scripts/d11_derep_supercluster_markers.py [depends[0]] [depends[1]] [targets[0]]",
    depends = [marker_dir + "supercluster_markers.fasta", marker_dir + "SC_clu.tsv"],
    targets = [marker_dir + "supercluster_markers_dereplicated.fasta"])

    workflow.add_task("cat [depends[0]] [depends[1]] > [targets[0]]",
    depends = [marker_dir + "clusterrep_markers.fasta", marker_dir + "supercluster_markers_dereplicated.fasta"],
    targets = [marker_dir + "all_markers.fasta"])

    workflow.add_task("bowtie2-build -f [depends[0]] [args[0]]",
    depends = [marker_dir + "all_markers.fasta"],
    args = [marker_dir + "allmarkers_bt2"],
    targets = [marker_dir + "allmarkers_bt2.1.bt2"])

    workflow.add_task("python scripts/d12_write_idmap.py [depends[0]] [depends[1]] [depends[2]] [targets[0]]",
    depends = [marker_dir + "clusterrep_marker_lengths.txt", marker_dir + "supercluster_marker_lengths.txt", marker_dir + "all_markers.fasta"],
    targets = [marker_info_dir + "idmap1.txt"])

    workflow.add_task("python scripts/d13_cluster_taxonomy.py [depends[0]] [depends[1]] [depends[2]] [depends[3]] [targets[0]]",
    depends = ["/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlavav02/baqlava/database_generation/input/reference_files/ICTV_ref2.txt", cluster_dir + "final_clusters.txt", final_database_dir + "BAQLaVa_nucleotidedb_reference_file.txt", marker_info_dir + "idmap1.txt"],
    targets = [marker_info_dir + "marker_taxonomy.txt"])

    workflow.add_task("python scripts/d14_DNA_to_AA.py [depends[0]] [depends[1]] [targets[0]] [targets[1]] [targets[2]]", 
    depends = [marker_dir + "all_markers.fasta", marker_info_dir + "idmap1.txt"], 
    targets = [marker_dir + "all_markers_translated.fasta", marker_info_dir + "idmap2.txt", marker_info_dir + "translated_markers_conversion.txt"])

    workflow.add_task("diamond makedb --in [depends[0]] -d [args[0]]",
    depends = [marker_dir + "all_markers_translated.fasta"],
    args = [marker_dir + "BAQLaVa.V0.2.201901b"],
    targets = [marker_dir + "BAQLaVa.V0.2.201901b.dmnd"])

    workflow.add_task("python scripts/d15_format_species.py [depends[0]] [targets[0]]",
    depends = [marker_info_dir + "marker_taxonomy.txt"],
    targets = [marker_info_dir + "species_conversion.txt"])

    #workflow.add_task("python scripts/d15_make_demo_databases.py [depends[0]] [depends[1]] [depends[2]] [targets[0]] [targets[1]]",
    #depends = [marker_info_dir + "marker_taxonomy.txt", marker_dir + "all_markers.fasta", marker_dir + "all_markers_translated.fasta"],
    #targets = [demo_dir + "BAQLaVa.V0.2.demo.nucleotide.fa", demo_dir + "BAQLaVa.V0.2.201901b.demo.translated.fa"])

    #workflow.add_task("bowtie2-build -f [depends[0]] [args[0]]",
    #depends = [demo_dir + "BAQLaVa.V0.2.demo.nucleotide.fa"],
    #args = [demo_dir + "BAQLaVa.V0.2.demo.nucleotide"],
    #targets = [demo_dir + "BAQLaVa.V0.2.demo.nucleotide.1.bt2"])

    #workflow.add_task("diamond makedb --in [depends[0]] -d [args[0]]",
    #depends = [demo_dir + "BAQLaVa.V0.2.201901b.demo.translated.fa"],
    #args = [demo_dir + "BAQLaVa.V0.2.201901b.demo.translated"],
    #targets = [demo_dir + "BAQLaVa.V0.2.201901b.demo.translated.dmnd"])

    #workflow.add_task("humann --input [args[0]] --output [args[1]] --bypass-nucleotide-index --nucleotide-database [args[2]] --id-mapping [depends[0]] --bypass-translated-search",
    #depends = [marker_info_dir + "idmap.txt"],
    #args = ["/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava_V0_sandbox/20230419_artillumina/ICTV_sample_2.fq", marker_info_dir + "/output/", marker_dir],
    #targets = [marker_info_dir + "/output/ICTV_sample_2_genefamilies.tsv"])

    workflow.go()
    
if __name__ == '__main__':
    main()



