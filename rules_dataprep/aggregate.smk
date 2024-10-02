# dirty implementation because I don't know how to connect the odb rules to the taxon rules
# but since ODB download is pretty fast compared to other parts of the workflow, it is safe to assume
# that everything will work out ok. If not, this rule will fail, and the workflow can be restarted
# after ODB download has finished.
rule aggregate_results:
    input: 
        rnaseq_done = "data/checkpoints_dataprep/{taxon}_B13_cleanup_rnaseq.done",
        genome_done = "data/checkpoints_dataprep/{taxon}_A10_delete_ncbi_readme.done",
        annotated_tbl_path = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl",
        blank_tbl_path = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl",
        annot_done = "data/checkpoints_dataprep/{taxon}_A13_gtf.done"
    output:
        table = "data/checkpoints_dataprep/{taxon}_C01_data.csv"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        input_csv = config['INPUT']['input_csv'],
        odb_dir = config['BRAKER']['orthodb_path']
    wildcard_constraints:
        taxon="[^_]+"
    run:
        # read the input_csv to get the odb partition for each taxon
        in_csv = pd.read_csv(params.input_csv, header=None, sep=' ', names=['taxa', 'odb_partition', 'busco_db'])
        # store busco_db
        busco_db = in_csv[in_csv['taxa'] == params.taxon]['busco_db'].values[0]
        # find out what odb parition applies to the currently given taxon
        odb_partition = in_csv[in_csv['taxa'] == params.taxon]['odb_partition'].values[0]
        # build path to odb file, they sit in data/params.odbdir/odb_partition.fa
        odb_file = params.odb_dir + "/" + odb_partition + ".fa"
        # check if the odb file exists, use try except and die if not
        try:
            with open(odb_file):
                pass
        except IOError:
            raise Exception("OrthoDB file not found for taxon: " + odb_file)
        # read the species tables
        df_anno = pd.read_csv(input.annotated_tbl_path, sep="\t")
        df_blank = pd.read_csv(input.blank_tbl_path, sep="\t")
        # this df contains all species in a taxon
        df = pd.concat([df_anno, df_blank])
        # generate an empty dataframe called out_data with the columns species, taxon, accession, odb_file, rnaseq_file, genome_file, legacy_prot_file, annotation_file
        out_data = pd.DataFrame(columns=['species', 'taxon', 'accession', 'genome_file', 'odb_file', 'rnaseq_file', 'legacy_prot_file', 'annotation_file', 'busco_lineage'])
        # iterate over the df
        print(df)
        for index, row in df.iterrows():
            # get the species, replace spaces by underscores
            species = row['species'].replace(" ", "_")
            # get the accession
            accession = row['accession']
            # build a base path for the species
            base_path = "data/species/" + species + "/"
            # retrieve genome file location
            genome_file = base_path + "genome/genome.fa"
            try:
                with open(genome_file):
                    pass
            except IOError:
                raise Exception("Genome file not found for species: " + species)
            # check whether an rnaseq.bam file exists
            rnaseq_file_path = base_path + "hisat2/" + species + ".bam"
            rnaseq_file = ""
            if os.path.isfile(rnaseq_file_path):
                rnaseq_file = rnaseq_file_path
            # check whether legacy proteins exist for the species
            legacy_prot_file_path = base_path + "prot_legacy/proteins.fa"
            legacy_prot_file = ""
            if os.path.isfile(legacy_prot_file_path):
                legacy_prot_file = legacy_prot_file_path
            # check wheter a reference annotation exists for the species
            annotation_file_path = base_path + "annot/annot.gtf"
            annotation_file = ""
            if os.path.isfile(annotation_file_path):
                annotation_file = annotation_file_path
            # append the data to the out_data dataframe
            new_row = pd.DataFrame([{'species': species, 'taxon': params.taxon, 'accession': accession, 'genome_file': genome_file, 'odb_file': odb_file, 'rnaseq_file': rnaseq_file, 'legacy_prot_file': legacy_prot_file, 'annotation_file': annotation_file, 'busco_lineage': busco_db}])
            out_data = pd.concat([out_data, new_row], ignore_index=True)
        # write the out_data dataframe to a csv file
        out_data.to_csv(output.table, index=False)
            
