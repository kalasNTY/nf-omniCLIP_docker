nextflow.enable.dsl = 2

// Define parameters directly
params.metadata = "./metadata.csv"
params.db_file = "/home/julian/DATA/OC/GRCh38_113.db"
params.bg_f = "/home/julian/DATA/OC/bg_dat"
params.clip_f = "/home/julian/DATA/OC/clip_dat"
params.genome_dir = "/home/julian/DATA/OC/hg38_noChr"
params.gff_file = "/home/julian/DATA/OC/Homo_sapiens.GRCh38.113.chromosome.3.gff3"
params.output_dir = "/home/julian/DATA/OC/Results"
params.max_it = 1
params.nb_cores = 1
params.omniclip_container = "kalasnty/omniclip_container:latest"
params.uid = 1000  // Replace with `id -u` output
params.gid = 1000  // Replace with `id -g` output

// Define all processes first

process PULL_CONTAINER {
    output:
    val true, emit: container_pulled

    script:
    """
    docker pull $params.omniclip_container
    """
}

process GENERATE_DB {
    input:
    path gff_file 
    path db_file

    output:
    path db_file, emit: db_ch

    script:
    """
    if [ ! -f $db_file ]; then
        docker run --rm -v ~/DATA/OC:/data --user ${params.uid}:${params.gid} $params.omniclip_container \
        generateDB --gff-file /data/${gff_file.name} --db-file /data/${db_file.name}
    fi
    """
}

process PARSING_BG {
    publishDir params.output_dir, mode: 'copy'

    input:
    path db_file
    path genome_dir
    path bg_file
    path bg_f

    output:
    path bg_f, emit: bg_dat_ch

    script:
    """
    if [ ! -f $bg_f ]; then
        docker run --rm -v ~/DATA/OC:/data --user ${params.uid}:${params.gid} $params.omniclip_container \
        parsingBG --db-file /data/${db_file.name} --genome-dir /data/${genome_dir.name} \
        ${bg_file.collect { "--bg-files /data/Input/${it.name}" }.join(" ")} \
        --out-file /data/${bg_f.name}
    fi
    """
}

process PARSING_CLIP {
    publishDir params.output_dir, mode: 'copy'

    input:
    path db_file
    path genome_dir
    path clip_file
    path clip_f

    output:
    path clip_f, emit: clip_dat_ch

    script:
    """
    if [ ! -f $clip_f ]; then
        docker run --rm -v ~/DATA/OC:/data --user ${params.uid}:${params.gid} $params.omniclip_container \
        parsingCLIP --db-file /data/${db_file.name} --genome-dir /data/${genome_dir.name} \
        ${clip_file.collect { "--clip-files /data/Input/${it.name}" }.join(" ")} \
        --out-file /data/${clip_f.name}
    fi
    """
}


process RUN_OMNICLIP {
    input:
    path db_file 
    path bg_dat_ch 
    path clip_dat_ch 
    path output_dir 
    val max_it 
    val nb_cores

    output:
    path output_dir

    script:
    """
    docker run --rm -v ~/DATA/OC:/data --user ${params.uid}:${params.gid} $params.omniclip_container \
    run_omniCLIP --db-file /data/${db_file.name} --bg-dat /data/${bg_dat_ch.name} --clip-dat /data/${clip_dat_ch.name} \
    --out-dir /data/${output_dir.name} --max-it $max_it --nb-cores $nb_cores
    """
}

workflow {
    def container_pulled = PULL_CONTAINER() 
    // Read metadata file and group files by Factor (bg or clip)
    def samples = Channel.fromPath(params.metadata)
        | splitCsv(header: true)
        | map { row -> tuple(row.Sample, row.Factor == "bg" ? file(row.File) : null, row.Factor == "clip" ? file(row.File) : null) }

    def db_ch = GENERATE_DB(params.gff_file, params.db_file)

    def bg_dat_ch = PARSING_BG(db_ch, params.genome_dir, samples.filter { it[1] != null }.map { it[1] }.collect(), params.bg_f)
    def clip_dat_ch = PARSING_CLIP(db_ch, params.genome_dir, samples.filter { it[2] != null }.map { it[2] }.collect(),params.clip_f)

    RUN_OMNICLIP(db_ch, bg_dat_ch, clip_dat_ch, params.output_dir, params.max_it, params.nb_cores)
}
