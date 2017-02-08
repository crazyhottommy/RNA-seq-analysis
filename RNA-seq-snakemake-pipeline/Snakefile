
import json
from os.path import join, basename, dirname

# globals ----------------------------

configfile: 'config.yml'
# Full path to an uncompressed FASTA file with all chromosome sequences.
CDNA = config['CDNA']

# Full path to a folder where intermediate output files will be created.
OUT_DIR = config['OUT_DIR']

FILES = json.load(open(config['SAMPLES_JSON']))

SAMPLES = sorted(FILES.keys())


# Rules ------------------------------

rule all:
	input:
		join(dirname(CDNA), 'salmon', basename(CDNA).rstrip(".fa")),
		[OUT_DIR + "/" + x for x in expand('{sample}/quant.sf', sample = SAMPLES)],
		[OUT_DIR + "/" + x for x in expand('{sample}/lib_format_counts.json', sample = SAMPLES)],
		'quant.tsv.gz'

rule salmon_index:
	input:
		cdna = CDNA
	output:
		index = join(dirname(CDNA), 'salmon', basename(CDNA).rstrip(".fa"))
	log:
		 'logs/salmon_index.log'
	shell:
		"""	
		salmon index -t {input} -i {output} &> {log}
		"""

rule salmon_quant:
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2'],
		index = rules.salmon_index.output.index
	output:
		join(OUT_DIR, '{sample}', 'quant.sf'),
		join(OUT_DIR, '{sample}', 'lib_format_counts.json')
	log:
		'logs/{sample}_salmons_quant.log'
	shell:
		"""
		salmon quant -p 4 -i {input.index} -l ISR -1 <(gunzip -c {input.r1}) -2 <(gunzip -c {input.r2}) -o OUT_DIR/{wildcards.sample} &> {log}'
		"""

rule collate_salmon:
    input:
        expand(join(OUT_DIR, '{sample}', 'quant.sf'), sample= SAMPLES)
    output:
        'quant.tsv.gz'
    run:
        import gzip

        b = lambda x: bytes(x, 'UTF8')

        # Create the output file.
        with gzip.open(output[0], 'wb') as out:

            # Print the header.
            header = open(input[0]).readline()
            out.write(b('sample\t' + header))

            for i in input:
                sample = basename(dirname(i))
                lines = open(i)
                # Skip the header in each file.
                lines.readline()
                for line in lines:
                    out.write(b(sample + '\t' + line))