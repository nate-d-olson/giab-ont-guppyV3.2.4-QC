## Sync Seq Summary files
rsync -rvm --delete \
	--include="*/" --include="*sequencing_summary.txt" --exclude="*" \
	sherlock:/scratch/groups/msalit/nanopore/processing/guppy-3.2.4-snakemake-pipe/ \
	data/guppy_V3.2.4/sequencing_summary/

## Sync BAM stats files
rsync -rvm --delete \
	--include="*/" --include="*.bam.stats.tsv.gz" --exclude="*" \
	sherlock:/scratch/groups/msalit/nanopore/processing/guppy-3.2.4-snakemake-pipe/ \
	data/guppy_V3.2.4/bam_stats/