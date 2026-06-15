/*
 * sab-v1.c - strand-allele-balance
 *
 * For each heterozygous site in a VCF/BCF file for a given sample,
 * finds sites covered by exactly two reads (above map quality threshold)
 * in a BAM file, then tallies allele counts by strand orientation.
 *
 * Output: 3 rows (++, +-, --) x 3 columns (ref/ref, ref/alt, alt/alt)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"

int VERSION = 1;

/* Count matrix: rows = strand combo (0=++, 1=+-, 2=--), cols = allele combo (0=RR,1=RA,2=AA) */
static long counts[3][3];

typedef struct {
    int32_t ref_allele_len;
    int32_t alt_allele_len;
    char ref_base;
    char alt_base;
} HetSite;

/* Data passed to the pileup callback */
typedef struct {
    int pos;           /* 0-based position */
    HetSite *site;
    int map_qual_cutoff;
    /* results */
    int n_reads;
    int strand[2];     /* 0=plus, 1=minus for each read */
    int allele[2];     /* 0=ref, 1=alt for each read, -1=unknown */
} PileupData;

static int pileup_callback(void *data, bam1_t *b)
{
    (void)data;
    (void)b;
    return 0;
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "strand-allele-balance version %d\n"
        "Usage: %s -v <vcf/bcf> -b <bam> -I <sample_id> [-m <mapq>]\n"
        "  -v  VCF or BCF file with genotype information\n"
        "  -b  BAM file with aligned sequence reads\n"
        "  -I  Sample identifier in the VCF/BCF file\n"
        "  -m  Map quality cutoff (default: 20)\n",
        VERSION, prog);
}

int main(int argc, char *argv[])
{
    char *vcf_file = NULL;
    char *bam_file = NULL;
    char *sample_id = NULL;
    int map_qual = 20;
    int opt;

    while ((opt = getopt(argc, argv, "v:b:m:I:h")) != -1) {
        switch (opt) {
        case 'v': vcf_file   = optarg; break;
        case 'b': bam_file   = optarg; break;
        case 'm': map_qual   = atoi(optarg); break;
        case 'I': sample_id  = optarg; break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }

    if (!vcf_file || !bam_file || !sample_id) {
        fprintf(stderr, "Error: -v, -b, and -I are required.\n");
        usage(argv[0]);
        return 1;
    }

    /* ------------------------------------------------------------------ */
    /* Open VCF/BCF (HTSlib auto-detects format)                           */
    /* ------------------------------------------------------------------ */
    htsFile *vcf_fp = hts_open(vcf_file, "r");
    if (!vcf_fp) {
        fprintf(stderr, "Error: cannot open VCF/BCF file '%s'\n", vcf_file);
        return 1;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(vcf_fp);
    if (!hdr) {
        fprintf(stderr, "Error: cannot read VCF/BCF header\n");
        return 1;
    }

    /* Find sample index */
    int sample_idx = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, sample_id);
    if (sample_idx < 0) {
        fprintf(stderr, "Error: sample '%s' not found in VCF/BCF\n", sample_id);
        return 1;
    }

    /* ------------------------------------------------------------------ */
    /* Open BAM and its index                                              */
    /* ------------------------------------------------------------------ */
    htsFile *bam_fp = hts_open(bam_file, "r");
    if (!bam_fp) {
        fprintf(stderr, "Error: cannot open BAM file '%s'\n", bam_file);
        return 1;
    }

    sam_hdr_t *sam_hdr = sam_hdr_read(bam_fp);
    if (!sam_hdr) {
        fprintf(stderr, "Error: cannot read BAM header\n");
        return 1;
    }

    hts_idx_t *bam_idx = sam_index_load(bam_fp, bam_file);
    if (!bam_idx) {
        fprintf(stderr, "Error: cannot load BAM index for '%s'\n"
                        "       (run 'samtools index %s' first)\n",
                        bam_file, bam_file);
        return 1;
    }

    /* ------------------------------------------------------------------ */
    /* Iterate over VCF/BCF records                                        */
    /* ------------------------------------------------------------------ */
    bcf1_t *rec = bcf_init();
    memset(counts, 0, sizeof(counts));

    while (bcf_read(vcf_fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        /* Only biallelic SNPs (ref and alt both single bases) */
        if (rec->n_allele != 2) continue;
        const char *ref_str = rec->d.allele[0];
        const char *alt_str = rec->d.allele[1];
        if (strlen(ref_str) != 1 || strlen(alt_str) != 1) continue;

        /* Skip if alt is not a real base (e.g. '.') */
        char ref_base = ref_str[0];
        char alt_base = alt_str[0];
        if (alt_base == '.' || alt_base == '*') continue;

        /* Get genotype for our sample */
        int32_t *gt_arr = NULL;
        int ngt = 0;
        int ret = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt);
        if (ret < 0) continue;

        int n_per_sample = ngt / bcf_hdr_nsamples(hdr);
        int32_t *gt = gt_arr + sample_idx * n_per_sample;

        /* Must be diploid, phased or unphased, heterozygous */
        if (n_per_sample < 2) { free(gt_arr); continue; }
        if (bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1])) { free(gt_arr); continue; }

        int a0 = bcf_gt_allele(gt[0]);
        int a1 = bcf_gt_allele(gt[1]);
        int is_het = (a0 != a1) && (a0 == 0 || a0 == 1) && (a1 == 0 || a1 == 1);
        free(gt_arr);
        if (!is_het) continue;

        /* ---------------------------------------------------------------- */
        /* Query BAM at this position                                        */
        /* ---------------------------------------------------------------- */
        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        hts_pos_t pos0 = rec->pos;  /* 0-based */

        int tid = sam_hdr_name2tid(sam_hdr, chrom);
        if (tid < 0) continue;

        hts_itr_t *itr = sam_itr_queryi(bam_idx, tid, (int)pos0, (int)pos0 + 1);
        if (!itr) continue;

        /* Collect reads that cover this position */
        typedef struct { int strand; int allele; } ReadInfo;
        ReadInfo reads[64];
        int nreads = 0;

        bam1_t *b = bam_init1();
        while (sam_itr_next(bam_fp, itr, b) >= 0) {
            /* Filter: unmapped, secondary, duplicate, QC fail */
            if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP | BAM_FQCFAIL))
                continue;
            /* Map quality filter */
            if (b->core.qual < map_qual) continue;

            /* Find the query base at pos0 using CIGAR */
            int32_t ref_pos = b->core.pos; /* current ref position (0-based) */
            uint32_t *cigar = bam_get_cigar(b);
            uint8_t  *seq   = bam_get_seq(b);
            int query_pos = 0;
            int found = 0;
            char query_base = 0;

            for (uint32_t ci = 0; ci < (uint32_t)b->core.n_cigar; ci++) {
                int op  = bam_cigar_op(cigar[ci]);
                int len = bam_cigar_oplen(cigar[ci]);

                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    for (int k = 0; k < len; k++) {
                        if (ref_pos == (int32_t)pos0) {
                            query_base = seq_nt16_str[bam_seqi(seq, query_pos)];
                            found = 1;
                            break;
                        }
                        ref_pos++;
                        query_pos++;
                    }
                } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                    ref_pos += len;
                } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
                    query_pos += len;
                } else if (op == BAM_CHARD_CLIP || op == BAM_CPAD) {
                    /* nothing */
                }
                if (found) break;
                if (ref_pos > (int32_t)pos0) break;
            }

            if (!found) continue;

            /* Determine allele */
            int allele = -1;
            if (query_base == ref_base) allele = 0;
            else if (query_base == alt_base) allele = 1;
            else continue; /* neither ref nor alt — skip */

            /* Strand: 0=plus(forward), 1=minus(reverse) */
            int strand = (b->core.flag & BAM_FREVERSE) ? 1 : 0;

            if (nreads < 64) {
                reads[nreads].strand = strand;
                reads[nreads].allele = allele;
                nreads++;
            }
        }
        bam_destroy1(b);
        hts_itr_destroy(itr);

        /* We only care about sites with exactly 2 qualifying reads */
        if (nreads != 2) continue;

        /* Determine strand combination */
        int s0 = reads[0].strand; /* 0=+, 1=- */
        int s1 = reads[1].strand;
        int row;
        if      (s0 == 0 && s1 == 0) row = 0; /* ++ */
        else if (s0 == 1 && s1 == 1) row = 2; /* -- */
        else                          row = 1; /* +- */

        /* Determine allele combination */
        int tot_alt = reads[0].allele + reads[1].allele; /* 0, 1, or 2 */
        counts[row][tot_alt]++;
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(vcf_fp);
    sam_hdr_destroy(sam_hdr);
    hts_idx_destroy(bam_idx);
    hts_close(bam_fp);

    /* ------------------------------------------------------------------ */
    /* Print results: 3 rows x 3 columns                                   */
    /* ------------------------------------------------------------------ */
    for (int r = 0; r < 3; r++) {
        printf("%ld %ld %ld\n", counts[r][0], counts[r][1], counts[r][2]);
    }

    return 0;
}
