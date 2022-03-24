# Full-length 16S rRNA gene amplicons from raw sewage spanning years across USA

Lou LaMartina, Angie Schmoldt, Ryan Newton


## Description

We generated full-length 16S rRNA gene sequences in wastewater treatment plant influent from a 5-year Milwaukee time series and USA-wide sample set. Samples were PCR-amplified with unique barcodes appended to 5' ends of the forward and reverse primers. Sequencing was done on a PacBio Sequel II in a single library and reads were demultiplexed using unique barcode combinations.


![image](https://github.com/loulanomics/Full16S_sewageDatabase/blob/main/Figures/dendrogram.png)


## Github contents
- <b>Sample metadata. </b> Includes location, lat/lon, sampling dates, and NCBI accession numbers.

- <b>ASV counts. </b> Number of ASV read counts by sample name.

- <b>ASV taxonomy. </b> Taxonomic assignments (from Kingdom to species) according to the [Silva v.138](https://www.arb-silva.de/documentation/release-138/) reference database.

- <b>FASTA. </b> Headers show ASV ID, taxonomic assignments, and read count in dataset. Both forward & reverse strands included.

- <b>Phyloseq object. </b> For easy exploration in R ([website](https://joey711.github.io/phyloseq/)).

- <b>Galaxy. </b> The [Galaxy online tool](https://usegalaxy.org) was used for the computationally expensive DADA2 steps. Here are analysis details and outputs.


### Curated vs. total

After processing in [DADA2](https://benjjneb.github.io/dada2/tutorial.html), we clustered ASVs to 99.5% similarity with [mothur](https://mothur.org/wiki/cluster/) to correct for sequencing errors potentially contributing to erroneous ASVs.


## Primers & barcodes

- 12 forward [barcodes](https://github.com/PacificBiosciences/Bioinformatics-Training/blob/master/barcoding/pacbio_384_barcodes.fasta) (0007_Forward - 0018_Forward)
- 4 reverse barcodes (reverse complements; 0004_Forward - 0007_Forward)
- 27F:AGRGTTYGATYMTGGCTCAG and 1492R:RGYTACCTTGTTACGACTT universal primer set


## Raw FASTQ availability

- NCBI Short Read Archive [PRJNA809416](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA809416)
- Corresponding V4 16S rRNA gene region archive [PRJNA597057](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA597057)



## FASTA format

Headers: [ASV ID] [kingdom] [phylum] [class] [order] [family] [genus] [species] [read count in dataset] : [read (1 or 2)]

For example...

```
>ASV0001__k_Bacteria__p_Proteobacteria__c_Gammaproteobacteria__o_Pseudomonadales__f_Moraxellaceae__g_Acinetobacter__g_johnsonii__count_8486:R2
GGTTAC...CAAACC
>ASV0001__k_Bacteria__p_Proteobacteria__c_Gammaproteobacteria__o_Pseudomonadales__f_Moraxellaceae__g_Acinetobacter__g_johnsonii__count_8486:R2
GGTTTG...GTAACC
>ASV0002__k_Bacteria__p_Firmicutes__c_Bacilli__o_Lactobacillales__f_Carnobacteriaceae__g_Trichococcus__g___count_6199:R2
GGTTAC...CAAACC
>ASV0002__k_Bacteria__p_Firmicutes__c_Bacilli__o_Lactobacillales__f_Carnobacteriaceae__g_Trichococcus__g___count_6199:R2
GGTTTG...GTAACC
```



## Sewage samples

- 22 wastewater treatment plants across USA ([LaMartina et al., 2021](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01038-5); [Newton et al., 2015](https://journals.asm.org/doi/10.1128/mBio.02574-14))
- 2 treatment plants from Milwaukee sampled monthly for five years ([LaMartina et al., 2021](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01038-5))





| Location | Date |
| :--- | :--- |
| Auburn, Alabama | 	8/14/12 |
| Austin, Texas | 	8/15/12 |
| Delano, Minnesota | 	1/16/13 |
| Denver, Colorado | 	1/23/13 |
| Great Falls, Montana | 	1/16/13 |
| Honolulu, Hawaii | 	9/7/12 |
| Iowa City, Iowa | 	1/15/13 |
| Johns Creek, Georgia | 	8/16/12 |
| Juneau, Alaska | 	1/23/13 |
| Kenedy, Texas | 	8/30/12 |
| Laramie, Wyoming | 	1/24/13 |
| Lincoln, Nebraska | 	1/23/13 |
| Madison, Wisconsin | 	1/27/13 |
| Marathon, Florida | 	8/21/12 |
| Memphis, Tennessee | 	8/15/12 |
| Portland, Oregon | 	1/16/13 |
| Stockton, California | 	8/14/12 |
| Vancouver, Washington | 	1/15/13 |
| West Palm Beach, Florida | 	8/8/12 |
| Whittier, California | 	8/21/12 |
| Woodmere, Ohio | 	1/17/13 |
| Yuma, Arizona | 	8/15/12 |
| Milwaukee, Wisconsin | 	4/7/16 |
| Milwaukee, Wisconsin | 	4/3/17 |
| Milwaukee, Wisconsin | 	8/3/16 |
| Milwaukee, Wisconsin | 	8/22/17 |
| Milwaukee, Wisconsin | 	12/7/16 |
| Milwaukee, Wisconsin | 	12/1/17 |
| Milwaukee, Wisconsin | 	2/8/16 |
| Milwaukee, Wisconsin | 	2/6/17 |
| Milwaukee, Wisconsin | 	1/6/16 |
| Milwaukee, Wisconsin | 	1/5/17 |
| Milwaukee, Wisconsin | 	7/18/16 |
| Milwaukee, Wisconsin | 	7/12/17 |
| Milwaukee, Wisconsin | 	6/8/16 |
| Milwaukee, Wisconsin | 	6/7/17 |
| Milwaukee, Wisconsin | 	3/2/16 |
| Milwaukee, Wisconsin | 	3/1/17 |
| Milwaukee, Wisconsin | 	5/2/16 |
| Milwaukee, Wisconsin | 	5/1/17 |
| Milwaukee, Wisconsin | 	11/3/16 |
| Milwaukee, Wisconsin | 	11/2/17 |
| Milwaukee, Wisconsin | 	10/5/16 |
| Milwaukee, Wisconsin | 	10/4/17 |
| Milwaukee, Wisconsin | 	9/21/16 |
| Milwaukee, Wisconsin | 	9/26/17 |







