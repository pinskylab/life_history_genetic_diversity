# **Microsatellite data entry**

05.24.2020

## **Do not enter data from:**

1. Captive, farmed, stocked, or any populations that aren't wild
2. EST-linked microsatellites
   * EST = Expressed Sequence Tag. We don't use EST-linked microsatellites as they are bits of genetic material copied from an organisms' transcriptome (mRNA) and is likely to be closely associated with genes that matter to the organisms' fitness and may be under natural selection. If the study uses EST-linked microsatellites, it will probably say so in the abstract. If it doesn't mention anything, the study is probably OK.
3. Loci that are not in Hardy-Weinberg equilibrium (they may be under selection)
4. Anything that isn't a microsatellite
   * Minisatellites, RFLPs, AFLPs, SNPs, etc. are NOT microsatellites
   * SSRs or STRs (Simple short repeats or short random repeats) ARE microsatellites
5. Species that are anadromous, catadromous, or freshwater (anything that is not marine)
6. Species that are not fish
   * Sharks and seahorses are fish
7. Microsatellites that have only one allele (are monomorphic)
8. Data on diversity averaged across a wide area (>3 degrees latitude or longitude), or where the site is only vaguely identified (accuracy worse than 3 degrees latitude or longitude)
9. Multiple samples from the same site (eg: in different years, if each year is reported separately)
10. A site that has already been reported. For example, a paper may report samples and data that were first analyzed and presented in another paper that you either have or will enter. If the data are the same, enter the data from the original paper.


## **Data entry instructions:**

In the spreadsheet for data entry, each line is the diversity for one microsatellite at one site in one species. If the paper only reports averages across multiple microsatellite loci, then the line will have that average.

If you have to choose between entering averaged data for many microsatellites at one site, or averaged data for many sites at one microsatellite, choose the average across microsatellites at one site. We care more about the spatial variation in diversity.

Sometimes, papers will report averaged microsatellite diversity at a site in the main text, but will break the diversity down by microsatellite per site in the supplemental material. Be sure to always check the supplemental material first (if it exists) before reporting any averages.

If a paper reports temporal replicates (i.e. same site is sampled in multiple years and reported separately), only report ONE time point per site. Pick the temporal replicate with the largest sample size. If this is a tie, pick the most recent sample. (ex: if 50 individuals were sampled from site A in 2012, and 25 were sampled from site A in 2014, report the diversity from the 2012 sample).

## **Data to enter:**

1. Spp
   * Scientific name. Check on Fishbase (www.fishbase.org) to make sure you have the most up-to-date scientific name, as it may differ from what is in the paper.
2. CommonName
   * Also check Fishbase for the most up-to-date/widely used common name.
3. Source
   * Enter in the form "Adams et al. (2010) Molecular Ecology 8:1000-1015,"
     * 8 is the volume and 1000-1015 are the page numbers
	   * If only one, two, or three authors, use a citation similar to Adams, Adams & Benner, or Adams, Benner & Hadley, respectively. More than three authors requires "et al."
4. PrimerNote
   * 1 for primer note, 0 for not
   * Some studies only exist to report the discovery of new microsatellites. These are called Primer Notes. Often they are labeled as such, but not always. They tend to be short and only evaluate a small number of individuals. Many of them are in journals such as Molecular Ecology Resources, Molecular Ecology NOtes, or Conservation Genetics Resources.
6. Country
   * Country where the sample was taken. If sampled in the middle of an ocean, record the appropirate ocean basin.
7. Lat_Deg
   * Degrees latitude. This can be decimal degrees. Use negative numbers for the Southern Hemisphere
   * If latitude coordinates are provided, use those. If they are not provided but the site name is specific enough that you can search for the coordinates on Google Maps, then do so. Right click and select "What's Here?" to get Google to give you the lat and long coordinates.
8. Lat_Min
   * Minutes latitude. Use this field only if you don't have decimal degrees for Lat_Deg.
9. Lat_Sec
   * Seconds latitude. Use this field only if you have neither decimal degrees nor decimal minutes for Lat_Deg and Lat_Min
10. Lon_Deg
    * Degrees longitude. Same rules as Lat_Deg.
11. Lon_Min
    * Minutes longitude. Same rules as Lat_Min.
12. Lon_Sec/Users/marialmalabag/Downloads/Marial_Papers/msat_entry_instructions.md
    * Seconds longitude. Same rules as Lat_Sec.
13. CollectionYear
    * Year(s) in which the samples were taken. Leave blank if not stated.
14. NumMarkers
    * Number of microsatellite markers whose data are entered on this line. Will often be 1. It will only be >1 if the paper only reports an average across multiple loci.
15. MarkerName
    * Name of the microsatellite marker, as listed in the paper.
    * Leave blank if there are multiple markers reported together on this line.
16. CrossSpp
    * Was the microsatellite originally developed in a different species? 1 if yes, 0 if no.
    * You may have to find the original reference for a microsatellite locus to learn whether it was developed in a different species.
    * If the line of data is for a number of markers averaged together and some were cross-species and some were not, input the proportion that were cross-speciese (e.g. 0.2 if 1 of 5 were cross-species.)
17. n
    * Number of individuals sampled. If <4, do not record the site.
    * If NumMarkers>1, then this might have to be an average across all of the markers reported on this line.
18. Repeat
    * The length of the microsatellite repeat. Will be 2-5 (typically).
    * You may have to find the original reference for a microsatellite locus to learn the repeat length.
    * A 2-nucleotide repeat is called a dinucleotide, 3 is called trinucleotide, 4 a tetranucleotide, etc. If the repeat is complex (has di's, tri's and tetra's, or some mix), then use the most common repeat motif length. For example, (GA)7(GACA)9(CT)6 would be 2, since 7+6 > 9.
    * If NumMarkers >1, then this will be an average across all the repeat sizes of those markers.
19. He
    * Expected heterozygosity (sometimes called gene diversity).
    * This will often be reported in a table. Sometimes it will be in the supplementary material for a paper.
20. Hese
    * Standard error of the He measurement.
    * Only enter this if NumMarkers >1 and the paper reports it.
    * If the paper reports standard deviation, you can convert to standard error by dividing by sqrt(n).
22. Notes
    * Anything of interest
    * If you pulled the latitude and longitude coordinates off of Google Maps, record so here by typing "geo coordinates pulled from place name"


## **Keep notes!**

Keep a dated logbook (a text file is good) where you make notes.

Every day you record a note, use the following format:

**DD.MM.YYY**

**Notes:**
 * Source A: Enter notes here (separated by ;)
 * Source B: Enter notes here

Things to note: 
 * If you did not include some microsatellite markers because they were out of HWE or because they were monomorphic, record that here (ex: MicrosatA (sample site it was excluded from) excluded bc not in HWE)
 * If you had to exclude a paper (did not pull any data from it) record why you did so
 * If you had to exclude some sites (either because they were temporal replicates, previously recorded, were from aquaculture sites, etc.) record so here

