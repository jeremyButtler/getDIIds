# Use:

Detects DI (defective influenza) reads in sequenced
  samples.

getDIIds (segment length reads only) and findDIFrag both
  use different methods to identify DI reads. So, they
  will likely not agree.

# License:

This coded is dual licensed under public domain and MIT
  license. Public domain is the default license, unless
  you, your institution/company, your government, or 
  any one else has a problem with it. Then it defaults to
  MIT.

# Install:

Requires some kind of C compiler. Default is `cc`.

## non-static build (any unix system)

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/getDIIdsSrc
make core
sudo make install
```

## static build (non-MAC)

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/getDIIdsSrc
make staticCore
sudo make install
```

## local install

During the install step the path is searched for the
  install location. If it is not found, then the installs
  bash.rc file is updated to export the old path plus the
  new install location. So, you can do a local install
  by changing the install location. This is likely a bad
  idea, but also reduces the knowledge needed to do
  a local install.

```
make install PREFIX=/path/to/location
```

# Core programs:

## getDIIds

### getDIIds Overview and usage

getDIIds uses 12 to 13 conserved nucleotides at 5' and 3'
  ends to identify segments. DI reads are called if the
  region between the primers has 50 or fewer nucleotides
  than the expected segment length.

Problem could come from id sights identify wrong segments.
  So, chimeric reads will be a nightmare. Also the method
  of calling DI reads is rough.

Run by `getDIIds -fq reads.fastq > out.tsv`. You can
  get a help message with `getDIIds -h`.

This outputs a tsv (tab deliminated file) with the read
  ids that were kept.

### getDIIds Output

- Columns:
  - Column 1 has read id
  - Column 2 has segment read is from
  - Column 3 has classification (vRNA/diRNA/mvRNA)
  - column 4 is the mapped length
    - distance between primer ends
  - column 5 is the read length
  - column 6 is the segments expected length
  - column 7 is
    - both for both primers supporting
    - for for forward primer support only
    - rev for forward primer support only
  - remaining columns are primers stats

### How getDIIds works

1. Find universal 5' and 3' primer sites in reads
   - reads missing one of these sites are discarded
   - uses the modified (updated/bug fixed) memwater
     from my alnSeq repository
   - or uses kmerFind code (uses waterman code) from
     tbSpoligo, which is a module from freezeTB
2. The segment the read is from is found by comparing the
   bases at the end of the primer sites
   - Reads that have no segment id are discarded
   - Reads that have differing segment ids are discarded
     (unless `-diff` is used).
   - Reads that segments that could only be identified by
     one site are kept (unless `-no-part` is used)
   - Reads were both segments agree are kept
3. The length of the mapping region between primers
   (end of reverse primer - start of forward primer) is
   found
   - reads having a mapping length less then 85% of the
     read are discarded
   - reads with a mapping length that is 10% larger (110%)
     than the expected segment length are discarded
4. Kept read ids are printed out

## findDIFrag

### findDIFrag overview

findDIFrag aligns reads to a set of user provided
  references using a waterman and then searches
  alignments to find number of DIs.

Problems include that this may not always pick the true
  reference. That the low gap extension penalty may
  sometimes create false deletions or may add the last
  few bases to the end. And that it is really slow (this
  one I can probably make better with some time).

Run by `findDIFrag -fq reads.fastq -ref refs.fasta -out-sam out.sam > out.tsv`.
  You can get a help message with `findIDFrag -h`.

The default output is a tsv file, but you can also output
  a sam file (includes tsv file) with -out-sam. The sam
  file can be passed to stdout by using `-out file.tsv`
  and `-sam-out -`.

The sam file does not have mapping qualities, but does
  include the alignment score as the `AS` (AS:i:) entry
  (number 12).

### findDIFrag output

- tsv columns:
  - Column 1 has read id
  - Column 2 has reference name
  - Column 3 has classification (diRNA/vRNA)
  - Column 4 has number of detected DI events
  - Column 5 has the direction (forward/reverse)
  - Column 6 has the read length
  - Column 7 has the aligned length (number reference
    bases covered)
  - Columns 8 and 9 have the mapping coordinates
  - Column 10 has the reference length
  - Column 11 has the alignment score
  - Column 12 has the maximum possible alignment score
    for the read
  - Columns 11 to 15 has the number of matches, SNPs,
    insertions, deletions, and masked bases
  - Column 16 has the number of shared kmers with the
    reference (kmer size is in column name)
  - Column 17 and 18 have the mean and median Q-scores
    - this is often over inflated (not sure if I am off or
      if it really is overinflated)

### How findDIFrag works

1. Counts number of kmers shared between references and
   read
   - The reference (and direction) sharing the most kmers
     with the read is kept
2. Checks to see if (shared kmers) / (reads total kmers)
   is over 40% (otherwise assumes no alignment)
3. Kept reads are mapped to best reference using a
   Waterman Smith alignment with a gap extension penalty
   of -0.05 (match = 5, snp = -4, gap open = -10)
4. Reads with alignment scores under 50% of maximum score
   for a read are discarded
5. The alignment is extracted from the directional matrix
   for the waterman as sam file entry
6. The sam entry cigar is scanned for large deletions (DI)
   (>= 20) that are at least 12 nucleotides away from
   either end
7. The number of detected DI events and other alignment
   information is printed out (+ sam file if requested)

## getDICoords

### getDICoords overview

getDiCoords is the DI detection part of findDIFrag, only
    the output prints the coordinates 

One problem here is that this can not detect missing
  reads. So, make sure your read mapper can find DI reads.

Run by `getDICoords -sam reads.sam > out.tsv`. You can get
  a help message with `getDICoords -h`.

getDICoords outputs a tsv with the coordinates of detected
  DI events. Each row represents one detected DI event,
  so there can be multiple rows per read.

### getDICoords output

- Columns:
  - Column 1 has the read id
  - Column 2 has the reference id
  - Column 3 has the number of DI events (rows) for this
    read
  - Column 4 is the start of one DI event
  - Column 5 is the end of one DI event

### How getDICoords works:

Basically step 6 of findDIFrag, with coordinate recording.

1. Read entry from sam file
2. Entry cigar is scanned for large deletions (DI)
   (>= 20) that are at least 12 nucleotides away from
   either end
3. The start and end of all DI events are recorded
4. Print read id, number DI events, and coordinates
   for each DI event

## samClust

Currently thinking about (no were near working). Goal is a
  clustering program for reads.

# Updates:

- 2024-07-16:
  - Added in a binning program to bin reads in a sam file
    by reference. This will be my first step in consensus
    building.

- 2024-07-15:
  - various small changes to tsv outputs (program
    versions not changed)
  - findDIFrag (only update):
    - changed references from being hardcoded to being
      user provided
      - providing multiple references fixed the problem
        with HA segments (reference was to distant)
    - improved speed (still slow) by filtering out reads
      with few kmers shared with references (less
      alignments)
    - added in options to control min score (as percent),
      min shared kmers (as percent), and kmer length

- 2024-07-12:
  - fixed the HA and NS miscalls do to not realizing both
    had the same reverse primer sequence in getDIIds
  - added findDIFrag for fragments (not very good, but
    then everything here is so/so at best)
  - added getDICoords to scan sam files for DI events and
    their coordinates
  - allowed the aligners to be run separately and added in
    a waterman

# Side programs:

## Installing everything (side programs)

### non-static build (any unix system)

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/getDIIdsSrc
make all
sudo make allInstall
```

### static build (non-MAC)

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/getDIIdsSrc
make staticAll
sudo make allInstall
```

## primFind (primer scanning)

### primFind overview

primFind is the generalized primer mapping step used
  in getDIIds.
It uses the primer search methods of getDIIds, but allows
  more than one pair of primers.

primFind works for the flu omni primers (not best), but
  really is step up for 18 to 25 nucleotide long primers
  (event then is likely not very good). It is setup for
  mapping very short sequences to very big sequences.

Use primFind
   with `primFind -fa primers.fa -fq reads.fq > out.tsv`.
   You can get the help message with `primFind -h`.

### primFind Output

- Output:
  - Column 1 is read id
  - Column 2 is primer id (name)
  - Column 3 is the aligned length (primer end to end)
  - Column 4 is the forward primers direction
  - Column 5 is the forward starting position
  - Column 6 is the forward ending position
  - Columns 7 to 9 are columns 4 to 6 for reverse primer
  - Column 8 is the forward alignment score
  - Column 9 is the maximum possible forward score
  - Column 10 is the number of times the primer mapped
    to the sequence (was above min score)
  - Column 11 is first foward primer mapped base
  - Column 12 is last foward primer mapped base
  - Columns 13 to 17 are columns 8 to 12 for reverse
    primer

### How primFind works

1. Find 5mer counts in primer
   - the kmer size is important for speed. 3mers are
     really slow (slower than doing a Waterman).
2. Split read into (maximum primer length) windows
3. Get 5mer counts in two windows
4. Find shared 5mers between both windows
5. If enough 5mers are shared, do Waterman alignment
   (memwater)
   - Else move on
6. If score is high enough, then keep coordinates and
   score if better than current best score (only one
   alignment kept)
7. Move down one window (one maximum primer length) and
   repeat steps 1 to 6

## samRefBin (read bining)

samRefBin is a program designed to split (bin) the reads
  in a sam file by the mapped reference. It will likely
  not be very usefull, since their are tools that already
  do this for bam files (bamtools).

samRefBin is mainly here to provided support for my other
  programs. Mainly the consensus building step I hope to
  make.

You can use samRefBin
  by `samRefBin -sam file.sam -prefix out`. 

- Output:
  - A sam file for each reference that has the reads that
    mapped to the sam file reference.

## alignment side programs:

These are nothing special, just standard waterman
  aligners. They do leverage branch less operations to
  reduce losses from branch prediction messing up
  maximizing scores. However, they are still much slower
  than the striped waterman programs. These will also do
  a kmer counting step to find if the read is reverse
  complement or forward.

For the alignment programs the default matrix is set up
  for long read DI mapping. The scores are adjusted by
  divding by 100.

You can change the scoring matrix or what counts as a
  match by providing a .txt file. See scoring-matrix.txt
  or match-matrix.txt for an example.
  
All these aligners only keep track of one row of scores
  at a time. So, most of the memory cost is from the
  directional matrix.

- alnwater:
  - does a waterman alignment
    - slow when compared to striped watermans
  - output is currently sam file only
- alnMemwater:
  - find the score, starting coordinate, and ending
    coordiante of highest scoring alignment.
    - this also converts the scoring to a single row. So,
      is very memory light, but is a little slower than
      alnwater (keeping track of starts slows things down)
  - TODO: find why alnMemwater and alnwater disagree on
    coordinates for larger alignments

## can not compile:

These are libraries I have built up to make my life
  easier. These may change between repositories as I
  become more used to coding and develop my style more.

- generalLib:
  - general libraries shared between programs
- generalAln:
  - general libraries shared between pairwise aligners

# Credits

- The idea of making a DI detection program was from
  Eric Bortz.
- The TB crew (Tara Ness and Bryce Inman) for being part
  of the TB program. This is the freezeTB project that
  contributed the source code for the primer search step.
- Family for being there.
