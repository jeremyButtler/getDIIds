# Use:

Detects DI (defective influenza) reads in long read
  sequenced samples.

This requires full length amplicons, so the old ONT rapid
  kits that fragmented DNA would not work. Also, short
  reads will not work.

Also, the IDs are based off the three nucleotides at the
  end of segments. So, there is some error.

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

# Updates:

- 2024-07-12:
  - fixed the HA and NS miscalls do to not realizing both
    had the same reverse primer sequence in getDIIds
  - added findDIFrag for fragments (not very good, but
    then everything here is so/so at best)
  - added getDICoords to scan sam files for DI events and
    their coordinates
  - allowed the aligners to be run separately and added in
    a waterman
    
# Programs:

## flu programs:

- getDIIds:
  - uses 12 to 13 conserved nucleotides at 5' and 3'
    ends to identify segments. DI reads are called if
    less than 50 nucleotides of the expected length
  - the problem here is when id sights identify the
    segment. also the method of calling DI reads is
    rough.
- findDIFrag:
  - aligns reads to reference using a waterman and then
    searches alignments to find number of DIs (50 or more
    nucleotide deletions).
  - this counts the number of shared 7mers between the
    references and read to find the best alignment
    direction and reference
  - the number of DIs are  found by looking for 20 or more
    nucleotide deletions that are at least 12 nucleotides
    away from the ends
  - default output is tsv file with number of detected
    DIs per read an alignment stats, but it can also
    output a sam file
  - problem is this may not always pick the true reference
  - problem the low gap extension penalty may sometimes
    create false deletions or may add the last few bases
    to the end
  - problem; this is really slow
- getDiCoords:
  - this is the DI detection part of findDIFrag, only
    the output prints the coordinates
  - outputs one row per DI event
    - each row has the read id, reference id, number DI
      events for the read, and the coordinates of one DI
      event
    - one read with multiple DI events will have multiple
      rows 

## alignment side programs:

These are nothing special, just standard waterman
  aligners. They do leverage branch less operations to
  reduce losses from branch prediction messing up
  maximizing scores. However, they are still much slower
  than the striped waterman programs. These will also do
  a kmer counting step to find if the read is reverse
  complement or forward.

For the alignment programs the default matrix is set up
  for long read mapping. The scores are adjusted by
  divding by 100.

You can change the scoring matrix or what counts as a
  match by providing a .txt file. See scoring-matrix.txt
  or match-matrix.txt for an example.
  
All these aligners only keep track of one row of scores
  at a time. So, most of the memory cost is from the
  directional matrix.

TODO: add in an alignment ouput option

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

- generalLib:
  - general libraries I have built up to support this code
- generalAln:
  - general libraries shared between pairwise aligners

# Installing everything

## non-static build (any unix system)

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/getDIIdsSrc
make all
sudo make allInstall
```

## static build (non-MAC)

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/getDIIdsSrc
make staticAll
sudo make allInstall
```

# Run:

You can run using `getDIIds -fq reads.fastq > ids.tsv`.
  For the help message do `getDIIds -h`.

# How works:

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

# Output:

This outputs a tsv (tab deliminated file) with the read
  ids that were kept.

- Columns:
  - Column 1 has read id
  - Column 2 has segment read is from
    - disagreements is "forward-segment/reverse-segment"
  - Column 3 has classification (vRNA/diRNA/mvRNA)
    - disagreements is "forward-class/reverse-class"
  - column 4 is the mapped length
    - distance between primer ends
  - column 5 is the read length
  - column 6 is the segments expected length
    - disagreements is "forward-length/reverse-length"
  - column 7 is
    - both for both primers supporting
    - for for forward primer support only
    - rev for forward primer support only
    - diff if primer positions disagree
  - remaining columns are primers stats

# Side program (primFind):

primFind is the generalized primer mapping step used
  in getDIIds.
It uses the primer search methods of getDIIds, but allows
  more than one pair of primers.

## Install primFind

### Static Linux build

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/getDIIdsSrc
make buildPrim
sudo make install
```

### Mac/non static Linux

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/getDIIdsSrc
make buildPrimMac
sudo make install
```

# Credits

- The idea of making a DI detection program was from
  Eric Bortz.
- The TB crew (Tara Ness and Bryce Inman) for being part
  of the TB program. This is the freezeTB project that
  contributed the source code for the primer search step.
- Family for being there.
