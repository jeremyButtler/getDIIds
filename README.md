# Use:

Detects DI (defective influenza) reads in long read
  sequenced samples.

This requires full length amplicons, so the old ONT rapid
  kits that fragmented DNA would not work. Also, short
  reads will not work.

# Install:

Requires some kind of C compiler. Default is `cc`.

## Static Linux build

```
git clone https://github.com/jeremybuttler/getDIIds
cd getDIIds/src
make
sudo make install
```

## Mac/non static Linux

```
cd getDIIds/src
git clone https://github.com/jeremybuttler/getDIIds
make mac
sudo make install
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
cd getDIIds/src
make buildPrim
sudo make install
```

### Mac/non static Linux

```
cd getDIIds/src
git clone https://github.com/jeremybuttler/getDIIds
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
