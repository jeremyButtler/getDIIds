# Use:

Detects DI (defective influenza) reads in long read
  sequenced samples.

This requires full length amplicons, so the old ONT rapid
  kits that fragmented DNA would not work. Also, short
  reads will not work.

# Install:

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
  of the TB program. This is the frezeTB project that
  contributed the source code for the primer
  searching step.
- Family for being there.
