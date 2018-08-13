#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use Carp;
use List::MoreUtils qw( zip );

=pod

    Input to this script is a table specifying the taxonomic group to which
    each OTU has been assigned. This assignment will have been done using
    USEARCH ('UCLUST' mode), typically via the Qiime script,
    parallel_assign_taxonomy_uclust.py . A typical location for such a file
    would be:

    pick_denovo_otus/uclust_assigned_taxonomy/seqs_rep_set_tax_assignments.txt

    Content of such a file is like so:

denovo11405	k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__S24-7; g__; s__	1.00	3
denovo10979	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__ 	1.00	3
denovo11151	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Dehalobacteriaceae; g__Dehalobacterium; s__	1.00	3
denovo10647	k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Rikenellaceae; g__; s__	1.00	3
denovo10646	k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__S24-7; g__; s__	1.00	3
denovo10645	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Ruminococcus]; s__gnavus   1.00    3
denovo10644	k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__S24-7; g__; s__	1.00	3
denovo10643	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Dorea; s__formicigenerans   1.00    3
denovo10642	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__	0.67	3
denovo10641	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Ruminococcus]; s__gnavus   1.00    3
denovo10640	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__	0.67	3
denovo11401	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Ruminococcus; s__   0.67    3
denovo10649	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__	0.67	3
denovo10648	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Oscillospira; s__   1.00    3
denovo11399	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Ruminococcus]; s__gnavus   1.00    3
denovo11235	k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__S24-7; g__; s__	1.00	3
denovo10908	k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Salinicoccus; s__  	1.00	3
denovo10909	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__ 	1.00	3
denovo8664	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae	0.67	3
denovo9048	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales	1.00	3

...
...
etc.

Note that there are 4 columns, tab-delimited.

Note that most taxonomies extend to species (s__), even if in some cases
they are anonymous species. Contrast that with the last two lines, which
extend only as far as family and order respectively.

It's not completely clear what the last 2 columns represent, but my
interpretation is that they are respectively:
the "fraction of database hits that ... have a specific taxonomic assignment"
the number of database hits.

That's because all the 3rd-column values are either 1.00 or 0.67, and the
manual for the Qiime script parallel_assign_taxonomy_uclust.py

	http://qiime.org/scripts/parallel_assign_taxonomy_uclust.html

says:

--min_consensus_fraction
    Minimum fraction of database hits that must have a specific taxonomic assignment to assign that taxonomy to a query [default: 0.51]

The number of database hits in each case is apparent in the USEARCH output
(USEARCH would have been run by file that Qiime script); note that in a
typical Qiime pipeline with default behaviour, USEARCH would have first been
used to cluster the reads into OTUs. Then, of interest here, USEARCH would
have been used to compare a representative sequence from each OTU with a
reference database of sequences with taxonomic assignments (in this task,
the reference sequences are in fact the 'library seeds' and the OTU
representative reads are the 'hits').

The USEARCH output from this comparison is in a 'log' file in the same dir
as the assigned-taxonomy file (produced from it). E.g. in this case the
log file would be:

    pick_denovo_otus/uclust_assigned_taxonomy/seqs_rep_set_tax_assignments.log

It seems that some kind of limit of a hit (ref OTU seq) appearing no more than
3 times is in place?

E.g. take this OTU:

$ grep 'denovo1 ' seqs_rep_set_tax_assignments.log 
H	30306	252	96.8	+	0	0	522I252M619I	denovo1 Ms23_6129132	276046
H	30376	252	96.8	+	0	0	522I252M615I	denovo1 Ms23_6129132	275627
H	30336	252	96.8	+	0	0	521I252M618I	denovo1 Ms23_6129132	275882

- many other OTUs appear 3 times, but incidences of 1 and 2 also occur; but no
more than 3. So, each time they appear, the context is a hit to a different
lib seed.

It's maybe counterintuitive that the reference sequences are the 'seeds' not
the hits (because otherwise, i.e. if the OTU rep seqs were the seeds and the
reference database sequences the hits, then an upper limit of 3 on the hits
returned from each query (seed) could be imposed. But here, it appears that
each OTU can appear as a hit no more than 3 times.

That's indeed the case, as the beginning of the log file stats the arguments
used, which is e.g.:

UclustConsensusTaxonAssigner parameters:
id_to_taxonomy_filepath:/usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
id_to_taxonomy_fp:/usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
max_accepts:3
min_consensus_fraction:0.51
reference_sequences_fp:/usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta
similarity:0.9
unassignable_label:Unassigned
Result path: /tmp/assign-taxDbGdTc

- note max_accepts (and min_consensus_fraction).

=cut

my $assignment_file = shift or croak;

my @level_names = qw( kingdom phylum class order family genus species );
my @level_keys  = map { substr $_, 0, 1 } @level_names;
my %name_of_level = zip @level_keys, @level_names;
$name_of_level{b} = q{binomial};

my %level_frq_all; # keys are single-letter codes for taxonomic levels
               # values are frequencies of lines (OTUs) which have
               # a taxonomy specified at this level (or deeper); this
               # includes cases which are specified as anonymous (e.g.
               # taxon at that level is blank).

my %level_frq_explicit; # similar to %level_frq_all, except that blank
               # taxa are treated as unspecified, so are NOT included

my %cumul_taxon_level_frq; # keys are single-letter codes for taxonomic levels
               # values are refs to hashes, whose keys are 'cumulative'
               # (whole or partial) taxonomy strings, e.g.
               # 'k__Bacteria'
               # 'k__Bacteria; p__Firmicutes'
               # 'k__Bacteria; p__Firmicutes; c__Clostridia'
               # 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales'
               # 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__'
               # 'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__'
               # - so does include those with anonymous taxa at one or more
               # levels; values of these hashes are frequencies

my %single_taxon_level_frq; # similar to %cumul_taxon_level_frq, except
              # the keys to the 'inner' hashes are single taxon names, without
              # the 'p__' etc prefix, and blanks are replaced with 'ANONYMOUS',
              # e.g. equivalent to above:
              # 'Bacteria'
              # 'Firmicutes'
              # 'Clostridia'
              # 'Clostridiales'
              # 'Lachnospiraceae'
              # 'ANONYMOUS'
              # 'ANONYMOUS'
              #
              # However, it's also of utility to have stats for individual
              # species; so there is an additional 'level' implemented, with
              # the code 'b' for 'binomial'; this concatenates the g__ and
              # s__ strings with intervening whitespace. So some will be
              # e.g. 'Ruminococcus ANONYMOUS', others 'ANONYMOUS ANONYMOUS'.
              #
              # My assumption is that the taxonomies come from the reference
              # database (the 'seed' sequences); so, that's the origin of the
              # blank/anonymous species, genera etc. Whereas the depth of the
              # taxonomic specification may come from the comparison process?
              # I.e. if the seed-hit relationships mean that no inference
              # can be made beyond family, then genus and species won't be
              # blank, but missing entirely.

open my $fh, '<', $assignment_file;

my $n_otus = 0;

while (defined (my $line = <$fh>) ) {

    next if $line !~ m{ \S }xms;

    chomp $line;

    my ($otu_id, $taxonomy, $prpn_matches_with_taxonomy, $n_matches) =
        split m{ \t }xms, $line;

    pos $taxonomy = 0;

    my $partial_taxonomy = '';

    my %this_taxonomy; # just for this OTU; keys single-letter codes; values
                       # are single-word taxon names

    while ($taxonomy =~ m{ \G ( ([a-z]) [_][_] ([^\s;]*) ) (?: [;] \s* | \z ) }gcxms) {
        my ($taxon_string, $level, $taxon) = ($1, $2, $3);
        ###print qq{'$taxon_string', '$level', '$taxon'\n};

        my $taxon_key = $taxon;
        $taxon_key = 'ANONYMOUS' if !$taxon_key;

        $level_frq_all{$level}++;
        $level_frq_explicit{$level}++ if $taxon;

        $partial_taxonomy .= '; ' if $partial_taxonomy;
        $partial_taxonomy .= $taxon_string;

        $cumul_taxon_level_frq{$level}{$partial_taxonomy}++;

        $single_taxon_level_frq{$level}{$taxon_key}++;

        if (lc $level eq 's') {
            my $binomial = qq{$this_taxonomy{g} $taxon_key};
            $single_taxon_level_frq{b}{$binomial}++;
        }

        $this_taxonomy{$level} = $taxon_key;

    }

    $n_otus++;
}

close $fh;

$level_frq_all{b}      = $level_frq_all{s};
$level_frq_explicit{b} = $level_frq_explicit{s};

print qq{OTUS:\t$n_otus\n};

for my $level (@level_keys) {

    printf qq{SUMMARY_ALL     : %6d of %6d OTUs ( %5.1f%s ) are specified to at least %s (including anonymous assignments)\n},
        $level_frq_all{$level}       , $n_otus,
        100 * $level_frq_all{$level} / $n_otus, q{%},
        $name_of_level{$level};
}

for my $level (@level_keys) {

    printf qq{SUMMARY_EXPLICIT: %6d of %6d OTUs ( %5.1f%s of all OTUs, and %5.1f%s of %7s-specified OTUs) are explicitly specified to at least %s (excluding anonymous assignments)\n},
        $level_frq_explicit{$level}       , $n_otus,
        100 * $level_frq_explicit{$level} / $n_otus, q{%},
        100 * $level_frq_explicit{$level} / $level_frq_all{$level}, q{%},
        $name_of_level{$level}, $name_of_level{$level};
}

display_taxon_frq(\%single_taxon_level_frq, 'SINGLE_TAXON_FREQUENCY', '', 'b');
display_taxon_frq(\%cumul_taxon_level_frq, 'TAXONOMY_FREQUENCY');

exit 0;

sub display_taxon_frq {

    my $taxon_frq_href = shift;
    my $taxonomy_type  = shift;
    my $width          = shift || 10;
    my @extra_levels   = @_;

    for my $level (@level_keys, @extra_levels) {

        my %frq_of = %{$taxon_frq_href->{$level}};

        my @taxa = sort { $frq_of{$b} <=> $frq_of{$a} } (keys %frq_of);

        for my $taxon (@taxa) {

            printf qq{%s: %7s %}.$width.qq{s %5d OTUs; %6.2f%s of all OTUs; }.
                   qq{%6.2f%s of OTUs assigned to at least this level; }.
                   qq{%6.2f%s of OTUs explicitly assigned to at least this level\n},
                $taxonomy_type, $name_of_level{$level}, $taxon,
                $taxon_frq_href->{$level}{$taxon},
                100 * $taxon_frq_href->{$level}{$taxon} / $n_otus,                     q{%},
                100 * $taxon_frq_href->{$level}{$taxon} / $level_frq_all{$level},      q{%},
                100 * $taxon_frq_href->{$level}{$taxon} / $level_frq_explicit{$level}, q{%},

        }
        $width += 5;

    }

}

=pod

my $width = 10;

for my $level (@level_keys, 'b') {

    my %frq_of = %{$single_taxon_level_frq{$level}};

    my @taxa = sort { $frq_of{$b} <=> $frq_of{$a} } (keys %frq_of);

    for my $taxon (@taxa) {

        printf qq{TAXONOMY_SINGLE_FREQUENCY: %7s %}.$width.qq{s %5d OTUs; %6.2f%s of all OTUs; }.
               qq{%6.2f%s of OTUs assigned to at least this level; }.
               qq{%6.2f%s of OTUs explicitly assigned to at least this level\n},
            $name_of_level{$level}, $taxon,
            $single_taxon_level_frq{$level}{$taxon},
            100 * $single_taxon_level_frq{$level}{$taxon} / $n_otus,                     q{%},
            100 * $single_taxon_level_frq{$level}{$taxon} / $level_frq_all{$level},      q{%},
            100 * $single_taxon_level_frq{$level}{$taxon} / $level_frq_explicit{$level}, q{%},

    }
    $width += 5;

}

