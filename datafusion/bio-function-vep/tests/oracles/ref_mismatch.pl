use strict;
use warnings;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;

# Minimal coordinate/cache objects. The codon, display_codon and alternate-CDS
# methods under test are the unmodified Ensembl methods installed in Docker.
{ package ProbeAttribute;
  sub code { $_[0]->{code} }
  sub value { $_[0]->{value} }
}
{ package ProbeTranscript;
  sub strand { $_[0]->{strand} // 1 }
  sub get_all_Attributes { $_[0]->{attrs} }
}
{ package ProbeVariant;
  sub strand { 1 }
}
{ package ProbeOverlap;
  sub transcript { $_[0]->{transcript} }
  sub variation_feature { bless {}, 'ProbeVariant' }
  sub translation_start { 2 }
  sub translation_end { int(($_[0]->cds_end + 2) / 3) }
  sub cds_start { $_[0]->{start} // 5 }
  sub cds_end { $_[0]->{end} // 5 }
  sub codon_position { ($_[0]->cds_start - 1) % 3 + 1 }
  sub _translateable_seq { 'ATGGCTGAATGA' }
  sub _codon_table { 1 }
  sub _seq_edits { [] }
  sub _three_prime_utr { undef }
  sub get_reference_TranscriptVariationAllele { $_[0]->{reference} }
}

sub genomic {
  my ($seq, $strand) = @_;
  return $seq if $strand == 1;
  $seq = reverse $seq;
  $seq =~ tr/ACGT/TGCA/;
  return $seq;
}

for my $strand (1, -1) {
  for my $case ([5, 'A', 'T'], [5, 'A', 'C'], [4, 'AA', 'TT'],
                [4, 'A', '-'], [4, 'AAA', '-'], [4, 'AA', 'T'],
                [4, 'AAA', 'TTTTTT'], [4, 'AAAA', 'T']) {
    my ($start, $ref_seq, $alt_seq) = @$case;
    my $tr = bless {attrs => [], strand => $strand}, 'ProbeTranscript';
    my $tv = bless {transcript => $tr, start => $start,
                   end => $start + length($ref_seq) - 1}, 'ProbeOverlap';
    my $ref = bless {base_variation_feature_overlap => $tv, is_reference => 1,
      variation_feature_seq => genomic($ref_seq, $strand), feature_seq => $ref_seq},
      'Bio::EnsEMBL::Variation::TranscriptVariationAllele';
    $tv->{reference} = $ref;
    my $alt = bless {base_variation_feature_overlap => $tv, is_reference => 0,
      variation_feature_seq => genomic($alt_seq, $strand), feature_seq => $alt_seq},
      'Bio::EnsEMBL::Variation::TranscriptVariationAllele';
    print join("\t", $strand, $start, $ref_seq, $alt_seq,
      $alt->display_codon_allele_string, $alt->pep_allele_string), "\n";
  }
}

for my $case (
  ['ordinary', undef, undef, 'GCT'],
  ['bam_state_only', undef, undef, 'GCT'],
  ['rna_edit', '_rna_edit', '5 5 A', 'GAT'],
  ['polya_edit', '_rna_edit', '13 12 AAAAAAAAAAAAAAAAAAAA', 'GAT'],
) {
  my ($name, $code, $value, $expected) = @$case;
  my @attrs = defined($code) ? (bless({code => $code, value => $value}, 'ProbeAttribute')) : ();
  my $tr = bless {attrs => \@attrs}, 'ProbeTranscript';
  $tr->{bam_edit_status} = 'ok' unless $name eq 'ordinary';
  my $tv = bless {transcript => $tr}, 'ProbeOverlap';
  my $ref = bless {
    base_variation_feature_overlap => $tv,
    is_reference => 1,
    variation_feature_seq => 'A',
    feature_seq => 'A',
  }, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele';
  $tv->{reference} = $ref;
  my $alt = bless {
    base_variation_feature_overlap => $tv,
    is_reference => 0,
    variation_feature_seq => 'T',
    feature_seq => 'T',
  }, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele';
  my $ref_codon = $ref->codon;
  my $alt_codon = $alt->codon;
  die "$name: expected $expected/GTT, got $ref_codon/$alt_codon\n"
    unless $ref_codon eq $expected && $alt_codon eq 'GTT';
  print join("\t", $name, $ref_codon, $alt_codon,
    $ref->display_codon, $alt->display_codon), "\n";
}
