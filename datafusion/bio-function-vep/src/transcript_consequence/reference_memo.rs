//! Memoize only immutable transcript + translation views. The contig owns the
//! shared slots; hydrated/owned/translation-only overrides get buffer-local slots.
use super::*;
use crate::cache_budget::{BudgetLease, CacheBudget};
use std::sync::atomic::{AtomicUsize, Ordering};

pub(crate) struct ReferenceMemo {
    entry: OnceLock<ReferenceEntry>,
    budget: Arc<CacheBudget>,
    counters: Arc<ReferenceCounters>,
    #[cfg(test)]
    pub(super) derivations: AtomicUsize,
    #[cfg(test)]
    pub(super) local_codons: AtomicUsize,
}

#[derive(Default)]
struct ReferenceCounters {
    enabled: bool,
    builds: AtomicUsize,
    hits: AtomicUsize,
    declined: AtomicUsize,
    local_variants: AtomicUsize,
    local_codons: AtomicUsize,
    fallbacks: AtomicUsize,
}

pub(super) enum ReferenceEntry {
    Missing,
    Declined,
    Ready(CodingReference),
}

pub(super) struct CodingReference {
    pub(super) cds: String,
    pub(super) aas: Option<Arc<[char]>>,
    pub(super) first_stop: Option<usize>,
    pub(super) local_substitution_safe: bool,
    _lease: BudgetLease,
}

impl ReferenceMemo {
    pub(super) fn new(budget: Arc<CacheBudget>) -> Arc<Self> {
        Arc::new(Self {
            entry: OnceLock::new(),
            budget,
            counters: Arc::new(ReferenceCounters::default()),
            #[cfg(test)]
            derivations: AtomicUsize::new(0),
            #[cfg(test)]
            local_codons: AtomicUsize::new(0),
        })
    }

    pub(super) fn get(
        &self,
        tx: &TranscriptFeature,
        translation: &TranslationFeature,
    ) -> &ReferenceEntry {
        if self.counters.enabled && matches!(self.entry.get(), Some(ReferenceEntry::Ready(_))) {
            self.counters.hits.fetch_add(1, Ordering::Relaxed);
        }
        self.entry.get_or_init(|| {
            if self.budget.limit() == 0 {
                return ReferenceEntry::Declined;
            }
            #[cfg(test)]
            self.derivations.fetch_add(1, Ordering::Relaxed);
            if self.counters.enabled {
                self.counters.builds.fetch_add(1, Ordering::Relaxed);
            }
            let Some(cds) = reference_translateable_seq_for_vep(tx, Some(translation)) else {
                return ReferenceEntry::Missing;
            };
            // Failed RefSeq references may use a longer cached peptide. Reserve
            // before deriving it, and retain the charge for the slot's lifetime.
            let peptide_len =
                (cds.len() / 3).max(translation.translation_seq.as_ref().map_or(0, String::len));
            let Some(bytes) = peptide_len
                .checked_mul(std::mem::size_of::<char>())
                .and_then(|n| n.checked_add(cds.capacity()))
                .and_then(|n| {
                    n.checked_add(
                        std::mem::size_of::<CodingReference>() + 2 * std::mem::size_of::<usize>(),
                    )
                })
            else {
                return ReferenceEntry::Declined;
            };
            let Some(lease) = self.budget.admit(bytes) else {
                if self.counters.enabled {
                    self.counters.declined.fetch_add(1, Ordering::Relaxed);
                }
                return ReferenceEntry::Declined;
            };
            let aas: Option<Arc<[char]>> =
                reference_aas_for_consequence(tx, Some(translation), &cds).map(Arc::from);
            let first_stop = aas
                .as_deref()
                .and_then(|aas| aas.iter().position(|&aa| aa == '*'));
            let local_substitution_safe = !bam_edit_failed(tx)
                && !tx.has_rna_edit
                && tx.refseq_edits.is_empty()
                && tx.translation_seq_edits.is_empty()
                && cds.len().is_multiple_of(3)
                && cds.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'));
            ReferenceEntry::Ready(CodingReference {
                cds,
                aas,
                first_stop,
                local_substitution_safe,
                _lease: lease,
            })
        })
    }

    pub(super) fn note_fallback(&self) {
        if self.counters.enabled {
            self.counters.fallbacks.fetch_add(1, Ordering::Relaxed);
        }
    }
}

pub(crate) struct SharedReferenceMemos {
    // Pointers identify the immutable objects, not just their stable IDs.
    slots: HashMap<usize, (usize, Arc<ReferenceMemo>)>,
    budget: Arc<CacheBudget>,
    counters: Arc<ReferenceCounters>,
}

impl SharedReferenceMemos {
    pub(crate) fn new(
        transcripts: &[TranscriptFeature],
        translations: &[TranslationFeature],
        budget: Arc<CacheBudget>,
    ) -> Self {
        let by_tx: HashMap<_, _> = translations
            .iter()
            .map(|tl| (tl.transcript_id.as_str(), tl))
            .collect();
        let counters = Arc::new(ReferenceCounters {
            enabled: std::env::var_os("VEP_ENGINE_PROFILE").is_some_and(|value| value == "1"),
            ..Default::default()
        });
        let slots = transcripts
            .iter()
            .map(|tx| {
                let translation = by_tx
                    .get(tx.transcript_id.as_str())
                    .map_or(0, |tl| *tl as *const TranslationFeature as usize);
                (
                    tx as *const TranscriptFeature as usize,
                    (translation, Self::make_slot(&budget, &counters)),
                )
            })
            .collect();
        Self {
            slots,
            budget,
            counters,
        }
    }

    fn make_slot(
        budget: &Arc<CacheBudget>,
        counters: &Arc<ReferenceCounters>,
    ) -> Arc<ReferenceMemo> {
        Arc::new(ReferenceMemo {
            entry: OnceLock::new(),
            budget: Arc::clone(budget),
            counters: Arc::clone(counters),
            #[cfg(test)]
            derivations: AtomicUsize::new(0),
            #[cfg(test)]
            local_codons: AtomicUsize::new(0),
        })
    }

    pub(super) fn slot(
        &self,
        tx: &TranscriptFeature,
        translation: Option<&TranslationFeature>,
    ) -> Arc<ReferenceMemo> {
        let translation = translation.map_or(0, |tl| tl as *const TranslationFeature as usize);
        if let Some((base_translation, slot)) =
            self.slots.get(&(tx as *const TranscriptFeature as usize))
            && *base_translation == translation
        {
            return Arc::clone(slot);
        }
        Self::make_slot(&self.budget, &self.counters)
    }
}

impl Drop for SharedReferenceMemos {
    fn drop(&mut self) {
        if self.counters.enabled {
            let c = &self.counters;
            eprintln!(
                "[VEP_REFERENCE_CACHE] builds={} hits={} admission_failures={} local_variants={} translated_local_codons={} fallback_variants={} execution_payload_peak_bytes={} base_slots={}",
                c.builds.load(Ordering::Relaxed),
                c.hits.load(Ordering::Relaxed),
                c.declined.load(Ordering::Relaxed),
                c.local_variants.load(Ordering::Relaxed),
                c.local_codons.load(Ordering::Relaxed),
                c.fallbacks.load(Ordering::Relaxed),
                self.budget.peak(),
                self.slots.len()
            );
        }
    }
}

/// Same-length, unambiguous codons only. Peptide assembly still copies the
/// protein, but translation work is proportional to the changed codon window.
pub(super) fn substitute_peptide_window(
    table: CodonTable,
    reference: &[char],
    mutated: &[u8],
    start: usize,
    end: usize,
    memo: &ReferenceMemo,
) -> Option<Vec<char>> {
    let first = start / 3;
    let last = end / 3;
    if first > last || last >= reference.len() || reference.len() * 3 != mutated.len() {
        return None;
    }
    let mut alternate = reference.to_vec();
    if memo.counters.enabled {
        memo.counters.local_variants.fetch_add(1, Ordering::Relaxed);
        memo.counters
            .local_codons
            .fetch_add(last - first + 1, Ordering::Relaxed);
    }
    for ci in first..=last {
        let codon: [u8; 3] = mutated.get(ci * 3..ci * 3 + 3)?.try_into().ok()?;
        alternate[ci] = translate_codon_bytes_in(table, codon)?;
        #[cfg(test)]
        memo.local_codons.fetch_add(1, Ordering::Relaxed);
    }
    Some(alternate)
}
